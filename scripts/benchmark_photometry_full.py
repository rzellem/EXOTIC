"""Time the full reduction with explicit native/reference photometry.

Science inputs and settings are passed unchanged to profile_full_reduction.
Outputs include every aperture-grid evaluation and final sampler inputs/draws
for numerical comparisons independent of adaptive pipeline decisions.
"""
import argparse
import copy
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--backend', choices=('cpp', 'python'), required=True)
    parser.add_argument('--likelihood-backend', choices=('auto', 'cpp', 'python'), default='auto')
    parser.add_argument('--preparation-backend', choices=('cpp', 'python'), default='cpp')
    parser.add_argument('--cprofile', action='store_true',
                        help='Include per-call profiling overhead for diagnostic profiles')
    parser.add_argument('--inits', type=Path, required=True)
    parser.add_argument('--science-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    import numpy as np
    import exotic.exotic as runtime
    from exotic.api import elca
    from exotic.api import hybrid_comparison
    from exotic.api.hybrid_comparison import fit_summary, json_value
    from scripts import profile_full_reduction as profiling

    runtime._PHOTOMETRY_BACKEND = args.backend
    runtime.bad_pixel_processing.PREPARATION_BACKEND = args.preparation_backend
    runtime._prepare_photometry_backend()
    grids, fits, scouts, fit_objects = [], [], [], []
    original_grid = runtime.compute_star_aperture_grid
    original_init = elca.lc_fitter.__init__
    original_scout = runtime.prepare_ultranest_bounds_with_lm_boundary_scout
    original_hybrid = hybrid_comparison.run_hybrid4_fit

    def save_trace():
        trace = dict(photometry_backend=args.backend, claret_backend='cpp',
                     orbital_replay_backend='cpp', likelihood_backend=args.likelihood_backend,
                     image_preparation_backend=args.preparation_backend,
                     plotting_uncertainty_reuse=True, fits=fits, scouts=scouts)
        (args.output/'bound_handoff_trace.json').write_text(
            json.dumps(json_value(trace), indent=2, allow_nan=False), encoding='utf-8')

    def grid(*a, **kw):
        result = original_grid(*a, **kw)
        grids.append(dict(star=int(a[1]), x=float(a[2]), y=float(a[3]),
            radii=np.asarray(a[4]), annuli=np.asarray(a[5]),
            flux=np.asarray(result[0]), background=np.asarray(result[1]),
            noise=result[2] if len(result)>2 else None))
        return result

    def initialize(self, *a, **kw):
        kw['likelihood_backend'] = args.likelihood_backend
        initial = copy.deepcopy(a[4] if len(a)>4 else kw['prior'])
        mode = kw.get('mode', a[7] if len(a)>7 else 'ns')
        original_init(self, *a, **kw)
        if mode != 'ns':
            return
        index = len(fits)
        fits.append(dict(index=index, points=len(self.time), initial_prior=initial,
                         actual_bounds=copy.deepcopy(self.bounds),
                         sample_bounds=copy.deepcopy(self.sample_bounds), summary=fit_summary(self)))
        fit_objects.append(self)
        checkpoint = dict(time=self.time, data=self.data, dataerr=self.dataerr, airmass=self.airmass,
            prior=initial, bounds=self.bounds, settings={k:v for k,v in kw.items() if k in (
                'fixed_flux_baseline', 'fixed_parameter_errors', 'duration_prior',
                'exposure_times_seconds', 'baseline_fit_mask', 'ultranest_min_num_live_points',
                'use_impactparameter_rather_than_inclination_to_fit', 'likelihood_backend')})
        (args.output/f'fit_input_{index}.json').write_text(
            json.dumps(json_value(checkpoint), indent=2), encoding='utf-8')
        weighted = (self.results or {}).get('weighted_samples', {})
        np.savez(args.output/f'fit_posterior_{index}.npz', points=weighted.get('points', []),
                 weights=weighted.get('weights', []), logl=weighted.get('logl', []),
                 sampled_keys=np.asarray(self.sampled_keys))
        save_trace()

    def hybrid(*a, **kw):
        result = original_hybrid(*a, **kw)
        record_hybrid(result)
        return result

    def record_hybrid(result):
        index = len(fits)
        fit_objects.append(result)
        fits.append(dict(index=index, points=len(result.time), initial_prior=copy.deepcopy(result.prior),
                         actual_bounds=copy.deepcopy(result.bounds),
                         sample_bounds=copy.deepcopy(result.sample_bounds), summary=fit_summary(result)))
        exposure = getattr(result, 'exposure_times_days', None)
        settings = dict(
            likelihood_backend=args.likelihood_backend,
            fixed_flux_baseline=result.fixed_flux_baseline,
            fixed_parameter_errors=result.fixed_parameter_errors, duration_prior=result.duration_prior,
            baseline_fit_mask=result.baseline_fit_mask,
            use_impactparameter_rather_than_inclination_to_fit=result.use_impactparameter_rather_than_inclination_to_fit,
            exposure_times_seconds=None if exposure is None else np.asarray(exposure)*86400.,
            exposure_smearing_supersample=result.exposure_smearing_supersample,
            exposure_smearing_change_tolerance=result.exposure_smearing_change_tolerance,
        )
        checkpoint = dict(time=result.time, data=result.data, dataerr=result.dataerr, airmass=result.airmass,
                          prior=result.prior, bounds=result.bounds, settings=settings)
        (args.output/f'fit_input_{index}.json').write_text(
            json.dumps(json_value(checkpoint), indent=2), encoding='utf-8')
        save_trace()

    def scout(*a, **kw):
        initial_prior, initial_bounds = copy.deepcopy(a[4]), copy.deepcopy(a[5])
        result = original_scout(*a, **kw)
        scouts.append(dict(points=len(a[0]), initial_prior=initial_prior, initial_bounds=initial_bounds,
                           returned_bounds=copy.deepcopy(result[1]), payload=copy.deepcopy(result[2])))
        save_trace()
        return result

    runtime.compute_star_aperture_grid = grid
    runtime.prepare_ultranest_bounds_with_lm_boundary_scout = scout
    elca.lc_fitter.__init__ = initialize
    hybrid_comparison.run_hybrid4_fit = hybrid
    old_argv = sys.argv
    sys.argv = ['profile_full_reduction.py', '--inits', str(args.inits),
                '--science-dir', str(args.science_dir), '--output', str(args.output)]
    if not args.cprofile:
        sys.argv.append('--no-cprofile')
    save_trace()
    try:
        result = profiling.main()
    finally:
        sys.argv = old_argv
        runtime.compute_star_aperture_grid = original_grid
        runtime.prepare_ultranest_bounds_with_lm_boundary_scout = original_scout
        elca.lc_fitter.__init__ = original_init
        hybrid_comparison.run_hybrid4_fit = original_hybrid
        # Legacy output_files imports can load api.hybrid_comparison separately
        # from exotic.api.hybrid_comparison. The final fit owns that secondary
        # result regardless of which import path ran it.
        for obj in tuple(fit_objects):
            secondary = getattr(obj, 'hybrid4_fit', None)
            if secondary is not None and all(secondary is not known for known in fit_objects):
                record_hybrid(secondary)
        # Retained samplers can be extended after construction. Save the final
        # draws/state for every candidate and the secondary hybrid4 fit.
        for record, obj in zip(fits, fit_objects):
            record.update(actual_bounds=copy.deepcopy(obj.bounds),
                          sample_bounds=copy.deepcopy(obj.sample_bounds), summary=fit_summary(obj))
            weighted = (obj.results or {}).get('weighted_samples', {})
            np.savez(args.output/f"fit_posterior_{record['index']}.npz",
                     points=weighted.get('points', []), weights=weighted.get('weights', []),
                     logl=weighted.get('logl', []), sampled_keys=np.asarray(obj.sampled_keys))
        arrays = {}
        metadata = []
        for i, record in enumerate(grids):
            row = {k:v for k,v in record.items() if k in ('star', 'x', 'y')}
            for key in ('radii', 'annuli', 'flux', 'background'):
                arrays[f'{i}_{key}'] = record[key]
            for key, value in (record['noise'] or {}).items():
                arrays[f'{i}_noise_{key}'] = value
            metadata.append(row)
        np.savez_compressed(args.output/'aperture_grids.npz', **arrays)
        (args.output/'aperture_grid_index.json').write_text(json.dumps(metadata), encoding='utf-8')
        save_trace()
    return result


if __name__ == '__main__':
    raise SystemExit(main())
