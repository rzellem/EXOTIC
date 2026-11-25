"""Refit saved FinalLightCurve CSVs with Claret and hybrid4 on identical data.

Run from this checkout. Each --case supplies an existing inits.json and a
FinalLightCurve CSV. All outputs remain in the requested local directory.
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import re
import sys
from time import perf_counter

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from exotic.api.elca import lc_fitter, transit_duration
from exotic.api.hybrid_comparison import comparison_for_output, json_value
from exotic.api.ld import LimbDarkening
from exotic.exotic import (
    build_initial_transit_bounds, build_search_restriction_prior_from_planet_dict,
    build_single_transit_duration_prior, enrich_search_restriction_prior_with_rprs_data_uncertainty,
    estimate_ephemeris_tmid_and_bounds, nonlinear_ld,
)


PLANET_FIELDS = {
    'pName': 'Planet Name', 'sName': 'Host Star Name',
    'pPer': 'Orbital Period (days)', 'pPerUnc': 'Orbital Period Uncertainty',
    'rprs': 'Ratio of Planet to Stellar Radius (Rp/Rs)',
    'rprsUnc': 'Ratio of Planet to Stellar Radius (Rp/Rs) Uncertainty',
    'aRs': 'Ratio of Distance to Stellar Radius (a/Rs)',
    'aRsUnc': 'Ratio of Distance to Stellar Radius (a/Rs) Uncertainty',
    'inc': 'Orbital Inclination (deg)', 'incUnc': 'Orbital Inclination (deg) Uncertainty',
    'ecc': 'Orbital Eccentricity (0 if null)', 'omega': 'Argument of Periastron (deg)',
    'midTUnc': 'Mid-Transit Time Uncertainty',
    'teff': 'Star Effective Temperature (K)',
    'teffUncPos': 'Star Effective Temperature (+) Uncertainty',
    'teffUncNeg': 'Star Effective Temperature (-) Uncertainty',
    'met': 'Star Metallicity ([FE/H])',
    'metUncPos': 'Star Metallicity (+) Uncertainty',
    'metUncNeg': 'Star Metallicity (-) Uncertainty',
    'logg': 'Star Surface Gravity (log(g))',
    'loggUncPos': 'Star Surface Gravity (+) Uncertainty',
    'loggUncNeg': 'Star Surface Gravity (-) Uncertainty',
}


def run_case(inits_path, photometry_path, output):
    config = json.loads(inits_path.read_text(encoding='utf-8-sig'))
    inputs = config['planetary_parameters']
    planet = {key: inputs.get(label, inputs.get(key)) for key, label in PLANET_FIELDS.items()}
    optional = config.get('optional_info', {})
    user = config.get('user_info', {})
    # Read only the existing reduced photometry, not the science-image stack.
    values = np.loadtxt(photometry_path, comments='#', delimiter=',', usecols=(0, 1, 2, 3, 4, 5))
    times, phase, flux, error, _, airmass = values.T
    # The saved phase and BJD_TDB retain the original fit's epoch convention.
    # No legacy UTC ephemeris is silently reinterpreted as TDB here.
    planet['midT'] = float(np.median(times - phase * planet['pPer']))
    dates = re.findall(r'\d{4}-\d{2}-\d{2}', photometry_path.name)
    date = dates[-1] if dates else 'saved-lightcurve'
    info = {
        'save': str(output), 'date': date,
        'filter': user.get('Filter Name (aavso.org/filters)', optional.get('Filter Name (aavso.org/filters)')),
        'wl_min': optional.get('Filter Minimum Wavelength (nm)'),
        'wl_max': optional.get('Filter Maximum Wavelength (nm)'),
    }
    ld = LimbDarkening(planet)
    nonlinear_ld(ld, info, non_interactive_run=True)
    prior = {
        'rprs': planet['rprs'], 'ars': planet['aRs'], 'per': planet['pPer'],
        'inc': planet['inc'], 'ecc': planet['ecc'], 'omega': planet['omega'],
        'tmid': planet['midT'], 'a0': 1.0, 'a2': 0.0,
        **{f'u{i}': coefficient[0] for i, coefficient in enumerate((ld.ld0, ld.ld1, ld.ld2, ld.ld3))},
    }
    timing = estimate_ephemeris_tmid_and_bounds(
        times, prior['tmid'], prior['per'], planet['midTUnc'], planet['pPerUnc'],
        expected_duration=transit_duration(prior), sigma_multiplier=35.0,
    )
    prior['tmid'] = timing['tmid']
    restrictions = enrich_search_restriction_prior_with_rprs_data_uncertainty(
        build_search_restriction_prior_from_planet_dict(planet), times, flux, error, prior,
        context_label='saved photometry Claret/hybrid4 comparison',
    )
    bounds = build_initial_transit_bounds(
        prior, timing['bounds'], ars_unc=planet['aRsUnc'],
        rprs_data_uncertainty=restrictions.get('rprs_data_uncertainty'),
        search_restriction_prior=restrictions,
    )
    exposure = optional.get('Exposure Time (s)')
    live_points = optional.get('minimum number of live points for ultranest')
    print(f'\nFitting {planet["pName"]}: {len(times)} retained points, {info["filter"]} band', flush=True)
    print('Shared bounds:', bounds, flush=True)
    started = perf_counter()
    fit = lc_fitter(times, flux, error, airmass, prior, bounds, mode='ns',
                    exposure_times_seconds=exposure,
                    duration_prior=build_single_transit_duration_prior(planet),
                    ultranest_min_num_live_points=live_points)
    fit.fit_runtime_seconds = perf_counter() - started
    report = comparison_for_output(fit, info, planet['pName'], date)
    if report is None or report.get('status') != 'completed':
        raise RuntimeError(f'Comparison did not complete: {report}')
    report['validation_inputs'] = {
        'inits': str(inits_path.resolve()), 'photometry': str(photometry_path.resolve()),
        'photometry_kind': 'Previously selected, detrended FinalLightCurve CSV; no new clipping or detrending',
        'timing_initialization': 'Epoch reconstructed from the saved BJD_TDB and orbital-phase columns',
        'exposure_seconds': exposure,
        'exposure_note': 'Recorded constant exposure from inits' if exposure is not None else
                         'No exposure time recorded in inits; both replay fits use instantaneous fluxes',
    }
    # Add the replay provenance to the full sidecar without duplicating profiles
    # in the compact aggregate report.
    sidecar = output / report['files']['json']
    full_report = json.loads(sidecar.read_text(encoding='utf-8'))
    full_report['validation_inputs'] = report['validation_inputs']
    sidecar.write_text(json.dumps(json_value(full_report), indent=2, allow_nan=False), encoding='utf-8')
    print('Completed:', sidecar, flush=True)
    print('Posterior changes:', report['differences'], flush=True)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', nargs=2, action='append', required=True,
                        metavar=('INITS_JSON', 'FINAL_LIGHTCURVE_CSV'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    results = []
    for index, (inits, photometry) in enumerate(args.case):
        np.random.seed(20261120 + index)
        try:
            result = run_case(Path(inits), Path(photometry), args.output)
        except Exception as exc:
            logging.exception('Case failed: %s', photometry)
            result = {'status': 'failed', 'inits': inits, 'photometry': photometry,
                      'error_type': type(exc).__name__, 'reason': str(exc)}
        results.append(result)
        (args.output / 'comparison_summary.json').write_text(
            json.dumps(json_value(results), indent=2, allow_nan=False), encoding='utf-8')
    return 1 if any(row.get('status') != 'completed' for row in results) else 0


if __name__ == '__main__':
    raise SystemExit(main())
