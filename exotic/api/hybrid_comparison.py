"""Paired Claret/hybrid4 inference on one selected photometry series."""

import copy
import csv
import json
import logging
from pathlib import Path

import numpy as np

from .hybrid_ld import HYBRID4_KEYS, LIMB_DARKENING_LAW_KEY

log = logging.getLogger(__name__)


def json_value(value):
    """Portable JSON: unavailable numeric diagnostics are null, never NaN."""
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_value(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, Path):
        return str(value)
    return value


def fit_summary(fit):
    """Keep full posterior intervals separate from EXOTIC's reported errors."""
    try:
        from ..transit_depth import observable_depth_percent
    except ImportError:
        from transit_depth import observable_depth_percent
    results = getattr(fit, 'results', None) or {}
    parameters = dict(fit.parameters)
    keys = list(getattr(fit, 'sampled_keys', fit.bounds))
    posterior = copy.deepcopy(getattr(fit, 'posterior_summaries', {}) or {})
    sampled = results.get('posterior', {})
    maximum = results.get('maximum_likelihood', {}).get('point')
    physical_ml = (
        fit._physical_values_from_sample_point(maximum) if maximum is not None else parameters
    )
    residual = np.asarray(fit.data) - np.asarray(fit.model)
    return json_value({
        'law': parameters.get(LIMB_DARKENING_LAW_KEY, 'claret'),
        'backend': (parameters.get('_hybrid4_backend', 'numpy')
                    if parameters.get(LIMB_DARKENING_LAW_KEY) == 'hybrid4'
                    else parameters.get('_claret_backend', 'cpp')),
        'sampled_parameters': keys,
        'bounds': fit.bounds,
        'posterior': posterior,
        'raw_sampled_posterior_mean': dict(zip(keys, sampled.get('mean', []))),
        'raw_sampled_posterior_stdev': dict(zip(keys, sampled.get('stdev', []))),
        'raw_maximum_likelihood_parameters': physical_ml,
        'reported_parameters': parameters,
        'reported_errors': getattr(fit, 'errors', {}),
        'posterior_error_adjustments': getattr(fit, 'ultranest_error_fallbacks', {}),
        'radius_prior_used_for_reporting': getattr(fit, 'rprs_prior_fallback_applied', False),
        'chi2_of_reported_model': float(np.sum((residual / fit.dataerr) ** 2)),
        'residual_rms_ppm': float(np.std(residual) / np.median(fit.data) * 1e6),
        'observable_depth_percent': observable_depth_percent(parameters)[0],
        'duration_days': getattr(fit, 'duration_expected', None),
        'logz': results.get('logz'),
        'logzerr': results.get('logzerr'),
        'ncall': results.get('ncall'),
        'sampling_budget': results.get('exotic_call_limit', {}).get('automatic_budget'),
        'sampling_call_control': results.get('exotic_call_limit'),
        'sampling_call_limit_reached': results.get('exotic_call_limit', {}).get('reached'),
        'runtime_seconds': getattr(fit, 'fit_runtime_seconds', None),
        'backend_setup_seconds': getattr(fit, 'backend_setup_seconds', None),
        'likelihood_backend': getattr(fit, 'likelihood_backend', None),
    })


def run_hybrid4_fit(reference, atmosphere, backend='cpp'):
    """Reuse the final fit inputs and priors; change only the intensity law.

    No clipping, comparison-star search, OOT selection, uncertainty rescaling,
    prior expansion, or warm start is repeated for the secondary fit.
    """
    from .elca import lc_fitter
    if not isinstance(atmosphere, dict):
        raise ValueError('The existing Claret fit has no shared stellar-atmosphere profile.')
    if atmosphere.get('status') == 'failed':
        raise ValueError(f'Atmosphere profile fitting failed: {atmosphere.get("reason")}')
    coefficients = np.asarray(atmosphere['coefficients'], dtype=float)
    if coefficients.shape != (4,) or not np.all(np.isfinite(coefficients)):
        raise ValueError('Four finite atmosphere-derived hybrid4 coefficients are required.')
    if any(key in reference.bounds for key in ('u0', 'u1', 'u2', 'u3')):
        raise ValueError('Freely sampled Claret coefficients need explicit hybrid4 coefficient priors; '
                         'Claret bounds do not define priors for the hybrid4 shape weights.')
    prior = copy.deepcopy(reference.prior)
    prior[LIMB_DARKENING_LAW_KEY] = 'hybrid4'
    prior['_hybrid4_backend'] = backend
    prior.update(zip(HYBRID4_KEYS, coefficients))
    exposures = getattr(reference, 'exposure_times_days', None)
    if exposures is not None:
        exposures = np.asarray(exposures, dtype=float) * 86400.0
    fit = lc_fitter(
        np.array(reference.time, copy=True), np.array(reference.data, copy=True),
        np.array(reference.dataerr, copy=True), np.array(reference.airmass, copy=True),
        prior, copy.deepcopy(reference.bounds), mode=None,
        neighbors=getattr(reference, 'neighbors', 200),
        jd_times=copy.deepcopy(getattr(reference, 'jd_times', None)),
        verbose=getattr(reference, 'verbose', True),
        use_impactparameter_rather_than_inclination_to_fit=getattr(
            reference, 'use_impactparameter_rather_than_inclination_to_fit', True),
        duration_prior=copy.deepcopy(getattr(reference, 'duration_prior', None)),
        baseline_fit_mask=copy.deepcopy(getattr(reference, 'baseline_fit_mask', None)),
        fixed_parameter_errors=copy.deepcopy(getattr(reference, 'fixed_parameter_errors', {})),
        fixed_flux_baseline=getattr(reference, 'fixed_flux_baseline', False),
        ultranest_min_num_live_points=getattr(reference, 'ultranest_min_num_live_points', None),
        exposure_times_seconds=exposures,
        exposure_smearing_supersample=getattr(reference, 'exposure_smearing_supersample', 7),
        exposure_smearing_change_tolerance=getattr(reference, 'exposure_smearing_change_tolerance', 1e-5),
        likelihood_backend=getattr(reference, 'likelihood_backend_requested', 'auto'),
    )
    fit.mode = 'ns'
    fit.max_ncalls = reference.max_ncalls
    fit.airmass_reference = reference._get_airmass_reference()
    fit.plot_time_range = copy.deepcopy(getattr(reference, 'plot_time_range', None))
    # Both fits hold the atmosphere coefficients fixed. Their profile-fit
    # covariance is reported separately, not treated as four independent errors.
    fit.hybrid4_atmosphere = copy.deepcopy(atmosphere)
    log.warning('Running hybrid4 (%s) UltraNest comparison on the same %d photometry points.', backend, len(fit.time))
    fit.fit_nested()
    return fit


def comparison_payload(reference, hybrid, atmosphere):
    claret = fit_summary(reference)
    alternative = fit_summary(hybrid)
    deltas = {}
    for key in sorted(set(claret['posterior']) & set(alternative['posterior'])):
        a, b = claret['posterior'][key], alternative['posterior'][key]
        if a.get('median') is not None and b.get('median') is not None:
            delta = b['median'] - a['median']
            deltas[key] = {'posterior_median_difference': delta}
            if key == 'tmid':
                deltas[key]['posterior_median_difference_seconds'] = delta * 86400.0
    if claret['logz'] is not None and alternative['logz'] is not None:
        deltas['logz_hybrid4_minus_claret'] = alternative['logz'] - claret['logz']
    return json_value({
        'status': 'completed',
        'sampler': 'UltraNest for both models',
        'comparison': 'hybrid4 minus Claret',
        'photometry': {
            'point_count': len(reference.time),
            'same_time_flux_uncertainty_airmass_and_exposures': True,
            'baseline_fit_mask': getattr(reference, 'baseline_fit_mask', None),
        },
        'shared_fit_settings': {
            'bounds': reference.bounds,
            'duration_prior': getattr(reference, 'duration_prior', None),
            'fixed_flux_baseline': getattr(reference, 'fixed_flux_baseline', False),
            'impact_parameter_sampling': getattr(reference, 'use_impactparameter_rather_than_inclination_to_fit', True),
            'minimum_live_points': getattr(reference, 'ultranest_min_num_live_points', None),
            'max_ncalls': reference.max_ncalls,
            'exposure_supersample': getattr(reference, 'exposure_smearing_supersample', None),
            'exposure_change_tolerance': getattr(reference, 'exposure_smearing_change_tolerance', None),
        },
        'stellar_atmosphere': atmosphere,
        'claret': claret,
        'hybrid4': alternative,
        'differences': deltas,
        'interpretation': (
            'Transit posterior intervals are conditional on each law\'s fixed atmosphere-derived '
            'coefficients. Atmosphere coefficient covariance is a separate diagnostic. '
            'The existing Claret result and its reporting adjustments are retained. '
            'The secondary fit uses the final Claret bounds and baseline treatment without repeating '
            'photometry selection or changing masks. Log-evidence differences use the existing '
            'likelihood, including any duration penalty and profiled baseline; they compare these '
            'conditional models, not atmosphere-marginalised models.'
        ),
    })


def save_comparison(reference, hybrid, atmosphere, directory, stem):
    """Write both parameter sets, identical inputs, paired curves and a plot."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    payload = comparison_payload(reference, hybrid, atmosphere)
    json_path = directory / f'{stem}.json'
    csv_path = directory / f'{stem}.csv'
    png_path = directory / f'{stem}.png'
    exposures = getattr(reference, 'exposure_times_days', None)
    exposures = np.zeros(len(reference.time)) if exposures is None else np.asarray(exposures) * 86400.0
    with csv_path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['BJD_TDB', 'Shared Flux', 'Shared Uncertainty', 'Airmass', 'Exposure Seconds',
                         'Claret Transit', 'Claret Reported Model', 'Claret Residual',
                         'Hybrid4 Transit', 'Hybrid4 Model', 'Hybrid4 Residual'])
        # Python float repr round-trips float64 even when another module has
        # changed NumPy's scalar print precision. This matters for BJD timing.
        writer.writerows(
            [float(value) for value in row]
            for row in zip(reference.time, reference.data, reference.dataerr, reference.airmass,
                           exposures, reference.transit, reference.model,
                           np.asarray(reference.data) - reference.model, hybrid.transit,
                           hybrid.model, np.asarray(hybrid.data) - hybrid.model)
        )
    plot_comparison(reference, hybrid, png_path, stem.replace('_', ' '))
    # Relative to this JSON: bulk archiving deletes the temporary workspace.
    payload['files'] = {'json': json_path.name, 'photometry_and_models_csv': csv_path.name,
                        'comparison_plot': png_path.name}
    json_path.write_text(json.dumps(payload, indent=2, allow_nan=False), encoding='utf-8')
    return payload


def plot_comparison(reference, hybrid, png_path, title, labels=None):
    """Evaluate continuous models across gaps instead of joining sparse samples."""
    import matplotlib.pyplot as plt
    epoch = float(reference.parameters['tmid'])
    hours = (np.asarray(reference.time) - epoch) * 24.0
    grid = np.linspace(np.min(reference.time), np.max(reference.time), 1000)
    grid_hours = (grid - epoch) * 24.0
    curves = []
    if labels is None:
        labels = ('Claret', 'hybrid4 (' + hybrid.parameters.get('_hybrid4_backend', 'numpy') + ')')
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(10, 8),
                             gridspec_kw={'height_ratios': [2, 1, 1]})
    axes[0].errorbar(hours, reference.data, yerr=reference.dataerr, fmt='.',
                     color='0.55', alpha=0.6, label='Shared photometry', zorder=1)
    for fit, label, color in ((reference, labels[0] + ' / UltraNest', '#2166ac'),
                              (hybrid, labels[1] + ' / UltraNest', '#d6604d')):
        curve = fit._transit_model(grid, fit.parameters) * fit._build_systematics_model_at(fit.parameters, grid)
        curves.append(curve)
        axes[0].plot(grid_hours, curve, color=color, label=label)
        axes[1].plot(hours, (np.asarray(fit.data) - fit.model) * 1e6, '.',
                     color=color, alpha=0.65, label=label)
    axes[0].set_ylabel('Flux supplied to both fits')
    axes[0].legend(loc='best')
    axes[1].set_ylabel('Residual [ppm]')
    axes[1].axhline(0, color='0.7', linewidth=0.7)
    axes[2].plot(grid_hours, (curves[1] - curves[0]) * 1e6,
                 color='#6a3d9a')
    axes[2].set_ylabel(f'{labels[1]} - {labels[0]}\nmodel [ppm]')
    axes[2].set_xlabel(f'Hours from {labels[0]} Tmid ({epoch:.8f} BJD_TDB)')
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(png_path, dpi=160)
    plt.close(fig)


def comparison_for_output(fit, info, planet_name, date):
    """Run once for a selected full reduction and expose actual failures."""
    if getattr(fit, 'ns_type', None) != 'ultranest' or getattr(fit, 'mode', None) != 'ns':
        return None
    if getattr(fit, 'stellar_variability_only', False):
        return None
    cached = getattr(fit, 'hybrid4_comparison', None)
    if cached is not None:
        return cached
    try:
        from ..utils import filename_date_token, safe_output_filename
    except ImportError:
        from utils import filename_date_token, safe_output_filename
    stem = Path(safe_output_filename('LimbDarkeningComparison', planet_name,
                                     filename_date_token(date), extension='json')).stem
    atmosphere = info.get('hybrid4_limb_darkening')
    if not isinstance(atmosphere, dict):
        result = {'status': 'unavailable', 'reason': 'No shared stellar-atmosphere profiles were '
                  'generated for this Claret fit; hybrid4 cannot use the same profiles.'}
        log.warning('hybrid4 comparison unavailable: %s', result['reason'])
    else:
        try:
            hybrid = run_hybrid4_fit(fit, atmosphere)
            result = save_comparison(fit, hybrid, atmosphere, info['save'], stem)
            # The companion JSON holds the atmosphere arrays. Keep the main
            # FinalParams file concise while reporting both numeric fit sets.
            result = {key: value for key, value in result.items() if key != 'stellar_atmosphere'}
            fit.hybrid4_fit = hybrid
        except Exception as exc:
            result = {'status': 'failed', 'error_type': type(exc).__name__, 'reason': str(exc)}
            log.exception('hybrid4 comparison failed: %s', exc)
            directory = Path(info['save'])
            directory.mkdir(parents=True, exist_ok=True)
            (directory / f'{stem}.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    fit.hybrid4_comparison = result
    return result
