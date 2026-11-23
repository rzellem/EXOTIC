"""Small, data-only LM geometry profiles; never choose production parameters."""
import copy
from pathlib import Path
from time import perf_counter

import numpy as np


def profile_geometry_likelihood(fit, published_geometry=None, bounds=None):
    if __package__:
        from .elca import impact_parameter_from_inclination
    else:
        from elca import impact_parameter_from_inclination

    started = perf_counter()
    fit = copy.copy(fit)
    fit.bounds = copy.deepcopy(fit.bounds if bounds is None else bounds)
    bounds = copy.deepcopy(fit.bounds)
    reference = ({**fit.prior, **getattr(fit, 'lm_initial_parameters', {})}
                 if published_geometry is None else dict(published_geometry))
    seed = dict(fit.parameters)
    original_inc_bounds = bounds.pop('inc', None)
    original_b_bounds = bounds.pop('b', None)
    result = {
        'available': False, 'duration_prior_applied': False,
        'diagnostic_only': True, 'point_count': len(fit.time),
        'reference_geometry_source': ('supplied input ars/inc' if published_geometry is not None
                                      else 'LM input values; no separate published reference supplied'),
        'nuisance_bounds': {key: list(map(float, value)) for key, value in bounds.items()},
        'original_inclination_bounds': list(map(float, original_inc_bounds)) if original_inc_bounds is not None else None,
        'note': ('Data-only local LM fits; same observations, errors, exposure integration and baseline. '
                 'Conditional on the recorded nuisance search bounds and fixed parameters. '
                 'The coarse profile is not a posterior or a confidence interval; failures and bound contacts are reported.'),
    }
    fixed_context = ('per', 'rprs', 'ars', 'ecc', 'omega', 'u0', 'u1', 'u2', 'u3')
    result['fixed_parameters'] = {key: float(seed[key]) for key in fixed_context
                                  if key in seed and key not in bounds}
    result['baseline_mode'] = ('nearest_neighbors' if np.ndim(fit.airmass) == 2 else
                               'free_coefficients' if fit._has_free_flux_baseline() else
                               'fixed' if fit._uses_fixed_flux_baseline() else 'analytically_profiled')
    upper_rprs = float(bounds.get('rprs', [seed['rprs'], seed['rprs']])[1])
    # Full transiting-b domain permitted by the radius search, as in the nested
    # impact-parameter fit. Inclination bounds do not truncate the fixed-b scan.
    upper_b = min(1.0 + upper_rprs, fit._get_impact_parameter_scale_upper_bound(seed))
    if original_b_bounds is not None:
        lower_b, upper_b = map(float, original_b_bounds)
    else:
        lower_b = 0.0
    result['impact_parameter_bounds'] = [lower_b, upper_b]
    seed['b'] = float(impact_parameter_from_inclination(seed, seed['inc']))
    reference_b = None
    if all(reference.get(key) is not None and np.isfinite(reference[key]) for key in ('ars', 'inc')):
        reference_b = float(impact_parameter_from_inclination({**seed, **{
            key: reference[key] for key in ('ars', 'inc')
        }}, reference['inc']))

    def solve(local_bounds, initial, *, fixed_b=None, fixed_geometry=None):
        tick = perf_counter()
        local = copy.copy(fit)  # observational arrays remain shared and read-only
        local.prior = dict(initial)
        local.bounds = copy.deepcopy(local_bounds)
        local.duration_prior = None
        local.fixed_parameter_errors = {}
        local.mode = 'lm'
        local.lm_model_evaluations = 0
        fixed_geometry = fixed_geometry or {}
        local.prior.update(fixed_geometry)
        for key, (lo, hi) in local.bounds.items():
            local.prior[key] = float(np.clip(local.prior[key], lo, hi))
        entry = {'available': False, 'fixed_b': fixed_b, 'fixed_geometry': dict(fixed_geometry),
                 'fitted_parameters': list(local.bounds)}
        try:
            local.fit_LM(fixed_impact_parameter=fixed_b, create_fit_variables=False,
                         allow_unbounded_retry=False)
            parameter_keys = set(local.bounds) | {
                'rprs', 'ars', 'inc', 'b', 'tmid', 'a0', 'a1', 'a2', 'ecc', 'omega', 'per',
            }
            parameters = {key: float(local.parameters[key]) for key in parameter_keys
                          if key in local.parameters and np.isfinite(local.parameters[key])}
            parameters['b'] = float(impact_parameter_from_inclination(parameters, parameters['inc']))
            entry.update({
                'available': bool(np.isfinite(local.lm_data_chi2)),
                'chi2': local.lm_data_chi2, 'parameters': parameters,
                'converged': local.lm_optimizer_success, 'message': local.lm_optimizer_message,
                'active_bounds': dict(local.lm_active_bounds),
                'model_evaluations': local.lm_model_evaluations,
            })
        except Exception as exc:
            entry['error'] = f'{type(exc).__name__}: {exc}'
            entry['model_evaluations'] = local.lm_model_evaluations
        entry['elapsed_seconds'] = perf_counter() - tick
        return entry

    free_bounds = {**bounds, 'b': [lower_b, upper_b]}
    free = solve(free_bounds, seed)
    result['free_geometry'] = free
    profile_seed = {**seed, **free.get('parameters', {})}
    grid = np.linspace(lower_b, upper_b, 5).tolist()
    grid.append(float(profile_seed['b']))
    if reference_b is not None and lower_b <= reference_b <= upper_b:
        grid.append(reference_b)
    points = [solve(bounds, profile_seed, fixed_b=float(b)) for b in sorted(set(grid))]
    result['profile'] = points

    if reference_b is None:
        published = {'available': False, 'error': 'Supplied geometry lacks finite ars/inc values.'}
    else:
        fixed_geometry = {key: float(reference[key]) for key in ('ars', 'inc')}
        published_bounds = {key: value for key, value in bounds.items() if key not in fixed_geometry}
        published = solve(published_bounds, profile_seed, fixed_geometry=fixed_geometry)
    result['published_geometry'] = published
    known = [entry['chi2'] for entry in [free, *points] if entry.get('available')]
    if known:
        best = float(min(known))
        result['available'] = True
        result['minimum_explored_chi2'] = best
        for entry in [free, *points, published]:
            if entry.get('available'):
                entry['delta_chi2'] = float(entry['chi2'] - best)
                if free.get('available'):
                    entry['tmid_shift_from_free_seconds'] = float(
                        (entry['parameters']['tmid'] - free['parameters']['tmid']) * 86400.0
                    )
        result['profile_improves_on_free_fit'] = bool(free.get('available') and best < free['chi2'])
    result['fit_count'] = 1 + len(points) + int(reference_b is not None)
    result['model_evaluations'] = sum(entry.get('model_evaluations', 0) for entry in [free, *points, published])
    result['elapsed_seconds'] = perf_counter() - started
    return result


def save_geometry_profile_plot(summary, path):
    """Plot the explored likelihood and timing shifts without confidence cutoffs."""
    import matplotlib.pyplot as plt

    points = summary.get('profile', [])
    valid = [entry for entry in points if entry.get('available')]
    if not valid:
        return None
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True, layout='constrained')
    try:
        # Axis transitions only; no fit is accepted/rejected at these values.
        axes[0].set_yscale('symlog', linthresh=1.0)
        axes[1].set_yscale('symlog', linthresh=1.0)
        b = [entry['fixed_b'] for entry in valid]
        axes[0].plot(b, [entry['delta_chi2'] for entry in valid], 'o-', label='Fixed-b profile; nuisance terms refitted')
        timing = [entry.get('tmid_shift_from_free_seconds', np.nan) for entry in valid]
        axes[1].plot(b, timing, 'o-', label='Fixed-b profile')
        for key, label, marker in [('free_geometry', 'Free geometry', '*'),
                                   ('published_geometry', 'Supplied i and a/R* fixed', 'D')]:
            entry = summary.get(key, {})
            if entry.get('available'):
                x = entry['parameters']['b']
                axes[0].scatter(x, entry['delta_chi2'], marker=marker, s=80, label=label, zorder=4)
                axes[1].scatter(x, entry.get('tmid_shift_from_free_seconds', np.nan), marker=marker, s=65, label=label, zorder=4)
        axes[0].set_ylabel('Delta chi-square')
        axes[1].set_ylabel('Tmid shift from free fit (s)')
        all_entries = [*valid, summary.get('free_geometry', {}), summary.get('published_geometry', {})]
        if all(entry.get('delta_chi2', 0) >= 0 for entry in all_entries):
            axes[0].set_ylim(bottom=-0.2)  # visual padding for markers at zero
        axes[1].set_xlabel('Impact parameter b')
        axes[1].axhline(0, color='0.6', linewidth=0.7)
        axes[0].legend(fontsize=8)
        for ax in axes:
            ax.grid(alpha=0.2)
        failures = sum(not entry.get('available') for entry in points)
        unconverged = sum(entry.get('available') and not entry.get('converged', True) for entry in points)
        bound_contacts = sum(bool(entry.get('active_bounds')) for entry in points)
        fig.suptitle('LM geometry diagnostic: duration prior excluded')
        axes[0].set_title(f'Coarse profile; {failures} failed, {unconverged} not converged, {bound_contacts} with bound contacts.\n'
                          'Current nuisance bounds retained; vertical axes use a symmetric log scale.', fontsize=9)
        published = summary.get('published_geometry', {})
        shift = published.get('tmid_shift_from_free_seconds')
        if shift is not None:
            axes[1].set_title(f'Supplied geometry minus free geometry: Tmid {shift:+.3f} s', fontsize=10)
        fig.savefig(path, dpi=150)
    finally:
        plt.close(fig)
    return path
