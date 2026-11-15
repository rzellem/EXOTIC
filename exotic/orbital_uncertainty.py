"""Post-fit orbital sensitivity, without inventing a posterior from a limit.

These are envelopes of transformations of a *conditional* transit posterior.
When the original likelihood is available, temporary importance weights also
measure conditional posterior sensitivity. No optimiser or sampler is run and
no fitted parameter, error, or stored posterior weight is changed.
"""

import numpy as np

try:
    from posterior import summarize_posterior, fit_posterior_sample_values
except ImportError:
    from .posterior import summarize_posterior, fit_posterior_sample_values


def _conditional_orbit_draws(fit, ecc, omega, values, weights):
    """Evaluate likelihood/prior ratios on existing joint draws, without fitting.

    The normalisation of the existing geometry-dependent uniform b prior is
    included. Extra prior support cannot be created by importance weighting and
    is reported as a diagnostic, not hidden by an acceptance threshold.
    """
    context = fit.orbital_likelihood_context
    n = len(weights)
    parameters = context['base_physical']
    sampled_keys = list(fit.sampled_keys)
    logl, prior_density, support_expansion = np.full(n, -np.inf), np.zeros(n), np.zeros(n)
    physical_values = {key: np.asarray(value).copy() for key, value in values.items()}
    for index in range(n):
        physical = dict(parameters)
        physical.update({key: float(values[key][index]) for key in sampled_keys})
        old_ecc, old_omega = physical['ecc'], physical['omega']
        physical.update(ecc=float(ecc), omega=float(omega))
        if 'b' in sampled_keys:
            old_f = (1 - old_ecc ** 2) / (1 + old_ecc * np.sin(np.radians(old_omega)))
            new_f = (1 - ecc ** 2) / (1 + ecc * np.sin(np.radians(omega)))
            old_upper = min(1 + physical['rprs'], physical['ars'] * old_f)
            new_upper = min(1 + physical['rprs'], physical['ars'] * new_f)
            support_expansion[index] = max(0., 1 - old_upper / new_upper)
            cosine = physical['b'] / (physical['ars'] * new_f)
            if not 0 <= cosine <= 1 or physical['b'] > new_upper:
                continue
            physical['inc'] = float(np.degrees(np.arccos(cosine)))
            physical_values['inc'][index] = physical['inc']
            prior_density[index] = np.log(old_upper / new_upper)
        elif 'inc' in sampled_keys:
            new_f = (1 - ecc ** 2) / (1 + ecc * np.sin(np.radians(omega)))
            physical_values['b'][index] = physical['ars'] * new_f * np.cos(np.radians(physical['inc']))
        model = np.asarray(fit._transit_model(context['time'], physical), dtype=float)
        model *= np.exp(physical.get('a2', 0.) * context['centered_airmass'])
        baseline_key = context.get('free_flux_baseline_key')
        if baseline_key is not None:
            model *= physical[baseline_key]
        elif context['uses_fixed_flux_baseline']:
            model *= context['fixed_flux_baseline_value']
        elif context['has_free_flux_baseline']:
            from exotic.api.elca import get_flux_baseline
            model *= get_flux_baseline(physical)
        else:
            mask = context['baseline_static_mask'] & np.isfinite(model) & (model != 0)
            denom = np.sum(context['baseline_weights'][mask] * model[mask] ** 2)
            # Preserve the original likelihood's profiled-baseline solver,
            # including its existing empty/degenerate-model behaviour.
            if not np.any(mask):
                baseline = 1.
            elif not np.isfinite(denom) or denom <= 0:
                ratios = context['data'][mask] / model[mask]
                ratios = ratios[np.isfinite(ratios)]
                baseline = float(np.median(ratios)) if ratios.size else 1.
            else:
                baseline = np.sum(context['baseline_weights'][mask] * context['data'][mask] * model[mask]) / denom
            model *= baseline if np.isfinite(baseline) else 1.
        residual = (context['data'] - model) * context['inverse_dataerr']
        score = -.5 * np.sum(residual ** 2)
        if context['duration_prior_valid']:
            from exotic.api.elca import transit_duration
            duration = transit_duration(physical)
            if not np.isfinite(duration) or duration <= 0:
                continue
            score -= .5 * (np.log(duration / context['expected_duration']) / context['sigma_log_duration']) ** 2
        if np.isfinite(score):
            logl[index] = score
    return physical_values, logl + prior_density, support_expansion


def orbital_posterior_sensitivity(fit, record, region):
    """Separate conditional random uncertainty from orbit-assumption changes.

    A deterministic grid is a sensitivity scan, not an eccentricity PDF. Only
    explicitly supplied joint RV posterior draws define an orbital random
    component. Their weights remain fixed (a modular external-orbit analysis).
    """
    values, weights, source = fit_posterior_sample_values(fit)
    unavailable = {'status': 'unavailable', 'fitted_values_changed': False}
    if not values or not hasattr(fit, 'orbital_likelihood_context'):
        return dict(unavailable, reason='Saved joint posterior and original likelihood context are required; geometry envelopes remain separately available.')
    if any(key in fit.sampled_keys for key in ('ecc', 'omega')):
        return dict(unavailable, reason='Orbit is already sampled jointly; do not double-count external orbital uncertainty.')
    weights = np.ones(len(next(iter(values.values())))) if weights is None else np.asarray(weights)
    valid = np.isfinite(weights) & (weights > 0)
    for array in values.values():
        valid &= np.isfinite(array)
    if not valid.any():
        return dict(unavailable, reason='There are no finite joint draws with positive posterior weight.')
    values, weights = {key: array[valid] for key, array in values.items()}, weights[valid]
    weights = weights / weights.sum()
    params = fit.orbital_likelihood_context['base_physical']
    e0, w0 = float(params['ecc']), float(params['omega'])
    _, baseline_logl, _ = _conditional_orbit_draws(fit, e0, w0, values, weights)
    cache = {}

    def evaluate(ecc, omega):
        orbit = (float(ecc), float(omega % 360) if ecc else 0.)
        if orbit not in cache:
            draws, score, expansion = _conditional_orbit_draws(fit, *orbit, values, weights)
            log_weights = np.log(weights) + score - baseline_logl
            finite = np.isfinite(log_weights)
            if not finite.any():
                cache[orbit] = {'status': 'no_likelihood_overlap', 'ecc': orbit[0], 'omega_planet_degrees': orbit[1]}
            else:
                normalized = np.zeros(len(weights))
                normalized[finite] = np.exp(log_weights[finite] - np.max(log_weights[finite]))
                normalized /= normalized.sum()
                summaries = {key: summarize_posterior(draws[key], normalized)
                             for key in ('tmid', 'rprs', 'ars', 'b', 'inc') if key in draws}
                means = {}
                for key in summaries:
                    anchor = float(draws[key][finite][0])
                    means[key] = float(anchor + np.sum(normalized[finite] * (draws[key][finite] - anchor)))
                cache[orbit] = {
                    'status': 'evaluated', 'ecc': orbit[0], 'omega_planet_degrees': orbit[1],
                    'effective_sample_count': float(1 / np.sum(normalized ** 2)),
                    'largest_normalized_weight': float(normalized.max()),
                    'maximum_conditional_b_prior_support_expansion_fraction': float(expansion.max()),
                    'random_conditional_posterior_intervals': summaries, 'posterior_means': means,
                }
        return cache[orbit]

    limit = (record.get('ecc') or {}).get('kind') == 'upper_limit'
    reference = evaluate(0., 0.) if limit else evaluate(e0, w0)
    result = {
        'status': 'evaluated', 'method': 'saved_joint_posterior_likelihood_and_prior_ratio_importance_sensitivity',
        'fitted_values_changed': False, 'stored_posterior_weights_changed': False,
        'posterior_sample_source': source, 'reference': reference,
        'description': 'Random intervals are conditional on each orbit. Changes of posterior means/medians '
                       'are orbital-assumption sensitivity, not new best-fit values. Grid extrema are explored '
                       'extrema, not certified global bounds or standard deviations. Importance weighting cannot '
                       'discover absent posterior support; ESS, largest weight and b-prior support expansion '
                       'are reported without a rejection cutoff. Baseline, exposure integration, duration prior '
                       'and the existing b-prior normalisation are retained. Independent sources are alternatives.',
    }
    stored_logl = ((getattr(fit, 'results', {}) or {}).get('weighted_samples') or {}).get('logl')
    if source == 'weighted_samples' and stored_logl is not None and np.asarray(stored_logl).shape == valid.shape:
        difference = baseline_logl - np.asarray(stored_logl)[valid]
        result['original_likelihood_replay_max_absolute_difference'] = float(np.max(np.abs(difference)))
    joint = record.get('joint_posterior')
    if joint:
        eccentricities, angles, orbital_weights = _joint_orbit_arrays(joint)
        scenarios = [evaluate(e, w) for e, w in zip(eccentricities, angles)]
        result['orbital_random_uncertainty'] = {'status': 'joint_external_posterior_supplied',
                                               'external_weights_updated_by_transit': False}
    else:
        eccentricities = np.linspace(*region['ecc_range'], 3 if limit else 5)
        scenarios = []
        coarse_orbits = set()
        for interval in region['omega_ranges_planet_degrees']:
            angles = np.arange(16) * 22.5 if interval is None else np.linspace(*interval, 9)
            for ei, e in enumerate(eccentricities):
                for wi, w in enumerate(angles):
                    scenario = evaluate(e, w)
                    if scenario not in scenarios:
                        scenarios.append(scenario)
                    if ei % 2 == 0 and wi % 2 == 0:
                        coarse_orbits.add((scenario['ecc'], scenario['omega_planet_degrees']))
        result['grid'] = {'ecc_nodes': len(eccentricities), 'omega_nodes_per_interval': len(angles),
                          'distinct_orbits': len(scenarios), 'probability_weights_assigned': False}
        result['orbital_random_uncertainty'] = {
            'status': 'not_identifiable_from_limit_or_marginal_errors', 'sigma': None,
            'reason': 'An upper limit does not define a PDF. Marginal asymmetric errors do not supply '
                      'the joint e/omega distribution or covariance. No uniform or Gaussian orbital prior is invented.',
        }
    evaluated = [scenario for scenario in scenarios if scenario['status'] == 'evaluated']
    result['unavailable_orbit_count'] = len(scenarios) - len(evaluated)
    result['scenarios'] = scenarios
    result['parameters'] = {}
    if reference['status'] != 'evaluated' or not evaluated:
        result['status'] = 'no_reference_or_scenario_overlap'
        return result
    for key, random_reference in reference['random_conditional_posterior_intervals'].items():
        scale = 86400. if key == 'tmid' else 1.
        unit = 'seconds' if key == 'tmid' else 'degrees' if key == 'inc' else 'dimensionless'
        ref_mean, ref_median = reference['posterior_means'][key], random_reference['median']
        shifts = np.array([(s['posterior_means'][key] - ref_mean) * scale for s in evaluated])
        median_shifts = [(s['random_conditional_posterior_intervals'][key]['median'] - ref_median) * scale for s in evaluated]
        sigmas = [s['random_conditional_posterior_intervals'][key]['stdev'] * scale for s in evaluated]
        detail = {
            'units': unit,
            'random_conditional_sigma_at_reference': random_reference['stdev'] * scale,
            'random_conditional_error_minus_at_reference': random_reference['error_minus'] * scale,
            'random_conditional_error_plus_at_reference': random_reference['error_plus'] * scale,
            'random_conditional_sigma_range_over_orbits': [min(sigmas), max(sigmas)],
            'orbital_assumption_mean_shift_range': [float(shifts.min()), float(shifts.max())],
            'orbital_assumption_median_shift_range': [min(median_shifts), max(median_shifts)],
            'maximum_absolute_explored_mean_shift': float(np.max(np.abs(shifts))),
        }
        if joint and len(evaluated) == len(scenarios):
            means = np.array([s['posterior_means'][key] for s in scenarios])
            anchor = means[0]
            mixed_mean = anchor + np.sum(orbital_weights * (means - anchor))
            within = np.sum(orbital_weights * np.array(sigmas) ** 2)
            between = np.sum(orbital_weights * ((means - mixed_mean) * scale) ** 2)
            detail.update(random_transit_sigma_mixture=float(np.sqrt(within)),
                          random_external_orbit_sigma=float(np.sqrt(between)),
                          total_random_sigma=float(np.sqrt(within + between)),
                          systematic_reference_to_mixture_mean_shift=float((mixed_mean - ref_mean) * scale))
        elif not joint:
            coarse = [s for s in evaluated if (s['ecc'], s['omega_planet_degrees']) in coarse_orbits]
            if coarse:
                coarse_shift = max(abs(s['posterior_means'][key] - ref_mean) * scale for s in coarse)
                detail['maximum_mean_shift_change_on_grid_refinement'] = detail['maximum_absolute_explored_mean_shift'] - coarse_shift
        result['parameters'][key] = detail
    return result


def _number(value):
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def _joint_orbit_arrays(joint):
    eccentricities = np.asarray(joint['ecc'], dtype=float)
    angles = np.asarray(joint['omega_planet_degrees'], dtype=float)
    weights = np.array(joint.get('weights', np.ones(eccentricities.size)), dtype=float, copy=True)
    if (eccentricities.ndim != 1 or eccentricities.shape != angles.shape or eccentricities.shape != weights.shape
            or not np.isfinite(weights).all() or np.any(weights < 0) or not weights.sum() > 0):
        raise ValueError('Joint orbital posterior requires paired finite elliptic ecc/planet omega draws and nonnegative weights.')
    positive = weights > 0
    eccentricities, angles, weights = eccentricities[positive], angles[positive], weights[positive]
    if (not np.isfinite(eccentricities).all() or not np.isfinite(angles).all()
            or np.any(eccentricities < 0) or np.any(eccentricities >= 1)):
        raise ValueError('Positive-weight joint orbital samples must have finite planet omega and 0 <= eccentricity < 1.')
    return eccentricities, angles, weights / weights.sum()


def archive_orbital_records(rows):
    """Retain errors, censoring and provenance from the SAME archive row.

    NASA does not specify the confidence of a limit or standardise the omega
    convention. Neither can be inferred from the limit flag or null errors.
    """
    records = []
    for row in rows:
        ecc = _number(row.get('pl_orbeccen'))
        omega = _number(row.get('pl_orblper'))
        if ecc is None and omega is None:
            continue
        record = {'source': row.get('pl_refname'), 'publication_date': row.get('pl_pubdate'),
                  'omega_convention': 'unspecified'}
        for key, column in [('ecc', 'pl_orbeccen'), ('omega', 'pl_orblper')]:
            flag = _number(row.get(column + 'lim'))
            record[key] = {
                'value': _number(row.get(column)),
                'error_plus': _number(row.get(column + 'err1')),
                'error_minus': _number(row.get(column + 'err2')),
                'kind': {1: 'upper_limit', -1: 'lower_limit', 0: 'measurement'}.get(flag, 'unspecified'),
                'confidence_probability': None,
            }
        records.append(record)
    return records


def nominal_eccentricity(fixed_value, records):
    """Use circular as an explicit nominal assumption for an upper limit.

    The first row with an eccentricity is the archive's latest relevant source.
    This applies when constructing a new model, never to a saved posterior.
    """
    for record in records or []:
        ecc = record.get('ecc') or {}
        if _number(ecc.get('value')) is not None:
            return 0. if ecc.get('kind') == 'upper_limit' else fixed_value
    return fixed_value


def _sine_range(omega):
    """Sine extrema over a continuous angle interval, including wraparound."""
    lo, hi = omega
    if hi - lo >= 360:
        return -1., 1.
    angles = [lo, hi]
    first = int(np.ceil((lo - 90) / 180))
    last = int(np.floor((hi - 90) / 180))
    angles.extend(90 + 180 * np.arange(first, last + 1))
    values = np.sin(np.radians(angles))
    return float(np.min(values)), float(np.max(values))


def orbital_factor_bounds(ecc_range, omega_range=None):
    """Exact extrema of f=(1-e^2)/(1+e sin w), g=(1+e sin w)/sqrt(1-e^2).

    Omitting omega covers every orientation; it does not assign a uniform
    probability to angles. Eccentricity must describe a bound elliptic orbit.
    """
    lo, hi = map(float, ecc_range)
    if not 0 <= lo <= hi < 1:
        raise ValueError('Eccentricity bounds must satisfy 0 <= lower <= upper < 1.')
    sine_bounds = (-1., 1.) if omega_range is None else _sine_range(omega_range)
    f_values, g_values = [], []
    for sine in sine_bounds:
        f_candidates = [lo, hi]
        g_candidates = [lo, hi]
        if sine < 0:
            stationary_f = -sine / (1 + np.sqrt(max(0., 1 - sine ** 2)))
            if lo <= stationary_f <= hi:
                f_candidates.append(stationary_f)
            if lo <= -sine <= hi:
                g_candidates.append(-sine)
        for ecc in f_candidates:
            f_values.append((1 - ecc ** 2) / (1 + ecc * sine))
        for ecc in g_candidates:
            g_values.append((1 + ecc * sine) / np.sqrt(1 - ecc ** 2))
    return {'separation_factor': [float(min(f_values)), float(max(f_values))],
            'speed_factor': [float(min(g_values)), float(max(g_values))]}


def _constraint_ranges(record):
    if record.get('joint_posterior'):
        ecc, omega, weights = _joint_orbit_arrays(record['joint_posterior'])
        ecc, omega = ecc[weights > 0], omega[weights > 0]
        e_range = [float(ecc.min()), float(ecc.max())]
        w_range = [float(omega.min()), float(omega.max())]
        return {'ecc_range': e_range, 'omega_ranges_planet_degrees': [w_range],
                'region_kind': 'conservative_marginal_draw_box_not_a_joint_confidence_region',
                **orbital_factor_bounds(e_range, w_range)}, None
    ecc = record.get('ecc') or {}
    value = _number(ecc.get('value'))
    if value is None:
        return None, 'No eccentricity constraint in this source.'
    if ecc.get('kind') == 'upper_limit':
        e_range = [0., value]
    elif ecc.get('kind') == 'measurement':
        ep, em = _number(ecc.get('error_plus')), _number(ecc.get('error_minus'))
        if ep is None or em is None:
            return None, 'Quoted value has no errors; it does not specify orbital uncertainty.'
        # The lower endpoint is intersected with the physical e >= 0 domain.
        e_range = [max(0., value - abs(em)), value + abs(ep)]
    else:
        return None, 'This source does not supply a finite eccentricity uncertainty region.'
    if not 0 <= e_range[0] <= e_range[1] < 1:
        return None, 'Quoted region reaches an unbound orbit; finite elliptic extrema are unavailable.'
    omega = record.get('omega') or {}
    ov = _number(omega.get('value'))
    op, om = _number(omega.get('error_plus')), _number(omega.get('error_minus'))
    if (omega.get('kind') == 'measurement' and ov is not None
            and op is not None and om is not None):
        intervals = [[ov - abs(om), ov + abs(op)]]
        convention = record.get('omega_convention', 'unspecified')
        if convention == 'star':
            intervals = [[a + 180, b + 180] for a, b in intervals]
        elif convention != 'planet':
            # A convention ambiguity is represented, never silently resolved.
            intervals += [[a + 180, b + 180] for a, b in intervals]
    else:
        intervals = [None]
    factors = [orbital_factor_bounds(e_range, interval) for interval in intervals]
    return {'ecc_range': e_range, 'omega_ranges_planet_degrees': intervals,
            **{key: [min(f[key][0] for f in factors), max(f[key][1] for f in factors)]
               for key in ('separation_factor', 'speed_factor')}}, None


def _geometry_samples(fit):
    """Use posterior weights, never equal weights for nested-sampler dead points."""
    results = getattr(fit, 'results', {}) or {}
    weighted = results.get('weighted_samples') or {}
    points, weights = weighted.get('points'), weighted.get('weights')
    source = 'weighted_samples'
    if points is None or weights is None:
        points, weights, source = results.get('samples'), None, 'equal_weight_posterior_samples'
    if points is None:
        return None
    points = np.asarray(points, dtype=float)
    keys = list(getattr(fit, 'sampled_keys', []))
    if points.ndim != 2 or not points.shape[0] or points.shape[1] < len(keys):
        return None
    if weights is None:
        weights = np.ones(points.shape[0])
    weights = np.asarray(weights, dtype=float)
    if weights.shape != (points.shape[0],) or not np.isfinite(weights).all() or np.any(weights < 0):
        return None
    if not np.sum(weights) > 0:
        return None
    params = getattr(fit, 'parameters', {})
    def values(key):
        if key in keys:
            return points[:, keys.index(key)]
        return np.full(points.shape[0], params.get(key, np.nan), dtype=float)
    ars, ecc, omega = values('ars'), values('ecc'), values('omega')
    f = (1 - ecc ** 2) / (1 + ecc * np.sin(np.radians(omega)))
    b = values('b') if 'b' in keys else ars * f * np.cos(np.radians(values('inc')))
    return {'b': b, 'ars': ars, 'ecc': ecc, 'omega': omega,
            'weights': weights, 'source': source}


def orbital_uncertainty_budget(fit, records):
    """Report external-orbit envelopes alongside, never replacing, fit errors.

    Joint transit b/a/R* correlations and weights are retained. Different
    publications are alternatives, since their RV data can overlap. A quoted
    region is not hard support: its outside tail is explicitly left unknown.
    """
    report = {
        'method': 'post_fit_orbital_uncertainty_components',
        'fitted_values_changed': False,
        'description': (
            'Envelopes cover the quoted orbital region and all unspecified orientations. '
            'They do not assume a probability density inside an upper limit. '
            'Endpoint quantiles describe the saved conditional transit posterior, not a joint credible interval. '
            'Probability outside a quoted region is not set to zero. '
            'Independent publications are reported separately, not multiplied. '
            'a/R* and density factors are small-angle duration-preserving sensitivity approximations; '
            'Rp/R*, impact parameter and transit-time inference are not marginalised over new orbital values.'
            ' When the original likelihood is present, a separate importance-sensitivity scan quantifies '
            'their conditional posterior changes without refitting. Random posterior intervals, orbital '
            'measurement uncertainty and assumption sensitivity are separate; an envelope is not added in quadrature.'
        ),
        'constraints': [],
    }
    values, random_weights, _ = fit_posterior_sample_values(fit)
    report['random_conditional_transit_uncertainty'] = {
        key: summarize_posterior(array, random_weights) for key, array in values.items()
        if key in ('tmid', 'rprs', 'ars', 'b', 'inc')
    }
    report['uncertainty_components_description'] = (
        'Random: the weighted transit posterior conditional on the fitted orbit and other model assumptions. '
        'Orbital random: between-orbit variance only when a joint external orbital posterior is supplied. '
        'Systematic/assumption sensitivity: shifts across the documented orbital region, relative to the circular '
        'assumption for an upper limit or the saved orbit for a measurement. Other systematics are not included. '
        'Different source constraints and nested confidence regions are not independent uncertainty terms.'
    )
    geometry = _geometry_samples(fit)
    report['posterior_sample_source'] = geometry['source'] if geometry is not None else 'unavailable'
    for record in records or []:
        entry = {'constraint': record, 'best_eccentricity_from_constraint': None}
        region, reason = _constraint_ranges(record)
        if region is None:
            entry.update(status='uncertainty_unspecified', reason=reason)
            report['constraints'].append(entry)
            continue
        entry.update(status='envelope', **region,
                     outside_region_probability='unknown unless documented by the source')
        is_limit = (record.get('ecc') or {}).get('kind') == 'upper_limit'
        entry['nominal_orbit'] = 'circular_assumed' if is_limit else 'published_measurement'
        entry['quoted_confidence_probability'] = (record.get('ecc') or {}).get('confidence_probability')
        if is_limit:
            entry['nominal_eccentricity'] = 0.
            entry['nominal_argument_of_periastron_degrees'] = None
            entry['best_value_description'] = 'Circular is an assumed nominal model, not an estimate derived from the upper limit.'
        if geometry is not None:
            b, ars, weights = geometry['b'], geometry['ars'], geometry['weights']
            f_lo, f_hi = region['separation_factor']
            with np.errstate(invalid='ignore', divide='ignore'):
                cos_lo, cos_hi = b / (ars * f_lo), b / (ars * f_hi)
                valid = (np.isfinite(b) & np.isfinite(ars) & (ars > 0) & (b >= 0)
                         & (cos_hi <= 1))
                full_region = valid & (cos_lo <= 1)
                # Intersect the region with the existence of a real inclination.
                # At f=b/ars the exact boundary is i=0; retain such samples and
                # explicitly report the weight with only partial compatibility.
                inc_lo = np.where(full_region, np.degrees(np.arccos(cos_lo)), 0.)
                inc_hi = np.degrees(np.arccos(cos_hi))
            entry['fully_transformable_posterior_weight'] = float(weights[full_region].sum() / weights.sum())
            entry['partly_transformable_posterior_weight'] = float(weights[valid & ~full_region].sum() / weights.sum())
            entry['untransformable_posterior_weight'] = float(weights[~valid].sum() / weights.sum())
            entry['inclination_endpoint_summaries_degrees'] = {
                'minimum': summarize_posterior(inc_lo[valid], weights[valid]),
                'maximum': summarize_posterior(inc_hi[valid], weights[valid]),
            }
            if is_limit:
                with np.errstate(invalid='ignore'):
                    circular_inc = np.degrees(np.arccos(b / ars))
                deviation = np.maximum(np.abs(inc_lo - circular_inc), np.abs(inc_hi - circular_inc))
                entry['maximum_inclination_departure_from_circular_per_sample_degrees'] = summarize_posterior(deviation[valid], weights[valid])
                entry['circular_reference_inclination_degrees'] = summarize_posterior(circular_inc, weights)
            g0 = (1 + geometry['ecc'] * np.sin(np.radians(geometry['omega']))) / np.sqrt(1 - geometry['ecc'] ** 2)
            g_lo, g_hi = region['speed_factor']
            entry['duration_preserving_ars_from_saved_fit_endpoint_summaries_approximate'] = {
                'minimum': summarize_posterior(ars * g0 / g_hi, weights),
                'maximum': summarize_posterior(ars * g0 / g_lo, weights),
            }
        params = getattr(fit, 'parameters', {})
        e0, w0 = _number(params.get('ecc')), _number(params.get('omega'))
        if e0 is not None and w0 is not None and 0 <= e0 < 1:
            g0 = (1 + e0 * np.sin(np.radians(w0))) / np.sqrt(1 - e0 ** 2)
            g_lo, g_hi = region['speed_factor']
            # Upper-limit deviations are referenced to the explicitly assumed
            # circular model, even if a historical fit fixed e at the limit.
            reference_g = 1. if is_limit else g0
            entry['sensitivity_reference_eccentricity'] = 0. if is_limit else e0
            entry['saved_fit_eccentricity'] = e0
            ratios = [reference_g / g_hi, reference_g / g_lo]
            entry['duration_preserving_ars_factor_approximate'] = ratios
            entry['duration_preserving_stellar_density_factor_approximate'] = [r ** 3 for r in ratios]
            entry['fixed_chord_duration_factor_approximate'] = ratios
            entry['maximum_absolute_fractional_departure_approximate'] = {
                'duration_preserving_ars': float(max(abs(r - 1) for r in ratios)),
                'duration_preserving_stellar_density': float(max(abs(r ** 3 - 1) for r in ratios)),
                'fixed_chord_duration': float(max(abs(r - 1) for r in ratios)),
            }
            if is_limit:
                ars0 = _number(params.get('ars'))
                b0 = _number(params.get('b'))
                if b0 is None and ars0 is not None and _number(params.get('inc')) is not None:
                    f0 = (1 - e0 ** 2) / (1 + e0 * np.sin(np.radians(w0)))
                    b0 = ars0 * f0 * np.cos(np.radians(params['inc']))
                if ars0 is not None and b0 is not None and ars0 > 0 and 0 <= b0 / ars0 <= 1:
                    f_lo, f_hi = region['separation_factor']
                    if b0 / (ars0 * f_lo) <= 1:
                        i0 = np.degrees(np.arccos(b0 / ars0))
                        endpoints = np.degrees(np.arccos(b0 / (ars0 * np.array([f_lo, f_hi]))))
                        entry['model_geometry_circular_inclination_degrees'] = float(i0)
                        entry['model_geometry_maximum_inclination_departure_degrees'] = float(np.max(np.abs(endpoints - i0)))
        entry['posterior_sensitivity'] = orbital_posterior_sensitivity(fit, record, region)
        report['constraints'].append(entry)
    return report
