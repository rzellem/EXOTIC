"""Batch the original conditional-orbit likelihood on saved posterior draws."""
import ctypes
import numpy as np

from .limbdarkfast_cpp import _library


def conditional_orbit_draws(fit, ecc, omega, values, weights):
    from .elca import (_exposure_smearing_supersample, _exposure_smearing_change_tolerance,
                       MIN_EXPOSURE_SMEARING_SECONDS, EXPOSURE_SMEARING_TRANSIT_WINDOW_PADDING_FACTOR)
    context = fit.orbital_likelihood_context
    base = context['base_physical']
    keys = list(fit.sampled_keys)
    n = len(weights)
    law = base.get('_limb_darkening_law', 'claret')
    coefficient_keys = [f'{"h" if law == "hybrid4" else "u"}{i}' for i in range(4)]
    columns = ['rprs','per','ars','ecc','inc','omega','tmid','a2',None,'b'] + coefficient_keys
    parameters = np.empty((n,14), dtype=np.float64)
    for column, key in enumerate(columns):
        if key is not None:
            parameters[:,column] = base.get(key, 0. if key == 'a2' else np.nan)
            if key in keys:
                parameters[:,column] = values[key]
    geometry = 1 if 'b' in keys else 2 if 'inc' in keys else 0
    if geometry == 1 and 'inc' in values:
        # Rejected draws retain the original derived inclinations.
        parameters[:,4] = values['inc']
    baseline_key = context.get('free_flux_baseline_key')
    mode = 1
    if baseline_key is not None:
        parameters[:,8] = values[baseline_key] if baseline_key in keys else base[baseline_key]
    elif context['uses_fixed_flux_baseline']:
        parameters[:,8] = context['fixed_flux_baseline_value']
    elif context['has_free_flux_baseline']:
        key = 'a0' if 'a0' in base else 'a1' if 'a1' in base else None
        parameters[:,8] = (values[key] if key in keys else base.get(key,1.)) if key else 1.
    else:
        parameters[:,8] = 1.
        mode = 0
    time = np.ascontiguousarray(context['time'], dtype=np.float64).reshape(-1)
    prepared = fit._values_with_exposure_smearing(base,time)
    exposure = np.empty(0,dtype=np.float64)
    raw = prepared.get('_exposure_time_days')
    if raw is not None:
        try:
            candidate = np.asarray(raw,dtype=np.float64)
            if candidate.shape == () or candidate.shape == time.shape:
                exposure = np.ascontiguousarray(candidate).reshape(-1)
        except (TypeError,ValueError):
            pass  # Existing exposure-input normalization, not a model fallback.
    arrays = [np.ascontiguousarray(context[key],dtype=np.float64).reshape(-1) for key in
              ['data','inverse_dataerr','centered_airmass','baseline_weights']]
    mask = np.ascontiguousarray(context['baseline_static_mask'],dtype=np.uint8).reshape(-1)
    scores, inclinations, impacts, expansion = [np.empty(n,dtype=np.float64) for _ in range(4)]
    error = ctypes.create_string_buffer(512)
    status = _library().exotic_orbital_replay(
        n,len(time),parameters.reshape(-1),time,exposure,len(exposure),*arrays,mask,
        0 if law == 'hybrid4' else 1,geometry,mode,float(ecc),float(omega),
        int(context['duration_prior_valid']),float(context.get('expected_duration') or 0.),
        float(context.get('sigma_log_duration') or 1.),
        _exposure_smearing_supersample(prepared),_exposure_smearing_change_tolerance(prepared),
        MIN_EXPOSURE_SMEARING_SECONDS/86400.,EXPOSURE_SMEARING_TRANSIT_WINDOW_PADDING_FACTOR,
        scores,inclinations,impacts,expansion,error,len(error),
    )
    if status:
        raise ValueError(error.value.decode('utf-8',errors='replace'))
    physical = {key:np.asarray(value).copy() for key,value in values.items()}
    if geometry == 1:
        physical['inc'][:] = inclinations
    elif geometry == 2:
        physical['b'][:] = impacts
    return physical,scores,expansion
