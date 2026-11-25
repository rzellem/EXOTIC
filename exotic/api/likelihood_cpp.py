"""One C++ call per UltraNest batch, with immutable per-fit scoring inputs."""
import ctypes

import numpy as np

from .limbdarkfast_cpp import _library


class NativeLogLikelihood:
    def __init__(self, fit, bound_keys, sampled_keys, impact_parameter_index,
                 observed_valid, fixed_airmass_scale):
        from .elca import (
            _exposure_smearing_supersample, _exposure_smearing_change_tolerance,
            MIN_EXPOSURE_SMEARING_SECONDS, EXPOSURE_SMEARING_TRANSIT_WINDOW_PADDING_FACTOR,
            get_flux_baseline,
        )
        context = fit.orbital_likelihood_context
        base = context['base_physical']
        self.dimension = len(sampled_keys)
        self.time = np.ascontiguousarray(context['time'], dtype=np.float64)
        law = base.get('_limb_darkening_law', 'claret')
        prefix = 'h' if law == 'hybrid4' else 'u'
        columns = ['rprs', 'per', 'ars', 'ecc', 'inc', 'omega', 'tmid', 'a2',
                   None, 'b'] + [f'{prefix}{i}' for i in range(4)]
        def number(key):
            try:
                return float(base.get(key, 0. if key in ('ecc', 'omega', 'a2') else np.nan))
            except (TypeError, ValueError):
                return np.nan
        self.base = np.array([number(key) if key else 1. for key in columns])
        mapping = np.full(14, -1, dtype=np.int32)
        for index, (key, _sampled_key) in enumerate(zip(bound_keys, sampled_keys)):
            if index == impact_parameter_index:
                mapping[9] = index
            elif key in columns:
                mapping[columns.index(key)] = index
        baseline_key = context['free_flux_baseline_key']
        baseline_mode = 1
        if baseline_key is not None:
            mapping[8] = sampled_keys.index(baseline_key)
        elif context['uses_fixed_flux_baseline']:
            self.base[8] = context['fixed_flux_baseline_value']
        elif context['has_free_flux_baseline']:
            self.base[8] = get_flux_baseline(base)
        else:
            baseline_mode = 0
        self.mapping = mapping
        if 'a2' in sampled_keys:
            mapping[7] = sampled_keys.index('a2')
        airmass_mode = 2 if 'a2' in sampled_keys else 1 if fixed_airmass_scale is not None else 0
        fixed_airmass = np.ascontiguousarray(
            fixed_airmass_scale if fixed_airmass_scale is not None else [], dtype=np.float64,
        )
        prepared = fit._values_with_exposure_smearing(base, self.time)
        exposure = np.empty(0, dtype=np.float64)
        raw = prepared.get('_exposure_time_days')
        if raw is not None:
            try:
                candidate = np.asarray(raw, dtype=np.float64)
                if candidate.shape == () or candidate.shape == self.time.shape:
                    exposure = np.ascontiguousarray(candidate).reshape(-1)
            except (TypeError, ValueError):
                pass  # Same exposure normalization as the existing transit wrapper.
        arrays = [np.ascontiguousarray(context[key], dtype=np.float64).reshape(-1)
                  for key in ('data', 'inverse_dataerr', 'centered_airmass', 'baseline_weights')]
        mask = np.ascontiguousarray(context['baseline_static_mask'], dtype=np.uint8).reshape(-1)
        self.arguments = (
            len(self.time), self.time, exposure, len(exposure), self.base, self.mapping,
            arrays[0], arrays[1], arrays[2], fixed_airmass, arrays[3], mask,
            0 if law == 'hybrid4' else 1, int(impact_parameter_index is not None),
            baseline_mode, airmass_mode, int(observed_valid), int(context['duration_prior_valid']),
            float(context['expected_duration']), float(context['sigma_log_duration']),
            _exposure_smearing_supersample(prepared), _exposure_smearing_change_tolerance(prepared),
            MIN_EXPOSURE_SMEARING_SECONDS / 86400., EXPOSURE_SMEARING_TRANSIT_WINDOW_PADDING_FACTOR,
        )
        # Resolve/compile before constructing UltraNest. Loading errors must not
        # be turned into invalid trial-point penalties or a Python fallback.
        self.function = _library().exotic_loglikelihood

    def __call__(self, points):
        points = np.asarray(points, dtype=np.float64)
        scalar = points.ndim == 1
        if scalar:
            points = points.reshape(1, -1)
        if points.ndim != 2 or points.shape[1] != self.dimension:
            raise ValueError('Likelihood points must have one column per sampled parameter.')
        points = np.ascontiguousarray(points)
        scores = np.empty(len(points), dtype=np.float64)
        error = ctypes.create_string_buffer(512)
        status = self.function(
            len(points), self.dimension, points.reshape(-1), *self.arguments,
            scores, error, len(error),
        )
        if status:
            raise RuntimeError(error.value.decode('utf-8', errors='replace'))
        return float(scores[0]) if scalar else scores
