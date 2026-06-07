# ########################################################################### #
#    Copyright (c) 2019-2020, California Institute of Technology.
#    All rights reserved.  Based on Government Sponsored Research under
#    contracts NNN12AA01C, NAS7-1407 and/or NAS7-03001.
#
#    Redistribution and use in source and binary forms, with or without
#    modification, are permitted provided that the following conditions
#    are met:
#      1. Redistributions of source code must retain the above copyright
#         notice, this list of conditions and the following disclaimer.
#      2. Redistributions in binary form must reproduce the above copyright
#         notice, this list of conditions and the following disclaimer in
#         the documentation and/or other materials provided with the
#         distribution.
#      3. Neither the name of the California Institute of
#         Technology (Caltech), its operating division the Jet Propulsion
#         Laboratory (JPL), the National Aeronautics and Space
#         Administration (NASA), nor the names of its contributors may be
#         used to endorse or promote products derived from this software
#         without specific prior written permission.
#
#    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE CALIFORNIA
#    INSTITUTE OF TECHNOLOGY BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
#    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# ########################################################################### #
#    EXOplanet Transit Interpretation Code (EXOTIC)
#    # NOTE: See companion file version.py for version info.
# ########################################################################### #
# Exoplanet light curve analysis
#
# Fit an exoplanet transit model to time series data.
# ########################################################################### #
from astropy.time import Time
import copy
from contextlib import redirect_stderr, redirect_stdout
import faulthandler
import io
from itertools import cycle, product
import math
import os
import sys
import bottleneck as bn
import matplotlib.pyplot as plt
import numpy as np
from scipy import spatial
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
from ultranest import ReactiveNestedSampler

try:
    from plotting import corner
except ImportError:
    from .plotting import corner

try:
    from ultranest_utils import run_reactive_sampler
except ImportError:
    from .ultranest_utils import run_reactive_sampler

BAD_LOG_LIKELIHOOD = -1.0e100
TRIANGLE_PLOT_EDGE_PEAK_FRACTION_MAX = 0.50
TRIANGLE_PLOT_EDGE_MIN_SAMPLE_COUNT = 30
TRIANGLE_PLOT_EDGE_EXPANSION_STEPS = 8
TRIANGLE_PLOT_FALLBACK_EXPANSION_BOUNDS = {
    'rprs': (0.0, 1.0),
}
TRANSIT_MODEL_UNCERTAINTY_KEYS = (
    'rprs', 'tmid', 'inc', 'ars', 'per', 'ecc', 'omega', 'u0', 'u1', 'u2', 'u3',
)
BASELINE_MODEL_UNCERTAINTY_KEYS = ('a0', 'a1', 'a2')
MODEL_UNCERTAINTY_POSTERIOR_SAMPLE_LIMIT = 2000
ULTRANEST_INFLATED_ERROR_REPLACEMENT_FACTOR = 3.0
ULTRANEST_LOCAL_UNCERTAINTY_MAX_DELTA_CHI2 = 9.0

def _pylightcurve_import_watchdog_seconds():
    try:
        return float(os.environ.get("EXOTIC_IMPORT_WATCHDOG_SECONDS", "120"))
    except (TypeError, ValueError):
        return 120.0


def _start_import_watchdog():
    timeout = _pylightcurve_import_watchdog_seconds()
    if timeout <= 0:
        return False

    try:
        if not faulthandler.is_enabled():
            faulthandler.enable(file=sys.__stdout__, all_threads=True)
        faulthandler.dump_traceback_later(timeout, repeat=True, file=sys.__stdout__)
        return True
    except Exception:
        return False


def _load_pylightcurve_transit():
    watchdog_started = _start_import_watchdog()
    try:
        with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
            from pylightcurve.models.exoplanet_lc import transit
        return transit
    finally:
        if watchdog_started:
            try:
                faulthandler.cancel_dump_traceback_later()
            except Exception:
                pass


pytransit = _load_pylightcurve_transit()


def weightedflux(flux, gw, nearest):
    return np.sum(flux[nearest] * gw, axis=-1)


def gaussian_weights(X, w=1, neighbors=50, feature_scale=1000):
    Xm = (X - np.median(X, 0)) * w
    kdtree = spatial.cKDTree(Xm * feature_scale)
    nearest = np.zeros((X.shape[0], neighbors))
    gw = np.zeros((X.shape[0], neighbors), dtype=float)
    for point in range(X.shape[0]):
        ind = kdtree.query(kdtree.data[point], neighbors + 1)[1][1:]
        dX = Xm[ind] - Xm[point]
        Xstd = np.std(dX, 0)
        gX = np.exp(-dX ** 2 / (2 * Xstd ** 2))
        gwX = np.product(gX, 1)
        gw[point, :] = gwX / gwX.sum()
        nearest[point, :] = ind
    gw[np.isnan(gw)] = 0.01
    return gw, nearest.astype(int)


def transit(times, values):
    model = pytransit([values['u0'], values['u1'], values['u2'], values['u3']],
                      values['rprs'], values['per'], values['ars'],
                      values['ecc'], values['inc'], values['omega'],
                      values['tmid'], times, method='claret', precision=3)
    return model


def impact_parameter_scale(values):
    ecc = values.get('ecc', 0.0)
    omega = np.deg2rad(values.get('omega', 0.0))
    denom = 1.0 + ecc * np.sin(omega)
    if np.any(np.isclose(denom, 0.0)):
        denom = np.where(np.isclose(denom, 0.0), np.finfo(float).eps, denom)
    return values['ars'] * (1.0 - ecc ** 2) / denom


def impact_parameter_from_inclination(values, inclination):
    return impact_parameter_scale(values) * np.cos(np.deg2rad(inclination))


def inclination_from_impact_parameter(values, impact_parameter):
    scale = impact_parameter_scale(values)
    if np.any(np.isclose(scale, 0.0)):
        scale = np.where(np.isclose(scale, 0.0), np.finfo(float).eps, scale)
    cosi = np.clip(np.asarray(impact_parameter, dtype=float) / scale, -1.0, 1.0)
    return np.rad2deg(np.arccos(cosi))


def grazing_impact_parameter(values):
    try:
        rprs = float(values['rprs'])
    except (KeyError, TypeError, ValueError):
        return np.nan

    if not np.isfinite(rprs) or rprs < 0:
        return np.nan
    return 1.0 + rprs


def transit_duration(values):
    try:
        period = float(values['per'])
        rprs = float(values['rprs'])
        ars = float(values['ars'])
        inc = float(values['inc'])
    except (KeyError, TypeError, ValueError):
        return np.nan

    if (
        not np.isfinite(period) or period <= 0
        or not np.isfinite(rprs) or rprs < 0
        or not np.isfinite(ars) or ars <= 0
        or not np.isfinite(inc)
    ):
        return np.nan

    ecc = values.get('ecc', 0.0)
    omega = np.deg2rad(values.get('omega', 0.0))
    sin_inc = np.sin(np.deg2rad(inc))
    if not np.isfinite(sin_inc) or sin_inc <= 0:
        return np.nan

    impact_scale = ars * (1.0 - ecc ** 2) / max(np.finfo(float).eps, 1.0 + ecc * np.sin(omega))
    impact_parameter = impact_scale * np.cos(np.deg2rad(inc))
    chord_sq = (1.0 + rprs) ** 2 - impact_parameter ** 2
    if not np.isfinite(chord_sq) or chord_sq <= 0 or not np.isfinite(impact_scale) or impact_scale <= 0:
        return np.nan

    argument = np.sqrt(chord_sq) / (impact_scale * sin_inc)
    argument = float(np.clip(argument, -1.0, 1.0))
    duration = (period / np.pi) * np.arcsin(argument)
    return float(duration) if np.isfinite(duration) and duration > 0 else np.nan


def get_phase(times, per, tmid):
    return (times - tmid + 0.25 * per) / per % 1 - 0.25


def normalize_time_range(time_range):
    if time_range is None:
        return None

    values = np.asarray(time_range, dtype=float).reshape(-1)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return None

    return float(np.min(finite)), float(np.max(finite))


def get_plot_phase(times, per, tmid, reference_times=None):
    times = np.asarray(times, dtype=float)
    if not np.isfinite(per) or per == 0:
        return times * np.nan

    raw_phase = (times - tmid) / per

    reference_range = normalize_time_range(reference_times)
    if reference_range is None:
        finite_phase = raw_phase[np.isfinite(raw_phase)]
        if finite_phase.size == 0:
            return raw_phase
        reference_epoch = float(np.rint(0.5 * (np.min(finite_phase) + np.max(finite_phase))))
    else:
        ref_phase = (np.asarray(reference_range, dtype=float) - tmid) / per
        reference_epoch = float(np.rint(np.mean(ref_phase)))

    return raw_phase - reference_epoch


def fallback_flux_baseline():
    return 1.0


def has_explicit_flux_baseline(bounds):
    return any(key in bounds for key in ('a0', 'a1'))


def get_flux_baseline(values, fallback=1.0):
    if 'a0' in values:
        return values['a0']
    if 'a1' in values:
        return values['a1']
    return fallback


def get_airmass_reference(airmass):
    airmass = np.asarray(airmass, dtype=float)
    finite_airmass = airmass[np.isfinite(airmass)]
    if finite_airmass.size == 0:
        return 0.0
    return float(np.nanmean(finite_airmass))


def center_airmass(airmass, reference=None):
    airmass = np.asarray(airmass, dtype=float)
    if reference is None:
        reference = get_airmass_reference(airmass)
    return airmass - float(reference)


def airmass_trend(a2, airmass, reference=None):
    return np.exp(np.asarray(a2, dtype=float) * center_airmass(airmass, reference=reference))


def airmass_trend_grid(a2_values, airmass, reference=None):
    centered = center_airmass(airmass, reference=reference)
    return np.exp(np.outer(np.asarray(a2_values, dtype=float), centered))


def normalized_optional_fit_mask(mask, shape):
    if mask is None:
        return None
    fit_mask = np.asarray(mask, dtype=bool)
    if fit_mask.shape != tuple(shape):
        return None
    if not np.any(fit_mask):
        return None
    return fit_mask


def solve_flux_baseline(model, data, dataerr=None, mask=None):
    model = np.asarray(model, dtype=float)
    data = np.asarray(data, dtype=float)
    weights = np.ones(model.shape, dtype=float)
    fit_mask = normalized_optional_fit_mask(mask, model.shape)

    if dataerr is not None:
        dataerr = np.asarray(dataerr, dtype=float)
        weights = np.zeros(model.shape, dtype=float)
        valid_err = np.isfinite(dataerr) & (dataerr > 0)
        weights[valid_err] = 1.0 / (dataerr[valid_err] ** 2)

    mask = np.isfinite(model) & np.isfinite(data) & (model != 0)
    if fit_mask is not None:
        mask &= fit_mask
    if dataerr is not None:
        mask &= np.isfinite(weights) & (weights > 0)

    if not np.any(mask):
        return fallback_flux_baseline()

    masked_model = model[mask]
    masked_data = data[mask]
    masked_weights = weights[mask]
    denom = np.sum(masked_weights * masked_model ** 2)

    if not np.isfinite(denom) or denom <= 0:
        ratio = masked_data / masked_model
        ratio = ratio[np.isfinite(ratio)]
        if ratio.size == 0:
            return fallback_flux_baseline()
        baseline = np.nanmedian(ratio)
        return baseline if np.isfinite(baseline) else fallback_flux_baseline()

    baseline = np.sum(masked_weights * masked_data * masked_model) / denom
    return baseline if np.isfinite(baseline) else fallback_flux_baseline()


def solve_flux_baseline_uncertainty(model, dataerr, mask=None):
    if dataerr is None:
        return 0.0
    model = np.asarray(model, dtype=float)
    dataerr = np.asarray(dataerr, dtype=float)
    fit_mask = normalized_optional_fit_mask(mask, model.shape)
    mask = np.isfinite(model) & np.isfinite(dataerr) & (dataerr > 0)
    if fit_mask is not None:
        mask &= fit_mask
    if not np.any(mask):
        return 0.0
    denom = np.sum((model[mask] / dataerr[mask]) ** 2)
    if not np.isfinite(denom) or denom <= 0:
        return 0.0
    return (1.0 / denom) ** 0.5


def mc_a1(m_a2, sig_a2, transit, airmass, data, dataerr=None, n=10000, mask=None):
    n = int(n)
    a2 = np.random.normal(m_a2, sig_a2, n)
    reference = get_airmass_reference(airmass)
    transit = np.asarray(transit, dtype=float)
    data = np.asarray(data, dtype=float)
    airmass = np.asarray(airmass, dtype=float)
    centered_airmass = center_airmass(airmass, reference=reference)
    weights = np.ones(transit.shape[0], dtype=float)
    fit_mask = normalized_optional_fit_mask(mask, transit.shape)

    if dataerr is not None:
        dataerr = np.asarray(dataerr, dtype=float)
        weights = np.zeros(transit.shape[0], dtype=float)
        valid_err = np.isfinite(dataerr) & (dataerr > 0)
        weights[valid_err] = 1.0 / (dataerr[valid_err] ** 2)

    mask = np.isfinite(data) & np.isfinite(transit) & np.isfinite(centered_airmass)
    if fit_mask is not None:
        mask &= fit_mask
    if dataerr is not None:
        mask &= np.isfinite(weights) & (weights > 0)

    if not np.any(mask):
        return fallback_flux_baseline(), 0.0

    masked_transit = transit[mask]
    masked_airmass = centered_airmass[mask]
    masked_data = data[mask]
    masked_weights = weights[mask]
    numer = np.empty(n, dtype=float)
    denom = np.empty(n, dtype=float)
    chunk_size = 1024

    for start in range(0, n, chunk_size):
        stop = min(start + chunk_size, n)
        masked_model = masked_transit * np.exp(np.outer(a2[start:stop], masked_airmass))
        numer[start:stop] = np.sum(masked_weights * masked_data * masked_model, axis=1)
        denom[start:stop] = np.sum(masked_weights * masked_model ** 2, axis=1)
    valid = np.isfinite(numer) & np.isfinite(denom) & (denom > 0)

    if not np.any(valid):
        best_model = transit * airmass_trend(m_a2, airmass, reference=reference)
        baseline = solve_flux_baseline(best_model, data, dataerr, mask=fit_mask)
        return baseline, solve_flux_baseline_uncertainty(best_model, dataerr, mask=fit_mask)

    baselines = numer[valid] / denom[valid]
    baseline = float(np.nanmean(baselines))
    baseline_unc = float(np.nanstd(baselines))

    if baseline_unc == 0.0:
        best_model = transit * airmass_trend(m_a2, airmass, reference=reference)
        baseline_unc = solve_flux_baseline_uncertainty(best_model, dataerr, mask=fit_mask)

    return baseline, baseline_unc


def round_to_2(*args):
    x = args[0]
    if len(args) == 1:
        y = args[0]
    else:
        y = args[1]
    if np.floor(y) >= 1.:
        roundval = 2
    else:
        try:
            roundval = -int(np.floor(np.log10(abs(y)))) + 1
        except:
            roundval = 1
    return round(x, roundval)


# average data into bins of dt from start to finish
def time_bin(time, flux, dt=1. / (60 * 24)):
    bins = int(np.floor((max(time) - min(time)) / dt))
    bflux = np.zeros(bins)
    btime = np.zeros(bins)
    bstds = np.zeros(bins)
    for i in range(bins):
        mask = (time >= (min(time) + i * dt)) & (time < (min(time) + (i + 1) * dt))
        if mask.sum() > 0:
            bflux[i] = bn.nanmean(flux[mask])
            btime[i] = bn.nanmean(time[mask])
            bstds[i] = bn.nanstd(flux[mask]) / (mask.sum() ** 0.5)
    zmask = (bflux == 0) | (btime == 0) | np.isnan(bflux) | np.isnan(btime)
    return btime[~zmask], bflux[~zmask], bstds[~zmask]


# Function that bins an array
def binner(arr, n, err=''):
    if len(err) == 0:
        ecks = np.pad(arr.astype(float), (0, ((n - arr.size % n) % n)), mode='constant',
                      constant_values=np.NaN).reshape(-1, n)
        arr = bn.nanmean(ecks, axis=1)
        return arr
    else:
        ecks = np.pad(arr.astype(float), (0, ((n - arr.size % n) % n)), mode='constant',
                      constant_values=np.NaN).reshape(-1, n)
        why = np.pad(err.astype(float), (0, ((n - err.size % n) % n)), mode='constant', constant_values=np.NaN).reshape(
            -1, n)
        weights = 1. / (why ** 2.)
        # Calculate the weighted average
        arr = bn.nansum(ecks * weights, axis=1) / bn.nansum(weights, axis=1)
        err = np.array([np.sqrt(1. / bn.nansum(1. / (np.array(i) ** 2.))) for i in why])
        return arr, err


class lc_fitter(object):

    def __init__(
        self,
        time,
        data,
        dataerr,
        airmass,
        prior,
        bounds,
        neighbors=200,
        mode='ns',
        jd_times=None,
        verbose=True,
        use_impactparameter_rather_than_inclination_to_fit=True,
        duration_prior=None,
        keep_ultranest_sampler=False,
        baseline_fit_mask=None,
        fixed_parameter_errors=None,
        fixed_flux_baseline=False,
        ultranest_min_num_live_points=None,
    ):
        self.time = time
        self.data = data
        self.dataerr = dataerr
        self.airmass = airmass
        self.airmass_reference = get_airmass_reference(airmass)
        self.prior = prior
        self.bounds = bounds
        self.max_ncalls = 2e5
        self.verbose = verbose
        self.jd_times = jd_times
        self.mode = mode
        self.neighbors = neighbors
        self.use_impactparameter_rather_than_inclination_to_fit = use_impactparameter_rather_than_inclination_to_fit
        self.duration_prior = copy.deepcopy(duration_prior) if isinstance(duration_prior, dict) else None
        self.keep_ultranest_sampler = bool(keep_ultranest_sampler)
        self.baseline_fit_mask = self._coerce_baseline_fit_mask(baseline_fit_mask)
        self.fixed_parameter_errors = (
            copy.deepcopy(fixed_parameter_errors)
            if isinstance(fixed_parameter_errors, dict)
            else {}
        )
        self.fixed_flux_baseline = bool(fixed_flux_baseline)
        self.ultranest_min_num_live_points = ultranest_min_num_live_points
        self._ultranest_resume_context = None
        self.results = None
        self.sampled_keys = list(bounds.keys())
        self.sample_bounds = copy.deepcopy(bounds)
        self.impact_parameter_sampled_directly = False
        self.sample_parameters = {}
        self.sample_errors = {}
        self.sample_quantiles = {}
        self.nested_fit_fallback = False
        self.nested_fit_failure_reason = None
        if self.mode == "lm":
            self.fit_LM()
        elif self.mode == "ns":
            try:
                self.fit_nested()
            except np.linalg.LinAlgError as exc:
                self.nested_fit_fallback = True
                self.nested_fit_failure_reason = f"{type(exc).__name__}: {exc}"
                self.ns_type = 'lm'
                self.mode = "lm"
                if self.verbose:
                    print(
                        "WARNING: Nested light curve fitting failed with a linear algebra error; "
                        "falling back to least-squares fit."
                    )
                    print(f"  Reason: {self.nested_fit_failure_reason}")
                self.fit_LM()

    def _validate_flux_baseline_keys(self):
        free_flux_keys = [key for key in self.bounds if key in ('a0', 'a1')]
        if len(free_flux_keys) > 1:
            raise ValueError("Use only one of 'a0' or 'a1' as a free baseline parameter.")

    def _has_free_flux_baseline(self):
        return has_explicit_flux_baseline(getattr(self, 'bounds', {}))

    def _uses_fixed_flux_baseline(self):
        return bool(getattr(self, 'fixed_flux_baseline', False))

    def _uses_analytic_flux_baseline(self):
        return (
            np.ndim(getattr(self, 'airmass', np.array([]))) != 2
            and not self._has_free_flux_baseline()
            and not self._uses_fixed_flux_baseline()
            and hasattr(self, 'time')
            and hasattr(self, 'data')
            and hasattr(self, 'dataerr')
        )

    def _set_flux_baseline(self, value, error=0.0):
        self.parameters['a0'] = value
        self.errors['a0'] = error
        self.parameters['a1'] = value
        self.errors['a1'] = error

    def _values_with_analytic_flux_baseline(self, values):
        values = copy.deepcopy(values)
        if not self._uses_analytic_flux_baseline():
            return values

        try:
            model = transit(self.time, values)
            model = np.asarray(model, dtype=float) * airmass_trend(
                values.get('a2', 0),
                self.airmass,
                reference=self._get_airmass_reference(),
            )
            flux_scale = solve_flux_baseline(
                model,
                self.data,
                self.dataerr,
                mask=self._get_baseline_fit_mask(),
            )
        except Exception:
            return values

        values['a0'] = flux_scale
        values['a1'] = flux_scale
        return values

    def _coerce_baseline_fit_mask(self, baseline_fit_mask):
        fit_mask = normalized_optional_fit_mask(baseline_fit_mask, np.asarray(self.time).shape)
        if fit_mask is None:
            return None
        return fit_mask

    def _get_baseline_fit_mask(self):
        return getattr(self, 'baseline_fit_mask', None)

    def _apply_fixed_parameter_errors(self):
        fixed_errors = getattr(self, 'fixed_parameter_errors', None)
        if not isinstance(fixed_errors, dict):
            return
        if not hasattr(self, 'parameters') or not isinstance(self.parameters, dict):
            return
        if not hasattr(self, 'errors') or not isinstance(self.errors, dict):
            self.errors = {}
        if not hasattr(self, 'quantiles') or not isinstance(self.quantiles, dict):
            self.quantiles = {}

        for key, error in fixed_errors.items():
            if key not in self.parameters or key in self.errors:
                continue
            try:
                error = float(error)
            except (TypeError, ValueError):
                continue
            if not np.isfinite(error) or error < 0:
                continue
            self.errors[key] = error
            self.quantiles[key] = [-error, error]

    def _get_airmass_reference(self):
        if hasattr(self, 'airmass_reference'):
            return self.airmass_reference
        return get_airmass_reference(self.airmass)

    def _get_plot_time_range(self):
        plot_time_range = normalize_time_range(getattr(self, 'plot_time_range', None))
        if plot_time_range is not None:
            return plot_time_range
        return normalize_time_range(self.time)

    def _update_plot_geometry(self):
        plot_time_range = self._get_plot_time_range()
        self.phase = get_plot_phase(self.time, self.parameters['per'], self.parameters['tmid'], plot_time_range)

        if plot_time_range is None:
            self.time_upsample = np.linspace(min(self.time), max(self.time), 1000)
        else:
            self.time_upsample = np.linspace(plot_time_range[0], plot_time_range[1], 1000)

        self.transit_upsample = transit(self.time_upsample, self.parameters)
        self.phase_upsample = get_plot_phase(
            self.time_upsample,
            self.parameters['per'],
            self.parameters['tmid'],
            plot_time_range,
        )

    def _build_systematics_model(self, values):
        return get_flux_baseline(values) * airmass_trend(
            values.get('a2', 0),
            self.airmass,
            reference=self._get_airmass_reference(),
        )

    def _build_systematics_model_at(self, values, times=None):
        if times is None:
            return self._build_systematics_model(values)

        times = np.asarray(times, dtype=float)
        if np.ndim(self.airmass) == 2:
            return np.full(times.shape, get_flux_baseline(values), dtype=float)

        source_times = np.asarray(self.time, dtype=float)
        source_airmass = np.asarray(self.airmass, dtype=float)
        finite = np.isfinite(source_times) & np.isfinite(source_airmass)
        if times.shape == source_times.shape and np.allclose(times, source_times, rtol=0.0, atol=0.0):
            airmass_values = source_airmass
        elif np.count_nonzero(finite) >= 2:
            order = np.argsort(source_times[finite])
            airmass_values = np.interp(
                times,
                source_times[finite][order],
                source_airmass[finite][order],
                left=source_airmass[finite][order][0],
                right=source_airmass[finite][order][-1],
            )
        elif np.count_nonzero(finite) == 1:
            airmass_values = np.full(times.shape, source_airmass[finite][0], dtype=float)
        else:
            airmass_values = np.zeros(times.shape, dtype=float)

        return get_flux_baseline(values) * airmass_trend(
            values.get('a2', 0),
            airmass_values,
            reference=self._get_airmass_reference(),
        )

    def _get_perturbed_transit_parameter_value(self, key, value):
        try:
            value = float(value)
        except (TypeError, ValueError):
            return np.nan

        if key in ('rprs', 'ars', 'per', 'a0', 'a1'):
            return max(value, np.finfo(float).eps)
        if key == 'ecc':
            return float(np.clip(value, 0.0, 0.999999))
        if key == 'inc':
            return float(np.clip(value, 0.0, 180.0))
        return value

    def _normalized_model_for_plot_times(self, times, values):
        values = self._values_with_analytic_flux_baseline(values)
        model = np.asarray(transit(times, values), dtype=float)
        if np.ndim(self.airmass) == 2:
            return model

        try:
            sample_systematics = self._build_systematics_model_at(values, times)
            best_systematics = self._build_systematics_model_at(self.parameters, times)
        except Exception:
            return model

        with np.errstate(divide='ignore', invalid='ignore'):
            normalized_model = model * sample_systematics / best_systematics
        if normalized_model.shape != model.shape or not np.any(np.isfinite(normalized_model)):
            return model
        return normalized_model

    def _posterior_model_uncertainty(self, times, sigma=1.0):
        if getattr(self, 'results', None) is None:
            return None

        try:
            sample_points, sample_logl, sample_weights = self._get_triangle_plot_samples()
        except Exception:
            return None

        sample_points = np.asarray(sample_points, dtype=float)
        if sample_points.ndim != 2 or sample_points.shape[0] < 2:
            return None

        finite_rows = np.all(np.isfinite(sample_points), axis=1)
        if sample_logl is not None:
            sample_logl = np.asarray(sample_logl, dtype=float)
            if sample_logl.shape[0] == sample_points.shape[0]:
                finite_rows &= np.isfinite(sample_logl)

        if np.count_nonzero(finite_rows) < 2:
            return None

        row_indices = np.flatnonzero(finite_rows)
        if row_indices.size > MODEL_UNCERTAINTY_POSTERIOR_SAMPLE_LIMIT:
            if sample_weights is not None:
                weights_array = np.asarray(sample_weights, dtype=float)
                if weights_array.shape[0] == sample_points.shape[0]:
                    row_weights = np.where(np.isfinite(weights_array[row_indices]), weights_array[row_indices], 0.0)
                    order = np.argsort(row_weights)[-MODEL_UNCERTAINTY_POSTERIOR_SAMPLE_LIMIT:]
                    row_indices = row_indices[np.sort(order)]
                else:
                    row_indices = row_indices[
                        np.linspace(0, row_indices.size - 1, MODEL_UNCERTAINTY_POSTERIOR_SAMPLE_LIMIT).astype(int)
                    ]
            else:
                row_indices = row_indices[
                    np.linspace(0, row_indices.size - 1, MODEL_UNCERTAINTY_POSTERIOR_SAMPLE_LIMIT).astype(int)
                ]

        selected_points = sample_points[row_indices]
        selected_weights = None
        if sample_weights is not None:
            weights_array = np.asarray(sample_weights, dtype=float)
            if weights_array.shape[0] == sample_points.shape[0]:
                selected_weights = weights_array[row_indices]
                selected_weights = np.where(np.isfinite(selected_weights) & (selected_weights >= 0), selected_weights, 0.0)
                if np.sum(selected_weights) <= 0:
                    selected_weights = None

        bound_keys = list(self.bounds.keys())
        sampled_keys = getattr(self, 'sampled_keys', None)
        if sampled_keys is None:
            sampled_keys = self._get_sampled_keys(bound_keys)
        models = []
        for point in selected_points:
            try:
                values = copy.deepcopy(self.parameters)
                values.update(self._physical_values_from_sample_point(point, bound_keys, sampled_keys))
                model = self._normalized_model_for_plot_times(times, values)
            except Exception:
                continue
            if model.shape == times.shape and np.all(np.isfinite(model)):
                models.append(model)

        if len(models) < 2:
            return None

        model_grid = np.asarray(models, dtype=float)
        if selected_weights is not None and selected_weights.shape[0] != model_grid.shape[0]:
            selected_weights = None

        try:
            sigma = float(sigma)
        except (TypeError, ValueError):
            sigma = 1.0
        if not np.isfinite(sigma) or sigma <= 0:
            sigma = 1.0
        coverage = math.erf(sigma / np.sqrt(2.0))
        q_lower = 0.5 * (1.0 - coverage)
        q_upper = 1.0 - q_lower

        if selected_weights is None:
            lower, median, upper = np.nanpercentile(
                model_grid,
                [100.0 * q_lower, 50.0, 100.0 * q_upper],
                axis=0,
            )
        else:
            lower = np.array([
                self._weighted_quantiles(model_grid[:, i], [q_lower], weights=selected_weights)[0]
                for i in range(model_grid.shape[1])
            ])
            median = np.array([
                self._weighted_quantiles(model_grid[:, i], [0.5], weights=selected_weights)[0]
                for i in range(model_grid.shape[1])
            ])
            upper = np.array([
                self._weighted_quantiles(model_grid[:, i], [q_upper], weights=selected_weights)[0]
                for i in range(model_grid.shape[1])
            ])

        try:
            best_model = self._normalized_model_for_plot_times(times, self.parameters)
        except Exception:
            best_model = median

        lower_width = median - lower
        upper_width = upper - median
        lower = best_model - np.maximum(lower_width, 0.0)
        upper = best_model + np.maximum(upper_width, 0.0)

        finite = np.isfinite(lower) & np.isfinite(upper) & (lower <= upper)
        if not np.any(finite):
            return None
        return lower, upper

    def transit_model_uncertainty(self, times=None, sigma=1.0):
        if times is None:
            times = getattr(self, 'time_upsample', self.time)
        times = np.asarray(times, dtype=float)
        if times.size == 0:
            return None

        try:
            model = self._normalized_model_for_plot_times(times, self.parameters)
        except Exception:
            return None

        posterior_envelope = self._posterior_model_uncertainty(times, sigma=sigma)
        if posterior_envelope is not None:
            return posterior_envelope

        sigma = float(sigma)
        variance = np.zeros_like(model, dtype=float)
        uncertainty_keys = list(TRANSIT_MODEL_UNCERTAINTY_KEYS)
        for key in BASELINE_MODEL_UNCERTAINTY_KEYS:
            if key == 'a1' and 'a0' in self.parameters:
                continue
            uncertainty_keys.append(key)

        for key in uncertainty_keys:
            if key not in self.parameters:
                continue
            error = self.errors.get(key)
            try:
                center = float(self.parameters[key])
                error = float(error)
            except (TypeError, ValueError):
                continue
            if not np.isfinite(center) or not np.isfinite(error) or error <= 0:
                continue

            lower_value = self._get_perturbed_transit_parameter_value(key, center - error)
            upper_value = self._get_perturbed_transit_parameter_value(key, center + error)
            if (
                not np.isfinite(lower_value)
                or not np.isfinite(upper_value)
                or np.isclose(lower_value, upper_value)
            ):
                continue

            lower_parameters = copy.deepcopy(self.parameters)
            upper_parameters = copy.deepcopy(self.parameters)
            lower_parameters[key] = lower_value
            upper_parameters[key] = upper_value
            try:
                lower_model = self._normalized_model_for_plot_times(times, lower_parameters)
                upper_model = self._normalized_model_for_plot_times(times, upper_parameters)
            except Exception:
                continue

            derivative = (upper_model - lower_model) / (upper_value - lower_value)
            contribution = derivative * error * sigma
            finite = np.isfinite(contribution)
            variance[finite] += contribution[finite] ** 2

        model_uncertainty = np.sqrt(variance)
        if not np.any(np.isfinite(model_uncertainty) & (model_uncertainty > 0)):
            return None
        return model - model_uncertainty, model + model_uncertainty

    def _posterior_baseline_model_uncertainty(self, times, sigma=1.0):
        if getattr(self, 'results', None) is None or np.ndim(getattr(self, 'airmass', np.array([]))) == 2:
            return None

        try:
            sample_points, sample_logl, sample_weights = self._get_triangle_plot_samples()
        except Exception:
            return None

        sample_points = np.asarray(sample_points, dtype=float)
        if sample_points.ndim != 2 or sample_points.shape[0] < 2:
            return None

        finite_rows = np.all(np.isfinite(sample_points), axis=1)
        if sample_logl is not None:
            sample_logl = np.asarray(sample_logl, dtype=float)
            if sample_logl.shape[0] == sample_points.shape[0]:
                finite_rows &= np.isfinite(sample_logl)

        if np.count_nonzero(finite_rows) < 2:
            return None

        row_indices = np.flatnonzero(finite_rows)
        if row_indices.size > MODEL_UNCERTAINTY_POSTERIOR_SAMPLE_LIMIT:
            if sample_weights is not None:
                weights_array = np.asarray(sample_weights, dtype=float)
                if weights_array.shape[0] == sample_points.shape[0]:
                    row_weights = np.where(np.isfinite(weights_array[row_indices]), weights_array[row_indices], 0.0)
                    order = np.argsort(row_weights)[-MODEL_UNCERTAINTY_POSTERIOR_SAMPLE_LIMIT:]
                    row_indices = row_indices[np.sort(order)]
                else:
                    row_indices = row_indices[
                        np.linspace(0, row_indices.size - 1, MODEL_UNCERTAINTY_POSTERIOR_SAMPLE_LIMIT).astype(int)
                    ]
            else:
                row_indices = row_indices[
                    np.linspace(0, row_indices.size - 1, MODEL_UNCERTAINTY_POSTERIOR_SAMPLE_LIMIT).astype(int)
                ]

        selected_points = sample_points[row_indices]
        selected_weights = None
        if sample_weights is not None:
            weights_array = np.asarray(sample_weights, dtype=float)
            if weights_array.shape[0] == sample_points.shape[0]:
                selected_weights = weights_array[row_indices]
                selected_weights = np.where(
                    np.isfinite(selected_weights) & (selected_weights >= 0),
                    selected_weights,
                    0.0,
                )
                if np.sum(selected_weights) <= 0:
                    selected_weights = None

        try:
            best_parameters = self._values_with_analytic_flux_baseline(self.parameters)
            best_systematics = self._build_systematics_model_at(best_parameters, times)
        except Exception:
            return None

        best_systematics = np.asarray(best_systematics, dtype=float)
        if best_systematics.shape != times.shape or not np.all(np.isfinite(best_systematics)):
            return None

        bound_keys = list(getattr(self, 'bounds', {}).keys())
        sampled_keys = getattr(self, 'sampled_keys', None)
        if sampled_keys is None:
            sampled_keys = self._get_sampled_keys(bound_keys)

        ratios = []
        ratio_weights = []
        for index, point in enumerate(selected_points):
            try:
                values = copy.deepcopy(self.parameters)
                values.update(self._physical_values_from_sample_point(point, bound_keys, sampled_keys))
                values = self._values_with_analytic_flux_baseline(values)
                sample_systematics = self._build_systematics_model_at(values, times)
                with np.errstate(divide='ignore', invalid='ignore'):
                    ratio = np.asarray(sample_systematics, dtype=float) / best_systematics
            except Exception:
                continue
            if ratio.shape != times.shape or not np.all(np.isfinite(ratio)):
                continue
            ratios.append(ratio)
            if selected_weights is not None:
                ratio_weights.append(selected_weights[index])

        if len(ratios) < 2:
            return None

        ratio_grid = np.asarray(ratios, dtype=float)
        if selected_weights is not None:
            selected_weights = np.asarray(ratio_weights, dtype=float)
            if selected_weights.shape[0] != ratio_grid.shape[0] or np.sum(selected_weights) <= 0:
                selected_weights = None

        try:
            sigma = float(sigma)
        except (TypeError, ValueError):
            sigma = 1.0
        if not np.isfinite(sigma) or sigma <= 0:
            sigma = 1.0
        coverage = math.erf(sigma / np.sqrt(2.0))
        q_lower = 0.5 * (1.0 - coverage)
        q_upper = 1.0 - q_lower

        if selected_weights is None:
            lower, median, upper = np.nanpercentile(
                ratio_grid,
                [100.0 * q_lower, 50.0, 100.0 * q_upper],
                axis=0,
            )
        else:
            lower = np.array([
                self._weighted_quantiles(ratio_grid[:, i], [q_lower], weights=selected_weights)[0]
                for i in range(ratio_grid.shape[1])
            ])
            median = np.array([
                self._weighted_quantiles(ratio_grid[:, i], [0.5], weights=selected_weights)[0]
                for i in range(ratio_grid.shape[1])
            ])
            upper = np.array([
                self._weighted_quantiles(ratio_grid[:, i], [q_upper], weights=selected_weights)[0]
                for i in range(ratio_grid.shape[1])
            ])

        lower_width = median - lower
        upper_width = upper - median
        lower = 1.0 - np.maximum(lower_width, 0.0)
        upper = 1.0 + np.maximum(upper_width, 0.0)
        finite = np.isfinite(lower) & np.isfinite(upper) & (lower <= upper)
        if not np.any(finite):
            return None
        return lower, upper

    def baseline_model_uncertainty(self, times=None, sigma=1.0):
        if times is None:
            times = getattr(self, 'time_upsample', self.time)
        times = np.asarray(times, dtype=float)
        if times.size == 0 or np.ndim(getattr(self, 'airmass', np.array([]))) == 2:
            return None

        posterior_envelope = self._posterior_baseline_model_uncertainty(times, sigma=sigma)
        if posterior_envelope is not None:
            return posterior_envelope

        try:
            best_parameters = self._values_with_analytic_flux_baseline(self.parameters)
            best_systematics = self._build_systematics_model_at(best_parameters, times)
        except Exception:
            return None

        if (
            np.asarray(best_systematics).shape != times.shape
            or not np.any(np.isfinite(best_systematics))
        ):
            return None

        try:
            sigma = float(sigma)
        except (TypeError, ValueError):
            sigma = 1.0
        if not np.isfinite(sigma) or sigma <= 0:
            sigma = 1.0

        variance = np.zeros_like(times, dtype=float)
        uses_analytic_flux_baseline = self._uses_analytic_flux_baseline()
        a2_error = self.errors.get('a2')
        try:
            a2_error = float(a2_error)
        except (TypeError, ValueError):
            a2_error = 0.0
        for key in ('a0', 'a1', 'a2'):
            if key == 'a1' and 'a0' in self.parameters:
                continue
            if key in ('a0', 'a1') and uses_analytic_flux_baseline and a2_error > 0:
                continue
            if key not in self.parameters:
                continue
            error = self.errors.get(key)
            try:
                center = float(self.parameters[key])
                error = float(error)
            except (TypeError, ValueError):
                continue
            if not np.isfinite(center) or not np.isfinite(error) or error <= 0:
                continue

            lower_value = self._get_perturbed_transit_parameter_value(key, center - error)
            upper_value = self._get_perturbed_transit_parameter_value(key, center + error)
            if (
                not np.isfinite(lower_value)
                or not np.isfinite(upper_value)
                or np.isclose(lower_value, upper_value)
            ):
                continue

            lower_parameters = copy.deepcopy(self.parameters)
            upper_parameters = copy.deepcopy(self.parameters)
            lower_parameters[key] = lower_value
            upper_parameters[key] = upper_value
            lower_parameters = self._values_with_analytic_flux_baseline(lower_parameters)
            upper_parameters = self._values_with_analytic_flux_baseline(upper_parameters)
            try:
                lower_systematics = self._build_systematics_model_at(lower_parameters, times)
                upper_systematics = self._build_systematics_model_at(upper_parameters, times)
            except Exception:
                continue

            with np.errstate(divide='ignore', invalid='ignore'):
                lower_ratio = lower_systematics / best_systematics
                upper_ratio = upper_systematics / best_systematics
            derivative = (upper_ratio - lower_ratio) / (upper_value - lower_value)
            contribution = derivative * error * sigma
            finite = np.isfinite(contribution)
            variance[finite] += contribution[finite] ** 2

        baseline_uncertainty = np.sqrt(variance)
        if not np.any(np.isfinite(baseline_uncertainty) & (baseline_uncertainty > 0)):
            return None
        return 1.0 - baseline_uncertainty, 1.0 + baseline_uncertainty

    def _plot_transit_model_uncertainty(self, ax, x_values, times, sort_index, label=None):
        envelope = self.transit_model_uncertainty(times)
        if envelope is None:
            return None

        lower, upper = envelope
        x_values = np.asarray(x_values, dtype=float)
        sort_index = np.asarray(sort_index, dtype=int)
        x_sorted = x_values[sort_index]
        lower_sorted = np.asarray(lower, dtype=float)[sort_index]
        upper_sorted = np.asarray(upper, dtype=float)[sort_index]
        band = ax.fill_between(
            x_sorted,
            lower_sorted,
            upper_sorted,
            color='red',
            alpha=0.16,
            linewidth=0,
            zorder=2.5,
            label=label,
        )
        ax.plot(
            x_sorted,
            lower_sorted,
            color='red',
            linestyle='--',
            linewidth=0.9,
            alpha=0.72,
            zorder=3.4,
        )
        ax.plot(
            x_sorted,
            upper_sorted,
            color='red',
            linestyle='--',
            linewidth=0.9,
            alpha=0.72,
            zorder=3.4,
        )
        return band

    def _plot_baseline_model_uncertainty(self, ax, x_values, times, sort_index, label=None):
        envelope = self.baseline_model_uncertainty(times)
        if envelope is None:
            return None

        lower, upper = envelope
        x_values = np.asarray(x_values, dtype=float)
        sort_index = np.asarray(sort_index, dtype=int)
        x_sorted = x_values[sort_index]
        lower_sorted = np.asarray(lower, dtype=float)[sort_index]
        upper_sorted = np.asarray(upper, dtype=float)[sort_index]
        band = ax.fill_between(
            x_sorted,
            lower_sorted,
            upper_sorted,
            color='gold',
            alpha=0.32,
            linewidth=0,
            zorder=2.1,
            label=label,
        )
        ax.plot(
            x_sorted,
            lower_sorted,
            color='gold',
            linestyle='--',
            linewidth=0.9,
            alpha=0.78,
            zorder=3.2,
        )
        ax.plot(
            x_sorted,
            upper_sorted,
            color='gold',
            linestyle='--',
            linewidth=0.9,
            alpha=0.78,
            zorder=3.2,
        )
        return band

    def _uses_internal_impact_parameter(self):
        return (
            self.use_impactparameter_rather_than_inclination_to_fit
            and self.mode == "ns"
            and 'inc' in self.bounds
            and 'b' not in self.bounds
        )

    def _get_sampled_keys(self, bound_keys=None):
        bound_keys = list(self.bounds.keys()) if bound_keys is None else list(bound_keys)
        if not self._uses_internal_impact_parameter():
            return bound_keys
        return ['b' if key == 'inc' else key for key in bound_keys]

    def _get_impact_parameter_scale_upper_bound(self, values):
        values = dict(values)
        scale_keys = ('ars', 'ecc', 'omega')
        endpoint_sets = []
        for key in scale_keys:
            if key in self.bounds:
                endpoints = np.asarray(self.bounds[key], dtype=float).reshape(-1)[:2]
            else:
                endpoints = np.asarray([values.get(key, 0.0)], dtype=float)
            finite_endpoints = [float(value) for value in endpoints if np.isfinite(value)]
            if not finite_endpoints:
                return np.nan
            endpoint_sets.append((key, finite_endpoints))

        scales = []
        for candidate_values in product(*[endpoints for _, endpoints in endpoint_sets]):
            candidate = dict(values)
            for key, value in zip([key for key, _ in endpoint_sets], candidate_values):
                candidate[key] = value
            try:
                scale = float(impact_parameter_scale(candidate))
            except (KeyError, TypeError, ValueError):
                continue
            if np.isfinite(scale) and scale > 0:
                scales.append(scale)

        return float(max(scales)) if scales else np.nan

    def _get_impact_parameter_sampling_bounds(self, values=None, use_search_bounds=False):
        values = self.prior if values is None else values
        if use_search_bounds and 'rprs' in self.bounds:
            values = dict(values)
            rprs_bounds = np.asarray(self.bounds['rprs'], dtype=float).reshape(-1)[:2]
            finite_rprs = rprs_bounds[np.isfinite(rprs_bounds) & (rprs_bounds >= 0)]
            if finite_rprs.size > 0:
                values['rprs'] = float(np.max(finite_rprs))

        grazing_upper = grazing_impact_parameter(values)
        if use_search_bounds:
            scale_upper = self._get_impact_parameter_scale_upper_bound(values)
        else:
            try:
                scale_upper = float(impact_parameter_scale(values))
            except (KeyError, TypeError, ValueError):
                scale_upper = np.nan

        upper_candidates = [
            float(value)
            for value in (grazing_upper, scale_upper)
            if np.isfinite(value) and value > 0
        ]
        upper = min(upper_candidates) if upper_candidates else 1.0
        return [0.0, float(max(0.0, upper))]

    def _get_impact_parameter_upper_bounds_for_sample_points(self, sample_points, bound_keys):
        sample_points = np.atleast_2d(np.asarray(sample_points, dtype=float))
        bound_index = {key: index for index, key in enumerate(bound_keys)}
        sample_count = sample_points.shape[0]

        def values_for(key, default):
            if key in bound_index:
                return sample_points[:, bound_index[key]]
            value = np.asarray(self.prior.get(key, default), dtype=float)
            if value.shape == ():
                return np.full(sample_count, float(value), dtype=float)
            return np.broadcast_to(value, (sample_count,)).astype(float)

        rprs = values_for('rprs', np.nan)
        grazing_upper = np.where(np.isfinite(rprs) & (rprs >= 0), 1.0 + rprs, np.nan)

        ars = values_for('ars', np.nan)
        ecc = values_for('ecc', 0.0)
        omega = np.deg2rad(values_for('omega', 0.0))
        denom = 1.0 + ecc * np.sin(omega)
        denom = np.where(np.isclose(denom, 0.0), np.finfo(float).eps, denom)
        scale_upper = ars * (1.0 - ecc ** 2) / denom

        valid_grazing = np.isfinite(grazing_upper) & (grazing_upper > 0)
        valid_scale = np.isfinite(scale_upper) & (scale_upper > 0)
        upper = np.full(sample_count, 1.0, dtype=float)

        both_valid = valid_grazing & valid_scale
        upper[both_valid] = np.minimum(grazing_upper[both_valid], scale_upper[both_valid])

        grazing_only = valid_grazing & ~valid_scale
        upper[grazing_only] = grazing_upper[grazing_only]

        scale_only = valid_scale & ~valid_grazing
        upper[scale_only] = scale_upper[scale_only]

        return np.maximum(0.0, upper)

    def _get_sample_bounds(self, bound_keys=None, values=None):
        bound_keys = list(self.bounds.keys()) if bound_keys is None else list(bound_keys)
        sampled_keys = self._get_sampled_keys(bound_keys)
        values = self.prior if values is None else values
        sample_bounds = {}
        for key, sampled_key in zip(bound_keys, sampled_keys):
            if key == 'inc' and sampled_key == 'b':
                sample_bounds[sampled_key] = self._get_impact_parameter_sampling_bounds(
                    values,
                    use_search_bounds=True,
                )
            else:
                sample_bounds[sampled_key] = list(self.bounds[key])
        return sample_bounds

    def _sample_point_from_unit_cube(self, upars, bound_keys=None):
        bound_keys = list(self.bounds.keys()) if bound_keys is None else list(bound_keys)
        upars_array = np.asarray(upars, dtype=float)
        boundarray = np.array([self.bounds[k] for k in bound_keys], dtype=float)
        lower_bounds = boundarray[:, 0]
        bound_widths = boundarray[:, 1] - lower_bounds
        uses_internal_impact_parameter = self._uses_internal_impact_parameter()
        inc_indices = [i for i, key in enumerate(bound_keys) if key == 'inc']

        sample_point = lower_bounds + bound_widths * upars_array
        if not uses_internal_impact_parameter or not inc_indices:
            return sample_point

        if len(inc_indices) == 1:
            inc_index = inc_indices[0]
            if upars_array.ndim == 2:
                upper_bounds = self._get_impact_parameter_upper_bounds_for_sample_points(
                    sample_point,
                    bound_keys,
                )
                sample_point[:, inc_index] = upper_bounds * upars_array[:, inc_index]
                return sample_point

            upper_bound = self._get_impact_parameter_upper_bounds_for_sample_points(
                sample_point.reshape(1, -1),
                bound_keys,
            )[0]
            sample_point[inc_index] = upper_bound * upars_array[inc_index]
            return sample_point

        if upars_array.ndim == 2:
            for row_index, row_sample_point in enumerate(sample_point):
                physical = dict(self.prior)
                for i, key in enumerate(bound_keys):
                    if i not in inc_indices:
                        physical[key] = row_sample_point[i]
                for i in inc_indices:
                    b_lower, b_upper = self._get_impact_parameter_sampling_bounds(physical)
                    row_sample_point[i] = b_lower + (b_upper - b_lower) * upars_array[row_index, i]
            return sample_point

        physical = dict(self.prior)
        for i, key in enumerate(bound_keys):
            if i not in inc_indices:
                physical[key] = sample_point[i]
        for i in inc_indices:
            b_lower, b_upper = self._get_impact_parameter_sampling_bounds(physical)
            sample_point[i] = b_lower + (b_upper - b_lower) * upars_array[i]

        return sample_point

    def _physical_values_from_sample_point(self, sample_point, bound_keys=None, sampled_keys=None):
        bound_keys = list(self.bounds.keys()) if bound_keys is None else list(bound_keys)
        sampled_keys = self._get_sampled_keys(bound_keys) if sampled_keys is None else list(sampled_keys)
        physical = dict(self.prior)
        impact_parameter = None

        for value, bound_key, sampled_key in zip(sample_point, bound_keys, sampled_keys):
            if sampled_key == 'b' and bound_key == 'inc':
                impact_parameter = value
                continue
            physical[bound_key] = value

        if impact_parameter is not None:
            physical['b'] = impact_parameter
            physical['inc'] = float(inclination_from_impact_parameter(physical, impact_parameter))

        return physical

    def _summarize_derived_parameter(self, samples, point_estimate):
        samples = np.asarray(samples, dtype=float)
        center = float(point_estimate)
        std = float(np.nanstd(samples))
        lower = float(np.nanpercentile(samples, 16))
        upper = float(np.nanpercentile(samples, 84))
        return center, std, [lower - center, upper - center]

    def _get_ultranest_weighted_sample_arrays(self):
        try:
            weighted_samples = self.results['weighted_samples']
            points = np.asarray(weighted_samples['points'], dtype=float)
            logl = np.asarray(weighted_samples['logl'], dtype=float)
        except Exception:
            return None, None

        if points.ndim != 2 or points.shape[0] == 0:
            return None, None
        if logl.shape[0] != points.shape[0]:
            return None, None
        return points, logl

    def _loglike_neighborhood_uncertainty(self, parameter_index, center, minimum_count=8, points=None, logl=None):
        if points is None or logl is None:
            points, logl = self._get_ultranest_weighted_sample_arrays()
        if points is None or parameter_index >= points.shape[1]:
            return None

        values = np.asarray(points[:, parameter_index], dtype=float)
        finite = np.isfinite(values) & np.isfinite(logl)
        if np.count_nonzero(finite) < 2:
            return None

        finite_values = values[finite]
        finite_logl = logl[finite]
        max_logl = float(np.nanmax(finite_logl))
        if not np.isfinite(max_logl):
            return None

        selected_values = None
        selected_delta = np.inf
        for delta_chi2 in (1.0, 4.0, 9.0, 16.0, 25.0, np.inf):
            if np.isfinite(delta_chi2):
                mask = 2.0 * (max_logl - finite_logl) <= delta_chi2
            else:
                mask = np.ones(finite_logl.shape, dtype=bool)
            if np.count_nonzero(mask) >= minimum_count or delta_chi2 == np.inf:
                selected_values = finite_values[mask]
                selected_delta = delta_chi2
                break

        if selected_values is None or selected_values.size < 2:
            return None

        lower, upper = np.nanpercentile(selected_values, [15.8655, 84.1345])
        std = float(np.nanstd(selected_values))
        half_width = float(0.5 * (upper - lower))
        candidates = [value for value in (std, half_width) if np.isfinite(value) and value > 0]
        if not candidates:
            return None

        error = float(max(candidates))
        return {
            'error': error,
            'quantiles': [float(lower), float(upper)],
            'sample_count': int(selected_values.size),
            'delta_chi2': float(selected_delta),
        }

    def _ultranest_error_needs_sample_fallback(self, parameter_index, center, reported_error, points=None):
        if points is None:
            points, _ = self._get_ultranest_weighted_sample_arrays()
        if points is None or parameter_index >= points.shape[1]:
            return False

        values = np.asarray(points[:, parameter_index], dtype=float)
        finite_values = values[np.isfinite(values)]
        if finite_values.size < 2:
            return False

        sample_scale = float(np.nanstd(finite_values))
        if not np.isfinite(sample_scale) or sample_scale <= 0:
            return False

        try:
            reported_error = float(reported_error)
        except (TypeError, ValueError):
            return True

        if not np.isfinite(reported_error) or reported_error <= 0:
            return True

        absolute_floor = max(abs(float(center)) * 1e-12, np.finfo(float).eps)
        return reported_error <= absolute_floor or reported_error < sample_scale * 1e-6

    def _ultranest_error_is_inflated_relative_to_local_fit(self, reported_error, local_uncertainty):
        if not isinstance(local_uncertainty, dict):
            return False

        try:
            reported_error = float(reported_error)
            local_error = float(local_uncertainty.get('error', np.nan))
            delta_chi2 = float(local_uncertainty.get('delta_chi2', np.inf))
        except (TypeError, ValueError):
            return False

        if (
            not np.isfinite(reported_error)
            or reported_error <= 0
            or not np.isfinite(local_error)
            or local_error <= 0
        ):
            return False
        if not np.isfinite(delta_chi2) or delta_chi2 > ULTRANEST_LOCAL_UNCERTAINTY_MAX_DELTA_CHI2:
            return False

        return reported_error > local_error * ULTRANEST_INFLATED_ERROR_REPLACEMENT_FACTOR

    def _get_plot_range(self, key):
        sample_parameters = getattr(self, 'sample_parameters', {})
        sample_errors = getattr(self, 'sample_errors', {})
        sample_bounds = getattr(self, 'sample_bounds', getattr(self, 'bounds', {}))
        center = sample_parameters[key] if key in sample_parameters else self.parameters[key]
        error = sample_errors[key] if key in sample_errors else self.errors[key]

        if isinstance(sample_bounds, dict) and key in sample_bounds:
            try:
                lower, upper = [
                    float(value) for value in np.asarray(sample_bounds[key], dtype=float).reshape(-1)[:2]
                ]
            except (TypeError, ValueError, IndexError):
                lower = np.nan
                upper = np.nan
            if np.isfinite(lower) and np.isfinite(upper) and lower < upper:
                return [lower, upper]

        lower = center - 5 * error
        upper = center + 5 * error
        if np.isfinite(lower) and np.isfinite(upper) and lower < upper:
            return [lower, upper]

        pad = error if np.isfinite(error) and error > 0 else max(abs(center) * 1e-6, 1e-6)
        return [center - pad, center + pad]

    def _expand_plot_range_for_sample_cloud(
        self,
        key,
        plot_range,
        sample_values,
        center,
        required_visible_fraction=1.0,
    ):
        sample_values = np.asarray(sample_values, dtype=float)
        finite_values = sample_values[np.isfinite(sample_values)]
        if finite_values.size < 2:
            return plot_range

        lower, upper = [float(value) for value in plot_range]
        in_range = (finite_values >= lower) & (finite_values <= upper)
        visible_fraction = np.count_nonzero(in_range) / float(finite_values.size)
        if visible_fraction >= required_visible_fraction:
            return plot_range

        new_lower = min(float(np.nanmin(finite_values)), float(center))
        new_upper = max(float(np.nanmax(finite_values)), float(center))
        padding = 0.05 * (new_upper - new_lower)
        if not np.isfinite(padding) or padding <= 0:
            padding = max(abs(float(center)) * 1e-6, 1e-6)
        new_lower -= padding
        new_upper += padding

        sample_bounds = getattr(self, 'sample_bounds', self.bounds)
        if key in sample_bounds:
            bound_lower, bound_upper = sample_bounds[key]
            new_lower = max(new_lower, float(bound_lower))
            new_upper = min(new_upper, float(bound_upper))

        if not np.isfinite(new_lower) or not np.isfinite(new_upper) or new_lower >= new_upper:
            return plot_range
        return [float(new_lower), float(new_upper)]

    def _histogram_edge_peak_fractions(self, sample_values, plot_range, bins, weights=None):
        sample_values = np.asarray(sample_values, dtype=float)
        finite_mask = np.isfinite(sample_values)
        finite_weights = None

        if weights is not None:
            weights = np.asarray(weights, dtype=float)
            if weights.shape == sample_values.shape:
                finite_mask &= np.isfinite(weights) & (weights >= 0)
                finite_weights = weights[finite_mask]
                finite_weight_sum = np.sum(finite_weights)
                if (
                    finite_weights.size == 0
                    or not np.isfinite(finite_weight_sum)
                    or finite_weight_sum <= 0
                ):
                    finite_weights = None

        finite_values = sample_values[finite_mask]
        if finite_values.size < 2:
            return np.nan, np.nan, np.nan

        try:
            lower, upper = [float(value) for value in np.asarray(plot_range, dtype=float).reshape(-1)[:2]]
        except (TypeError, ValueError, IndexError):
            return np.nan, np.nan, np.nan

        if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
            return np.nan, np.nan, np.nan

        bins = max(1, int(bins))
        counts, _ = np.histogram(
            finite_values,
            bins=bins,
            range=(lower, upper),
            weights=finite_weights,
        )
        counts = np.asarray(counts, dtype=float)
        if counts.size == 0:
            return np.nan, np.nan, np.nan

        peak = float(np.nanmax(counts))
        if not np.isfinite(peak) or peak <= 0:
            return np.nan, np.nan, np.nan

        lower_fraction = float(counts[0] / peak)
        upper_fraction = float(counts[-1] / peak)
        return lower_fraction, upper_fraction, peak

    def _get_plot_range_expansion_bounds(self, key, sample_values, center):
        sample_bounds = getattr(self, 'sample_bounds', getattr(self, 'bounds', {}))
        if isinstance(sample_bounds, dict) and key in sample_bounds:
            try:
                bound_lower, bound_upper = [
                    float(value) for value in np.asarray(sample_bounds[key], dtype=float).reshape(-1)[:2]
                ]
            except (TypeError, ValueError, IndexError):
                bound_lower = np.nan
                bound_upper = np.nan

            if np.isfinite(bound_lower) and np.isfinite(bound_upper) and bound_lower < bound_upper:
                fallback_bounds = TRIANGLE_PLOT_FALLBACK_EXPANSION_BOUNDS.get(key)
                if fallback_bounds is not None:
                    fallback_lower, fallback_upper = fallback_bounds
                    if (
                        np.isfinite(fallback_lower)
                        and np.isfinite(fallback_upper)
                        and fallback_lower < fallback_upper
                    ):
                        return [
                            float(min(bound_lower, fallback_lower)),
                            float(max(bound_upper, fallback_upper)),
                        ]
                return [bound_lower, bound_upper]

        fallback_bounds = TRIANGLE_PLOT_FALLBACK_EXPANSION_BOUNDS.get(key)
        if fallback_bounds is not None:
            fallback_lower, fallback_upper = fallback_bounds
            if np.isfinite(fallback_lower) and np.isfinite(fallback_upper) and fallback_lower < fallback_upper:
                return [float(fallback_lower), float(fallback_upper)]

        sample_values = np.asarray(sample_values, dtype=float)
        finite_values = sample_values[np.isfinite(sample_values)]
        try:
            center = float(center)
        except (TypeError, ValueError):
            center = np.nan
        if np.isfinite(center):
            finite_values = np.concatenate([finite_values, [center]])
        if finite_values.size < 2:
            return None

        bound_lower = float(np.nanmin(finite_values))
        bound_upper = float(np.nanmax(finite_values))
        width = bound_upper - bound_lower
        if not np.isfinite(width) or width <= 0:
            padding = max(abs(float(center)) * 1e-6 if np.isfinite(center) else 0.0, 1e-6)
        else:
            padding = 0.05 * width
        return [bound_lower - padding, bound_upper + padding]

    def _expand_plot_range_for_histogram_edge_dropoff(
        self,
        key,
        plot_range,
        sample_values,
        center,
        bins=None,
        weights=None,
        max_edge_peak_fraction=TRIANGLE_PLOT_EDGE_PEAK_FRACTION_MAX,
        minimum_count=TRIANGLE_PLOT_EDGE_MIN_SAMPLE_COUNT,
        max_steps=TRIANGLE_PLOT_EDGE_EXPANSION_STEPS,
    ):
        sample_values = np.asarray(sample_values, dtype=float)
        finite_values = sample_values[np.isfinite(sample_values)]
        if finite_values.size < int(minimum_count):
            return plot_range

        try:
            lower, upper = [float(value) for value in np.asarray(plot_range, dtype=float).reshape(-1)[:2]]
        except (TypeError, ValueError, IndexError):
            return plot_range

        if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
            return plot_range

        expansion_bounds = self._get_plot_range_expansion_bounds(key, finite_values, center)
        if expansion_bounds is None:
            return plot_range

        bound_lower, bound_upper = expansion_bounds
        if not np.isfinite(bound_lower) or not np.isfinite(bound_upper) or bound_lower >= bound_upper:
            return plot_range

        lower = max(lower, bound_lower)
        upper = min(upper, bound_upper)
        if lower >= upper:
            return plot_range

        if bins is None:
            bins = int(np.clip(np.sqrt(finite_values.size), 10, 80))
        bins = max(1, int(bins))
        epsilon = max((bound_upper - bound_lower) * 1e-12, np.finfo(float).eps)

        for _ in range(max(0, int(max_steps)) + 1):
            lower_fraction, upper_fraction, _ = self._histogram_edge_peak_fractions(
                sample_values,
                [lower, upper],
                bins,
                weights=weights,
            )
            if not np.isfinite(lower_fraction) or not np.isfinite(upper_fraction):
                return [float(lower), float(upper)]

            needs_lower = lower_fraction >= max_edge_peak_fraction
            needs_upper = upper_fraction >= max_edge_peak_fraction
            if not needs_lower and not needs_upper:
                return [float(lower), float(upper)]

            width = upper - lower
            if not np.isfinite(width) or width <= 0:
                return [float(lower), float(upper)]

            new_lower = lower
            new_upper = upper
            if needs_lower and lower > bound_lower + epsilon:
                new_lower = max(bound_lower, lower - width)
            if needs_upper and upper < bound_upper - epsilon:
                new_upper = min(bound_upper, upper + width)

            if new_lower == lower and new_upper == upper:
                return [float(lower), float(upper)]

            lower, upper = new_lower, new_upper

        return [float(lower), float(upper)]

    def _get_mirrored_geometry_sample_cloud_range(self, sample_values, center, percentile_padding=0.5):
        sample_values = np.asarray(sample_values, dtype=float)
        finite_values = sample_values[np.isfinite(sample_values)]
        if finite_values.size < 2:
            return None

        try:
            center = float(center)
        except (TypeError, ValueError):
            center = np.nan
        if not np.isfinite(center):
            return None

        percentile_padding = float(percentile_padding)
        percentile_padding = min(max(percentile_padding, 0.0), 49.0)
        q_lower, q_upper = np.nanpercentile(
            finite_values,
            [percentile_padding, 100.0 - percentile_padding],
        )
        plot_lower = min(float(q_lower), center)
        plot_upper = max(float(q_upper), center)
        width = plot_upper - plot_lower
        if not np.isfinite(width) or width <= 0:
            return None

        padding = 0.05 * width
        plot_lower -= padding
        plot_upper += padding
        max_distance = float(np.nanmax(np.abs([plot_lower - center, plot_upper - center])))
        if not np.isfinite(max_distance) or max_distance <= 0:
            return None
        return [-max_distance, max_distance]

    def _get_mirrored_geometry_full_range(
        self,
        key,
        sample_values,
        center,
        sample_weights=None,
    ):
        try:
            center = float(center)
        except (TypeError, ValueError):
            center = np.nan
        if not np.isfinite(center):
            return None

        plot_lower, plot_upper = self._get_plot_range(key)
        plot_lower, plot_upper = self._expand_plot_range_for_sample_cloud(
            key,
            [plot_lower, plot_upper],
            sample_values,
            center,
            required_visible_fraction=1.0,
        )
        plot_bins = int(max(1, np.sqrt(np.asarray(sample_values).size)))
        plot_lower, plot_upper = self._expand_plot_range_for_histogram_edge_dropoff(
            key,
            [plot_lower, plot_upper],
            sample_values,
            center,
            bins=plot_bins,
            weights=sample_weights,
        )

        max_distance = float(np.nanmax(np.abs([plot_lower - center, plot_upper - center])))
        if not np.isfinite(max_distance) or max_distance <= 0:
            finite_offsets = np.asarray(sample_values, dtype=float) - center
            finite_offsets = finite_offsets[np.isfinite(finite_offsets)]
            if finite_offsets.size > 0:
                max_distance = float(np.nanmax(np.abs(finite_offsets)))
        if not np.isfinite(max_distance) or max_distance <= 0:
            max_distance = max(abs(center) * 1e-6, 1e-6)
        return [-max_distance, max_distance]

    def _get_triangle_plot_samples(self):
        if self.ns_type == 'ultranest':
            weighted_samples = self.results['weighted_samples']
            points = np.asarray(weighted_samples['points'], dtype=float)
            logl = np.asarray(weighted_samples['logl'], dtype=float)
            weights = self._get_triangle_plot_sample_weights(
                weighted_samples.get('weights'),
                points.shape[0],
            )
            return points, logl, weights

        raise RuntimeError("Triangle plots require an UltraNest nested-sampling result.")

    def _get_triangle_plot_sample_weights(self, weights, sample_count):
        if weights is None:
            return None

        weights = np.asarray(weights, dtype=float)
        if weights.ndim != 1 or weights.shape[0] != sample_count:
            return None

        finite = np.isfinite(weights) & (weights >= 0)
        if not np.all(finite):
            weights = np.where(finite, weights, 0.0)

        if not np.isfinite(np.sum(weights)) or np.sum(weights) <= 0:
            return None
        return weights

    def get_parameter_posterior_samples(self, key):
        try:
            sample_points, _, _ = self._get_triangle_plot_samples()
        except Exception:
            return np.array([], dtype=float)

        sample_points = np.asarray(sample_points, dtype=float)
        if sample_points.ndim != 2 or sample_points.shape[0] == 0:
            return np.array([], dtype=float)

        sampled_keys = list(getattr(self, 'sampled_keys', self._get_sampled_keys()))
        if key in sampled_keys:
            key_index = sampled_keys.index(key)
            if key_index < sample_points.shape[1]:
                return np.asarray(sample_points[:, key_index], dtype=float)

        bound_keys = list(self.bounds.keys())
        physical_samples = [
            self._physical_values_from_sample_point(point, bound_keys, sampled_keys).get(key, np.nan)
            for point in sample_points
        ]
        return np.asarray(physical_samples, dtype=float)

    def _estimate_histogram_mode(self, samples, bounds=None, bins=None, weights=None):
        samples = np.asarray(samples, dtype=float)
        if weights is None:
            finite_mask = np.isfinite(samples)
            finite_weights = None
        else:
            weights = np.asarray(weights, dtype=float)
            if weights.shape != samples.shape:
                finite_mask = np.isfinite(samples)
                finite_weights = None
            else:
                finite_mask = np.isfinite(samples) & np.isfinite(weights) & (weights >= 0)
                finite_weights = weights[finite_mask]
                if finite_weights.size == 0 or np.sum(finite_weights) <= 0:
                    finite_weights = None

        finite_samples = samples[finite_mask]
        if finite_samples.size == 0:
            return np.nan, np.nan
        if finite_samples.size == 1:
            return float(finite_samples[0]), np.nan

        if bounds is None:
            lower = float(np.nanmin(finite_samples))
            upper = float(np.nanmax(finite_samples))
        else:
            lower, upper = np.asarray(bounds, dtype=float).reshape(-1)[:2]
            if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
                lower = float(np.nanmin(finite_samples))
                upper = float(np.nanmax(finite_samples))

        if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
            return float(np.nanmedian(finite_samples)), np.nan

        if bins is None:
            bins = int(np.clip(np.sqrt(finite_samples.size), 10, 80))
        bins = max(1, int(bins))

        counts, edges = np.histogram(finite_samples, bins=bins, range=(lower, upper), weights=finite_weights)
        if counts.size == 0:
            return float(np.nanmedian(finite_samples)), np.nan
        if not np.any(counts > 0):
            return float(np.nanmedian(finite_samples)), np.nan

        mode_index = int(np.argmax(counts))
        mode = float(0.5 * (edges[mode_index] + edges[mode_index + 1]))
        bin_width = float(edges[1] - edges[0]) if edges.size > 1 else np.nan
        return mode, bin_width

    def _format_triangle_plot_parameter_title(self, value, error):
        try:
            value = float(value)
        except (TypeError, ValueError):
            return "n/a"
        try:
            error = float(error)
        except (TypeError, ValueError):
            error = np.nan
        if not np.isfinite(value):
            return "n/a"
        if not np.isfinite(error) or error < 0:
            return str(round_to_2(value))
        return f"{round_to_2(value, error)} +/- {round_to_2(error)}"

    def _weighted_quantiles(self, values, quantiles, weights=None):
        values = np.asarray(values, dtype=float)
        quantiles = np.asarray(quantiles, dtype=float)
        finite_mask = np.isfinite(values)

        finite_weights = None
        if weights is not None:
            weights = np.asarray(weights, dtype=float)
            if weights.shape == values.shape:
                finite_mask &= np.isfinite(weights) & (weights >= 0)
                finite_weights = weights[finite_mask]
                if finite_weights.size == 0 or np.sum(finite_weights) <= 0:
                    finite_weights = None

        finite_values = values[finite_mask]
        if finite_values.size == 0:
            return np.full(quantiles.shape, np.nan, dtype=float)
        if finite_weights is None:
            return np.nanpercentile(finite_values, 100.0 * quantiles)

        order = np.argsort(finite_values)
        sorted_values = finite_values[order]
        sorted_weights = finite_weights[order]
        cumulative = np.cumsum(sorted_weights)
        total = cumulative[-1]
        if not np.isfinite(total) or total <= 0:
            return np.nanpercentile(finite_values, 100.0 * quantiles)

        cumulative = (cumulative - 0.5 * sorted_weights) / total
        cumulative = np.clip(cumulative, 0.0, 1.0)
        return np.interp(quantiles, cumulative, sorted_values)

    def _triangle_plot_display_estimate(
        self,
        sample_values,
        fallback_center,
        fallback_error,
        plot_range=None,
        weights=None,
        min_informative_peak_ratio=1.5,
        bins=None,
        force_histogram_mode=False,
    ):
        sample_values = np.asarray(sample_values, dtype=float)
        finite_mask = np.isfinite(sample_values)
        finite_weights = None
        if weights is not None:
            weights = np.asarray(weights, dtype=float)
            if weights.shape == sample_values.shape:
                finite_mask &= np.isfinite(weights) & (weights >= 0)
                finite_weights = weights[finite_mask]
                if finite_weights.size == 0 or np.sum(finite_weights) <= 0:
                    finite_weights = None

        finite_values = sample_values[finite_mask]
        if finite_values.size < 2:
            return fallback_center, fallback_error

        q16, q50, q84 = self._weighted_quantiles(
            sample_values,
            [0.158655, 0.5, 0.841345],
            weights=weights,
        )
        estimate = q50
        if plot_range is None:
            bounds = [float(np.nanmin(finite_values)), float(np.nanmax(finite_values))]
        else:
            try:
                bounds = [float(value) for value in np.asarray(plot_range, dtype=float).reshape(-1)[:2]]
            except (TypeError, ValueError, IndexError):
                bounds = [float(np.nanmin(finite_values)), float(np.nanmax(finite_values))]

        if np.all(np.isfinite(bounds)) and bounds[0] < bounds[1]:
            if bins is None:
                bins = int(np.clip(np.sqrt(finite_values.size), 10, 80))
            else:
                bins = max(1, int(bins))
            counts, edges = np.histogram(
                finite_values,
                bins=max(1, bins),
                range=bounds,
                weights=finite_weights,
            )
            positive_counts = counts[counts > 0]
            if positive_counts.size > 0:
                peak = float(np.nanmax(positive_counts))
                typical = float(np.nanmedian(positive_counts))
                total = float(np.nansum(positive_counts))
                informative_peak = (
                    np.isfinite(peak)
                    and np.isfinite(typical)
                    and np.isfinite(total)
                    and total > 0
                    and (
                        force_histogram_mode
                        or (
                            typical > 0
                            and peak >= min_informative_peak_ratio * typical
                            and peak >= 0.05 * total
                        )
                    )
                )
                if informative_peak:
                    mode_index = int(np.argmax(counts))
                    estimate = float(0.5 * (edges[mode_index] + edges[mode_index + 1]))

        try:
            fallback_center = float(fallback_center)
        except (TypeError, ValueError):
            fallback_center = np.nan
        if not np.isfinite(estimate):
            estimate = fallback_center

        spread_candidates = [
            abs(float(q84) - float(estimate)) if np.isfinite(q84) and np.isfinite(estimate) else np.nan,
            abs(float(estimate) - float(q16)) if np.isfinite(q16) and np.isfinite(estimate) else np.nan,
            0.5 * abs(float(q84) - float(q16)) if np.isfinite(q16) and np.isfinite(q84) else np.nan,
        ]
        try:
            fallback_error = float(fallback_error)
        except (TypeError, ValueError):
            fallback_error = np.nan

        finite_spreads = [value for value in spread_candidates if np.isfinite(value) and value >= 0]
        if finite_spreads and max(finite_spreads) > 0:
            error = float(max(finite_spreads))
        elif np.isfinite(fallback_error) and fallback_error > 0:
            error = fallback_error
        else:
            error = np.nan

        return float(estimate), error

    def _visible_triangle_plot_values(self, values, plot_range, weights=None):
        values = np.asarray(values, dtype=float)
        try:
            lower, upper = [float(value) for value in np.asarray(plot_range, dtype=float).reshape(-1)[:2]]
        except (TypeError, ValueError, IndexError):
            finite_mask = np.isfinite(values)
            return values[finite_mask], None

        finite_mask = np.isfinite(values)
        if np.isfinite(lower) and np.isfinite(upper) and lower < upper:
            finite_mask &= (values >= lower) & (values <= upper)

        visible_weights = None
        if weights is not None:
            weights = np.asarray(weights, dtype=float)
            if weights.shape == values.shape:
                visible_weights = weights[finite_mask]
        return values[finite_mask], visible_weights

    def get_parameter_posterior_recenter_diagnostics(self, key, sigma_scale=5.0, bins=None):
        diagnostics = {
            'key': key,
            'clipped': False,
            'edge': None,
            'mode': np.nan,
            'std': np.nan,
            'full_std': np.nan,
            'bounds': None,
            'original_bounds': None,
            'sample_size': 0,
            'peak_height': np.nan,
            'lower_edge_height': np.nan,
            'upper_edge_height': np.nan,
            'lower_edge_peak_fraction': np.nan,
            'upper_edge_peak_fraction': np.nan,
            'reason': None,
        }

        bounds = getattr(self, 'sample_bounds', {}).get(key, self.bounds.get(key))
        if bounds is None:
            diagnostics['reason'] = "parameter bounds are unavailable."
            return diagnostics

        try:
            lower_bound, upper_bound = np.asarray(bounds, dtype=float).reshape(-1)[:2]
        except (TypeError, ValueError, IndexError):
            diagnostics['reason'] = "parameter bounds are malformed."
            return diagnostics

        diagnostics['original_bounds'] = [float(lower_bound), float(upper_bound)]
        if not np.isfinite(lower_bound) or not np.isfinite(upper_bound) or lower_bound >= upper_bound:
            diagnostics['reason'] = "parameter bounds are not finite."
            return diagnostics

        samples = self.get_parameter_posterior_samples(key)
        finite_samples = np.asarray(samples, dtype=float)
        finite_samples = finite_samples[np.isfinite(finite_samples)]
        diagnostics['sample_size'] = int(finite_samples.size)
        if finite_samples.size < 8:
            diagnostics['reason'] = "too few posterior samples are available."
            return diagnostics

        mode, bin_width = self._estimate_histogram_mode(finite_samples, bounds=(lower_bound, upper_bound), bins=bins)
        full_std = float(np.nanstd(finite_samples))
        diagnostics['mode'] = mode
        diagnostics['full_std'] = full_std

        if not np.isfinite(mode):
            diagnostics['reason'] = "posterior mode could not be estimated."
            return diagnostics

        q05, q16, q50, q84, q95 = np.nanpercentile(finite_samples, [5, 16, 50, 84, 95])
        width = float(upper_bound - lower_bound)
        histogram_bins = int(np.clip(np.sqrt(finite_samples.size), 10, 80)) if bins is None else max(1, int(bins))
        histogram_counts, _ = np.histogram(
            finite_samples,
            bins=histogram_bins,
            range=(lower_bound, upper_bound),
        )
        histogram_counts = np.asarray(histogram_counts, dtype=float)
        peak_height = float(np.nanmax(histogram_counts)) if histogram_counts.size else np.nan
        lower_edge_height = float(histogram_counts[0]) if histogram_counts.size else np.nan
        upper_edge_height = float(histogram_counts[-1]) if histogram_counts.size else np.nan
        if np.isfinite(peak_height) and peak_height > 0:
            lower_edge_peak_fraction = float(lower_edge_height / peak_height)
            upper_edge_peak_fraction = float(upper_edge_height / peak_height)
        else:
            lower_edge_peak_fraction = np.nan
            upper_edge_peak_fraction = np.nan

        diagnostics['peak_height'] = peak_height
        diagnostics['lower_edge_height'] = lower_edge_height
        diagnostics['upper_edge_height'] = upper_edge_height
        diagnostics['lower_edge_peak_fraction'] = lower_edge_peak_fraction
        diagnostics['upper_edge_peak_fraction'] = upper_edge_peak_fraction

        scale_floor = max(
            2.0 * bin_width if np.isfinite(bin_width) and bin_width > 0 else 0.0,
            0.01 * width,
            np.finfo(float).eps,
        )
        tail_gap_threshold = max(0.5 * full_std if np.isfinite(full_std) and full_std > 0 else 0.0, scale_floor)
        mode_gap_threshold = max(
            1.0 * full_std if np.isfinite(full_std) and full_std > 0 else 0.0,
            3.0 * bin_width if np.isfinite(bin_width) and bin_width > 0 else 0.0,
            0.05 * width,
            np.finfo(float).eps,
        )

        upper_gap_q95 = float(upper_bound - q95)
        lower_gap_q05 = float(q05 - lower_bound)
        upper_gap_mode = float(upper_bound - mode)
        lower_gap_mode = float(mode - lower_bound)

        upper_clipped = upper_gap_q95 <= tail_gap_threshold and upper_gap_mode <= mode_gap_threshold
        lower_clipped = lower_gap_q05 <= tail_gap_threshold and lower_gap_mode <= mode_gap_threshold
        edge_peak_fraction_floor = 0.20
        rejected_edges = []

        if upper_clipped and np.isfinite(upper_edge_peak_fraction) and upper_edge_peak_fraction < edge_peak_fraction_floor:
            upper_clipped = False
            rejected_edges.append(
                f"upper edge histogram height is only {upper_edge_peak_fraction:.3f} of the posterior peak"
            )
        if lower_clipped and np.isfinite(lower_edge_peak_fraction) and lower_edge_peak_fraction < edge_peak_fraction_floor:
            lower_clipped = False
            rejected_edges.append(
                f"lower edge histogram height is only {lower_edge_peak_fraction:.3f} of the posterior peak"
            )

        if upper_clipped and lower_clipped:
            clipped_edge = 'upper' if upper_gap_mode <= lower_gap_mode else 'lower'
        elif upper_clipped:
            clipped_edge = 'upper'
        elif lower_clipped:
            clipped_edge = 'lower'
        else:
            diagnostics['bounds'] = [float(lower_bound), float(upper_bound)]
            if rejected_edges:
                diagnostics['reason'] = (
                    "posterior reaches a search bound, but "
                    + " and ".join(rejected_edges)
                    + ", so it is not treated as truncated."
                )
            else:
                diagnostics['reason'] = "posterior support is comfortably inside the sampled bounds."
            return diagnostics

        diagnostics['clipped'] = True
        diagnostics['edge'] = clipped_edge

        if clipped_edge == 'upper':
            side_distances = mode - finite_samples[finite_samples <= mode]
        else:
            side_distances = finite_samples[finite_samples >= mode] - mode

        side_distances = np.asarray(side_distances, dtype=float)
        side_distances = side_distances[np.isfinite(side_distances)]
        side_distances = side_distances[side_distances >= 0]

        if side_distances.size >= 2:
            mirrored = np.concatenate([side_distances, -side_distances])
            estimated_std = float(np.nanstd(mirrored))
        else:
            estimated_std = full_std

        min_std = max(
            bin_width if np.isfinite(bin_width) and bin_width > 0 else 0.0,
            width * 1e-3,
            np.finfo(float).eps,
        )
        if not np.isfinite(estimated_std) or estimated_std <= 0:
            estimated_std = full_std
        if not np.isfinite(estimated_std) or estimated_std <= 0:
            estimated_std = min_std
        estimated_std = float(max(estimated_std, min_std))
        diagnostics['std'] = estimated_std

        radius = float(max(sigma_scale * estimated_std, min_std))
        new_lower = float(mode - radius)
        new_upper = float(mode + radius)
        if lower_bound >= 0:
            new_lower = max(0.0, new_lower)
        diagnostics['bounds'] = [new_lower, new_upper]
        diagnostics['reason'] = (
            f"posterior peaks against the {clipped_edge} search bound "
            f"(mode={mode:.6g}, sigma={estimated_std:.6g})."
        )
        diagnostics['q16'] = float(q16)
        diagnostics['q50'] = float(q50)
        diagnostics['q84'] = float(q84)
        diagnostics['q05'] = float(q05)
        diagnostics['q95'] = float(q95)
        return diagnostics

    def _get_triangle_plot_display_spec(
        self,
        sampled_keys,
        sample_parameters,
        sample_errors,
        sample_points,
        sample_weights=None,
    ):
        if 'b' in sampled_keys:
            key = 'b'
            label = r'Impact parameter $b$'
            mirror = False
        elif 'inc' in sampled_keys:
            key = 'inc'
            label = r'$\Delta i$'
            mirror = True
        else:
            return None

        geometry_index = sampled_keys.index(key)
        center = float(sample_parameters.get(key, self.parameters.get(key, 0.0)))
        sample_values = np.asarray(sample_points[:, geometry_index], dtype=float)
        magnitude_samples = np.abs(sample_values - center)
        error = float(sample_errors.get(key, np.nanstd(magnitude_samples)))
        if not np.isfinite(error) or error <= 0:
            error = float(np.nanstd(magnitude_samples))

        if mirror:
            display_range = self._get_mirrored_geometry_full_range(
                key,
                sample_values,
                center,
                sample_weights=sample_weights,
            )
        else:
            display_range = None

        if display_range is None:
            plot_lower, plot_upper = self._get_plot_range(key)
            plot_lower, plot_upper = self._expand_plot_range_for_sample_cloud(
                key,
                [plot_lower, plot_upper],
                sample_values,
                center,
            )
            plot_bins = int(max(1, np.sqrt(sample_points.shape[0])))
            plot_lower, plot_upper = self._expand_plot_range_for_histogram_edge_dropoff(
                key,
                [plot_lower, plot_upper],
                sample_values,
                center,
                bins=plot_bins,
                weights=sample_weights,
            )
            if mirror:
                max_distance = float(np.nanmax(np.abs([plot_lower - center, plot_upper - center])))
                if not np.isfinite(max_distance) or max_distance <= 0:
                    max_distance = float(np.nanmax(magnitude_samples))
                if not np.isfinite(max_distance) or max_distance <= 0:
                    max_distance = max(abs(center) * 1e-6, 1e-6)
                display_range = [-max_distance, max_distance]
            else:
                display_range = [float(plot_lower), float(plot_upper)]
        return {
            'key': key,
            'index': geometry_index,
            'label': label,
            'mirror': mirror,
            'center': center,
            'mask_center': 0.0 if mirror else center,
            'mask_error': error,
            'magnitude_samples': magnitude_samples,
            'range': display_range,
            'truth': 0.0 if mirror else center,
            'reference_lines': self._get_triangle_plot_geometry_reference_lines(
                key,
                center,
                sample_parameters,
            ),
        }

    def _get_triangle_plot_geometry_reference_lines(self, key, center, sample_parameters):
        if key != 'b':
            return []

        try:
            center = float(center)
        except (TypeError, ValueError):
            center = np.nan
        if not np.isfinite(center):
            return []

        rprs = sample_parameters.get('rprs')
        if rprs is None:
            rprs = getattr(self, 'parameters', {}).get(
                'rprs',
                getattr(self, 'prior', {}).get('rprs', np.nan),
            )
        try:
            rprs = float(rprs)
        except (TypeError, ValueError):
            rprs = np.nan

        reference_lines = [
            {
                'value': 1.0,
                'color': '#707070',
                'linestyle': ':',
                'linewidth': 0.9,
                'alpha': 0.9,
            },
        ]
        if np.isfinite(rprs) and rprs >= 0:
            reference_lines.append(
                {
                    'value': 1.0 + rprs,
                    'color': '#a35d00',
                    'linestyle': '-.',
                    'linewidth': 0.9,
                    'alpha': 0.9,
                }
            )
        return reference_lines

    def _get_triangle_plot_geometry_reference_offsets(self, display_spec):
        reference_lines = display_spec.get('reference_lines', []) if isinstance(display_spec, dict) else []
        if not reference_lines:
            return []

        try:
            center = float(display_spec['center'])
        except (KeyError, TypeError, ValueError):
            return []
        if not np.isfinite(center):
            return []

        offsets = []
        for reference in reference_lines:
            try:
                value = float(reference['value'])
            except (KeyError, TypeError, ValueError):
                continue
            if not np.isfinite(value):
                continue

            distance = abs(value - center)
            if not np.isfinite(distance):
                continue
            reference_offsets = [0.0] if distance <= np.finfo(float).eps else [-distance, distance]
            for offset in reference_offsets:
                offsets.append({
                    'offset': float(offset),
                    'color': reference.get('color', '#707070'),
                    'linestyle': reference.get('linestyle', ':'),
                    'linewidth': reference.get('linewidth', 0.9),
                    'alpha': reference.get('alpha', 0.9),
                })
        return offsets

    def _draw_triangle_plot_geometry_reference_lines(self, ax, display_spec, axis='x', limits=None):
        if ax is None:
            return

        if limits is None:
            limits = ax.get_xlim() if axis == 'x' else ax.get_ylim()
        lower, upper = np.sort(np.asarray(limits, dtype=float).reshape(-1)[:2])
        if not display_spec.get('mirror', True):
            for reference in display_spec.get('reference_lines', []):
                try:
                    value = float(reference['value'])
                except (KeyError, TypeError, ValueError):
                    continue
                if not np.isfinite(value) or value < lower or value > upper:
                    continue
                line_kwargs = {
                    'color': reference.get('color', '#707070'),
                    'linestyle': reference.get('linestyle', ':'),
                    'linewidth': reference.get('linewidth', 0.9),
                    'alpha': reference.get('alpha', 0.9),
                    'zorder': 2,
                }
                if axis == 'y':
                    ax.axhline(value, **line_kwargs)
                else:
                    ax.axvline(value, **line_kwargs)
            return

        if lower <= 0.0 <= upper:
            center_kwargs = {
                'color': '#4682b4',
                'linestyle': '--',
                'linewidth': 0.9,
                'alpha': 0.85,
                'zorder': 2,
            }
            if axis == 'y':
                ax.axhline(0.0, **center_kwargs)
            else:
                ax.axvline(0.0, **center_kwargs)
        for reference in self._get_triangle_plot_geometry_reference_offsets(display_spec):
            offset = reference['offset']
            if offset < lower or offset > upper:
                continue
            line_kwargs = {
                'color': reference['color'],
                'linestyle': reference['linestyle'],
                'linewidth': reference['linewidth'],
                'alpha': reference['alpha'],
                'zorder': 2,
            }
            if axis == 'y':
                ax.axhline(offset, **line_kwargs)
            else:
                ax.axvline(offset, **line_kwargs)

    def _get_triangle_plot_geometry_overlay(self, display_spec, sample_points):
        if display_spec is None:
            return None
        if not display_spec.get('mirror', True):
            return None

        geometry_index = display_spec['index']
        center = display_spec['center']
        offsets = np.asarray(sample_points[:, geometry_index], dtype=float) - center
        left_offsets = offsets[offsets <= 0]
        right_offsets = offsets[offsets >= 0]

        def mirrored_offsets(branch_offsets):
            if branch_offsets.size == 0:
                return np.array([], dtype=float)
            return np.concatenate([branch_offsets, -branch_offsets])

        return {
            'index': geometry_index,
            'left_count': left_offsets.size,
            'right_count': right_offsets.size,
            'left_mirrored': mirrored_offsets(left_offsets),
            'right_mirrored': mirrored_offsets(right_offsets),
        }

    def _format_triangle_plot_geometry_value(self, value, error, suffix=''):
        if value is None or not np.isfinite(value):
            return f"n/a{suffix}"
        if error is None or not np.isfinite(error) or error < 0:
            return f"{round_to_2(value)}{suffix}"
        return f"{round_to_2(value, error)} +/- {round_to_2(error)}{suffix}"

    def _get_triangle_plot_geometry_summary(self, sampled_keys, sample_points):
        bound_keys = list(self.bounds.keys())
        sample_points = np.asarray(sample_points, dtype=float)
        sample_parameters = getattr(self, 'sample_parameters', {})
        sample_errors = getattr(self, 'sample_errors', {})

        physical_samples = [
            self._physical_values_from_sample_point(point, bound_keys, sampled_keys)
            for point in sample_points
        ]
        inc_samples = np.array([sample['inc'] for sample in physical_samples], dtype=float)
        b_samples = np.array([
            sample.get('b', impact_parameter_from_inclination(sample, sample['inc']))
            for sample in physical_samples
        ], dtype=float)

        inc_center = float(self.parameters.get('inc', np.nanmedian(inc_samples)))
        inc_error = float(self.errors.get('inc', np.nanstd(inc_samples)))
        if 'b' in sample_parameters:
            b_center = float(sample_parameters['b'])
        elif 'b' in self.parameters:
            b_center = float(self.parameters['b'])
        else:
            b_center = float(np.nanmedian(b_samples))
        b_error = float(sample_errors.get('b', self.errors.get('b', np.nanstd(b_samples))))

        title = (
            f"b={self._format_triangle_plot_geometry_value(b_center, b_error)}\n"
            f"i={self._format_triangle_plot_geometry_value(inc_center, inc_error, ' deg')}"
        )
        return {
            'b_center': b_center,
            'b_error': b_error,
            'inc_center': inc_center,
            'inc_error': inc_error,
            'title': title,
        }

    def _smooth_triangle_plot_counts(self, counts):
        counts = np.asarray(counts, dtype=float)
        if counts.size <= 1 or not np.any(counts > 0):
            return counts

        sigma_bins = max(1.0, counts.size / 18.0)
        radius = max(1, int(np.ceil(3 * sigma_bins)))
        grid = np.arange(-radius, radius + 1, dtype=float)
        kernel = np.exp(-0.5 * (grid / sigma_bins) ** 2)
        kernel /= np.sum(kernel)
        return np.convolve(counts, kernel, mode='same')

    def _build_triangle_plot_geometry_curves(self, geometry_overlay, hist_range, bins_1d):
        hist_range = np.sort(np.asarray(hist_range, dtype=float))
        bins_1d = max(1, int(bins_1d))
        edges = np.linspace(hist_range[0], hist_range[1], bins_1d + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        half_bin = 0.5 * (edges[1] - edges[0]) if edges.size > 1 else 0.0

        def branch_curve(samples):
            samples = np.asarray(samples, dtype=float)
            counts, _ = np.histogram(samples, bins=edges)
            support = np.zeros_like(counts, dtype=bool)
            if samples.size > 0:
                support = np.abs(centers) <= (np.max(np.abs(samples)) + half_bin)
            return counts.astype(float), support

        left_curve, left_support = branch_curve(geometry_overlay['left_mirrored'])
        right_curve, right_support = branch_curve(geometry_overlay['right_mirrored'])

        left_count = int(geometry_overlay.get('left_count', 0))
        right_count = int(geometry_overlay.get('right_count', 0))
        max_count = max(left_count, right_count)
        min_count = min(left_count, right_count)

        if max_count == 0:
            main_curve = np.zeros_like(centers, dtype=float)
        elif min_count == 0 or (min_count / max_count) < 0.35:
            main_curve = left_curve if left_count >= right_count else right_curve
        else:
            stacked = np.vstack([
                np.where(left_support, left_curve, np.nan),
                np.where(right_support, right_curve, np.nan),
            ])
            valid_counts = np.sum(np.isfinite(stacked), axis=0)
            summed = np.nansum(stacked, axis=0)
            main_curve = np.divide(
                summed,
                valid_counts,
                out=np.zeros_like(summed, dtype=float),
                where=valid_counts > 0,
            )

        main_curve = self._smooth_triangle_plot_counts(main_curve)

        return {
            'centers': centers,
            'left_curve': left_curve,
            'right_curve': right_curve,
            'main_curve': main_curve,
        }

    def _get_triangle_plot_payload(self):
        sampled_keys = getattr(self, 'sampled_keys', list(self.bounds.keys()))
        sample_parameters = getattr(self, 'sample_parameters', self.parameters)
        sample_errors = getattr(self, 'sample_errors', self.errors)
        sample_points, sample_logl, sample_weights = self._get_triangle_plot_samples()
        display_spec = self._get_triangle_plot_display_spec(
            sampled_keys,
            sample_parameters,
            sample_errors,
            sample_points,
            sample_weights=sample_weights,
        )
        geometry_overlay = self._get_triangle_plot_geometry_overlay(display_spec, sample_points)
        geometry_summary = self._get_triangle_plot_geometry_summary(sampled_keys, sample_points)

        display_points = np.array(sample_points, copy=True)
        display_logl = np.array(sample_logl, copy=True)
        display_weights = None if sample_weights is None else np.array(sample_weights, copy=True)
        mask_values = np.array(sample_points, copy=True)

        if display_spec is not None and display_spec.get('mirror', True):
            geometry_index = display_spec['index']
            positive_points = np.array(sample_points, copy=True)
            negative_points = np.array(sample_points, copy=True)
            positive_points[:, geometry_index] = display_spec['magnitude_samples']
            negative_points[:, geometry_index] = -display_spec['magnitude_samples']
            display_points = np.vstack([positive_points, negative_points])
            display_logl = np.concatenate([sample_logl, sample_logl])
            if sample_weights is not None:
                display_weights = np.concatenate([sample_weights, sample_weights])
            mask_values = np.array(display_points, copy=True)

        plot_bins = int(max(1, np.sqrt(display_points.shape[0])))

        flabels = {
            'rprs': r'R$_{p}$/R$_{s}$',
            'per': r'Period [day]',
            'tmid': r'T$_{mid}$',
            'ars': r'a/R$_{s}$',
            'inc': r'Inc. [deg]',
            'b': r'Impact parameter',
            'u1': r'u$_1$',
            'fpfs': r'F$_{p}$/F$_{s}$',
            'omega': r'$\omega$ [deg]',
            'mplanet': r'M$_{p}$ [M$_{\oplus}$]',
            'mstar': r'M$_{s}$ [M$_{\odot}$]',
            'ecc': r'$e$',
            'c0': r'$c_0$',
            'c1': r'$c_1$',
            'c2': r'$c_2$',
            'c3': r'$c_3$',
            'c4': r'$c_4$',
            'a0': r'$a_0$',
            'a1': r'$a_1$',
            'a2': r'$a_2$'
        }

        labels = []
        titles = []
        ranges = []
        mask_centers = []
        mask_errors = []
        truths = []

        for i, key in enumerate(sampled_keys):
            center = sample_parameters.get(key, self.parameters.get(key, 0.0))
            error = sample_errors.get(key, self.errors.get(key, 0.0))
            label = flabels.get(key, key)
            plot_range = self._get_plot_range(key)
            if sample_points.ndim == 2 and i < sample_points.shape[1]:
                plot_range = self._expand_plot_range_for_sample_cloud(
                    key,
                    plot_range,
                    sample_points[:, i],
                    center,
                )
                plot_range = self._expand_plot_range_for_histogram_edge_dropoff(
                    key,
                    plot_range,
                    sample_points[:, i],
                    center,
                    bins=plot_bins,
                    weights=sample_weights,
                )
                center, error = self._triangle_plot_display_estimate(
                    sample_points[:, i],
                    center,
                    error,
                    plot_range=plot_range,
                    weights=sample_weights,
                )
            title = self._format_triangle_plot_parameter_title(center, error)
            truth = center

            if display_spec is not None and key == display_spec['key']:
                label = display_spec['label']
                title = geometry_summary['title']
                plot_range = display_spec['range']
                center = display_spec['mask_center']
                error = display_spec['mask_error']
                truth = display_spec['truth']

            labels.append(label)
            titles.append(title)
            ranges.append(plot_range)
            mask_centers.append(center)
            mask_errors.append(error)
            try:
                truth = float(truth)
            except (TypeError, ValueError):
                truth = np.nan
            truths.append(truth if np.isfinite(truth) else None)

        return {
            'sampled_keys': sampled_keys,
            'display_points': display_points,
            'display_logl': display_logl,
            'display_weights': display_weights,
            'mask_values': mask_values,
            'display_spec': display_spec,
            'geometry_overlay': geometry_overlay,
            'geometry_summary': geometry_summary,
            'labels': labels,
            'titles': titles,
            'ranges': ranges,
            'mask_centers': mask_centers,
            'mask_errors': mask_errors,
            'truths': truths,
        }

    def _triangle_plot_sigma_window_ranges(self, payload, sigma):
        try:
            sigma = float(sigma)
        except (TypeError, ValueError):
            return payload['ranges']
        if not np.isfinite(sigma) or sigma <= 0:
            return payload['ranges']

        zoomed_ranges = []
        display_points = np.asarray(payload.get('display_points', []), dtype=float)
        for i, plot_range in enumerate(payload['ranges']):
            try:
                range_lower, range_upper = [
                    float(value) for value in np.asarray(plot_range, dtype=float).reshape(-1)[:2]
                ]
            except (TypeError, ValueError, IndexError):
                zoomed_ranges.append(plot_range)
                continue

            if not np.isfinite(range_lower) or not np.isfinite(range_upper) or range_lower >= range_upper:
                zoomed_ranges.append(plot_range)
                continue

            try:
                center = float(payload['mask_centers'][i])
                error = float(payload['mask_errors'][i])
            except (TypeError, ValueError, IndexError):
                zoomed_ranges.append(plot_range)
                continue

            if not np.isfinite(center) or not np.isfinite(error) or error <= 0:
                zoomed_ranges.append(plot_range)
                continue

            lower = max(range_lower, center - sigma * error)
            upper = min(range_upper, center + sigma * error)
            if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
                zoomed_ranges.append(plot_range)
                continue

            if display_points.ndim == 2 and i < display_points.shape[1]:
                values = display_points[:, i]
                finite_values = values[np.isfinite(values)]
                if finite_values.size and not np.any((finite_values >= lower) & (finite_values <= upper)):
                    zoomed_ranges.append(plot_range)
                    continue

            zoomed_ranges.append([float(lower), float(upper)])

        return zoomed_ranges

    def _recenter_triangle_plot_payload_for_visible_ranges(self, payload):
        display_points = np.asarray(payload.get('display_points', []), dtype=float)
        if display_points.ndim != 2 or display_points.shape[1] == 0:
            return payload

        updated = dict(payload)
        titles = list(payload.get('titles', []))
        truths = list(payload.get('truths', []))
        mask_centers = list(payload.get('mask_centers', []))
        mask_errors = list(payload.get('mask_errors', []))
        ranges = list(payload.get('ranges', []))
        sampled_keys = list(payload.get('sampled_keys', []))

        display_weights = payload.get('display_weights')
        if display_weights is not None:
            display_weights = np.asarray(display_weights, dtype=float)
            if display_weights.ndim != 1 or display_weights.shape[0] != display_points.shape[0]:
                display_weights = None

        plot_bins = int(max(1, np.sqrt(display_points.shape[0])))
        display_spec = payload.get('display_spec')
        geometry_summary = payload.get('geometry_summary') or {}

        for i, key in enumerate(sampled_keys):
            if i >= display_points.shape[1] or i >= len(ranges):
                continue
            if (
                display_spec is not None
                and key == display_spec.get('key')
                and display_spec.get('mirror', False)
            ):
                continue

            visible_values, visible_weights = self._visible_triangle_plot_values(
                display_points[:, i],
                ranges[i],
                weights=display_weights,
            )
            if visible_values.size < 2:
                continue

            fallback_center = truths[i] if i < len(truths) else np.nan
            if fallback_center is None or not np.isfinite(fallback_center):
                fallback_center = mask_centers[i] if i < len(mask_centers) else np.nan
            fallback_error = mask_errors[i] if i < len(mask_errors) else np.nan
            center, error = self._triangle_plot_display_estimate(
                visible_values,
                fallback_center,
                fallback_error,
                plot_range=ranges[i],
                weights=visible_weights,
                bins=plot_bins,
                force_histogram_mode=True,
            )
            if not np.isfinite(center):
                continue

            if i < len(truths):
                truths[i] = center
            if i < len(mask_centers):
                mask_centers[i] = center
            if i < len(mask_errors):
                mask_errors[i] = error
            if i < len(titles):
                if display_spec is not None and key == display_spec.get('key'):
                    inc_center = geometry_summary.get('inc_center')
                    inc_error = geometry_summary.get('inc_error')
                    titles[i] = (
                        f"b={self._format_triangle_plot_geometry_value(center, error)}\n"
                        f"i={self._format_triangle_plot_geometry_value(inc_center, inc_error, ' deg')}"
                    )
                else:
                    titles[i] = self._format_triangle_plot_parameter_title(center, error)

        updated['titles'] = titles
        updated['truths'] = truths
        updated['mask_centers'] = mask_centers
        updated['mask_errors'] = mask_errors
        return updated

    def _triangle_contour_levels(self, chi2, mask1, mask2, mask3):
        raw_levels = np.array([
            np.percentile(chi2[mask1], 95),
            np.percentile(chi2[mask2], 95),
            np.percentile(chi2[mask3], 95),
        ], dtype=float)
        finite_levels = np.sort(raw_levels[np.isfinite(raw_levels)])
        if finite_levels.size == 0:
            return []

        unique_levels = []
        min_spacing = max(
            np.finfo(float).eps,
            np.nanmax(np.abs(finite_levels)) * 1e-12,
        )
        for level in finite_levels:
            if not unique_levels or level > unique_levels[-1] + min_spacing:
                unique_levels.append(float(level))
        return unique_levels

    def _overlay_triangle_plot_geometry_histograms(self, fig, payload, title_kwargs=None, label_kwargs=None):
        if not hasattr(fig, 'axes'):
            return

        display_spec = payload.get('display_spec')
        if display_spec is None:
            return

        sampled_keys = payload['sampled_keys']
        if len(fig.axes) != len(sampled_keys) ** 2:
            return

        axes = np.array(fig.axes).reshape((len(sampled_keys), len(sampled_keys)))
        if not display_spec.get('mirror', True):
            geometry_index = display_spec['index']
            for row in range(len(sampled_keys)):
                for col in range(len(sampled_keys)):
                    panel = axes[row, col]
                    if row == geometry_index and col == geometry_index:
                        self._draw_triangle_plot_geometry_reference_lines(panel, display_spec, axis='x')
                    elif col == geometry_index and row > col:
                        self._draw_triangle_plot_geometry_reference_lines(panel, display_spec, axis='x')
                    elif row == geometry_index and col < row:
                        self._draw_triangle_plot_geometry_reference_lines(panel, display_spec, axis='y')
            return

        geometry_overlay = payload.get('geometry_overlay')
        if geometry_overlay is None:
            return

        geometry_index = geometry_overlay['index']
        ax = axes[geometry_index, geometry_index]
        hist_range = np.sort(payload['ranges'][geometry_index])
        bins_1d = int(max(25, np.round(np.sqrt(payload['display_points'].shape[0]) * 6)))
        curves = self._build_triangle_plot_geometry_curves(geometry_overlay, hist_range, bins_1d)

        title = payload['titles'][geometry_index]
        x_label = payload['labels'][geometry_index]
        branch_left_color = '#6f8fcf'
        branch_right_color = '#d79b9b'
        title_kwargs = {} if title_kwargs is None else dict(title_kwargs)
        label_kwargs = {} if label_kwargs is None else dict(label_kwargs)

        ax.cla()
        ax.plot(curves['centers'], curves['main_curve'], color='black', linewidth=1.5, zorder=4)
        ax.plot(curves['centers'], curves['left_curve'], color=branch_left_color, linestyle='--',
                linewidth=0.75, alpha=0.75, zorder=3)
        ax.plot(curves['centers'], curves['right_curve'], color=branch_right_color, linestyle='--',
                linewidth=0.75, alpha=0.75, zorder=3)
        self._draw_triangle_plot_geometry_reference_lines(
            ax,
            display_spec,
            axis='x',
            limits=hist_range,
        )
        ax.set_title(title, **title_kwargs)
        if 'fontsize' in title_kwargs:
            ax.title.set_fontsize(title_kwargs['fontsize'])
        ax.set_xlim(hist_range)

        max_y = max(
            np.max(curves['main_curve']) if curves['main_curve'].size > 0 else 0.0,
            np.max(curves['left_curve']) if curves['left_curve'].size > 0 else 0.0,
            np.max(curves['right_curve']) if curves['right_curve'].size > 0 else 0.0,
        )
        ax.set_ylim(0, 1.1 * max(max_y, 1e-6))
        ax.set_yticks([])

        if geometry_index < len(sampled_keys) - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel(x_label, **label_kwargs)

        for row in range(len(sampled_keys)):
            for col in range(len(sampled_keys)):
                if row == geometry_index and col == geometry_index:
                    continue
                panel = axes[row, col]
                if col == geometry_index and row > col:
                    self._draw_triangle_plot_geometry_reference_lines(panel, display_spec, axis='x')
                if row == geometry_index and col < row:
                    self._draw_triangle_plot_geometry_reference_lines(panel, display_spec, axis='y')

    def _adjust_triangle_plot_layout(self, fig):
        if not hasattr(fig, 'subplots_adjust'):
            return
        subplotpars = getattr(fig, 'subplotpars', None)
        if subplotpars is None:
            return

        fig.subplots_adjust(
            left=max(subplotpars.left, 0.08),
            bottom=max(subplotpars.bottom, 0.12),
            right=min(subplotpars.right, 0.97),
            top=min(subplotpars.top, 0.94),
            wspace=subplotpars.wspace,
            hspace=subplotpars.hspace,
        )

    def fit_LM(self):
        freekeys = list(self.bounds.keys())
        boundarray = np.array([self.bounds[k] for k in freekeys])
        self._validate_flux_baseline_keys()

        # trim data around predicted transit/eclipse time
        if np.ndim(self.airmass) == 2:
            print(f'Computing nearest neighbors and gaussian weights for {len(self.time)} npts...')
            self.gw, self.nearest = gaussian_weights(self.airmass, neighbors=self.neighbors)

        def lc2min_nneighbor(pars):
            for i in range(len(pars)):
                self.prior[freekeys[i]] = pars[i]
            lightcurve = transit(self.time, self.prior)
            detrended = self.data / lightcurve
            wf = weightedflux(detrended, self.gw, self.nearest)
            model = lightcurve * wf
            return ((self.data - model) / self.dataerr) ** 2

        def lc2min_airmass(pars):
            for i in range(len(pars)):
                self.prior[freekeys[i]] = pars[i]
            model = transit(self.time, self.prior)
            model *= airmass_trend(
                self.prior.get('a2', 0),
                self.airmass,
                reference=self._get_airmass_reference(),
            )
            if self._has_free_flux_baseline():
                model *= get_flux_baseline(self.prior)
            elif self._uses_fixed_flux_baseline():
                model *= get_flux_baseline(self.prior)
            else:
                model *= solve_flux_baseline(
                    model,
                    self.data,
                    self.dataerr,
                    mask=self._get_baseline_fit_mask(),
                )
            return ((self.data - model) / self.dataerr) ** 2

        try:
            if np.ndim(self.airmass) == 2:
                res = least_squares(lc2min_nneighbor, x0=[self.prior[k] for k in freekeys],
                                    bounds=[boundarray[:, 0], boundarray[:, 1]], jac='3-point', loss='linear')
            else:
                res = least_squares(lc2min_airmass, x0=[self.prior[k] for k in freekeys],
                                    bounds=[boundarray[:, 0], boundarray[:, 1]], jac='3-point', loss='linear')
        except Exception as e:
            print(f"{e} \nbounded light curve fitting failed...check priors "
                  "(e.g. estimated mid-transit time + orbital period)")

            for i, k in enumerate(freekeys):
                if not boundarray[i, 0] < self.prior[k] < boundarray[i, 1]:
                    print(f"bound: [{boundarray[i, 0]}, {boundarray[i, 1]}] prior: {self.prior[k]}")

            print("removing bounds and trying again...")

            if np.ndim(self.airmass) == 2:
                res = least_squares(lc2min_nneighbor, x0=[self.prior[k] for k in freekeys],
                                    method='lm', jac='3-point', loss='linear')
            else:
                res = least_squares(lc2min_airmass, x0=[self.prior[k] for k in freekeys],
                                    method='lm', jac='3-point', loss='linear')

        self.parameters = copy.deepcopy(self.prior)
        self.errors = {}
        self.quantiles = {}

        for i, k in enumerate(freekeys):
            self.parameters[k] = res.x[i]
            self.errors[k] = 0
            self.quantiles[k] = [0, 0]

        self.sampled_keys = list(freekeys)
        self.sample_bounds = copy.deepcopy(self.bounds)
        self.sample_parameters = {k: self.parameters[k] for k in self.sampled_keys}
        self.sample_errors = {k: self.errors[k] for k in self.sampled_keys}
        self.sample_quantiles = {k: self.quantiles[k] for k in self.sampled_keys}

        self.create_fit_variables()

    def create_fit_variables(self):
        self.transit = transit(self.time, self.parameters)
        self._apply_fixed_parameter_errors()
        self._update_plot_geometry()
        if np.ndim(self.airmass) != 2:
            if self._has_free_flux_baseline():
                flux_scale = get_flux_baseline(self.parameters)
                flux_scale_err = self.errors.get('a0', self.errors.get('a1', 0.0))
            elif self._uses_fixed_flux_baseline():
                flux_scale = get_flux_baseline(self.parameters)
                flux_scale_err = self.errors.get(
                    'a0',
                    self.errors.get('a1', self.fixed_parameter_errors.get('a0', 0.0)),
                )
            elif self.mode == "ns":
                flux_scale, flux_scale_err = mc_a1(
                    self.parameters.get('a2', 0),
                    self.errors.get('a2', 1e-6),
                    self.transit,
                    self.airmass,
                    self.data,
                    self.dataerr,
                    mask=self._get_baseline_fit_mask(),
                )
            else:
                systematics = self.transit * airmass_trend(
                    self.parameters.get('a2', 0),
                    self.airmass,
                    reference=self._get_airmass_reference(),
                )
                flux_scale = solve_flux_baseline(
                    systematics,
                    self.data,
                    self.dataerr,
                    mask=self._get_baseline_fit_mask(),
                )
                flux_scale_err = self.errors.get(
                    'a0',
                    self.errors.get(
                        'a1',
                        solve_flux_baseline_uncertainty(
                            systematics,
                            self.dataerr,
                            mask=self._get_baseline_fit_mask(),
                        ),
                    ),
                )
            self._set_flux_baseline(flux_scale, flux_scale_err)
        if np.ndim(self.airmass) == 2:
            detrended = self.data / self.transit
            self.wf = weightedflux(detrended, self.gw, self.nearest)
            self.model = self.transit * self.wf
            self.detrended = self.data / self.wf
            self.detrendederr = self.dataerr / self.wf
        else:
            self.airmass_model = self._build_systematics_model(self.parameters)
            self.model = self.transit * self.airmass_model
            self.detrended = self.data / self.airmass_model
            self.detrendederr = self.dataerr / self.airmass_model

        self.residuals = self.data - self.model
        self.res_stdev = np.std(self.residuals)/np.median(self.data)
        self.chi2 = np.sum(self.residuals ** 2 / self.dataerr ** 2)
        self.bic = len(self.bounds) * np.log(len(self.time)) - 2 * np.log(self.chi2)

        # compare fit chi2 to smoothed data chi2
        dt = np.diff(np.sort(self.time)).mean()
        si = np.argsort(self.time)
        try:
            self.sdata = savgol_filter(self.data[si], 1 + 2 * int(0.5 / 24 / dt), 2)
        except:
            self.sdata = np.ones(len(self.time))

        schi2 = np.sum((self.data[si] - self.sdata) ** 2 / self.dataerr[si] ** 2)
        self.quality = schi2 / self.chi2

        # measured duration
        tdur = (self.transit < 1).sum() * np.median(np.diff(np.sort(self.time)))

        # test for partial transit
        newdur = transit_duration(self.parameters)
        if not np.isfinite(newdur) or newdur <= 0:
            newtime = np.linspace(self.parameters['tmid'] - 0.2, self.parameters['tmid'] + 0.2, 10000)
            newtran = transit(newtime, self.parameters)
            masktran = newtran < 1
            newdur = np.diff(newtime).mean() * masktran.sum()

        self.duration_measured = tdur
        self.duration_expected = newdur

    def _finalize_ultranest_fit_results(self, bound_keys, sampled_keys, physical_from_sample_point):
        self.sample_parameters = {}
        self.sample_errors = {}
        self.sample_quantiles = {}
        self.errors = {}
        self.quantiles = {}
        self.parameters = copy.deepcopy(self.prior)

        ml_point = self.results['maximum_likelihood']['point']
        self.sample_bounds = self._get_sample_bounds(bound_keys, physical_from_sample_point(ml_point))
        self.ultranest_error_fallbacks = {}
        weighted_points, weighted_logl = self._get_ultranest_weighted_sample_arrays()

        for i, key in enumerate(sampled_keys):
            self.sample_parameters[key] = ml_point[i]
            reported_error = self.results['posterior']['stdev'][i]
            reported_quantiles = [
                self.results['posterior']['errlo'][i],
                self.results['posterior']['errup'][i]]
            if self._ultranest_error_needs_sample_fallback(
                i,
                ml_point[i],
                reported_error,
                points=weighted_points,
            ):
                fallback = self._loglike_neighborhood_uncertainty(
                    i,
                    ml_point[i],
                    points=weighted_points,
                    logl=weighted_logl,
                )
                if fallback is not None:
                    fallback['reason'] = 'degenerate_posterior_summary'
            else:
                fallback = None
                local_uncertainty = self._loglike_neighborhood_uncertainty(
                    i,
                    ml_point[i],
                    points=weighted_points,
                    logl=weighted_logl,
                )
                if self._ultranest_error_is_inflated_relative_to_local_fit(reported_error, local_uncertainty):
                    fallback = local_uncertainty
                    fallback['reported_error'] = float(reported_error)
                    fallback['reason'] = 'posterior_summary_inflated_relative_to_local_fit'
            if fallback is not None:
                self.sample_errors[key] = fallback['error']
                self.sample_quantiles[key] = fallback['quantiles']
                self.ultranest_error_fallbacks[key] = fallback
            else:
                self.sample_errors[key] = reported_error
                self.sample_quantiles[key] = reported_quantiles

        physical_ml = physical_from_sample_point(ml_point)
        self.parameters.update(physical_ml)

        for bound_key, sampled_key in zip(bound_keys, sampled_keys):
            if bound_key == 'inc' and sampled_key == 'b':
                continue
            self.errors[bound_key] = self.sample_errors[sampled_key]
            self.quantiles[bound_key] = self.sample_quantiles[sampled_key]

        if 'inc' in bound_keys and 'b' in sampled_keys and weighted_points is not None:
            bound_index = {key: index for index, key in enumerate(bound_keys)}

            def weighted_sample_values(key, default=0.0):
                index = bound_index.get(key)
                if index is not None and index < weighted_points.shape[1] and sampled_keys[index] != 'b':
                    return weighted_points[:, index]
                value = physical_ml.get(key, self.prior.get(key, default))
                return np.full(weighted_points.shape[0], float(value), dtype=float)

            b_index = sampled_keys.index('b')
            scale_values = {
                'ars': weighted_sample_values('ars', np.nan),
                'ecc': weighted_sample_values('ecc', 0.0),
                'omega': weighted_sample_values('omega', 0.0),
            }
            inc_samples = np.asarray(
                inclination_from_impact_parameter(scale_values, weighted_points[:, b_index]),
                dtype=float,
            )
            center, std, quantiles = self._summarize_derived_parameter(inc_samples, physical_ml['inc'])
            self.parameters['inc'] = center
            self.errors['inc'] = std
            self.quantiles['inc'] = quantiles
        self._apply_fixed_parameter_errors()

    def extend_ultranest_fit(self, min_num_live_points=None, max_ncalls=None):
        context = getattr(self, '_ultranest_resume_context', None)
        if getattr(self, 'ns_type', None) != 'ultranest' or not isinstance(context, dict):
            return False

        sampler = context.get('sampler')
        if sampler is None:
            return False

        run_kwargs = {"max_ncalls": int(max_ncalls if max_ncalls is not None else self.max_ncalls)}
        if min_num_live_points is not None:
            run_kwargs["min_num_live_points"] = int(min_num_live_points)

        self.results = run_reactive_sampler(
            sampler,
            run_kwargs=run_kwargs,
            verbose=self.verbose,
        )
        self._finalize_ultranest_fit_results(
            context['bound_keys'],
            context['sampled_keys'],
            context['physical_from_sample_point'],
        )
        self.create_fit_variables()
        return True

    def clear_ultranest_resume_state(self):
        self._ultranest_resume_context = None

    def fit_nested(self):
        bound_keys = list(self.bounds.keys())
        sampled_keys = self._get_sampled_keys(bound_keys)
        self._validate_flux_baseline_keys()
        self.sampled_keys = list(sampled_keys)
        self.sample_bounds = self._get_sample_bounds(bound_keys, self.prior)
        self.impact_parameter_sampled_directly = self._uses_internal_impact_parameter()

        if len(set(self.sampled_keys)) != len(self.sampled_keys):
            raise ValueError("Free-parameter labels must be unique after internal parameter transforms.")

        # alloc data for best fit + error
        self.sample_parameters = {}
        self.sample_errors = {}
        self.sample_quantiles = {}
        self.errors = {}
        self.quantiles = {}
        self.parameters = copy.deepcopy(self.prior)

        base_physical = dict(self.prior)
        direct_sample_assignments = [
            (bound_key, index)
            for index, (bound_key, sampled_key) in enumerate(zip(bound_keys, sampled_keys))
            if not (sampled_key == 'b' and bound_key == 'inc')
        ]
        impact_parameter_index = next(
            (
                index
                for index, (bound_key, sampled_key) in enumerate(zip(bound_keys, sampled_keys))
                if sampled_key == 'b' and bound_key == 'inc'
            ),
            None,
        )
        time = self.time
        data = np.asarray(self.data, dtype=float)
        dataerr = np.asarray(self.dataerr, dtype=float)
        data_shape = data.shape
        dataerr_shape_matches = dataerr.shape == data_shape
        finite_dataerr = np.isfinite(dataerr) & (dataerr > 0) if dataerr_shape_matches else False
        observed_values_valid = (
            dataerr_shape_matches
            and np.all(np.isfinite(data))
            and np.all(finite_dataerr)
        )
        inverse_dataerr = np.zeros(data_shape, dtype=float)
        baseline_weights = np.zeros(data_shape, dtype=float)
        baseline_static_mask = np.isfinite(data)
        if dataerr_shape_matches:
            inverse_dataerr[finite_dataerr] = 1.0 / dataerr[finite_dataerr]
            baseline_weights[finite_dataerr] = inverse_dataerr[finite_dataerr] ** 2
            baseline_static_mask &= np.isfinite(baseline_weights) & (baseline_weights > 0)
        else:
            baseline_static_mask &= False

        baseline_fit_mask = self._get_baseline_fit_mask()
        if baseline_fit_mask is not None:
            baseline_static_mask &= baseline_fit_mask

        centered_airmass = center_airmass(self.airmass, reference=self._get_airmass_reference())
        has_free_flux_baseline = self._has_free_flux_baseline()
        uses_fixed_flux_baseline = self._uses_fixed_flux_baseline()
        sampled_key_index = {key: index for index, key in enumerate(sampled_keys)}
        sampled_a2_index = sampled_key_index.get('a2')
        fixed_airmass_scale = None
        if sampled_a2_index is None:
            try:
                fixed_a2 = float(base_physical.get('a2', 0.0))
            except (TypeError, ValueError):
                fixed_a2 = np.nan
            if not np.isfinite(fixed_a2):
                observed_values_valid = False
            elif fixed_a2 != 0.0:
                fixed_airmass_scale = np.exp(fixed_a2 * centered_airmass)

        free_flux_baseline_index = next(
            (sampled_key_index[key] for key in ('a0', 'a1') if key in sampled_key_index),
            None,
        )
        fixed_flux_baseline_value = None
        if free_flux_baseline_index is None and uses_fixed_flux_baseline:
            fixed_flux_baseline_value = get_flux_baseline(base_physical)

        duration_prior = self.duration_prior if isinstance(self.duration_prior, dict) else None
        duration_prior_applied = bool(duration_prior and duration_prior.get('applied'))
        try:
            expected_duration = float(duration_prior.get('expected_duration', np.nan)) if duration_prior else np.nan
            sigma_log_duration = float(duration_prior.get('sigma_log_duration', np.nan)) if duration_prior else np.nan
        except (TypeError, ValueError):
            expected_duration = np.nan
            sigma_log_duration = np.nan
        duration_prior_valid = (
            duration_prior_applied
            and np.isfinite(expected_duration)
            and expected_duration > 0
            and np.isfinite(sigma_log_duration)
            and sigma_log_duration > 0
        )

        def solve_flux_baseline_for_model(model):
            mask = baseline_static_mask & np.isfinite(model) & (model != 0)
            if not np.any(mask):
                return fallback_flux_baseline()

            masked_model = model[mask]
            masked_data = data[mask]
            masked_weights = baseline_weights[mask]
            denom = np.sum(masked_weights * masked_model ** 2)

            if not np.isfinite(denom) or denom <= 0:
                ratio = masked_data / masked_model
                ratio = ratio[np.isfinite(ratio)]
                if ratio.size == 0:
                    return fallback_flux_baseline()
                baseline = np.nanmedian(ratio)
                return baseline if np.isfinite(baseline) else fallback_flux_baseline()

            baseline = np.sum(masked_weights * masked_data * masked_model) / denom
            return baseline if np.isfinite(baseline) else fallback_flux_baseline()

        def physical_from_sample_point(sample_point):
            physical = base_physical.copy()
            for bound_key, index in direct_sample_assignments:
                physical[bound_key] = sample_point[index]
            if impact_parameter_index is not None:
                impact_parameter = sample_point[impact_parameter_index]
                physical['b'] = impact_parameter
                physical['inc'] = float(inclination_from_impact_parameter(physical, impact_parameter))
            return physical

        def single_loglike(pars):
            if not observed_values_valid:
                return BAD_LOG_LIKELIHOOD

            physical = physical_from_sample_point(pars)
            duration_loglike = 0.0
            if duration_prior_valid:
                duration = transit_duration(physical)
                if not np.isfinite(duration) or duration <= 0:
                    return BAD_LOG_LIKELIHOOD
                duration_log_residual = np.log(duration / expected_duration)
                duration_loglike = -0.5 * (duration_log_residual / sigma_log_duration) ** 2
            try:
                model = np.asarray(transit(time, physical), dtype=float)
                if sampled_a2_index is not None:
                    model *= np.exp(float(pars[sampled_a2_index]) * centered_airmass)
                elif fixed_airmass_scale is not None:
                    model *= fixed_airmass_scale

                if free_flux_baseline_index is not None:
                    model *= pars[free_flux_baseline_index]
                elif fixed_flux_baseline_value is not None:
                    model *= fixed_flux_baseline_value
                elif has_free_flux_baseline:
                    model *= get_flux_baseline(physical)
                else:
                    model *= solve_flux_baseline_for_model(model)
            except Exception:
                return BAD_LOG_LIKELIHOOD

            if model.shape != data_shape or not np.all(np.isfinite(model)):
                return BAD_LOG_LIKELIHOOD

            residuals = (data - model) * inverse_dataerr
            chi2 = np.sum(residuals * residuals)
            logl = -0.5 * chi2 + duration_loglike
            return float(logl) if np.isfinite(logl) else BAD_LOG_LIKELIHOOD

        def loglike(pars):
            pars_array = np.asarray(pars, dtype=float)
            if pars_array.ndim == 2:
                return np.fromiter(
                    (single_loglike(row) for row in pars_array),
                    dtype=float,
                    count=pars_array.shape[0],
                )
            return single_loglike(pars_array)

        prior_boundarray = np.array([self.bounds[k] for k in bound_keys], dtype=float)
        prior_lower_bounds = prior_boundarray[:, 0]
        prior_bound_widths = prior_boundarray[:, 1] - prior_lower_bounds
        prior_bound_index = {key: index for index, key in enumerate(bound_keys)}
        prior_inc_index = prior_bound_index.get('inc')

        def prior_values_for(sample_points, key, default):
            if key in prior_bound_index:
                return sample_points[:, prior_bound_index[key]]
            value = np.asarray(base_physical.get(key, default), dtype=float)
            if value.shape == ():
                return np.full(sample_points.shape[0], float(value), dtype=float)
            return np.broadcast_to(value, (sample_points.shape[0],)).astype(float)

        def prior_impact_upper_bounds(sample_points):
            rprs = prior_values_for(sample_points, 'rprs', np.nan)
            grazing_upper = np.where(np.isfinite(rprs) & (rprs >= 0), 1.0 + rprs, np.nan)

            ars = prior_values_for(sample_points, 'ars', np.nan)
            ecc = prior_values_for(sample_points, 'ecc', 0.0)
            omega = np.deg2rad(prior_values_for(sample_points, 'omega', 0.0))
            denom = 1.0 + ecc * np.sin(omega)
            denom = np.where(np.isclose(denom, 0.0), np.finfo(float).eps, denom)
            scale_upper = ars * (1.0 - ecc ** 2) / denom

            valid_grazing = np.isfinite(grazing_upper) & (grazing_upper > 0)
            valid_scale = np.isfinite(scale_upper) & (scale_upper > 0)
            upper = np.full(sample_points.shape[0], 1.0, dtype=float)

            both_valid = valid_grazing & valid_scale
            upper[both_valid] = np.minimum(grazing_upper[both_valid], scale_upper[both_valid])
            upper[valid_grazing & ~valid_scale] = grazing_upper[valid_grazing & ~valid_scale]
            upper[valid_scale & ~valid_grazing] = scale_upper[valid_scale & ~valid_grazing]
            return np.maximum(0.0, upper)

        def prior_transform(upars):
            upars_array = np.asarray(upars, dtype=float)
            sample_points = prior_lower_bounds + prior_bound_widths * upars_array
            if not self.impact_parameter_sampled_directly or prior_inc_index is None:
                return sample_points

            if upars_array.ndim == 2:
                upper_bounds = prior_impact_upper_bounds(sample_points)
                sample_points[:, prior_inc_index] = upper_bounds * upars_array[:, prior_inc_index]
                return sample_points

            upper_bound = prior_impact_upper_bounds(sample_points.reshape(1, -1))[0]
            sample_points[prior_inc_index] = upper_bound * upars_array[prior_inc_index]
            return sample_points

        self.ns_type = 'ultranest'
        test = ReactiveNestedSampler(sampled_keys, loglike, prior_transform, vectorized=True)

        run_kwargs = {"max_ncalls": int(self.max_ncalls)}
        if self.ultranest_min_num_live_points is not None:
            run_kwargs["min_num_live_points"] = int(self.ultranest_min_num_live_points)

        self.results = run_reactive_sampler(
            test,
            run_kwargs=run_kwargs,
            verbose=self.verbose,
        )

        if self.keep_ultranest_sampler:
            self._ultranest_resume_context = {
                'sampler': test,
                'bound_keys': list(bound_keys),
                'sampled_keys': list(sampled_keys),
                'physical_from_sample_point': physical_from_sample_point,
            }
        else:
            self._ultranest_resume_context = None
        self._finalize_ultranest_fit_results(bound_keys, sampled_keys, physical_from_sample_point)

        if not self.sample_parameters:
            self.sample_parameters = {
                key: self.parameters.get(key, self.sample_parameters.get(key))
                for key in self.sampled_keys
            }
        for bound_key, sampled_key in zip(bound_keys, sampled_keys):
            if sampled_key not in self.sample_errors and bound_key in self.errors:
                self.sample_errors[sampled_key] = self.errors[bound_key]
            if sampled_key not in self.sample_quantiles and bound_key in self.quantiles:
                self.sample_quantiles[sampled_key] = self.quantiles[bound_key]

        # final model
        self.create_fit_variables()

    def plot_bestfit(
        self,
        title="",
        bin_dt=30. / (60 * 24),
        zoom=False,
        phase=True,
        show_flux_baseline_label=True,
        show_model_uncertainty=False,
        show_baseline_uncertainty=False,
    ):
        f = plt.figure(figsize=(9, 6))
        f.subplots_adjust(top=0.92, bottom=0.09, left=0.14, right=0.98, hspace=0)
        ax_lc = plt.subplot2grid((4, 5), (0, 0), colspan=5, rowspan=3)
        ax_res = plt.subplot2grid((4, 5), (3, 0), colspan=5, rowspan=1)
        axs = [ax_lc, ax_res]

        axs[0].set_title(title)
        axs[0].set_ylabel("Relative Flux", fontsize=14)
        axs[0].grid(True, ls='--')

        rprs2 = self.parameters['rprs'] ** 2
        rprs2err = 2 * self.parameters['rprs'] * self.errors['rprs']
        lclabel1 = r"Area ratio $(R_{p}/R_{s})^{2}$ = %s $\pm$ %s" % (
            str(round_to_2(rprs2, rprs2err)),
            str(round_to_2(rprs2err))
        )

        lclabel2 = r"$T_{mid}$ = %s $\pm$ %s BJD$_{TDB}$" % (
            str(round_to_2(self.parameters['tmid'], self.errors.get('tmid', 0))),
            str(round_to_2(self.errors.get('tmid', 0)))
        )

        lclabel = lclabel1 + "\n" + lclabel2
        if show_flux_baseline_label and 'a0' in self.parameters:
            lclabel3 = r"$a_0$ = %s $\pm$ %s" % (
                str(round_to_2(self.parameters['a0'], self.errors.get('a0', 0))),
                str(round_to_2(self.errors.get('a0', 0)))
            )
            lclabel += "\n" + lclabel3

        if zoom:
            axs[0].set_ylim([1 - 1.25 * self.parameters['rprs'] ** 2, 1 + 0.5 * self.parameters['rprs'] ** 2])
        else:
            if phase:
                axs[0].errorbar(self.phase, self.detrended, yerr=np.std(self.residuals) / np.median(self.data),
                                ls='none', marker='.', color='black', ecolor='0.72',
                                elinewidth=1.0, zorder=1, alpha=1.0)
            else:
                axs[0].errorbar(self.time, self.detrended, yerr=np.std(self.residuals) / np.median(self.data),
                                ls='none', marker='.', color='black', ecolor='0.72',
                                elinewidth=1.0, zorder=1, alpha=1.0)

        if phase:
            si = np.argsort(self.phase)
            bt2, br2, _ = time_bin(self.phase[si] * self.parameters['per'],
                                   self.residuals[si] / np.median(self.data) * 1e2, bin_dt)
            axs[1].plot(self.phase, self.residuals / np.median(self.data) * 1e2, 'k.', alpha=1.0,
                        label=r'$\sigma$ = {:.2f} %'.format(np.std(self.residuals / np.median(self.data) * 1e2)))
            axs[1].plot(bt2 / self.parameters['per'], br2, 'bs', alpha=1, zorder=2)
            axs[1].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
            axs[1].set_xlabel("Phase", fontsize=14)

            si = np.argsort(self.phase)
            bt2, bf2, bs = time_bin(self.phase[si] * self.parameters['per'], self.detrended[si], bin_dt)
            axs[0].errorbar(bt2 / self.parameters['per'], bf2, yerr=bs, alpha=1, zorder=2, color='blue', ls='none',
                            marker='s')
            # axs[0].plot(self.phase[si], self.transit[si], 'r-', zorder=3, label=lclabel)
            sii = np.argsort(self.phase_upsample)
            if show_baseline_uncertainty:
                self._plot_baseline_model_uncertainty(
                    axs[0],
                    self.phase_upsample,
                    self.time_upsample,
                    sii,
                    label='_nolegend_',
                )
            if show_model_uncertainty:
                self._plot_transit_model_uncertainty(
                    axs[0],
                    self.phase_upsample,
                    self.time_upsample,
                    sii,
                    label='_nolegend_',
                )
            axs[0].plot(self.phase_upsample[sii], self.transit_upsample[sii], 'r-', zorder=3, label=lclabel)
            axs[0].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
            axs[0].set_xlabel("Phase ", fontsize=14)
        else:
            bt, br, _ = time_bin(self.time, self.residuals / np.median(self.data) * 1e2, bin_dt)
            axs[1].plot(self.time, self.residuals / np.median(self.data) * 1e2, 'k.', alpha=1.0,
                        label=r'$\sigma$ = {:.2f} %'.format(np.std(self.residuals / np.median(self.data) * 1e2)))
            axs[1].plot(bt, br, 'bs', alpha=1, zorder=2, label=r'$\sigma$ = {:.2f} %'.format(np.std(br)))
            axs[1].set_xlim([min(self.time_upsample), max(self.time_upsample)])
            axs[1].set_xlabel("Time [day]", fontsize=14)

            bt, bf, bs = time_bin(self.time, self.detrended, bin_dt)
            si = np.argsort(self.time)
            sii = np.argsort(self.time_upsample)
            axs[0].errorbar(bt, bf, yerr=bs, alpha=1, zorder=2, color='blue', ls='none', marker='s')
            if show_baseline_uncertainty:
                self._plot_baseline_model_uncertainty(
                    axs[0],
                    self.time_upsample,
                    self.time_upsample,
                    sii,
                    label='_nolegend_',
                )
            if show_model_uncertainty:
                self._plot_transit_model_uncertainty(
                    axs[0],
                    self.time_upsample,
                    self.time_upsample,
                    sii,
                    label='_nolegend_',
                )
            axs[0].plot(self.time_upsample[sii], self.transit_upsample[sii], 'r-', zorder=3, label=lclabel)
            axs[0].set_xlim([min(self.time_upsample), max(self.time_upsample)])
            axs[0].set_xlabel("Time [day]", fontsize=14)

        axs[0].get_xaxis().set_visible(False)
        axs[1].legend(loc='best')
        axs[0].legend(loc='best')
        axs[1].set_ylabel("Residuals [%]", fontsize=14)
        axs[1].grid(True, ls='--', axis='y')
        return f, axs

    def plot_triangle(self, plot_title=None, zoom_sigma=None):
        payload = self._get_triangle_plot_payload()
        if zoom_sigma is not None:
            payload = dict(payload)
            payload['ranges'] = self._triangle_plot_sigma_window_ranges(payload, zoom_sigma)
            payload = self._recenter_triangle_plot_payload_for_visible_ranges(payload)

        chi2 = payload['display_logl'] * -2
        parameter_count = max(1, len(payload['sampled_keys']))
        fig_size = max(9.0, 2.35 * parameter_count)
        mask1 = np.ones(len(chi2), dtype=bool)
        mask2 = np.ones(len(chi2), dtype=bool)
        mask3 = np.ones(len(chi2), dtype=bool)

        for i, key in enumerate(payload['sampled_keys']):
            if key in ('a0', 'a1', 'a2'):
                continue

            center = payload['mask_centers'][i]
            error = payload['mask_errors'][i]
            if not np.isfinite(center) or not np.isfinite(error) or error <= 0:
                continue

            values = payload['mask_values'][:, i]
            mask3 = mask3 & (values > (center - 3 * error)) & (values < (center + 3 * error))
            mask1 = mask1 & (values > (center - error)) & (values < (center + error))
            mask2 = mask2 & (values > (center - 2 * error)) & (values < (center + 2 * error))

        if not np.any(mask1):
            mask1 = np.ones(len(chi2), dtype=bool)
        if not np.any(mask2):
            mask2 = np.ones(len(chi2), dtype=bool)
        if not np.any(mask3):
            mask3 = np.ones(len(chi2), dtype=bool)

        label_kwargs = {
            'labelpad': 10,
            'fontsize': 10,
        }
        title_kwargs = {
            'loc': 'left',
            'pad': 4,
            'fontsize': 11,
        }

        fig = corner(payload['display_points'],
                     labels=payload['labels'],
                     bins=int(np.sqrt(payload['display_points'].shape[0])),
                     range=payload['ranges'],
                      weights=payload['display_weights'],
                      plot_contours=True,
                      levels=self._triangle_contour_levels(chi2, mask1, mask2, mask3),
                      plot_density=False,
                      titles=payload['titles'],
                      truths=payload['truths'],
                      data_kwargs={
                          'c': chi2,
                          'vmin': np.percentile(chi2[mask3], 1),
                          'vmax': np.percentile(chi2[mask3], 95),
                          'cmap': 'viridis',
                          's': 1.6,
                          'alpha': 0.38,
                      },
                      label_kwargs=label_kwargs,
                      title_kwargs=title_kwargs,
                      hist_kwargs={
                          'color': 'black',
                      }
                      )
        if hasattr(fig, 'set_size_inches'):
            fig.set_size_inches(fig_size, fig_size, forward=True)
        if plot_title and hasattr(fig, 'suptitle'):
            fig.suptitle(plot_title, fontsize=13, y=0.99)
        self._adjust_triangle_plot_layout(fig)
        self._overlay_triangle_plot_geometry_histograms(
            fig,
            payload,
            title_kwargs=title_kwargs,
            label_kwargs=label_kwargs,
        )
        return fig

# simultaneously fit multiple data sets with global and local parameters
class glc_fitter(lc_fitter):
    # needed for lc_fitter
    ns_type = 'ultranest'

    def __init__(self, input_data, global_bounds, local_bounds, individual_fit=False, stdev_cutoff=0.03, verbose=False):
        # keys for input_data: time, flux, ferr, airmass, priors all numpy arrays
        self.lc_data = copy.deepcopy(input_data)
        self.global_bounds = global_bounds
        self.local_bounds = local_bounds
        self.individual_fit = individual_fit
        self.stdev_cutoff = stdev_cutoff
        self.verbose = verbose
        self.results = None
        self.fit_nested()

    def fit_nested(self):

        # create bound arrays for generating samples
        nobs = len(self.lc_data)
        gfreekeys = list(self.global_bounds.keys())

        # if isinstance(self.local_bounds, dict):
        #     lfreekeys = list(self.local_bounds.keys())
        #     boundarray = np.vstack([ [self.global_bounds[k] for k in gfreekeys], [self.local_bounds[k] for k in lfreekeys]*nobs ])
        # else:
        #     # if list type
        lfreekeys = []
        boundarray = [self.global_bounds[k] for k in gfreekeys]
        for i in range(nobs):
            lfreekeys.append(list(self.local_bounds[i].keys()))
            boundarray.extend([self.local_bounds[i][k] for k in lfreekeys[-1]])
        boundarray = np.array(boundarray)

        # fit individual light curves to constrain priors
        if self.individual_fit:
            for i in range(nobs):

                print(f"Fitting individual light curve {i+1}/{nobs}")
                try:
                    mybounds = dict(**self.local_bounds[i], **self.global_bounds)
                except:
                    mybounds = {}
                    for k in self.local_bounds[i]:
                        mybounds[k] = self.local_bounds[i][k]
                    for k in self.global_bounds:
                        mybounds[k] = self.global_bounds[k]
                if 'per' in mybounds: del(mybounds['per'])
                if 'inc' in mybounds and 'rprs' in mybounds: del(mybounds['inc'])
                if 'tmid' in mybounds:
                    # find the closet mid transit time to the last observation
                    phase = (self.lc_data[i]['time'][-1] - self.lc_data[i]['priors']['tmid']) / self.lc_data[i]['priors']['per']
                    nepochs = np.round(phase)
                    newtmid = self.lc_data[i]['priors']['tmid'] + nepochs * self.lc_data[i]['priors']['per']
                    err = np.diff(mybounds['tmid'])[0]/2.
                    mybounds['tmid'] = [newtmid - err, newtmid + err ]

                # fit individual light curve
                myfit = lc_fitter(
                    self.lc_data[i]['time'],
                    self.lc_data[i]['flux'],
                    self.lc_data[i]['ferr'],
                    self.lc_data[i]['airmass'],
                    self.lc_data[i]['priors'],
                    mybounds
                )

                # check stdev_cutoff and residuals
                if myfit.res_stdev > self.stdev_cutoff:
                    print(f"WARNING: Stdev of residuals is large! {myfit.res_stdev:.3f} > {self.stdev_cutoff:.3f}")
                    #raise ValueError(f"Stdev of residuals is too large, please remove id: {self.lc_data[i]['name']}

                # copy data over for individual fits
                self.lc_data[i]['individual'] = myfit.parameters.copy()
                self.lc_data[i]['individual_err'] = myfit.errors.copy()
                self.lc_data[i]['res_stdev'] = myfit.res_stdev
                self.lc_data[i]['quality'] = myfit.quality

                ti = sum([len(self.local_bounds[k]) for k in range(i)])
                # update local priors
                for j, key in enumerate(self.local_bounds[i].keys()):

                    boundarray[j+ti+len(gfreekeys),0] = myfit.parameters[key] - 5*myfit.errors[key]
                    boundarray[j+ti+len(gfreekeys),1] = myfit.parameters[key] + 5*myfit.errors[key]

                    if key == 'rprs':
                        boundarray[j+ti+len(gfreekeys),0] = max(0,myfit.parameters[key] - 5*myfit.errors[key])

                # print name and stdev of residuals
                mint = np.min(self.lc_data[i]['time'])
                maxt = np.max(self.lc_data[i]['time'])
                try:
                    print(f"{self.lc_data[i]['name']} & {Time(mint,format='jd').isot} & {Time(maxt,format='jd').isot} & {np.std(myfit.residuals)} & {len(self.lc_data[i]['time'])}")
                except:
                    print(f"{self.lc_data[i]['name']} & {mint} & {maxt} & {np.std(myfit.residuals)} & {len(self.lc_data[i]['time'])}")

                del(myfit)

        # transform unit cube to prior volume
        bounddiff = np.diff(boundarray,1).reshape(-1)
        def prior_transform(upars):
            return (boundarray[:,0] + bounddiff*upars)

        def loglike(pars):
            chi2 = 0

            # for each light curve
            for i in range(nobs):

                # global keys
                for j, key in enumerate(gfreekeys):
                    self.lc_data[i]['priors'][key] = pars[j]

                # local keys
                ti = sum([len(self.local_bounds[k]) for k in range(i)])
                for j, key in enumerate(lfreekeys[i]):
                    self.lc_data[i]['priors'][key] = pars[j+ti+len(gfreekeys)]

                # compute model
                model = transit(self.lc_data[i]['time'], self.lc_data[i]['priors'])
                model *= airmass_trend(
                    self.lc_data[i]['priors'].get('a2', 0),
                    self.lc_data[i]['airmass'],
                )
                if has_explicit_flux_baseline(self.global_bounds) or has_explicit_flux_baseline(self.local_bounds[i]):
                    model *= get_flux_baseline(self.lc_data[i]['priors'])
                else:
                    model *= solve_flux_baseline(model, self.lc_data[i]['flux'], self.lc_data[i]['ferr'])

                # add to chi2
                chi2 += np.sum( ((self.lc_data[i]['flux']-model)/self.lc_data[i]['ferr'])**2 )

            # maximization metric for nested sampling
            return -0.5*chi2

        freekeys = []+gfreekeys
        for n in range(nobs):
            for k in lfreekeys[n]:
                #clean_name = self.lc_data[n].get('name', n).replace(' ','_').replace('(','').replace(')','').replace('[','').replace(']','').replace('-','_').split('-')[0]
                freekeys.append(f"local_{k}_{n}")

        sampler = ReactiveNestedSampler(freekeys, loglike, prior_transform)
        self.results = run_reactive_sampler(
            sampler,
            run_kwargs={"max_ncalls": int(1e6)},
            verbose=self.verbose,
        )

        self.quantiles = {}
        self.errors = {}
        self.parameters = self.lc_data[0]['priors'].copy()

        for i, key in enumerate(freekeys):
            self.parameters[key] = self.results['maximum_likelihood']['point'][i]
            #self.errors[key] = self.results['posterior']['median'][i]

            self.errors[key] = self.results['posterior']['stdev'][i]
            self.quantiles[key] = [
                self.results['posterior']['errlo'][i],
                self.results['posterior']['errup'][i]]

        # create an average Rp/Rs if it is not in global keys
        # check if 'rprs' is in lfreekeys
        rprs_in_local = False
        for i in range(nobs):
            if 'rprs' in lfreekeys[i]:
                rprs_in_local = True
                break

        if rprs_in_local:
            local_rprs = [] # used for creating an average value
            local_rprs_err = []

        # loop over observations
        for n in range(nobs):
            self.lc_data[n]['errors'] = {}
            
            # set global parameters without overwriting everything
            for gk in gfreekeys:
                self.lc_data[n]['priors'][gk] = self.parameters[gk]
                self.lc_data[n]['errors'][gk] = self.errors[gk]

            # loop over local keys and save best fit values
            for k in lfreekeys[n]:

                # create key to get results
                pkey = f"local_{k}_{n}"
                
                # overwrite priors with best fit value
                self.lc_data[n]['priors'][k] = self.parameters[pkey]
                self.lc_data[n]['errors'][k] = self.errors[pkey]

                # update key for final bestfit plot if needed
                if k == 'rprs':
                    local_rprs.append(self.lc_data[n]['priors'][k])
                    local_rprs_err.append(self.lc_data[n]['errors'][k])

            # solve for the local baseline flux scale
            model = transit(self.lc_data[n]['time'], self.lc_data[n]['priors'])
            airmass = airmass_trend(
                self.lc_data[n]['priors'].get('a2', 0),
                self.lc_data[n]['airmass'],
            )
            if has_explicit_flux_baseline(self.global_bounds) or has_explicit_flux_baseline(self.local_bounds[n]):
                flux_scale = get_flux_baseline(self.lc_data[n]['priors'])
                flux_scale_err = self.lc_data[n]['errors'].get('a0', self.lc_data[n]['errors'].get('a1', 0))
            else:
                flux_scale = solve_flux_baseline(model * airmass, self.lc_data[n]['flux'], self.lc_data[n]['ferr'])
                flux_scale_err = solve_flux_baseline_uncertainty(model * airmass, self.lc_data[n]['ferr'])
            self.lc_data[n]['priors']['a0'] = flux_scale
            self.lc_data[n]['priors']['a1'] = flux_scale
            self.lc_data[n]['errors']['a0'] = flux_scale_err
            self.lc_data[n]['errors']['a1'] = flux_scale_err
            self.lc_data[n]['residuals'] = self.lc_data[n]['flux'] - model * airmass * flux_scale
            self.lc_data[n]['detrend'] = self.lc_data[n]['flux'] / (airmass * flux_scale)

            # phase
            plot_time_range = normalize_time_range(self.lc_data[n].get('plot_time_range'))
            if plot_time_range is None:
                plot_time_range = normalize_time_range(self.lc_data[n]['time'])
            self.lc_data[n]['plot_time_range'] = plot_time_range
            self.lc_data[n]['phase'] = get_plot_phase(
                self.lc_data[n]['time'],
                self.lc_data[n]['priors']['per'],
                self.lc_data[n]['priors']['tmid'],
                plot_time_range,
            )
            self.lc_data[n]['time_upsample'] = np.linspace(plot_time_range[0], plot_time_range[1], 1000)
            self.lc_data[n]['phase_upsample'] = get_plot_phase(
                self.lc_data[n]['time_upsample'],
                self.lc_data[n]['priors']['per'],
                self.lc_data[n]['priors']['tmid'],
                plot_time_range,
            )
            self.lc_data[n]['transit_upsample'] = transit(self.lc_data[n]['time_upsample'], self.lc_data[n]['priors'])

        # create an average value from all the local fits, used for plotting final best fit
        if rprs_in_local:
            self.parameters['rprs'] = np.mean(local_rprs)
            self.errors['rprs'] = np.std(local_rprs)


    def plot_bestfits(self):
        nrows = len(self.lc_data)//4+1
        # make sure there isn't an extra row
        if len(self.lc_data)%4 == 0:
            nrows -= 1

        fig,ax = plt.subplots(nrows, 4, figsize=(5+5*nrows, 5*nrows))

        # turn off all axes
        for i in range(nrows*4):
            ri = int(i/4)
            ci = i%4
            if ax.ndim == 1:
                ax[i].axis('off')
            else:
                ax[ri,ci].axis('off')

        # cycle the colors and markers
        markers = cycle(['o','v','^','<','>','s','*','h','H','D','d','P','X'])
        colors = cycle(['black','blue','green','orange','purple','grey','magenta','cyan','lime'])

        # plot observations
        for i in range(len(self.lc_data)):
            ri = int(i/4)
            ci = i%4
            ncolor = next(colors)
            nmarker = next(markers)

            model = transit(self.lc_data[i]['time'], self.lc_data[i]['priors'])
            airmass = airmass_trend(
                self.lc_data[i]['priors'].get('a2', 0),
                self.lc_data[i]['airmass'],
            )
            detrend = self.lc_data[i]['flux'] / (model * airmass)

            if ax.ndim == 1:
                ax[i].axis('on')
                ax[i].errorbar(self.lc_data[i]['time'], self.lc_data[i]['flux']/airmass/detrend.mean(), yerr=self.lc_data[i]['ferr']/airmass/detrend.mean(), 
                                ls='none', marker=nmarker, color=ncolor, alpha=0.5, zorder=1)
                
                ax[i].plot(self.lc_data[i]['time_upsample'], self.lc_data[i]['transit_upsample'], 'r-', zorder=2)
                ax[i].set_xlabel("Time [BJD]", fontsize=14)
                ax[i].set_ylabel("Relative Flux", fontsize=14)
                ax[i].set_title(f"{self.lc_data[i].get('name','')}", fontsize=16)
            else:
                ax[ri,ci].axis('on')
                ax[ri,ci].errorbar(self.lc_data[i]['time'], self.lc_data[i]['flux']/airmass/detrend.mean(), yerr=self.lc_data[i]['ferr']/airmass/detrend.mean(), 
                                   ls='none', marker=nmarker, color=ncolor, alpha=0.5, zorder=1)
                ax[ri,ci].plot(self.lc_data[i]['time_upsample'], self.lc_data[i]['transit_upsample'], 'r-', zorder=2)
                ax[ri,ci].set_xlabel("Time[BJD]", fontsize=14)
                ax[ri,ci].set_ylabel("Relative Flux", fontsize=14)
                ax[ri,ci].set_title(f"{self.lc_data[i].get('name','')}", fontsize=16)

        plt.tight_layout()
        return fig

    def plot_bestfit(self, title="", bin_dt=30./(60*24), alpha=0.05, ylim_sigma=5, phase_limits='median', show_legend=True, limit_legend=False, show_individual_fits=False):
        """
        Plot the best fit model and residuals

        Parameters
        ----------
        title : str
            Title for the plot

        bin_dt : float
            Bin size for plotting the residuals

        alpha : float
            Alpha value for plotting the data

        ylim_sigma : float
            Number of sigma to plot the residuals

        phase_limits : str
            'median' or 'all' to set the phase limits

        show_legend : bool
            Show the legend

        limit_legend : bool
            Limit the legend to 3 entries
        """
        f = plt.figure(figsize=(15,12))
        f.subplots_adjust(top=0.92,bottom=0.09,left=0.1,right=0.98, hspace=0)
        ax_lc = plt.subplot2grid((4,5), (0,0), colspan=5,rowspan=3)
        ax_res = plt.subplot2grid((4,5), (3,0), colspan=5, rowspan=1)
        axs = [ax_lc, ax_res]

        axs[0].set_title(title, fontsize=18)
        axs[0].set_ylabel("Relative Flux", fontsize=14)
        axs[0].grid(True,ls='--')

        try:
            rprs2 = self.parameters['rprs']**2
            rprs2err = 2*self.parameters['rprs']*self.errors['rprs']
        except:
            rprs2 = self.lc_data[0]['priors']['rprs']**2
            rprs2err = 2*self.lc_data[0]['priors']['rprs']*self.lc_data[0]['errors']['rprs']

        lclabel1 = r"Area ratio $(R_{p}/R_{s})^{2}$ = %s $\pm$ %s" %(
            str(round_to_2(rprs2, rprs2err)),
            str(round_to_2(rprs2err))
        )
        
        lclabel2 = r"$T_{mid}$ = %s $\pm$ %s BJD$_{TDB}$" %(
            str(round_to_2(self.parameters['tmid'], self.errors.get('tmid',0))),
            str(round_to_2(self.errors.get('tmid',0)))
        )

        lclabel = lclabel1 + "\n" + lclabel2
        minp = 1
        maxp = 0

        min_std = 1
        # cycle the colors and markers
        markers = cycle(['o','v','^','<','>','s','*','h','H','D','d','P','X'])
        colors = cycle(['black','blue','green','orange','purple','grey','magenta','cyan','lime'])

        alldata = {
            'time': [],
            'phase': [],
            'flux': [],
            'detrend': [],
            'ferr': [],
            'residuals': [],
        }

        for n in range(len(self.lc_data)):
            ncolor = next(colors)
            nmarker = next(markers)
            alldata['time'].extend(self.lc_data[n]['time'].tolist())
            alldata['phase'].extend(self.lc_data[n]['phase'].tolist())
            alldata['detrend'].extend(self.lc_data[n]['detrend'].tolist())
            alldata['flux'].extend(self.lc_data[n]['flux'].tolist())
            alldata['ferr'].extend(self.lc_data[n]['ferr'].tolist())
            alldata['residuals'].extend(self.lc_data[n]['residuals'].tolist())
            
            phase = self.lc_data[n]['phase']
            si = np.argsort(phase)
            #bt2, br2, _ = time_bin(phase[si]*self.parameters['per'], self.lc_data[n]['residuals'][si]/np.median(self.lc_data[n]['flux'])*1e2, bin_dt)

            # plot data
            axs[0].errorbar(phase, self.lc_data[n]['detrend'], yerr=np.std(self.lc_data[n]['residuals'])/np.median(self.lc_data[n]['flux']), 
                            ls='none', marker=nmarker, color=ncolor, zorder=1, alpha=alpha)

            # plot residuals
            axs[1].plot(phase, self.lc_data[n]['residuals']/np.median(self.lc_data[n]['flux'])*1e2, color=ncolor, marker=nmarker, ls='none',
                         alpha=0.2)

            # plot binned data
            bt2, bf2, bs = time_bin(phase[si]*self.lc_data[n]['priors']['per'], self.lc_data[n]['detrend'][si], bin_dt)

            if limit_legend:
                axs[0].errorbar(bt2/self.lc_data[n]['priors']['per'],bf2,yerr=bs,alpha=1,zorder=2,color=ncolor,ls='none',marker=nmarker)
            else:
                axs[0].errorbar(bt2/self.lc_data[n]['priors']['per'],bf2,yerr=bs,alpha=1,zorder=2,color=ncolor,ls='none',marker=nmarker,
                                label=r'{}: {:.2f} %'.format(self.lc_data[n].get('name',''),np.std(self.lc_data[n]['residuals']/np.median(self.lc_data[n]['flux'])*1e2)))

            # replace min and max for upsampled lc model
            minp = min(minp, min(self.lc_data[n]['phase_upsample']))
            maxp = max(maxp, max(self.lc_data[n]['phase_upsample']))
            min_std = min(min_std, np.std(self.lc_data[n]['residuals']/np.median(self.lc_data[n]['flux'])))

            # plot individual best fit models
            if show_individual_fits:
                axs[0].plot(self.lc_data[n]['phase_upsample'], self.lc_data[n]['transit_upsample'], color=ncolor, zorder=3, alpha=0.5)

        # create binned plot for all the data
        for k in alldata.keys():
            alldata[k] = np.array(alldata[k])
            
        phase = alldata['phase']
        si = np.argsort(phase)
        bt, br, _ = time_bin(phase[si]*self.parameters['per'], alldata['residuals'][si]/np.median(alldata['flux']), 2*bin_dt)
        bt, bf, bs = time_bin(phase[si]*self.parameters['per'], alldata['detrend'][si], 2*bin_dt)

        axs[0].errorbar(bt/self.parameters['per'],bf,yerr=bs,alpha=1,zorder=2,color='white',ls='none',marker='o',ms=15,
                        markeredgecolor='black',
                        ecolor='black',
                        label=r'Binned Data: {:.2f} %'.format(np.std(br)*1e2))

        axs[1].plot(bt/self.parameters['per'],br*1e2,color='white',ls='none',marker='o',ms=11,markeredgecolor='black')

        # best fit model
        self.phase_upsample = np.linspace(minp, maxp, 10000)
        self.time_upsample = self.parameters['tmid'] + self.phase_upsample * self.parameters['per']
        self.transit_upsample = transit(self.time_upsample, self.parameters)
        axs[0].plot(self.phase_upsample, self.transit_upsample, 'r-', zorder=3, label=lclabel, lw=3)

        # set up axes limits
        axs[0].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
        axs[0].set_xlabel("Phase ", fontsize=14)
        axs[0].set_ylim([1-self.parameters['rprs']**2-ylim_sigma*min_std, 1+ylim_sigma*min_std])
        axs[1].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
        axs[1].set_xlabel("Phase", fontsize=14)
        axs[1].set_ylim([-5*min_std*1e2, 5*min_std*1e2])

        # compute average min and max for all the data
        mins = []; maxs = []
        for n in range(len(self.lc_data)):
            mins.append(min(self.lc_data[n]['phase_upsample']))
            maxs.append(max(self.lc_data[n]['phase_upsample']))

        # set up phase limits
        if isinstance(phase_limits, str):
            if phase_limits == "minmax":
                axs[0].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
                axs[1].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
            elif phase_limits == "median":
                axs[0].set_xlim([np.median(mins), np.median(maxs)])
                axs[1].set_xlim([np.median(mins), np.median(maxs)])
            else:
                axs[0].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
                axs[1].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
        elif isinstance(phase_limits, list):
            axs[0].set_xlim([phase_limits[0], phase_limits[1]])
            axs[1].set_xlim([phase_limits[0], phase_limits[1]])
        elif isinstance(phase_limits, tuple):
            axs[0].set_xlim([phase_limits[0], phase_limits[1]])
            axs[1].set_xlim([phase_limits[0], phase_limits[1]])
        else:
            axs[0].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
            axs[1].set_xlim([min(self.phase_upsample), max(self.phase_upsample)])

        axs[0].get_xaxis().set_visible(False)
        axs[1].set_ylabel("Residuals [%]", fontsize=14)
        axs[1].grid(True,ls='--',axis='y')
    
        if show_legend:
            axs[0].legend(loc='best',ncol=len(self.lc_data)//7+1)
    
        return f,axs

    def plot_stack(self, title="", bin_dt=30./(60*24), dy=0.02):
        f, ax = plt.subplots(1,figsize=(9,12))
        
        ax.set_title(title)
        ax.set_ylabel("Relative Flux", fontsize=14)
        ax.grid(True,ls='--')

        rprs2 = self.parameters['rprs']**2
        rprs2err = 2*self.parameters['rprs']*self.errors['rprs']
        lclabel1 = r"Area ratio $(R_{p}/R_{s})^{2}$ = %s $\pm$ %s" %(
            str(round_to_2(rprs2, rprs2err)),
            str(round_to_2(rprs2err))
        )
        
        lclabel2 = r"$T_{mid}$ = %s $\pm$ %s BJD$_{TDB}$" %(
            str(round_to_2(self.parameters['tmid'], self.errors.get('tmid',0))),
            str(round_to_2(self.errors.get('tmid',0)))
        )

        lclabel = lclabel1 + "\n" + lclabel2
        minp = 1
        maxp = 0

        min_std = 1
        # cycle the colors and markers
        markers = cycle(['o','v','^','<','>','s','*','h','H','D','d','P','X'])
        colors = cycle(['black','blue','green','orange','purple','grey','magenta','cyan','lime'])
        for n in range(len(self.lc_data)):
            ncolor = next(colors)
            nmarker = next(markers)

            phase = self.lc_data[n]['phase']
            si = np.argsort(phase)
            bt2, br2, _ = time_bin(phase[si]*self.parameters['per'], self.lc_data[n]['residuals'][si]/np.median(self.lc_data[n]['flux'])*1e2, bin_dt)
            
            # plot data
            ax.errorbar(phase, self.lc_data[n]['detrend']-n*dy, yerr=np.std(self.lc_data[n]['residuals'])/np.median(self.lc_data[n]['flux']), 
                            ls='none', marker=nmarker, color=ncolor, zorder=1, alpha=0.25)
        
            # plot binned data
            bt2, bf2, bs = time_bin(phase[si]*self.lc_data[n]['priors']['per'], self.lc_data[n]['detrend'][si]-n*dy, bin_dt)
            ax.errorbar(bt2/self.lc_data[n]['priors']['per'],bf2,yerr=bs,alpha=1,zorder=2,color=ncolor,ls='none',marker=nmarker)

            # replace min and max for upsampled lc model
            minp = min(minp, min(self.lc_data[n]['phase_upsample']))
            maxp = max(maxp, max(self.lc_data[n]['phase_upsample']))
            min_std = min(min_std, np.std(self.lc_data[n]['residuals']/np.median(self.lc_data[n]['flux'])))

            # best fit model
            self.phase_upsample = np.linspace(minp, maxp, 10000)
            self.time_upsample = self.parameters['tmid'] + self.phase_upsample * self.parameters['per']
            self.transit_upsample = transit(self.time_upsample, self.parameters)
            ax.plot(self.phase_upsample, self.transit_upsample-n*dy, ls='-', color=ncolor, zorder=3, label=self.lc_data[n].get('name',''))

        ax.set_xlim([min(self.phase_upsample), max(self.phase_upsample)])
        ax.set_xlabel("Phase ", fontsize=14)
        ax.set_ylim([1-self.parameters['rprs']**2-5*min_std-n*dy, 1+5*min_std])
        ax.get_xaxis().set_visible(False)
        ax.legend(loc='best')
        return f,ax


if __name__ == "__main__":

    prior = {
        'rprs': 0.02,  # Rp/Rs
        'ars': 14.25,  # a/Rs
        'per': 3.33,  # Period [day]
        'inc': 88.5,  # Inclination [deg]
        'u0': 0, 'u1': 0, 'u2': 0, 'u3': 0,  # limb darkening (nonlinear)
        'ecc': 0.5,  # Eccentricity
        'omega': 120,  # Arg of periastron
        'tmid': 0.75,  # Time of mid transit [day],
        'a0': 50,  # Baseline flux normalization
        'a2': 0.,  # trend = a0 * np.exp(a2 * (airmass - mean(airmass)))

        'teff': 5000,
        'tefferr': 50,
        'met': 0,
        'meterr': 0,
        'logg': 3.89,
        'loggerr': 0.01
    }

    # example generating LD coefficients
    from pylightcurve import exotethys

    u0, u1, u2, u3 = exotethys(prior['logg'], prior['teff'], prior['met'], 'TESS', method='claret',
                               stellar_model='phoenix')

    prior['u0'], prior['u1'], prior['u2'], prior['u3'] = u0, u1, u2, u3

    time = np.linspace(0.7, 0.8, 1000)  # [day]

    # simulate extinction from airmass
    stime = time - time[0]
    alt = 90 * np.cos(4 * stime - np.pi / 6)
    # airmass = 1./np.cos(np.deg2rad(90-alt))
    airmass = np.zeros(time.shape[0])

    # GENERATE NOISY DATA
    data = transit(time, prior) * prior['a0'] * airmass_trend(prior['a2'], airmass)
    data += np.random.normal(0, prior['a0'] * 250e-6, len(time))
    dataerr = np.random.normal(300e-6, 50e-6, len(time)) + np.random.normal(300e-6, 50e-6, len(time))

    # add bounds for free parameters only
    mybounds = {
        'rprs': [0, 0.1],
        'tmid': [prior['tmid'] - 0.01, prior['tmid'] + 0.01],
        'ars': [13, 15],
        # 'a0': [0.95 * prior['a0'], 1.05 * prior['a0']],  # optional explicit baseline offset
        # 'a2': [0, 0.3] # uncomment if you want to fit for airmass
        # if a0 is omitted, the normalization is solved analytically during the fit
        # never list both 'a0' and 'a1' in bounds because they are the same scale term
    }

    myfit = lc_fitter(time, data, dataerr, airmass, prior, mybounds, mode='ns')

    for k in myfit.bounds.keys():
        print(f"{myfit.parameters[k]:.6f} +- {myfit.errors[k]}")

    fig, axs = myfit.plot_bestfit()
    plt.tight_layout()
    plt.show()

    fig = myfit.plot_triangle()
    plt.tight_layout()
    plt.show()
