"""Timing search regressions using a deterministic photometric likelihood.

The grid fitter integrates the likelihood directly, avoiding sampler noise when
checking that timing spread follows the light curve rather than its input box.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from exotic import exotic as ex


def gaussian_transit(times, center, depth=.018):
    return 1 - depth * np.exp(-.5 * ((np.asarray(times) - center) / .01) ** 2)


def grid_fitter(calls, depth=.018):
    def fit(times, flux, errors, airmass, prior, bounds, mode=None, **kwargs):
        calls.append({key: list(value) for key, value in bounds.items()})
        lower, upper = bounds['tmid']
        grid = np.linspace(lower, upper, 4001)
        models = gaussian_transit(np.asarray(times)[:, None], grid, depth)
        logl = -.5 * np.sum(((np.asarray(flux)[:, None] - models) / np.asarray(errors)[:, None]) ** 2, axis=0)
        weights = np.exp(logl - logl.max())
        weights /= weights.sum()
        mean = grid[0] + np.sum(weights * (grid - grid[0]))
        sigma = np.sqrt(np.sum(weights * (grid - mean) ** 2))
        result = SimpleNamespace(
            parameters={**prior, 'tmid': grid[np.argmax(logl)]}, errors={'tmid': sigma},
            sampled_keys=['tmid'], bounds=bounds, fixed_parameter_errors=kwargs.get('fixed_parameter_errors', {}),
        )
        result._get_triangle_plot_samples = lambda: (grid[:, None], None, weights)
        return result
    return fit


@pytest.mark.parametrize('initialization', ['disabled', 'failed'])
@pytest.mark.parametrize('partial', [False, True])
def test_precise_ephemeris_does_not_set_photometric_timing_or_uncertainty(monkeypatch, initialization, partial):
    t0, period = 2455706.82976910, 2.470613372
    predicted = t0 - 96 * period
    injected = predicted - 1.4 / 1440
    times = predicted + np.linspace(-.12, -.005 if partial else .12, 121)
    flux = gaussian_transit(times, injected)
    prior = {'tmid': predicted, 'per': period, 'rprs': np.sqrt(.018), 'ars': 12., 'inc': 89.}
    summary = ex.estimate_ephemeris_tmid_and_bounds(
        times, t0, period, 4.1e-6, 2.4e-8, expected_duration=.075, sigma_multiplier=35.,
    )
    bounds = summary['bounds']
    calls, fitted = [], []
    monkeypatch.setattr(ex, 'lc_fitter', grid_fitter(calls))
    if initialization == 'failed':
        class FailingBLS:
            def __init__(self, *args, **kwargs):
                pass

            def power(self, *args, **kwargs):
                raise RuntimeError('Injected initializer failure')
        monkeypatch.setattr(ex, 'BoxLeastSquares', FailingBLS)
        eebls = ex.estimate_tmid_and_bounds_with_eebls(times, flux, np.full(times.size, .003), prior, bounds)
        assert not eebls['applied']
        assert 'Injected initializer failure' in eebls['note']
        bounds = eebls['bounds']
    for noise in (.003, .006):
        fit = ex.run_nested_lightcurve_fit_with_rprs_posterior_retry(
            times, flux, np.full(times.size, noise), np.ones(times.size), prior, {'tmid': bounds},
            use_lm_boundary_scout=False,
        )
        fitted.append(fit)
        assert abs(fit.parameters['tmid'] - injected) * 86400 < 5
        assert fit.errors['tmid'] * 86400 > 30  # Original uniform box reported 3.1 seconds.
        assert not fit.tmid_bound_diagnostics['bound_influenced']
    if not partial:
        assert fitted[1].errors['tmid'] / fitted[0].errors['tmid'] == pytest.approx(2., rel=.08)


def test_posterior_refits_a_tiny_interval_until_shifted_transit_is_measured(monkeypatch):
    times = np.linspace(-.12, .12, 121)
    injected = -.004
    calls = []
    monkeypatch.setattr(ex, 'lc_fitter', grid_fitter(calls))
    fit = ex.run_nested_lightcurve_fit_with_rprs_posterior_retry(
        times, gaussian_transit(times, injected), np.full(times.size, .003), np.ones(times.size),
        {'tmid': 0., 'per': 2.47}, {'tmid': [-.000063, .000063]},
    )
    assert len(calls) > 1
    assert abs(fit.parameters['tmid'] - injected) * 86400 < 5
    assert fit.errors['tmid'] * 86400 > 30
    assert not fit.tmid_bound_diagnostics['bound_influenced']
    assert fit.tmid_posterior_refit_applied
    assert fit.tmid_bound_diagnostics['refit_history'][-1]['stdev_change_seconds'] > 0
    assert fit.posterior_refit_final_bounds == calls[-1]
    for previous, current in zip(calls, calls[1:]):
        assert current['tmid'][0] <= previous['tmid'][0]
        assert current['tmid'][1] >= previous['tmid'][1]


def test_flat_posterior_remains_reportable_at_existing_period_cap(monkeypatch):
    times = np.linspace(-.1, .1, 61)
    calls = []
    monkeypatch.setattr(ex, 'lc_fitter', grid_fitter(calls, depth=0.))
    fit = ex.run_nested_lightcurve_fit_with_rprs_posterior_retry(
        times, np.ones(times.size), np.full(times.size, .01), np.ones(times.size),
        {'tmid': 0., 'per': 1.}, {'tmid': [-.01, .01]},
    )
    assert len(calls) == 4
    assert fit.posterior_refit_final_bounds['tmid'] == [-.25, .25]
    assert fit.tmid_bound_diagnostics['fills_window']
    assert 'remains sensitive' in fit.tmid_bound_diagnostics['note']
    assert 'quarter-period' in fit.tmid_bound_diagnostics['note']
    assert fit.parameters['tmid'] == -.25  # Retain the fitted estimate; never substitute the ephemeris.


def test_noisy_timing_matches_a_broad_independent_likelihood_after_expansion(monkeypatch):
    times = np.linspace(-.12, .12, 121)
    injected = .004
    errors = np.full(times.size, .006)
    flux = gaussian_transit(times, injected) + np.random.default_rng(17).normal(0, errors)
    calls = []
    fitter = grid_fitter(calls)
    monkeypatch.setattr(ex, 'lc_fitter', fitter)
    fit = ex.run_nested_lightcurve_fit_with_rprs_posterior_retry(
        times, flux, errors, np.ones(times.size), {'tmid': 0., 'per': 2.47},
        {'tmid': [-.000063, .000063]},
    )
    reference = fitter(times, flux, errors, np.ones(times.size), {'tmid': 0.}, {'tmid': [-.06, .06]})
    assert abs(fit.parameters['tmid'] - reference.parameters['tmid']) * 86400 < 5
    assert fit.errors['tmid'] == pytest.approx(reference.errors['tmid'], rel=.03)
    assert abs(fit.parameters['tmid'] - injected) < 3 * fit.errors['tmid']


def test_timing_expansion_preserves_fixed_geometry_and_its_reported_uncertainty(monkeypatch):
    times = np.linspace(-.12, .12, 121)
    calls = []
    fitter = grid_fitter(calls)
    fixed_errors_seen = []

    def fixed_geometry_fitter(*args, **kwargs):
        fixed_errors_seen.append(kwargs.get('fixed_parameter_errors', {}).copy())
        fit = fitter(*args, **kwargs)
        fit.rprs_prior_fallback_applied = True
        fit.rprs_prior_fallback_note = 'Existing fixed-geometry choice.'
        fit.errors['rprs'] = .007
        return fit

    monkeypatch.setattr(ex, 'lc_fitter', fixed_geometry_fitter)
    fit = ex.run_nested_lightcurve_fit_with_rprs_posterior_retry(
        times, gaussian_transit(times, -.004), np.full(times.size, .003), np.ones(times.size),
        {'tmid': 0., 'per': 2.47, 'rprs': .134}, {'tmid': [-.000063, .000063]},
    )
    assert fit.tmid_posterior_refit_applied
    assert all('rprs' not in bounds for bounds in calls)
    assert all(errors['rprs'] == .007 for errors in fixed_errors_seen[1:])
    assert fit.rprs_prior_fallback_note == 'Existing fixed-geometry choice.'


def test_lm_timing_covariance_expands_directly_from_seconds_to_minutes(monkeypatch):
    def fit(times, flux, errors, airmass, prior, bounds, **kwargs):
        return SimpleNamespace(parameters={**prior, 'tmid': .0005}, errors={'tmid': 5.5 / 1440})

    monkeypatch.setattr(ex, 'lc_fitter', fit)
    _, bounds, payload = ex.prepare_ultranest_bounds_with_lm_boundary_scout(
        np.linspace(-.1, .1, 101), np.ones(101), np.full(101, .01), np.ones(101),
        {'tmid': 0., 'per': 2.47}, {'tmid': [-.000063, .000063]}, enabled=True,
    )
    assert bounds['tmid'][0] <= .0005 - 4 * 5.5 / 1440
    assert bounds['tmid'][1] >= .0005 + 4 * 5.5 / 1440
    assert payload['adjusted'] and payload['iteration_count'] <= 2


def test_timing_search_floor_obeys_existing_period_cap():
    result = ex.estimate_ephemeris_tmid_and_bounds(
        np.linspace(-.5, .5, 11), 0., .1, 1e-8, 1e-9, expected_duration=.08,
    )
    assert result['half_width'] == pytest.approx(.025)
