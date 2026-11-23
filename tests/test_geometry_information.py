"""Pre-UltraNest geometry decisions use information, not a point-count floor."""
from types import SimpleNamespace

import numpy as np
import pytest

from exotic.api.elca import lc_fitter, transit_duration


def model_at_solution(*, noise=0.001, inc=87.0, times=None, exposure_seconds=None):
    fit = object.__new__(lc_fitter)
    fit.mode = 'lm'
    fit.parameters = {
        'rprs': 0.1, 'ars': 10.0, 'inc': inc, 'per': 3.0, 'tmid': 0.0,
        'ecc': 0.0, 'omega': 90.0, 'u0': 0.0, 'u1': 0.0, 'u2': 0.0, 'u3': 0.0,
        'a0': 1.0, 'a2': 0.0,
    }
    fit.prior = dict(fit.parameters)
    fit.bounds = {'rprs': [0.05, 0.15], 'ars': [7, 13], 'inc': [80, 90], 'tmid': [-0.01, 0.01]}
    fit.time = np.linspace(-0.09, 0.09, 600) if times is None else np.asarray(times)
    fit.airmass = np.ones(fit.time.shape)
    fit.airmass_reference = 1.0
    fit.baseline_fit_mask = None
    fit.fixed_flux_baseline = False
    fit.exposure_times_days = None
    if exposure_seconds is not None:
        from exotic.api.elca import DEFAULT_EXPOSURE_SMEARING_SUPERSAMPLE, DEFAULT_EXPOSURE_SMEARING_CHANGE_TOLERANCE
        fit.exposure_times_days = np.full(fit.time.shape, exposure_seconds / 86400)
        fit.exposure_smearing_supersample = DEFAULT_EXPOSURE_SMEARING_SUPERSAMPLE
        fit.exposure_smearing_change_tolerance = DEFAULT_EXPOSURE_SMEARING_CHANGE_TOLERANCE
    fit.dataerr = np.full(fit.time.shape, noise)
    fit.data = fit._transit_model(fit.time, fit.parameters)
    return fit


def test_precision_tracks_noise_even_with_identical_point_counts_and_zero_residuals():
    precise = model_at_solution(noise=0.0002).estimate_geometry_information()
    noisy = model_at_solution(noise=0.004).estimate_geometry_information()
    assert precise['available'] and noisy['available']
    assert precise['ingress_points'] == noisy['ingress_points']
    assert 0 < precise['ingress_fractional_uncertainty'] < 0.1
    assert noisy['ingress_fractional_uncertainty'] > 0.1
    assert noisy['ingress_fractional_uncertainty'] / precise['ingress_fractional_uncertainty'] == pytest.approx(20)
    assert precise['noise_scale'] == 1.0


def test_central_transit_has_finite_information_in_b_squared():
    assessment = model_at_solution(noise=0.0002, inc=90).estimate_geometry_information()
    assert assessment['available']
    assert 0 < assessment['b_squared_uncertainty'] < 1
    assert 0 < assessment['ingress_fractional_uncertainty'] < 0.1


def test_two_precise_points_in_each_slope_can_remain_permitted():
    reference = model_at_solution()
    outer = transit_duration(reference.parameters)
    inner = transit_duration(reference.parameters, 1 - reference.parameters['rprs'])
    ingress_times = -0.5 * outer + (outer - inner) * np.array([0.15, 0.35])
    times = np.sort(np.concatenate([
        np.linspace(-0.09, -0.065, 12), ingress_times,
        np.linspace(-0.02, 0.02, 15), -ingress_times, np.linspace(0.065, 0.09, 12),
    ]))
    assessment = model_at_solution(noise=1e-6, times=times).estimate_geometry_information()
    assert assessment['ingress_points'] == assessment['egress_points'] == 2
    assert 0 < assessment['ingress_fractional_uncertainty'] < 0.1


def test_flat_bottom_only_is_unconstrained_instead_of_false_zero_covariance():
    assessment = model_at_solution(times=np.linspace(-0.015, 0.015, 200)).estimate_geometry_information()
    assert assessment['available']
    assert assessment['ingress_points'] == assessment['egress_points'] == 0
    assert assessment['b_squared_uncertainty'] is None or assessment['b_squared_uncertainty'] > 1
    assert assessment['ingress_fractional_uncertainty'] is None or assessment['ingress_fractional_uncertainty'] > 0.1


def test_singular_geometry_never_gets_false_zero_errors():
    fit = model_at_solution(times=np.linspace(-0.015, 0.015, 200))
    fit._transit_model = lambda times, values: np.full(np.shape(times), 1 - values['rprs'] ** 2)
    assessment = fit.estimate_geometry_information()
    assert assessment['available']
    assert assessment['b_squared_uncertainty'] is None
    assert assessment['ingress_fractional_uncertainty'] is None


def test_duration_prior_does_not_masquerade_as_data_information():
    fit = model_at_solution()
    before = fit.estimate_geometry_information()
    fit.duration_prior = {'applied': True, 'expected_duration': 0.1, 'sigma_log_duration': 1e-9}
    assert fit.estimate_geometry_information() == before


def test_long_integrations_reduce_ingress_information():
    short = model_at_solution(exposure_seconds=15).estimate_geometry_information()
    long = model_at_solution(exposure_seconds=900).estimate_geometry_information()
    assert long['ingress_fractional_uncertainty'] > short['ingress_fractional_uncertainty']


def test_free_ars_uncertainty_is_included():
    fit = model_at_solution()
    joint = fit.estimate_geometry_information()
    fit.bounds.pop('ars')
    conditional = fit.estimate_geometry_information()
    assert joint['ingress_fractional_uncertainty'] > conditional['ingress_fractional_uncertainty']


def test_estimate_from_an_actual_least_squares_fit():
    source = model_at_solution(noise=0.0003)
    rng = np.random.default_rng(2851)
    fit = lc_fitter(
        source.time, source.data + rng.normal(0, source.dataerr), source.dataerr, source.airmass,
        {**source.parameters, 'inc': 87.1, 'ars': 10.1}, source.bounds, mode='lm',
    )
    assessment = fit.estimate_geometry_information()
    assert assessment['available']
    assert 0 < assessment['ingress_fractional_uncertainty'] < 0.1
    assert 0 < assessment['b_squared_uncertainty'] < 1


@pytest.mark.parametrize('fractional_error, keep_free', [(0.05, True), (0.1, True), (0.3, False), (None, False)])
@pytest.mark.parametrize('forced', [False, True])
@pytest.mark.parametrize('sample_impact_parameter', [False, True])
def test_scout_decides_before_first_ultranest_and_restores_input_geometry(
    monkeypatch, fractional_error, keep_free, forced, sample_impact_parameter,
):
    import exotic.exotic as ex
    calls = []

    def fake_fitter(times, flux, unc, airmass, prior, bounds, mode=None, **kwargs):
        calls.append((mode, dict(prior), dict(bounds), kwargs))
        fit = SimpleNamespace(parameters={**prior, 'inc': 87.2} if mode == 'lm' else dict(prior), errors={})
        fit.estimate_geometry_information = lambda: {
            'available': True, 'ingress_fractional_uncertainty': fractional_error,
            'b_squared_uncertainty': 0.1 if fractional_error is not None else None,
            'ingress_points': 12, 'egress_points': 12,
        }
        fit.get_parameter_posterior_recenter_diagnostics = lambda key: None
        return fit

    monkeypatch.setattr(ex, 'lc_fitter', fake_fitter)
    prior = dict(model_at_solution().parameters, inc=88.0, inc_unc=0.2)
    bounds = {'rprs': [0.05, 0.15], 'ars': [7, 13], 'inc': [80, 90], 'tmid': [-0.01, 0.01]}
    fit = ex.run_nested_lightcurve_fit_with_rprs_posterior_retry(
        np.linspace(-0.09, 0.09, 100), np.ones(100), np.full(100, 0.001), np.ones(100),
        prior, bounds, use_lm_boundary_scout=True,
        force_fit_inclination=forced,
        use_impactparameter_rather_than_inclination_to_fit=sample_impact_parameter,
        max_rprs_retries=0, max_ars_retries=0, max_impact_parameter_retries=0,
    )
    assert calls[0][0] == 'lm' and calls[-1][0] == 'ns'
    assert ('inc' in calls[-1][2]) is (keep_free or forced)
    assert calls[-1][3]['use_impactparameter_rather_than_inclination_to_fit'] is sample_impact_parameter
    assert set(calls[-1][2]) >= {'rprs', 'ars', 'tmid'}
    assert ('inc' in bounds) and prior['inc'] == 88.0  # caller input is untouched
    if keep_free or forced:
        assert calls[-1][1]['inc'] == 87.2
        assert 'inc' not in calls[-1][3].get('fixed_parameter_errors', {})
    else:
        assert calls[-1][1]['inc'] == 88.0
        assert calls[-1][3]['fixed_parameter_errors']['inc'] == 0.2
        assert fit.pre_ultranest_geometry_information['conditional_on_supplied_geometry']
    assert fit.pre_ultranest_geometry_information['meets_target'] is keep_free
    assert fit.pre_ultranest_geometry_information['force_fit_inclination'] is forced
    assert fit.force_fit_inclination is forced
    # Later baseline/full-resolution fits can reassess the original free
    # geometry, rather than inheriting a decision from a coarse input series.
    assert 'inc' in ex.get_posterior_refit_final_bounds(fit, bounds)


@pytest.mark.parametrize('pre_points,post_points', [(5, 0), (0, 5), (0, 0)])
def test_force_keeps_geometry_free_with_partial_coverage_and_poor_precision(monkeypatch, pre_points, post_points):
    import exotic.exotic as ex
    calls = []

    def fake_fitter(times, flux, unc, airmass, prior, bounds, mode=None, **kwargs):
        calls.append((mode, dict(bounds), kwargs))
        fit = SimpleNamespace(parameters=dict(prior), errors={},
                              get_parameter_posterior_recenter_diagnostics=lambda key: None)
        fit.estimate_geometry_information = lambda: {
            'available': True, 'ingress_fractional_uncertainty': 0.8,
            'b_squared_uncertainty': 1.0, 'ingress_points': 1, 'egress_points': 0,
        }
        return fit

    monkeypatch.setattr(ex, 'lc_fitter', fake_fitter)
    source = model_at_solution()
    fit = ex.run_nested_lightcurve_fit_with_rprs_posterior_retry(
        source.time, source.data, np.full(source.time.shape, 0.05), source.airmass,
        source.parameters, source.bounds, force_fit_inclination=True, use_lm_boundary_scout=True,
        pre_ultranest_coverage_assessment={
            'valid': True, 'transit_fraction_observed': 0.2, 'in_transit_points': 2,
            'pre_ingress_points': pre_points, 'post_egress_points': post_points,
        },
    )
    assert calls[0][0] == 'lm' and calls[-1][0] == 'ns'
    assert all('inc' in bounds for _, bounds, _ in calls)
    assert all('inc' not in kwargs.get('fixed_parameter_errors', {}) for _, _, kwargs in calls)
    assert 'inc' not in fit.partial_transit_geometry_prior_assumption_fixed_parameters
    assert 'inc' in fit.partial_transit_geometry_prior_assumption_sampled_parameters
    assert fit.pre_ultranest_geometry_information['meets_target'] is False
    assert 'force_fit_inclination=true' in fit.partial_transit_geometry_prior_assumption_note
    assert 'remains free' in fit.pre_ultranest_geometry_information['note']


def test_force_restores_geometry_without_a_scout_or_existing_free_bound(monkeypatch):
    import exotic.exotic as ex
    calls = []

    def fake_fitter(*args, **kwargs):
        calls.append((args[5], kwargs))
        return SimpleNamespace(parameters=dict(args[4]), get_parameter_posterior_recenter_diagnostics=lambda key: None)

    monkeypatch.setattr(ex, 'lc_fitter', fake_fitter)
    source = model_at_solution()
    fit = ex.run_nested_lightcurve_fit_with_rprs_posterior_retry(
        source.time, source.data, source.dataerr, source.airmass, source.parameters,
        {'tmid': [-0.01, 0.01]}, fixed_parameter_errors={'inc': 0.2},
        force_fit_inclination=True, use_lm_boundary_scout=False,
    )
    assert len(calls) == 1
    assert 'inc' in calls[0][0]
    assert 'inc' not in calls[0][1].get('fixed_parameter_errors', {})
    assert fit.force_fit_inclination is True


def test_force_configuration_resets_and_reaches_workers(monkeypatch):
    import exotic.exotic as ex
    monkeypatch.setattr(ex, 'FORCE_FIT_INCLINATION', False)
    cleanup = []
    monkeypatch.setattr(ex, 'suppress_inherited_tk_cleanup_in_worker', lambda: cleanup.append(True))
    assert ex.configure_force_fit_inclination(True) is True
    assert ex.FORCE_FIT_INCLINATION is True
    assert ex.configure_force_fit_inclination(None) is False
    ex.initialize_lightcurve_candidate_worker(True)
    assert ex.FORCE_FIT_INCLINATION is True
    source = model_at_solution()
    _, bounds, _, _ = ex.apply_partial_transit_geometry_prior_assumption(
        source.parameters, source.bounds,
        {'valid': True, 'transit_fraction_observed': 0.2, 'in_transit_points': 2,
         'pre_ingress_points': 0, 'post_egress_points': 4},
    )
    assert 'inc' in bounds  # worker/default callers inherit the configured choice
    ex.initialize_lightcurve_candidate_worker(False)
    assert ex.FORCE_FIT_INCLINATION is False
    assert cleanup == [True, True]


def test_forced_partial_geometry_remains_free_during_posterior_retry(monkeypatch):
    import exotic.exotic as ex
    calls = []

    def fake_fitter(times, flux, unc, airmass, prior, bounds, **kwargs):
        calls.append((dict(bounds), kwargs))
        first = len(calls) == 1
        fit = SimpleNamespace(parameters=dict(prior),
                              sampled_keys=['b' if key == 'inc' else key for key in bounds])
        fit.get_parameter_posterior_recenter_diagnostics = lambda key: {
            'clipped': first and key == 'b', 'edge': 'upper', 'mode': 0.9,
            'std': 0.1, 'bounds': [0.0, 1.05], 'reason': 'interior' if not first else 'edge',
        }
        return fit

    monkeypatch.setattr(ex, 'lc_fitter', fake_fitter)
    source = model_at_solution()
    bounds = {**source.bounds, 'inc': [85.0, 90.0]}
    fit = ex.run_nested_lightcurve_fit_with_rprs_posterior_retry(
        source.time, source.data, source.dataerr, source.airmass, source.parameters, bounds,
        force_fit_inclination=True, use_lm_boundary_scout=False, max_impact_parameter_retries=1,
        pre_ultranest_coverage_assessment={
            'valid': True, 'success_label': 'low', 'transit_fraction_observed': 0.2,
            'in_transit_points': 2, 'pre_ingress_points': 0, 'post_egress_points': 4,
        },
    )
    assert len(calls) == 2
    assert all('inc' in call_bounds for call_bounds, _ in calls)
    assert all('inc' not in kwargs.get('fixed_parameter_errors', {}) for _, kwargs in calls)
    assert calls[-1][0]['inc'][0] < calls[0][0]['inc'][0]
    assert fit.b_posterior_refit_applied is True
    assert fit.force_fit_inclination is True


def test_failed_estimate_is_reported_without_suppressing_free_geometry(monkeypatch):
    import exotic.exotic as ex

    def failed_estimate():
        raise RuntimeError('model derivative failed')

    monkeypatch.setattr(ex, 'lc_fitter', lambda *args, **kwargs: SimpleNamespace(
        parameters=dict(args[4]), errors={}, estimate_geometry_information=failed_estimate,
    ))
    source = model_at_solution()
    prior, bounds, payload = ex.prepare_ultranest_bounds_with_lm_boundary_scout(
        source.time, source.data, source.dataerr, source.airmass, source.parameters,
        source.bounds, enabled=True,
    )
    assert 'inc' in bounds
    assert payload['geometry_information']['available'] is False
    assert 'model derivative failed' in payload['geometry_information']['note']
    assert 'remains free' in payload['geometry_information']['note']
