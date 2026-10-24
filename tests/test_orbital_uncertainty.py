from copy import deepcopy
from types import SimpleNamespace

import numpy as np
import pandas
import pytest

from exotic.api.nea import NASAExoplanetArchive
from exotic.orbital_uncertainty import (
    archive_orbital_records, nominal_eccentricity, orbital_factor_bounds,
    orbital_uncertainty_budget,
)
from exotic.output_files import posterior_reporting_metadata


def upper_limit(value=.025, probability=.9545):
    return {'source': 'published RV analysis', 'omega_convention': 'unspecified',
            'ecc': {'value': value, 'kind': 'upper_limit', 'confidence_probability': probability},
            'omega': {'value': None}}


def saved_fit():
    return SimpleNamespace(
        parameters={'ars': 10., 'b': .4, 'ecc': .018, 'omega': 97., 'inc': 87.5},
        errors={'inc': 1.2}, posterior_summaries={'inc': {'median': 87.5}},
        sampled_keys=['b', 'ars'],
        results={'weighted_samples': {'points': np.array([[.2, 8.], [.4, 10.], [.9, 14.]]),
                                      'weights': np.array([.1, .8, .1])}},
    )


@pytest.mark.parametrize('ecc', [0., .018, .025, .9])
def test_upper_limit_extrema_are_exact_and_include_every_orientation(ecc):
    factors = orbital_factor_bounds([0, ecc])
    assert factors['separation_factor'] == pytest.approx([1 - ecc, 1 + ecc])
    g = np.sqrt((1 + ecc) / (1 - ecc))
    assert factors['speed_factor'] == pytest.approx([1 / g, g])


@pytest.mark.parametrize('angle_range', [[350, 390], [60, 110], [200, 290], [-30, 30]])
def test_measured_region_extrema_agree_with_independent_dense_grid(angle_range):
    e = np.linspace(.01, .35, 1100)[:, None]
    sine = np.sin(np.radians(np.linspace(*angle_range, 1500)))[None, :]
    f = (1 - e ** 2) / (1 + e * sine)
    g = (1 + e * sine) / np.sqrt(1 - e ** 2)
    bounds = orbital_factor_bounds([.01, .35], angle_range)
    assert bounds['separation_factor'] == pytest.approx([f.min(), f.max()], abs=2e-7)
    assert bounds['speed_factor'] == pytest.approx([g.min(), g.max()], abs=2e-7)


def test_upper_limit_uses_circular_reference_without_changing_saved_fit():
    fit = saved_fit()
    original_params, original_errors = deepcopy(fit.parameters), deepcopy(fit.errors)
    original_points = fit.results['weighted_samples']['points'].copy()
    budget = orbital_uncertainty_budget(fit, [upper_limit()])
    entry = budget['constraints'][0]
    assert entry['nominal_orbit'] == 'circular_assumed'
    assert entry['nominal_eccentricity'] == 0
    assert entry['best_eccentricity_from_constraint'] is None
    assert entry['quoted_confidence_probability'] == .9545
    assert entry['saved_fit_eccentricity'] == .018
    assert entry['sensitivity_reference_eccentricity'] == 0
    assert entry['nominal_argument_of_periastron_degrees'] is None
    assert entry['maximum_absolute_fractional_departure_approximate']['duration_preserving_stellar_density'] == pytest.approx((1.025 / .975) ** 1.5 - 1)
    assert entry['duration_preserving_ars_factor_approximate'] == pytest.approx(
        [np.sqrt(.975 / 1.025), np.sqrt(1.025 / .975)])
    nominal = np.degrees(np.arccos(.4 / 10))
    worst = nominal - np.degrees(np.arccos(.4 / (10 * .975)))
    assert entry['model_geometry_maximum_inclination_departure_degrees'] == pytest.approx(worst)
    assert entry['fully_transformable_posterior_weight'] == 1
    assert fit.parameters == original_params
    assert fit.errors == original_errors
    np.testing.assert_array_equal(fit.results['weighted_samples']['points'], original_points)
    metadata = posterior_reporting_metadata(fit, budget)
    assert metadata['external_orbital_uncertainty'] == budget
    assert metadata['model_values'] == original_params


def test_joint_transit_correlations_and_weights_are_retained():
    fit = saved_fit()
    base = orbital_uncertainty_budget(fit, [upper_limit()])['constraints'][0]
    # Preserve both marginals but change their pairing. Geometry must change.
    fit.results['weighted_samples']['points'][:, 1] = [14., 10., 8.]
    swapped = orbital_uncertainty_budget(fit, [upper_limit()])['constraints'][0]
    assert base['circular_reference_inclination_degrees'] != swapped['circular_reference_inclination_degrees']
    fit.results['weighted_samples']['weights'] = [.99, .005, .005]
    weighted = orbital_uncertainty_budget(fit, [upper_limit()])['constraints'][0]
    assert weighted['circular_reference_inclination_degrees']['median'] > 89


def test_missing_limit_confidence_and_missing_errors_remain_unknown():
    record = upper_limit(probability=None)
    fixed = {'ecc': {'value': 0, 'kind': 'measurement'}, 'omega': {'value': 97}}
    budget = orbital_uncertainty_budget(saved_fit(), [record, fixed])
    assert budget['constraints'][0]['quoted_confidence_probability'] is None
    assert budget['constraints'][1]['status'] == 'uncertainty_unspecified'
    assert 'no errors' in budget['constraints'][1]['reason']
    assert len(budget['constraints']) == 2  # alternatives, never combined


def test_convention_is_explicit_and_star_omega_is_shifted_by_180():
    record = {'ecc': {'value': .1, 'kind': 'measurement', 'error_plus': .02, 'error_minus': -.01},
              'omega': {'value': 97, 'kind': 'measurement', 'error_plus': 44, 'error_minus': -20},
              'omega_convention': 'star'}
    entry = orbital_uncertainty_budget(saved_fit(), [record])['constraints'][0]
    assert entry['omega_ranges_planet_degrees'] == [[257, 321]]
    record['omega_convention'] = 'unspecified'
    unknown = orbital_uncertainty_budget(saved_fit(), [record])['constraints'][0]
    assert unknown['omega_ranges_planet_degrees'] == [[77, 141], [257, 321]]
    assert unknown['separation_factor'][0] < entry['separation_factor'][0]


def test_unweighted_dead_points_are_not_used_as_a_posterior():
    fit = saved_fit()
    del fit.results['weighted_samples']['weights']
    budget = orbital_uncertainty_budget(fit, [upper_limit()])
    assert budget['posterior_sample_source'] == 'unavailable'
    assert 'inclination_endpoint_summaries_degrees' not in budget['constraints'][0]
    fit.results['samples'] = np.array([[.4, 10.]])
    assert orbital_uncertainty_budget(fit, [upper_limit()])['posterior_sample_source'] == 'equal_weight_posterior_samples'


def test_archive_preserves_each_sources_errors_and_censoring(monkeypatch, tmp_path):
    nea = NASAExoplanetArchive('WASP-16 b')
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(nea, 'planet_names', lambda **kwargs: None)
    rows = pandas.DataFrame([
        {'pl_name': 'WASP-16 b', 'pl_pubdate': '2017', 'pl_refname': 'Bonomo',
         'pl_orbeccen': .018, 'pl_orbeccenerr1': np.nan, 'pl_orbeccenerr2': np.nan,
         'pl_orbeccenlim': 1, 'pl_orblper': np.nan, 'pl_orblpererr1': np.nan,
         'pl_orblpererr2': np.nan, 'pl_orblperlim': np.nan},
        {'pl_name': 'WASP-16 b', 'pl_pubdate': '2014', 'pl_refname': 'Knutson',
         'pl_orbeccen': .015, 'pl_orbeccenerr1': .012, 'pl_orbeccenerr2': -.011,
         'pl_orbeccenlim': 0, 'pl_orblper': 97, 'pl_orblpererr1': 44,
         'pl_orblpererr2': -20, 'pl_orblperlim': 0},
    ])
    queries = []
    def query(base, fields):
        queries.append(dict(fields))
        return rows.iloc[:1].copy() if 'default_flag' in fields['where'] else rows.copy()
    monkeypatch.setattr(nea, '_tap_query', query)
    nea._new_scrape()
    import json
    result = json.loads((tmp_path / 'eaConf.json').read_text())[0]
    records = result['orbital_constraints']
    assert records[0]['ecc']['kind'] == 'upper_limit'
    assert records[0]['ecc']['error_plus'] is None
    assert records[0]['omega']['value'] is None
    assert records[1]['omega']['error_plus'] == 44
    assert records[1]['ecc']['error_plus'] == .012
    assert result['pl_orbeccenerr1'] is None  # never borrowed from 2014
    assert 'pl_orbeccenlim' in queries[0]['select']
    assert nominal_eccentricity(result['pl_orbeccen'], records) == 0


def test_nominal_circular_assumption_does_not_replace_a_measured_eccentricity():
    record = archive_orbital_records([{'pl_orbeccen': .3, 'pl_orbeccenlim': 0}])
    assert nominal_eccentricity(.3, record) == .3
    assert nominal_eccentricity(.018, [upper_limit()]) == 0


def test_partial_geometry_compatibility_is_retained_and_reported():
    fit = saved_fit()
    fit.results['weighted_samples']['points'] = np.array([[.99, 1.], [.4, 10.], [1.04, 1.]])
    entry = orbital_uncertainty_budget(fit, [upper_limit()])['constraints'][0]
    assert entry['fully_transformable_posterior_weight'] == pytest.approx(.8)
    assert entry['partly_transformable_posterior_weight'] == pytest.approx(.1)
    assert entry['untransformable_posterior_weight'] == pytest.approx(.1)
    assert entry['inclination_endpoint_summaries_degrees']['minimum'] is not None


def test_archive_parameter_model_uses_zero_but_retains_published_upper_limit():
    nea = NASAExoplanetArchive('WASP-16 b')
    nea._get_params({'pl_name': 'WASP-16 b', 'hostname': 'WASP-16',
                     'pl_orbper': 3.1, 'pl_ratdor': 10., 'pl_orbincl': 88.,
                     'pl_ratror': .1, 'pl_ratrorerr1': .01, 'pl_ratrorerr2': -.01,
                     'pl_orbeccen': .025, 'pl_orbeccenlim': 1, 'pl_orblper': 97.})
    assert nea.pl_dict['ecc'] == 0
    assert nea.orbital_constraints[0]['ecc']['value'] == .025
    assert nea.orbital_constraints[0]['ecc']['confidence_probability'] is None


def timing_sensitivity_fit():
    """A known location likelihood, so orbit response has an exact answer."""
    times = np.linspace(-.65, .65, 501)
    weights = np.exp(-.5 * (times / .1) ** 2)
    parameters = {'tmid': 0., 'ecc': 0., 'omega': 0., 'a0': 1.}
    fit = SimpleNamespace(
        parameters=parameters, errors={}, sampled_keys=['tmid'],
        results={'weighted_samples': {'points': times[:, None], 'weights': weights}},
        _transit_model=lambda t, p: np.full(len(t), p['tmid'] + .2 * p['ecc'] * np.cos(np.radians(p['omega']))),
        orbital_likelihood_context={
            'base_physical': parameters.copy(), 'time': np.array([0.]), 'data': np.array([0.]),
            'inverse_dataerr': np.array([10.]), 'centered_airmass': np.array([0.]),
            'uses_fixed_flux_baseline': True, 'has_free_flux_baseline': False,
            'fixed_flux_baseline_value': 1., 'free_flux_baseline_key': None,
            'duration_prior_valid': False,
        },
    )
    return fit


def test_non_refit_timing_sensitivity_matches_exact_likelihood_and_preserves_fit():
    fit = timing_sensitivity_fit()
    before = deepcopy(fit.parameters)
    weights = fit.results['weighted_samples']['weights'].copy()
    sensitivity = orbital_uncertainty_budget(fit, [upper_limit()])['constraints'][0]['posterior_sensitivity']
    detail = sensitivity['parameters']['tmid']
    assert detail['orbital_assumption_mean_shift_range'] == pytest.approx([-432., 432.], abs=.001)
    assert detail['random_conditional_sigma_at_reference'] == pytest.approx(8640., abs=.001)
    assert detail['maximum_mean_shift_change_on_grid_refinement'] == pytest.approx(0, abs=.001)
    assert sensitivity['grid']['probability_weights_assigned'] is False
    assert sensitivity['orbital_random_uncertainty']['sigma'] is None
    assert fit.parameters == before
    np.testing.assert_array_equal(fit.results['weighted_samples']['weights'], weights)


def test_joint_external_draws_partition_random_variance_without_fitting_orbit():
    record = upper_limit()
    record['joint_posterior'] = {'ecc': [.025, .025], 'omega_planet_degrees': [0, 180], 'weights': [1, 1]}
    detail = orbital_uncertainty_budget(timing_sensitivity_fit(), [record])['constraints'][0]['posterior_sensitivity']['parameters']['tmid']
    assert detail['random_transit_sigma_mixture'] == pytest.approx(8640., abs=.001)
    assert detail['random_external_orbit_sigma'] == pytest.approx(432., abs=.001)
    assert detail['total_random_sigma'] == pytest.approx(np.hypot(8640, 432), abs=.001)
    assert detail['systematic_reference_to_mixture_mean_shift'] == pytest.approx(0, abs=.001)


def test_no_joint_distribution_is_invented_from_measured_marginal_errors():
    record = {'ecc': {'value': .02, 'kind': 'measurement', 'error_minus': -.01, 'error_plus': .01},
              'omega': {'value': 0, 'kind': 'measurement', 'error_minus': -20, 'error_plus': 30},
              'omega_convention': 'planet'}
    sensitivity = orbital_uncertainty_budget(timing_sensitivity_fit(), [record])['constraints'][0]['posterior_sensitivity']
    assert sensitivity['orbital_random_uncertainty']['sigma'] is None
    assert sensitivity['parameters']['tmid']['maximum_absolute_explored_mean_shift'] > 0


def test_joint_orbital_samples_need_no_marginal_error_prerequisite():
    record = {'source': 'full RV posterior',
              'joint_posterior': {'ecc': [.025, .025], 'omega_planet_degrees': [0., 180.]}}
    entry = orbital_uncertainty_budget(timing_sensitivity_fit(), [record])['constraints'][0]
    assert entry['posterior_sensitivity']['parameters']['tmid']['random_external_orbit_sigma'] == pytest.approx(432, abs=.001)
    assert entry['region_kind'] == 'conservative_marginal_draw_box_not_a_joint_confidence_region'


def test_zero_weight_orbital_samples_do_not_restrict_the_analysis():
    record = {'joint_posterior': {'ecc': [.025, .025, np.nan],
                                 'omega_planet_degrees': [0., 180., np.nan], 'weights': [1, 1, 0]}}
    sensitivity = orbital_uncertainty_budget(timing_sensitivity_fit(), [record])['constraints'][0]['posterior_sensitivity']
    assert len(sensitivity['scenarios']) == 2
    assert sensitivity['parameters']['tmid']['random_external_orbit_sigma'] == pytest.approx(432, abs=.001)
