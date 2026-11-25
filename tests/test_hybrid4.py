"""Model, shared-atmosphere and paired-input regression tests."""
import copy
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from exotic.api import elca, gael_ld, hybrid_comparison
from exotic.api.hybrid_ld import fit_hybrid4_profiles
from exotic.api.limbdarkfast import hybrid_flux, intensity
from exotic.transit_depth import observable_depth_percent


def prior():
    return dict(rprs=0.1, per=3.0, ars=10.0, inc=88.0, ecc=0.25, omega=120.0,
                tmid=2460000.0, u0=0.0, u1=0.0, u2=0.0, u3=0.4, a0=1.0, a2=0.0)


def atmosphere():
    mu = np.linspace(0.0, 1.0, 41)
    return fit_hybrid4_profiles(mu, 1 - 0.4 * (1 - mu ** 2), np.full(41, 0.01))[0]


def reference_fit():
    times = 2460000.0 + np.linspace(-0.1, 0.1, 50)
    flux = elca.transit(times, prior())
    fit = elca.lc_fitter(times, flux, np.full(50, 0.001), np.linspace(1.1, 1.6, 50),
                         prior(), {'rprs': [0.0, 0.5], 'tmid': [2459999.98, 2460000.02]},
                         mode=None, exposure_times_seconds=np.linspace(10.0, 120.0, 50),
                         fixed_flux_baseline=True, baseline_fit_mask=np.arange(50) % 2 == 0,
                         duration_prior={'applied': True, 'expected_duration': 0.1, 'sigma_log_duration': 0.1},
                         fixed_parameter_errors={'a0': 0.002}, ultranest_min_num_live_points=600)
    fit.mode = 'ns'
    fit.ns_type = 'ultranest'
    fit.max_ncalls = 123456
    return fit


def test_hybrid4_matches_independent_high_precision_reference_vectors():
    # Published 40-digit annulus-quadrature references, upstream revision
    # 120caecc; these values do not come from the model under test.
    weights = [0.3709047608728925, 0.00711438176464726, 0.3162293789763985, 0.01888876666845794]
    expected = [0.9875524096687214, 0.9875986911300181, 0.991591406808181, 1.0]
    np.testing.assert_allclose(hybrid_flux([0.0, 0.1, 0.9, 1.1], 0.1, weights), expected, rtol=0, atol=1e-12)


@pytest.mark.parametrize('ecc', [0.0, 0.25, 0.7])
@pytest.mark.parametrize('inc', [88.0, 84.5])
def test_model_dispatch_uses_same_orbit_and_front_back_convention(ecc, inc):
    p = {**prior(), 'ecc': ecc, 'inc': inc}
    times = p['tmid'] + np.linspace(-1.6, 1.6, 4001)
    alternative = {**p, '_limb_darkening_law': 'hybrid4', 'h0': 0.4, 'h1': 0, 'h2': 0, 'h3': 0}
    np.testing.assert_allclose(elca.transit(times, p), elca.transit(times, alternative),
                               rtol=0, atol=2e-7)  # PyLightcurve's existing precision=3


def test_hybrid4_uses_existing_variable_exposure_integration():
    p = {**prior(), '_limb_darkening_law': 'hybrid4', 'h0': 0.4, 'h1': 0.0, 'h2': 0.0, 'h3': 0.0}
    times = p['tmid'] + np.linspace(-0.06, 0.06, 50)
    exposures = np.linspace(30, 600, len(times)) / 86400.0
    integrated = elca.transit(times, {**p, '_exposure_time_days': exposures,
                                    '_exposure_smearing_change_tolerance': 0.0})
    offsets = (np.arange(7) + 0.5) / 7 - 0.5
    expected = elca._instantaneous_transit((times[:, None] + exposures[:, None] * offsets).ravel(), p)
    np.testing.assert_allclose(integrated, expected.reshape(-1, 7).mean(axis=1), atol=1e-14)
    assert np.max(np.abs(integrated - elca.transit(times, p))) > 1e-5


def test_depth_reporting_dispatches_to_hybrid4():
    p = {**prior(), '_limb_darkening_law': 'hybrid4', 'h0': 0.4, 'h1': 0.0, 'h2': 0.0, 'h3': 0.0}
    first = observable_depth_percent(p)[0]
    p['u0'], p['u1'], p['u2'], p['u3'] = 0.9, 0.0, 0.0, 0.0
    assert observable_depth_percent(p)[0] == first
    p['h0'] = 0.0
    assert observable_depth_percent(p)[0] != first


def test_profile_fit_allows_negative_weights_and_diagnostic_nonphysical_profile():
    mu = np.linspace(0, 1, 100)
    weights = np.array([0.5, -0.2, 0.3, 0.6])
    fitted = fit_hybrid4_profiles(mu, intensity(mu, weights, 'hybrid4'), np.full(100, 0.002))[0]
    np.testing.assert_allclose(fitted['coefficients'], weights, atol=1e-12)
    assert fitted['minimum_model_intensity'] < 0
    assert fitted['fit_mask'] == (mu > 0.1).tolist()
    assert np.asarray(fitted['coefficient_covariance']).shape == (4, 4)


def test_claret_and_hybrid4_use_one_atmosphere_generation(monkeypatch):
    mu = np.linspace(0, 1, 50)
    profile = 1 - 0.4 * (1 - mu ** 2)
    ps = SimpleNamespace(profile_mu=mu, profile_averages=profile[None, :],
                         profile_uncertainties=np.full((1, 50), 0.002))
    calls = []
    class Creator:
        def __init__(self, **kwargs):
            pass
        def create_profiles(self, nsamples):
            calls.append(nsamples)
            return ps
    def claret_fit(x, y, e, **kwargs):
        assert x is ps.profile_mu and y is ps.profile_averages and e is ps.profile_uncertainties
        return np.array([[0, 0, 0, 0.4]]), np.full((1, 4), 0.01)
    monkeypatch.setattr(gael_ld.ldtk, 'LDPSetCreator', Creator)
    monkeypatch.setattr(gael_ld, 'ldx', claret_fit)
    stellar = {'T*': 5000, 'T*_uperr': 50, 'T*_lowerr': -50, 'FEH*': 0,
               'FEH*_uperr': 0.1, 'FEH*_lowerr': -0.1, 'LOGG*': 4.5,
               'LOGG*_uperr': 0.1, 'LOGG*_lowerr': -0.1}
    result = gael_ld.createldgrid(np.array([0.55]), np.array([0.7]), stellar)
    assert calls == [10000]
    metadata = result['HYBRID4_PROFILES'][0]
    np.testing.assert_allclose(metadata['coefficients'], [0.4, 0, 0, 0], atol=1e-12)
    np.testing.assert_array_equal(metadata['intensity'], profile)
    assert metadata['claret_coefficients'] == [0, 0, 0, 0.4]
    assert metadata['wavelength_min_nm'] == 550


def test_secondary_fit_copies_inputs_and_settings_without_mutating_reference(monkeypatch):
    reference = reference_fit()
    snapshots = {key: copy.deepcopy(getattr(reference, key)) for key in
                 ('time', 'data', 'dataerr', 'airmass', 'prior', 'bounds')}
    seen = []
    monkeypatch.setattr(elca.lc_fitter, 'fit_nested', lambda fit: seen.append(fit))
    alternative = hybrid_comparison.run_hybrid4_fit(reference, atmosphere())
    assert seen == [alternative]
    for key in ('time', 'data', 'dataerr', 'airmass'):
        np.testing.assert_array_equal(getattr(reference, key), snapshots[key])
        np.testing.assert_array_equal(getattr(alternative, key), snapshots[key])
        assert not np.shares_memory(getattr(reference, key), getattr(alternative, key))
    assert reference.prior == snapshots['prior'] and reference.bounds == snapshots['bounds']
    assert alternative.bounds == reference.bounds
    for key in ('fixed_flux_baseline', 'duration_prior', 'ultranest_min_num_live_points',
                'max_ncalls', 'fixed_parameter_errors', 'exposure_smearing_supersample',
                'exposure_smearing_change_tolerance'):
        assert getattr(alternative, key) == getattr(reference, key)
    np.testing.assert_allclose(alternative.exposure_times_days, reference.exposure_times_days)
    np.testing.assert_array_equal(alternative.baseline_fit_mask, reference.baseline_fit_mask)
    assert alternative.prior['_limb_darkening_law'] == 'hybrid4'
    assert '_limb_darkening_law' not in reference.prior
    assert alternative.ultranest_warmstart_source is None


def test_secondary_preserves_fit_timer_separate_from_backend_setup(monkeypatch):
    def fit_nested(fit):
        fit.fit_runtime_seconds = 3.25
        fit.backend_setup_seconds = 2.0
    monkeypatch.setattr(elca.lc_fitter, 'fit_nested', fit_nested)
    fit = hybrid_comparison.run_hybrid4_fit(reference_fit(), atmosphere())
    assert fit.fit_runtime_seconds == 3.25
    assert fit.backend_setup_seconds == 2.0


def test_missing_profile_and_failed_comparison_are_explicit(tmp_path, monkeypatch):
    fit = reference_fit()
    result = hybrid_comparison.comparison_for_output(fit, {'save': str(tmp_path)}, 'Test b', '2026-11-20')
    assert result['status'] == 'unavailable'
    del fit.hybrid4_comparison
    def fail(*args):
        raise RuntimeError('demonstrated failure')
    monkeypatch.setattr(hybrid_comparison, 'run_hybrid4_fit', fail)
    result = hybrid_comparison.comparison_for_output(fit, {'save': str(tmp_path),
                'hybrid4_limb_darkening': atmosphere()}, 'Test b', '2026-11-20')
    assert result['status'] == 'failed' and result['reason'] == 'demonstrated failure'
    assert json.loads(next(tmp_path.glob('*.json')).read_text())['status'] == 'failed'


def test_sparse_atmosphere_profile_is_fitted_and_rank_is_reported():
    mu = np.array([0.3, 0.6, 1.0])
    result = fit_hybrid4_profiles(mu, 1 - 0.3 * (1 - mu ** 2), np.full(3, 0.01))[0]
    assert result['profile_design_rank'] < 4
    assert len(result['coefficients']) == 4


def test_real_secondary_atmosphere_failure_does_not_discard_claret(monkeypatch):
    from exotic.api import hybrid_ld
    mu = np.linspace(0, 1, 50)
    ps = SimpleNamespace(profile_mu=mu, profile_averages=np.ones((1, 50)),
                         profile_uncertainties=np.full((1, 50), 0.01))
    class Creator:
        def __init__(self, **kwargs):
            pass
        def create_profiles(self, **kwargs):
            return ps
    def fail(*args, **kwargs):
        raise np.linalg.LinAlgError('actual matrix failure')
    monkeypatch.setattr(gael_ld.ldtk, 'LDPSetCreator', Creator)
    monkeypatch.setattr(gael_ld, 'ldx', lambda *args, **kwargs: (np.ones((1, 4)), np.ones((1, 4))))
    monkeypatch.setattr(hybrid_ld, 'fit_hybrid4_profiles', fail)
    stellar = {'T*': 5000, 'T*_uperr': 50, 'T*_lowerr': -50, 'FEH*': 0,
               'FEH*_uperr': 0.1, 'FEH*_lowerr': -0.1, 'LOGG*': 4.5,
               'LOGG*_uperr': 0.1, 'LOGG*_lowerr': -0.1}
    result = gael_ld.createldgrid(np.array([0.55]), np.array([0.7]), stellar)
    np.testing.assert_array_equal(result['LD'], np.ones((4, 1)))
    assert result['HYBRID4_PROFILES'][0]['status'] == 'failed'
    assert result['HYBRID4_PROFILES'][0]['reason'] == 'actual matrix failure'


def test_reports_contain_both_models_and_the_same_photometry(tmp_path):
    reference = reference_fit()
    alternate = copy.deepcopy(reference)
    alternate.prior.update({'_limb_darkening_law': 'hybrid4', **dict(zip(('h0','h1','h2','h3'), [0.4,0,0,0]))})
    for fit in (reference, alternate):
        fit.parameters = copy.deepcopy(fit.prior)
        fit.errors = {'rprs': 0.001, 'tmid': 0.0001}
        fit.quantiles = {}
        fit.posterior_summaries = {'rprs': {'median': 0.1, 'lower': 0.099, 'upper': 0.101, 'stdev': 0.001}}
        fit.create_fit_variables()
    report = hybrid_comparison.save_comparison(reference, alternate, atmosphere(), tmp_path, 'comparison')
    assert report['claret']['law'] == 'claret' and report['hybrid4']['law'] == 'hybrid4'
    with (tmp_path/'comparison.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    np.testing.assert_array_equal([float(row['BJD_TDB']) for row in rows], reference.time)
    np.testing.assert_array_equal([float(row['Shared Flux']) for row in rows], reference.data)
    np.testing.assert_array_equal([float(row['Shared Uncertainty']) for row in rows], reference.dataerr)
    assert (tmp_path/'comparison.png').is_file()
    assert json.loads((tmp_path/'comparison.json').read_text())['status'] == 'completed'
    # The bulk runner removes its temporary workspace after copying outputs.
    import shutil
    destination = tmp_path.parent / (tmp_path.name + '_archived')
    shutil.move(str(tmp_path), destination)
    archived = json.loads((destination/'comparison.json').read_text())
    for relative_path in archived['files'].values():
        assert not Path(relative_path).is_absolute()
        assert (destination / relative_path).is_file()
