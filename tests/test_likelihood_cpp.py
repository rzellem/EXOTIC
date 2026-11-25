"""Native live scores must match the actual reference callback, without sampling."""
import copy

import numpy as np
import pytest

from exotic.api import elca


class Captured(BaseException):
    pass


def capture(monkeypatch, fit, backend):
    result = {}
    def sampler(keys, loglike, transform, **kwargs):
        result.update(keys=keys, loglike=loglike, transform=transform)
        raise Captured()
    monkeypatch.setattr(elca, 'ReactiveNestedSampler', sampler)
    fit.mode = 'ns'
    fit.likelihood_backend_requested = backend
    with pytest.raises(Captured):
        fit.fit_nested()
    return result


def fixture(law='claret', geometry='impact', baseline='profile', duration=True, count=51):
    time = 2460000. + np.linspace(-.13, .13, count)
    prior = dict(rprs=.12, per=3.12, ars=8.21, ecc=.018, inc=84., omega=97.,
                 tmid=2460000., a2=.01, a0=1.008,
                 u0=1.77248, u1=-2.96806, u2=3.44175, u3=-1.24617,
                 h0=.4, h1=.3, h2=.1, h3=.05,
                 _limb_darkening_law=law, _claret_backend='cpp', _hybrid4_backend='cpp')
    bounds = dict(rprs=[.03, .2], ars=[3., 20.], tmid=[2460000. - .13, 2460000. + .13], a2=[-.1, .1])
    if geometry in ('impact', 'inclination'):
        bounds['inc'] = [75., 90.]
    elif geometry == 'direct_b':
        bounds['b'] = [-2., 2.]
    if baseline in ('a0', 'a1'):
        bounds[baseline] = [.98, 1.02]
        if baseline == 'a1':
            prior['a1'] = prior.pop('a0')
    data = elca.transit(time, prior) + .001 * np.sin(np.arange(count))
    return elca.lc_fitter(
        time, data, np.full(count, .003), np.linspace(1., 2., count), prior, bounds,
        mode=None, verbose=False, fixed_flux_baseline=baseline == 'fixed',
        use_impactparameter_rather_than_inclination_to_fit=geometry != 'inclination',
        baseline_fit_mask=np.arange(count) % 3 != 0,
        duration_prior=dict(applied=duration, expected_duration=elca.transit_duration(prior),
                            sigma_log_duration=.15),
        exposure_times_seconds=np.linspace(15., 40., count),
    )


def compare(monkeypatch, fit):
    reference = capture(monkeypatch, fit, 'python')
    native = capture(monkeypatch, fit, 'cpp')
    assert reference['keys'] == native['keys']
    rng = np.random.default_rng(2910)
    cube = rng.uniform(.03, .97, (29, len(reference['keys'])))
    points = reference['transform'](cube)
    np.testing.assert_array_equal(points, native['transform'](cube))
    expected = reference['loglike'](points)
    actual = native['loglike'](points)
    np.testing.assert_allclose(actual, expected, rtol=3.e-12, atol=3.e-7)
    assert native['loglike'](points[:0]).shape == (0,)
    assert native['loglike'](points[0]) == actual[0]
    assert isinstance(native['loglike'](points[0]), float)
    np.testing.assert_allclose(native['loglike'](points[::-1]), actual[::-1], rtol=0, atol=0)
    return reference, native, points


@pytest.mark.parametrize('law', ['claret', 'hybrid4'])
@pytest.mark.parametrize('geometry', ['impact', 'inclination', 'fixed', 'direct_b'])
@pytest.mark.parametrize('baseline', ['profile', 'fixed', 'a0', 'a1'])
@pytest.mark.parametrize('duration', [False, True])
def test_live_likelihood_modes_match_reference(monkeypatch, law, geometry, baseline, duration):
    fit = fixture(law, geometry, baseline, duration)
    prior = copy.deepcopy(fit.prior)
    compare(monkeypatch, fit)
    assert prior == fit.prior


@pytest.mark.parametrize('count', [0, 1, 7, 8, 127, 128, 129, 240, 337])
def test_reduction_lengths_and_fixed_airmass(monkeypatch, count):
    fit = fixture(count=count)
    fit.bounds.pop('a2')
    compare(monkeypatch, fit)


@pytest.mark.parametrize('law', ['claret', 'hybrid4'])
def test_sampled_orbit_and_atmosphere_coefficients(monkeypatch, law):
    fit = fixture(law=law)
    fit.bounds.update(per=[3., 3.2], ecc=[0., .3], omega=[0., 360.], unused=[0., 1.])
    fit.bounds['h0' if law == 'hybrid4' else 'u0'] = [.3, .7]
    compare(monkeypatch, fit)


@pytest.mark.parametrize('field', ['data', 'dataerr', 'airmass'])
def test_invalid_observations_keep_existing_penalty(monkeypatch, field):
    fit = fixture()
    getattr(fit, field)[0] = np.nan
    reference, native, points = compare(monkeypatch, fit)
    np.testing.assert_array_equal(native['loglike'](points), np.full(len(points), elca.BAD_LOG_LIKELIHOOD))


def test_clipped_geometry_has_no_replay_rejection_or_prior(monkeypatch):
    fit = fixture(duration=False)
    reference, native, points = compare(monkeypatch, fit)
    i = native['keys'].index('b')
    points[:, i] = np.linspace(-30., 30., len(points))
    expected = reference['loglike'](points)
    actual = native['loglike'](points)
    np.testing.assert_allclose(actual, expected, rtol=3.e-12, atol=3.e-7)
    assert np.any(actual != elca.BAD_LOG_LIKELIHOOD)


def test_auto_keeps_custom_models_and_explicit_python_model(monkeypatch):
    fit = fixture()
    fit.prior['_claret_backend'] = 'python'
    assert capture(monkeypatch, fit, 'auto')
    assert fit.likelihood_backend == 'python'
    fit.prior['_claret_backend'] = 'cpp'
    monkeypatch.setattr(fit, '_transit_model', lambda times, values: np.ones_like(times))
    capture(monkeypatch, fit, 'auto')
    assert fit.likelihood_backend == 'python'


def test_custom_sample_label_for_airmass_matches_reference(monkeypatch):
    fit = fixture(geometry='fixed')
    fit.bounds = {('trend' if key == 'a2' else key): bounds for key, bounds in fit.bounds.items()}
    monkeypatch.setattr(fit, '_get_sampled_keys', lambda keys: ['rprs', 'ars', 'tmid', 'a2'])
    compare(monkeypatch, fit)


def test_short_sample_mapping_preserves_reference_scoring(monkeypatch):
    fit = fixture(geometry='fixed')
    monkeypatch.setattr(fit, '_get_sampled_keys', lambda keys: ['rprs', 'ars', 'tmid'])
    reference = capture(monkeypatch, fit, 'python')
    native = capture(monkeypatch, fit, 'cpp')
    # Exercise scoring independently of the box transform, which still owns
    # its original bounds. Unassigned physical values retain their priors.
    points = np.array([[.11, 8., 2460000.], [.13, 8.5, 2460000.001]])
    np.testing.assert_allclose(native['loglike'](points), reference['loglike'](points), rtol=3.e-12, atol=3.e-7)


def test_auto_uses_native_and_load_failures_propagate_before_sampler(monkeypatch):
    fit = fixture()
    capture(monkeypatch, fit, 'auto')
    assert fit.likelihood_backend == 'cpp'
    import exotic.api.likelihood_cpp as native
    def failed():
        raise RuntimeError('test native load failure')
    monkeypatch.setattr(native, '_library', failed)
    fit.likelihood_backend_requested = 'cpp'
    with pytest.raises(RuntimeError, match='test native load failure'):
        fit.fit_nested()
