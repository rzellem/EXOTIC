"""Compare the batched kernel with the original saved-draw likelihood loop."""
import copy

import numpy as np
import pytest

from exotic.api.elca import lc_fitter, transit, transit_duration
from exotic.orbital_uncertainty import _conditional_orbit_draws, _conditional_orbit_draws_python


def replay_fixture(law='claret', geometry='b', baseline='profile', duration=True, count=75):
    rng = np.random.default_rng(2910)
    time = 2460000. + np.linspace(-.13, .13, count)
    base = dict(rprs=.12, per=3.12, ars=8.21, ecc=.018, inc=84., omega=97.,
                tmid=2460000., a2=.01, a0=1.008, b=.8,
                u0=1.77248, u1=-2.96806, u2=3.44175, u3=-1.24617,
                h0=.4, h1=.3, h2=.1, h3=.05,
                _limb_darkening_law=law, _claret_backend='python', _hybrid4_backend='numpy')
    keys = ['rprs', 'ars', 'tmid', 'a2'] + ([geometry] if geometry else [])
    if baseline in ('a0', 'a1'):
        keys += [baseline]
        if baseline == 'a1':
            base['a1'] = base.pop('a0')
    fit = lc_fitter(time, np.ones(count), np.full(count,.003), np.linspace(1,2,count),
                    base, {key: [0,1] for key in keys}, mode=None,
                    exposure_times_seconds=np.linspace(15,40,count))
    fit.sampled_keys = keys
    data = fit._transit_model(time, base) + rng.normal(0,.003,count)
    fit.orbital_likelihood_context = dict(
        base_physical=copy.deepcopy(base), time=time.copy(), data=data,
        inverse_dataerr=np.full(count,1/.003), centered_airmass=np.linspace(-.5,.5,count),
        baseline_static_mask=np.arange(count)%3!=0, baseline_weights=np.full(count,1/.003**2),
        uses_fixed_flux_baseline=baseline=='fixed', has_free_flux_baseline=baseline in ('a0','a1','stored'),
        fixed_flux_baseline_value=1.008, free_flux_baseline_key=baseline if baseline in ('a0','a1') else None,
        duration_prior_valid=duration, expected_duration=transit_duration(base), sigma_log_duration=.15,
    )
    n = 29
    values = dict(rprs=rng.uniform(.08,.16,n), ars=rng.uniform(.7,12,n),
                  tmid=base['tmid']+rng.uniform(-.01,.01,n), a2=rng.normal(.01,.02,n),
                  b=rng.uniform(0,1.2,n), inc=rng.uniform(75,90,n))
    if baseline in ('a0','a1'):
        values[baseline] = rng.uniform(.99,1.02,n)
    weights = rng.uniform(.1,1,n)
    return fit, values, weights


@pytest.mark.parametrize('law',['claret','hybrid4'])
@pytest.mark.parametrize('geometry',['b','inc',None])
@pytest.mark.parametrize('baseline',['profile','fixed','a0','a1','stored'])
def test_saved_draw_likelihoods_and_geometry_match_python(law,geometry,baseline):
    fit,values,weights = replay_fixture(law,geometry,baseline)
    before = copy.deepcopy(values)
    original_weights = weights.copy()
    for ecc,omega in [(0.,0.),(.025,270.),(.7,90.)]:
        expected = _conditional_orbit_draws_python(fit,ecc,omega,values,weights)
        fit._orbital_replay_backend = 'cpp'
        actual = _conditional_orbit_draws(fit,ecc,omega,values,weights)
        for key in expected[0]:
            np.testing.assert_allclose(actual[0][key],expected[0][key],rtol=0,atol=3.e-14)
        np.testing.assert_allclose(actual[1],expected[1],rtol=2.e-12,atol=2.e-7)
        np.testing.assert_allclose(actual[2],expected[2],rtol=0,atol=5.e-16)
        valid = np.isfinite(expected[1])
        if valid.any():
            old = np.exp(expected[1][valid]-np.max(expected[1][valid]))*weights[valid]
            new = np.exp(actual[1][valid]-np.max(actual[1][valid]))*weights[valid]
            np.testing.assert_allclose(new/new.sum(),old/old.sum(),rtol=1.e-7,atol=1.e-8)
    for key in before:
        np.testing.assert_array_equal(values[key],before[key])
    np.testing.assert_array_equal(weights,original_weights)


@pytest.mark.parametrize('mask,weight',[(False,1.),(True,0.),(True,np.nan)])
def test_existing_degenerate_profiled_baseline_behavior(mask,weight):
    fit,values,weights = replay_fixture(duration=False)
    fit.orbital_likelihood_context['baseline_static_mask'][:] = mask
    fit.orbital_likelihood_context['baseline_weights'][:] = weight
    expected = _conditional_orbit_draws_python(fit,.018,97.,values,weights)
    fit._orbital_replay_backend = 'cpp'
    actual = _conditional_orbit_draws(fit,.018,97.,values,weights)
    np.testing.assert_allclose(actual[1],expected[1],rtol=2.e-12,atol=2.e-7)


def test_native_dispatch_and_full_numpy_pairwise_sum(monkeypatch):
    fit,values,weights = replay_fixture(count=337)
    expected = _conditional_orbit_draws_python(fit,.025,270.,values,weights)
    fit.orbital_likelihood_context['base_physical']['_claret_backend'] = 'cpp'
    def wrong_path(*args):
        raise AssertionError('native dispatch called scalar Python replay')
    monkeypatch.setattr('exotic.orbital_uncertainty._conditional_orbit_draws_python',wrong_path)
    actual = _conditional_orbit_draws(fit,.025,270.,values,weights)
    np.testing.assert_allclose(actual[1],expected[1],rtol=2.e-12,atol=2.e-7)
