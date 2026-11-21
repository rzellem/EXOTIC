"""Native orbit/integration parity against the retained Python model path."""
import numpy as np
import pytest

from exotic.api import elca, limbdarkfast_cpp


def parameters(**updates):
    return dict(rprs=.13, per=3.2, ars=9., ecc=0., inc=88., omega=120.,
                tmid=2460000., h0=.37, h1=.007, h2=.316, h3=.019,
                _limb_darkening_law='hybrid4', **updates)


@pytest.mark.parametrize('ecc', [0., .05, .4, .8])
@pytest.mark.parametrize('omega', [-90., 0., 90., 270.])
@pytest.mark.parametrize('inclination', [84., 89.5])
@pytest.mark.parametrize('exposure', [None, 60., 'variable'])
def test_native_matches_reference_through_ingress_egress_and_multiple_orbits(ecc, omega, inclination, exposure):
    p = parameters()
    p.update(ecc=ecc, omega=omega, inc=inclination)
    # Include negative epochs, full orbits and points behind the star.
    time = p['tmid'] + np.r_[np.linspace(-.15, .15, 301), np.linspace(-2.3*p['per'], 2.3*p['per'], 181)]
    if exposure is not None:
        seconds = np.resize([0., .5, 1., 1.01, 15., 60., 180., np.nan, -1., np.inf],time.size) if exposure == 'variable' else exposure
        p['_exposure_time_days'] = np.asarray(seconds)/86400.
    old = elca.transit(time, {**p, '_hybrid4_backend': 'cpp_flux'})
    new = elca.transit(time, {**p, '_hybrid4_backend': 'cpp'})
    np.testing.assert_allclose(new, old, rtol=0, atol=2.e-12)


@pytest.mark.parametrize('samples', [3, 4, 7, 16])
@pytest.mark.parametrize('tolerance', [0., 1.e-5, 1.])
def test_custom_exposure_controls_and_scalar_arrays(samples, tolerance):
    p = parameters()
    p.update(ecc=.25, _exposure_time_days=180./86400.,
             _exposure_smearing_supersample=samples,
             _exposure_smearing_change_tolerance=tolerance)
    for t in [2460000.+np.linspace(-.12,.12,60).reshape(6,10)[:,::2],
              np.asarray(2460000.), np.empty(0)]:
        old = elca.transit(t,{**p,'_hybrid4_backend':'cpp_flux'})
        new = elca.transit(t,{**p,'_hybrid4_backend':'cpp'})
        assert new.shape == np.asarray(t).shape
        np.testing.assert_allclose(new,old,rtol=0,atol=2.e-12)


@pytest.mark.parametrize('exposure', [[1.,2.], 'nonnumeric', None, np.nan, 1./86400.])
def test_existing_exposure_normalization_is_preserved(exposure):
    p = parameters()
    p['_exposure_time_days'] = exposure
    t = p['tmid']+np.linspace(-.1,.1,50)
    old = elca.transit(t,{**p,'_hybrid4_backend':'cpp_flux'})
    new = elca.transit(t,{**p,'_hybrid4_backend':'cpp'})
    np.testing.assert_allclose(new,old,rtol=0,atol=2.e-12)


def test_cpp_path_does_not_call_python_orbit_or_integration(monkeypatch):
    t = 2460000.+np.linspace(-.1,.1,50)
    p = parameters()
    p.update(_hybrid4_backend='cpp',_exposure_time_days=120./86400.)
    expected = elca.transit(t,p)
    def forbidden(*args,**kwargs):
        raise AssertionError('Python orbit/exposure path was used')
    for name in ('_instantaneous_transit','_exposure_time_days_for_times',
                 '_exposure_smearing_candidate_mask','transit_duration'):
        monkeypatch.setattr(elca,name,forbidden)
    np.testing.assert_array_equal(elca.transit(t,p),expected)


def test_native_errors_are_not_replaced_with_an_instantaneous_model():
    p = parameters()
    p.update(rprs=-.1,_hybrid4_backend='cpp',_exposure_time_days=120./86400.)
    with pytest.raises(ValueError, match='r >= 0'):
        elca.transit([2460000.],p)


def test_native_integration_averages_the_same_midpoint_samples():
    p = parameters()
    t = p['tmid']+np.linspace(-.08,.08,100)
    exposure=300./86400.
    settings={**p,'_hybrid4_backend':'cpp','_exposure_time_days':exposure,
              '_exposure_smearing_change_tolerance':0.,'_exposure_smearing_supersample':7}
    samples=t[:,None]+exposure*((np.arange(7)+.5)/7-.5)
    expected=elca.transit(samples,{**p,'_hybrid4_backend':'cpp_flux'}).mean(axis=1)
    np.testing.assert_allclose(elca.transit(t,settings),expected,rtol=0,atol=2.e-12)
