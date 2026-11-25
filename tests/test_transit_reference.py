"""Frozen PyLightcurve 4.0.4 vectors; no external package or databases needed."""
import json
from pathlib import Path

import numpy as np
import pytest

from exotic.api import transit_reference as reference
from exotic.api import joint_fitter as joint


VECTORS = json.loads((Path(__file__).parent / 'data' / 'pylightcurve_4_0_4_reference.json').read_text())


@pytest.mark.parametrize('case', VECTORS['calculations'], ids=lambda case: case['id'])
def test_original_numerical_results(case):
    args = [np.asarray(arg) if isinstance(arg, list) else arg for arg in case['args']]
    actual = getattr(reference, case['function'])(*args, **case.get('kwargs', {}))
    np.testing.assert_allclose(actual, case['expected'], rtol=3.e-12, atol=3.e-10)


def joint_inputs():
    """Two light curves, two RV instruments and transit/eclipse ephemerides."""
    prior = dict(rprs=.11, per=3., ars=9., ecc=.15, inc=88.5, omega=70., tmid=0.,
                 u0=1.77248, u1=-2.96806, u2=3.44175, u3=-1.24617,
                 rstar=1.1, mstar=1.1, mplanet=.8, rv_linear=.3, rv_quad=.01,
                 fpfs=.02, c1=.00002, c2=.00001, c3=.00003, c4=-.00001)
    prior['mu'] = prior['mplanet'] * joint.Mjup / (prior['mstar'] * joint.Msun)
    light_curves, velocities = [], []
    for i in range(2):
        times = np.linspace(-.1, .1, 29) + i * 3.
        flux = 1. - (.011 + .001 * i) * np.exp(-.5 * ((times - 3. * i) / .025) ** 2)
        flux += .0003 * np.sin(np.arange(len(times)))
        light_curves.append(dict(time=times, flux=flux, ferr=np.full(len(times), .001),
                                 airmass=np.ones(len(times)), priors=prior.copy()))
        rv_times = np.linspace(-2.1, 4.2, 17) + i * .04
        velocities.append(dict(time=rv_times, vel=20. + 60. * np.sin(2. * np.pi * rv_times / 3.),
                               velerr=np.full(len(rv_times), 4.), priors=prior.copy()))
    ephemeris = dict(tmid=np.array([[-3.0001, .0004], [3.0002, .0003]]),
                     emid=np.array([[1.61, .004], [4.611, .005]]),
                     historic=np.empty((0, 2)), noise=.0001,
                     prior=dict(tmid=[0., .001], per=[3., .001], emid=[1.61, .005]))
    bounds = dict(per=[2.99, 3.01], ecc=[.12, .18], omega=[60., 80.], inc=[87.5, 89.5],
                  tmid=[-.002, .002], rstar=[1.1, .04])
    return dict(lc_data=light_curves, local_lc_bounds=[dict(rprs=[.09, .13])] * 2,
                rv_data=velocities, local_rv_bounds=[{}, {}], global_bounds=bounds,
                ephemeris=ephemeris, individual_fit=False, verbose=False)


def capture_joint_callbacks():
    callbacks = {}

    class Captured(BaseException):
        pass

    def sampler(keys, loglike, transform, **kwargs):
        callbacks.update(keys=keys, loglike=loglike, transform=transform)
        raise Captured()

    original = joint.ReactiveNestedSampler
    joint.ReactiveNestedSampler = sampler
    try:
        with pytest.raises(Captured):
            joint.joint_fitter(**joint_inputs())
    finally:
        joint.ReactiveNestedSampler = original
    return callbacks


def test_joint_likelihood_and_prior_transform_match_original():
    callbacks = capture_joint_callbacks()
    expected = VECTORS['joint']
    assert callbacks['keys'] == expected['keys']
    points = [callbacks['transform'](np.asarray(cube)) for cube in expected['unit_cube']]
    np.testing.assert_allclose(points, expected['points'], rtol=3.e-13, atol=1.e-14)
    scores = [callbacks['loglike'](point) for point in points]
    np.testing.assert_allclose(scores, expected['loglike'], rtol=5.e-11, atol=2.e-7)


@pytest.mark.parametrize('model', ['transit', 'eclipse', 'phasecurve', 'rv_model'])
def test_joint_models_match_original(model):
    prior = joint_inputs()['lc_data'][0]['priors']
    case = VECTORS['joint_models'][model]
    np.testing.assert_allclose(getattr(joint, model)(np.asarray(case['times']), prior),
                               case['expected'], rtol=3.e-12, atol=3.e-10)
