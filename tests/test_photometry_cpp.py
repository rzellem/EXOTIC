"""Compare native photometry to the independently retained SciPy/Photutils path."""
import numpy as np
import pytest
from scipy.optimize import least_squares

from exotic import exotic as runtime
from exotic.api import photometry_cpp


@pytest.mark.parametrize('center', [False, True])
@pytest.mark.parametrize('position', [(40.2, 35.7), (.1, 4.4), (-40., 150.)])
@pytest.mark.parametrize('kind', ['normal', 'ties', 'nonfinite', 'bigendian', 'strided'])
def test_aperture_grid_matches_reference(monkeypatch, center, position, kind):
    rng = np.random.default_rng(567)
    data = rng.normal(100., 4., (77, 93))
    if kind == 'ties':
        data = np.round(data)
    if kind == 'nonfinite':
        data[34, 40], data[30, 45], data[10, 7] = np.nan, np.inf, -np.inf
    if kind == 'bigendian':
        data = data.astype('>f4')
    if kind == 'strided':
        data = data[::-1, ::-1]
    args = (data, 0, *position, [2.25, 5.3, 13., np.nan, -1.], [0., 3.7, 8., np.nan, -1.])
    kwargs = dict(fast_mode=center, sigma_hint=2.7, return_noise=True,
                  exposure_s=60., airmass=1.5, aperture_correction_factors=[1.02, 1., .99, 1., 1.],
                  noise_config={'gain_e_per_adu': 1.8, 'read_noise_electrons': 7.})
    monkeypatch.setattr(runtime, '_PHOTOMETRY_BACKEND', 'python')
    reference = runtime.compute_star_aperture_grid(*args, **kwargs)
    monkeypatch.setattr(runtime, '_PHOTOMETRY_BACKEND', 'cpp')
    actual = runtime.compute_star_aperture_grid(*args, **kwargs)
    for i in range(2):
        np.testing.assert_allclose(actual[i], reference[i], rtol=2e-13, atol=2e-10, equal_nan=True)
    for key in reference[2]:
        np.testing.assert_allclose(actual[2][key], reference[2][key], rtol=2e-12, atol=2e-10, equal_nan=True)


@pytest.mark.parametrize('center', [False, True])
@pytest.mark.parametrize('position', [(30.2, 25.1), (.2, 1.3), (-100., 0.)])
def test_standalone_aperture_and_sky_match(monkeypatch, center, position):
    data = np.random.default_rng(13).normal(100, 5, (60, 70))
    results = []
    for backend in ['python', 'cpp']:
        monkeypatch.setattr(runtime, '_PHOTOMETRY_BACKEND', backend)
        results.append((runtime.aperPhot(data, 0, *position, r=3.7, dr=4.2, fast_mode=center),
                        runtime.skybg_phot(data, 0, *position, r=8., dr=6., ptol=96., fast_mode=center)))
    for i in range(2):
        np.testing.assert_allclose(results[1][i], results[0][i], rtol=2e-13, atol=2e-10, equal_nan=True)


@pytest.mark.parametrize('method', ['2-point', '3-point'])
@pytest.mark.parametrize('widths,angle', [((1.2, 2.8), -.4), ((4., 3.), .32), ((2., 2.), 0.)])
def test_native_gaussian_fit_matches_scipy(method, widths, angle):
    y, x = np.mgrid[:31, :29]
    truth = [14.3, 15.2, 1300., *widths, angle, 72.]
    data = runtime.gaussian_psf(x, y, *truth) + np.random.default_rng(61).normal(0, 3, x.shape)
    initial = [14., 15., data.max()-data.min(), 1., 1., 0., data.min()]
    bounds = ([7., 8., 0., .5, .5, -np.pi/4, data.min()-1],
              [21., 22., 1e7, 20., 20., np.pi/4, data.max()+1])
    reference = least_squares(lambda p: (data-runtime.gaussian_psf(x, y, *p)).ravel(), initial,
                              bounds=bounds, jac=method, method='trf')
    actual = photometry_cpp.fit_gaussian(x, y, data, initial, bounds, method)
    # Near circular profiles can exchange widths/angle. Compare their actual
    # pixel models, integrated flux and residual cost instead of angle labels.
    np.testing.assert_allclose(runtime.gaussian_psf(x, y, *actual.x),
                               runtime.gaussian_psf(x, y, *reference.x), rtol=2e-7, atol=2e-5)
    assert actual.cost == pytest.approx(reference.cost, rel=1e-10)
    assert actual.x[2]*actual.x[3]*actual.x[4] == pytest.approx(
        reference.x[2]*reference.x[3]*reference.x[4], rel=1e-7)


def test_native_setup_failure_is_not_converted_to_a_moment_fallback(monkeypatch):
    def fail():
        raise RuntimeError('Cannot build native photometry')
    monkeypatch.setattr(photometry_cpp, 'prepare', fail)
    with pytest.raises(RuntimeError, match='Cannot build'):
        runtime.fit_centroid(np.ones((50, 50)), [25., 25.], 0)


def test_modern_flux_fit_is_stable_for_a_rotated_elliptical_star(monkeypatch):
    # The old circular starting profile selected different angle-bound minima
    # on this noisy elliptical star, changing integrated flux by about 3.4%.
    y, x = np.mgrid[:71, :73]
    data = runtime.gaussian_psf(x, y, 36.2, 35.3, 770., 1.86, .94, -.385, 700.)
    data += np.random.default_rng(4).normal(0., 12., data.shape)
    rows = []
    for backend in ['python', 'cpp']:
        monkeypatch.setattr(runtime, '_PHOTOMETRY_BACKEND', backend)
        centroid = runtime.fit_centroid(data, [36.2, 35.3], 0)
        for offset in [0., 1e-8, -1e-8]:
            seed = centroid.copy()
            seed[:2] += offset
            rows.append(runtime.fit_psf_photometry_flux_row(data, seed, 0))
    flux = [np.prod(row[2:5]) for row in rows]
    np.testing.assert_allclose(flux, flux[0], rtol=1e-7)
    for row in rows[1:]:
        np.testing.assert_allclose(runtime.gaussian_psf(x, y, *row),
                                   runtime.gaussian_psf(x, y, *rows[0]), rtol=1e-7, atol=2e-5)


@pytest.mark.parametrize('radius', [0., -1., np.nan])
def test_standalone_invalid_radius_keeps_photutils_error(monkeypatch, radius):
    for backend in ['python', 'cpp']:
        monkeypatch.setattr(runtime, '_PHOTOMETRY_BACKEND', backend)
        with pytest.raises(ValueError, match='positive scalar'):
            runtime.aperPhot(np.ones((30, 30)), 0, 15., 15., r=radius, dr=0.)


def test_standalone_infinite_circle_keeps_reference_empty_mask_behavior(monkeypatch):
    results = []
    for backend in ['python', 'cpp']:
        monkeypatch.setattr(runtime, '_PHOTOMETRY_BACKEND', backend)
        results.append(runtime.aperPhot(np.ones((30, 30)), 0, 15., 15., r=np.inf, dr=0.))
    assert results == [(0., 0.), (0., 0.)]


@pytest.mark.parametrize('nan_setting', ['percentile', 'sigma'])
def test_nonfinite_sky_statistics_keep_reference_nan_semantics(monkeypatch, nan_setting):
    data = np.random.default_rng(41).normal(100., 3., (60, 65))
    kwargs = {'ptol': np.nan} if nan_setting == 'percentile' else {}
    if nan_setting == 'sigma':
        monkeypatch.setattr(runtime, 'SKY_BACKGROUND_SIGMA_CLIP', np.nan)
    values = []
    for backend in ['python', 'cpp']:
        monkeypatch.setattr(runtime, '_PHOTOMETRY_BACKEND', backend)
        values.append(runtime.skybg_phot(data, 0, 32., 30., **kwargs))
    np.testing.assert_allclose(values[1], values[0], rtol=1e-14, equal_nan=True)


def test_custom_psf_keeps_the_requested_model(monkeypatch):
    def custom(x, y, *p):
        return runtime.gaussian_psf(x, y, *p)
    def fail():
        raise AssertionError('A custom PSF must not enter the Gaussian native kernel')
    monkeypatch.setattr(photometry_cpp, 'prepare', fail)
    y, x = np.mgrid[:55, :55]
    data = custom(x, y, 27., 27., 500., 2., 2.1, .1, 100.)
    result = runtime.fit_centroid(data, [27., 27.], 0, psf_function=custom)
    np.testing.assert_allclose(result[:2], [27., 27.], atol=1e-6)
