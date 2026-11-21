"""Exact parity of low-tail decisions, local medians and nearest-pixel ties."""
from fractions import Fraction

import numpy as np
import pytest
from astropy.io import fits

from exotic.api import bad_pixels as reference
from exotic.api import image_preparation_cpp as native


def images():
    rng = np.random.default_rng(71)
    return [
        np.array([]), np.array(3.), np.full((4, 5), np.nan),
        np.array([[np.nan, np.inf, -np.inf, 2.2], [-4.9, -4.1, -.9, .9]]),
        np.repeat([-17., 2., 2**54, -1e300, 1e300], [31, 31, 1, 1, 1]),
        rng.normal(-200, 4, (30, 41)), rng.normal(300, 8, (53, 61))[::-2, ::-2],
        np.array([0., -0., -1.5, -1., 2., 2., 3., 3.]),
    ]


@pytest.mark.parametrize('image', images(), ids=range(8))
@pytest.mark.parametrize('gap', [-1, 0, 1, 17, 31])
@pytest.mark.parametrize('population', [0., .5, 25., np.inf, np.nan])
def test_low_threshold_and_mask_match(image, gap, population):
    expected, threshold = reference._low_pixel_mask_python(image, gap, population)
    actual, actual_threshold = native.low_pixel_mask(image, gap, population)
    assert actual.shape == image.shape
    np.testing.assert_array_equal(actual, expected)
    assert actual_threshold == threshold


@pytest.mark.parametrize('dtype', [np.int16, np.uint16, np.float32, np.float64, '>f8'])
@pytest.mark.parametrize('shape', [(0, 0), (1, 1), (1, 17), (17, 1), (17, 23), (41, 43)])
@pytest.mark.parametrize('detect_low', [False, True])
def test_repair_preserves_values_and_original_masks(dtype, shape, detect_low):
    rng = np.random.default_rng(104)
    frame = rng.integers(100, 500, shape).astype(dtype)
    mask = rng.uniform(size=shape) < .14
    if min(shape) > 9:
        mask[3:9, 3:9] = True  # Exercise interior clusters/EDT distance ties.
        mask[0:3, 0:3] = True
    ref = {'mask': mask, 'detect_low': detect_low}
    before, original_mask = frame.copy(), mask.copy()
    try:
        expected = reference._repair_frame_python(frame, ref)
    except ValueError as error:
        with pytest.raises(ValueError, match=str(error)):
            native.repair_frame(frame, ref)
    else:
        result = native.repair_frame(frame, ref)
        np.testing.assert_array_equal(result, expected)
        assert result.dtype == np.float64
    np.testing.assert_array_equal(frame, before)
    np.testing.assert_array_equal(mask, original_mask)


def test_nonfinite_noncontiguous_and_fortran_repair():
    base = np.arange(25*30, dtype=float).reshape(25, 30)
    base[2:8, 2:8] = np.nan
    base[9, 0] = np.inf
    for frame in (base[::-1, ::-1], np.asfortranarray(base)):
        ref = {'mask': np.zeros(frame.shape, bool), 'detect_low': False}
        np.testing.assert_array_equal(native.repair_frame(frame, ref), reference._repair_frame_python(frame, ref))


def test_even_neighbour_median_keeps_existing_overflow_behavior():
    frame = np.full((3, 3), 1e308)
    mask = np.zeros((3, 3), bool)
    mask[1, 1] = True
    with np.errstate(over='ignore'):
        expected = reference._repair_frame_python(frame, {'mask': mask})
        actual = native.repair_frame(frame, {'mask': mask})
    np.testing.assert_array_equal(actual, expected)
    assert np.isinf(actual[1, 1])


def test_owned_frame_can_reuse_buffer_without_using_previous_repairs():
    frame = np.arange(40*45, dtype=float).reshape(40, 45).copy()
    mask = np.zeros(frame.shape, bool)
    mask[1:8, 1:8] = True
    ref = {'mask': mask, 'detect_low': True}
    expected = reference._repair_frame_python(frame, ref)
    actual = native.repair_frame(frame, ref, copy=False)
    assert actual is frame
    np.testing.assert_array_equal(actual, expected)


def test_custom_histogram_arguments_retain_python_semantics(monkeypatch):
    monkeypatch.setattr(reference, 'PREPARATION_BACKEND', 'cpp')
    frame = np.repeat([1., 2., 3.], 30)
    for gap, population in [(0, 'ignored'), (1, Fraction(25, 1)), (1, 2**2048)]:
        a, ta = reference.low_pixel_mask(frame, gap, population)
        b, tb = reference._low_pixel_mask_python(frame, gap, population)
        np.testing.assert_array_equal(a, b)
        assert ta == tb


def test_header_access_does_not_decode_image_and_uses_one_open(tmp_path, monkeypatch):
    import exotic.exotic as runtime
    path = tmp_path / 'compressed.fits.fz'
    fits.HDUList([fits.PrimaryHDU(), fits.CompImageHDU(np.ones((9, 13)),
                  header=fits.Header({'EXPTIME': 15.}))]).writeto(path)
    original = fits.open
    calls = []
    def opened(*args, **kwargs):
        calls.append(kwargs)
        return original(*args, **kwargs)
    monkeypatch.setattr(runtime.bad_pixel_processing, 'PREPARATION_BACKEND', 'cpp')
    monkeypatch.setattr(fits, 'open', opened)
    def no_pixels(hdu):
        raise AssertionError('header access decoded image pixels')
    monkeypatch.setattr(fits.CompImageHDU, 'data', property(no_pixels))
    header = runtime.read_reduction_image_header(path)
    assert header['NAXIS1'] == 13 and header['NAXIS2'] == 9 and header['EXPTIME'] == 15.
    assert len(calls) == 1 and calls[0]['lazy_load_hdus'] is True
