import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS, Sip

from exotic.api.bayer_binning import bin2x2, binned_header, repair_bayer_pixels


def test_flux_preserved_without_integer_overflow():
    image = np.full((6, 8), 60000, dtype=np.uint16)
    result = bin2x2(image)
    assert result.shape == (3, 4)
    assert result.dtype == np.float64
    assert np.all(result == 240000)
    assert result.sum() == image.sum()


def test_odd_edges_and_signed_calibrated_values():
    image = np.arange(35).reshape(5, 7) - 10
    result = bin2x2(image)
    assert result.shape == (2, 3)
    assert result[0, 0] == -24
    assert result.sum() == image[:4, :6].sum()
    with pytest.raises(ValueError):
        bin2x2(np.zeros((1, 10)))


def test_repairs_native_bayer_phase_before_sum():
    image = np.tile([[10., 100.], [1000., 10000.]], (4, 4))
    image[2, 2] = 1e6
    mask = np.zeros(image.shape, bool)
    mask[2, 2] = True
    repaired = repair_bayer_pixels(image, {'mask': mask, 'detect_low': False})
    assert repaired[2, 2] == 10
    assert np.all(bin2x2(repaired) == 11110)


@pytest.mark.parametrize('use_cd', [False, True])
@pytest.mark.parametrize('use_sip', [False, True])
def test_wcs_preserves_sky_at_superpixel_centres(use_cd, use_sip):
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [60.3, 48.8]
    wcs.wcs.crval = [125., -30.]
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    if use_cd:
        wcs.wcs.cd = np.array([[-0.0002, 0.00001], [0.00001, 0.0002]])
    else:
        wcs.wcs.cdelt = [-0.0002, 0.0002]
    if use_sip:
        wcs.wcs.ctype = ['RA---TAN-SIP', 'DEC--TAN-SIP']
        a = np.zeros((3, 3)); b = a.copy()
        a[2, 0] = 1e-5; b[0, 2] = -2e-5
        wcs.sip = Sip(a, b, None, None, wcs.wcs.crpix)
    header = wcs.to_header(relax=True)
    new_wcs = WCS(binned_header(header))
    points = np.array([[10., 12.], [24., 20.], [40., 35.]])
    np.testing.assert_allclose(new_wcs.all_pix2world(points, 0),
                               wcs.all_pix2world(points * 2 + .5, 0), atol=1e-10)


def test_header_units_and_input_preserved():
    header = fits.Header({'BAYERPAT': 'RGGB', 'PIXSCALE': 1.2, 'XBINNING': 1,
                          'YBINNING': 1, 'SATURATE': 65535, 'RDNOISE': 3., 'GAIN': 2.})
    result = binned_header(header)
    assert result['PIXSCALE'] == 2.4
    assert result['XBINNING'] == result['YBINNING'] == 2
    assert result['SATURATE'] == 262140
    assert result['RDNOISE'] == 6
    assert result['GAIN'] == 2
    assert 'BAYERPAT' not in result
    assert header['BAYERPAT'] == 'RGGB'


@pytest.mark.parametrize('pattern', ['RGGB', 'BGGR', 'GRBG', 'GBRG'])
def test_reduction_stages_calibration_then_native_repair_then_binning(tmp_path, pattern):
    from exotic.exotic import prepare_bin2x2_reduction, load_calibrated_reduction_frame

    image = np.tile([[10., 100.], [1000., 10000.]], (4, 4))
    image[2, 2] = 50000
    mask = np.zeros(image.shape, np.uint8); mask[2, 2] = 1
    source = tmp_path / 'raw.fits'; mask_path = tmp_path / 'mask.fits'
    # Bias + dark rate * exposure, then divide by the flat.
    fits.writeto(source, image * 2 + 5 + 3 * 10, fits.Header({'EXPTIME': 10, 'SATURATE': 1e6}))
    fits.writeto(mask_path, mask)
    info = dict(save=str(tmp_path), demosaic_fmt=pattern, demosaic_out='bin2x2',
                tar_coords=[4.5, 6.5], comp_stars=[[2.5, 2.5]], pixel_scale=1.2,
                pixel_bin='1x1', bad_pixel_map=str(mask_path),
                detect_bad_pixels_from_darks=False, detect_low_pixels_before_photometry=False)
    paths, configured = prepare_bin2x2_reduction(
        np.array([str(source)]), info, np.full(image.shape, 3.),
        np.full(image.shape, 5.), np.full(image.shape, 2.))
    assert np.all(fits.getdata(paths[0]) == 11110)
    assert fits.getdata(paths[0]).shape == (4, 4)
    assert fits.getdata(source)[2, 2] == 100035
    assert configured['tar_coords'] == [2., 3.]
    assert configured['comp_stars'] == [[1., 1.]]
    assert configured['pixel_scale'] == 2.4
    assert configured['pixel_bin'] == '2x2'
    assert configured['bad_pixel_map'] is None
    assert info['demosaic_out'] == 'bin2x2'
    empty = np.empty((0, 0))
    header, reloaded = load_calibrated_reduction_frame(paths[0], empty, empty, empty, None, None, None)
    assert header['EXOBIN'] == '2x2'
    np.testing.assert_array_equal(reloaded, fits.getdata(paths[0]))


def test_one_saturated_bayer_site_remains_rejectable(tmp_path):
    from exotic.exotic import prepare_bin2x2_reduction, saturation_value_from_header
    image = np.full((4, 4), 5.)
    image[0, 0] = 1000
    source = tmp_path / 'saturated.fits'
    fits.writeto(source, image, fits.Header({'EXPTIME': 1, 'SATURATE': 1000}))
    info = dict(save=str(tmp_path), demosaic_fmt='BGGR',
                detect_bad_pixels_from_darks=False, detect_low_pixels_before_photometry=False)
    empty = np.empty((0, 0))
    paths, configured = prepare_bin2x2_reduction([str(source)], info, empty, empty, empty)
    result = fits.getdata(paths[0]); header = fits.getheader(paths[0])
    assert saturation_value_from_header(header) == 4000
    assert result[0, 0] >= 4000
    assert result[1, 1] == 20
