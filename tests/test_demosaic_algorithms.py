import json
import pickle

import numpy as np
import pytest
from astropy.io import fits
import colour_demosaicing

from exotic.api.demosaicing import normalize_demosaic_algorithm
from exotic.inputs import Inputs
import exotic.exotic as reduction


@pytest.mark.parametrize('value, expected', [
    (None, 'bilinear'), ('', 'bilinear'), (' Bilinear ', 'bilinear'),
    ('Malvar2004', 'malvar2004'), ('MENON2007', 'menon2007'), ('DDFAPD', 'menon2007'),
])
def test_algorithm_names(value, expected):
    assert normalize_demosaic_algorithm(value) == expected


def test_invalid_algorithm_fails_explicitly():
    with pytest.raises(ValueError, match='Invalid Demosaic Algorithm'):
        normalize_demosaic_algorithm('menon2008')


@pytest.mark.parametrize('section', ['user_info', 'optional_info'])
@pytest.mark.parametrize('output', ['green', 'green_binned', 'blue_binned', 'red_binned'])
def test_init_file_accepts_algorithm_and_output_in_either_section(tmp_path, section, output):
    data = dict(user_info={}, optional_info={}, planetary_parameters={})
    data[section] = {'Demosaic Format': 'RGGB', 'Demosaic Output': output,
                     'Demosaic Algorithm': 'DDFAPD'}
    path = tmp_path / 'init.json'; path.write_text(json.dumps(data))
    inputs = Inputs(init_opt='y'); inputs.comp_params(path, {})
    assert inputs.info_dict['demosaic_algorithm'] == 'menon2007'
    assert inputs.info_dict['demosaic_fmt'] == 'RGGB'
    assert inputs.info_dict['demosaic_out'] == output


def test_init_default_and_optional_precedence(tmp_path):
    data = dict(user_info={'Demosaic Algorithm': 'malvar2004'},
                optional_info={'Demosaic Algorithm': 'menon2007'}, planetary_parameters={})
    path = tmp_path / 'init.json'; path.write_text(json.dumps(data))
    inputs = Inputs(init_opt='y'); inputs.comp_params(path, {})
    assert inputs.info_dict['demosaic_algorithm'] == 'menon2007'
    data['user_info'] = {}; data['optional_info'] = {}
    path.write_text(json.dumps(data))
    inputs = Inputs(init_opt='y'); inputs.comp_params(path, {})
    assert inputs.info_dict['demosaic_algorithm'] == 'bilinear'


@pytest.mark.parametrize('algorithm, function_name', [
    ('bilinear', 'demosaicing_CFA_Bayer_bilinear'),
    ('malvar2004', 'demosaicing_CFA_Bayer_Malvar2004'),
    ('menon2007', 'demosaicing_CFA_Bayer_Menon2007'),
])
@pytest.mark.parametrize('pattern', ['RGGB', 'BGGR', 'GRBG', 'GBRG'])
@pytest.mark.parametrize('output, weights', [('green', [0., 1., 0.]),
                                          ('gray', [.299, .587, .114]),
                                          ([1, 1, 0], [.5, .5, 0.])])
def test_reconstruction_and_worker_roundtrip(algorithm, function_name, pattern, output, weights):
    image = np.random.default_rng(51).uniform(-20., 300., (12, 14))
    # Exercise the same serialization used by worker initializer arguments.
    mix = pickle.loads(pickle.dumps(reduction.calculate_demosaic_mult(output, algorithm)))
    expected = getattr(colour_demosaicing, function_name)(image, pattern) @ np.asarray(weights)
    result = reduction.demosaic_img(image, pattern, output, mix, 1)
    assert result.shape == image.shape
    np.testing.assert_allclose(result, expected)


@pytest.mark.parametrize('algorithm', ['malvar2004', 'menon2007'])
def test_new_algorithms_keep_negative_float_values_from_integer_inputs(algorithm):
    image = np.zeros((12, 12), np.uint16); image[6, 6] = 60000
    if algorithm == 'menon2007':
        image = np.random.default_rng(5).integers(0, 65535, (12, 12), dtype=np.uint16)
    mix = reduction.calculate_demosaic_mult('green', algorithm)
    result = reduction.demosaic_img(image, 'RGGB', 'green', mix, 1)
    assert np.issubdtype(result.dtype, np.floating)
    assert np.min(result) < 0


def test_existing_bilinear_integer_behavior_is_preserved():
    image = np.random.default_rng(7).integers(0, 60000, (12, 12), dtype=np.uint16)
    mix = reduction.calculate_demosaic_mult('gray')
    expected = (colour_demosaicing.demosaicing_CFA_Bayer_bilinear(image, 'RGGB') @ mix).astype(image.dtype)
    result = reduction.demosaic_img(image, 'RGGB', 'gray', mix, 1)
    assert result.dtype == image.dtype
    np.testing.assert_array_equal(result, expected)


def test_default_output_is_green_and_bin2x2_bypasses_algorithm(monkeypatch):
    image = np.random.default_rng(11).uniform(10, 100, (12, 12))
    mix = reduction.calculate_demosaic_mult(None, 'menon2007')
    np.testing.assert_allclose(reduction.demosaic_img(image, 'RGGB', None, mix, 1),
                              colour_demosaicing.demosaicing_CFA_Bayer_Menon2007(image, 'RGGB')[:, :, 1])
    def fail(*args, **kwargs):
        raise AssertionError('bin2x2 must not reconstruct RGB')
    monkeypatch.setattr('exotic.api.demosaicing.reconstruct_bayer', fail)
    mix = reduction.calculate_demosaic_mult('bin2x2', 'menon2007')
    result = reduction.demosaic_img(image, 'RGGB', 'bin2x2', mix, 1)
    assert result.shape == (6, 6)
    assert result.sum() == pytest.approx(image.sum())


@pytest.mark.parametrize('algorithm', ['malvar2004', 'menon2007'])
def test_frame_and_pool_loaders_use_the_selected_algorithm(tmp_path, monkeypatch, algorithm):
    image = np.random.default_rng(2).uniform(10, 300, (12, 12))
    path = tmp_path / 'science.fits'
    fits.writeto(path, image + 10, fits.Header({'EXPTIME': 10.}))
    empty = np.empty((0, 0)); bias = np.full(image.shape, 10.)
    mix = reduction.calculate_demosaic_mult('green', algorithm)
    expected = reduction.demosaic_img(image, 'RGGB', 'green', mix, 1)
    _, serial = reduction.load_calibrated_reduction_frame(str(path), empty, bias, empty, 'RGGB', 'green', mix)
    context = dict(generalDark=empty, generalBias=bias, generalFlat=empty,
                   demosaic_fmt='RGGB', demosaic_out='green', demosaic_mult=pickle.loads(pickle.dumps(mix)),
                   bad_pixel_reference=None)
    monkeypatch.setattr(reduction, '_ALIGNMENT_POOL_CONTEXT', context)
    monkeypatch.setattr(reduction, '_BAD_PIXEL_PRECHECK_POOL_CONTEXT', context)
    _, alignment = reduction._load_alignment_worker_frame(str(path))
    precheck = reduction._load_bad_pixel_precheck_worker_frame(str(path))
    for result in (serial, alignment, precheck):
        np.testing.assert_allclose(result, expected)


@pytest.mark.parametrize('output', ['green_binned', 'blue_binned', 'red_binned'])
@pytest.mark.parametrize('algorithm', ['bilinear', 'malvar2004', 'menon2007'])
def test_channel_binning_never_interpolates(monkeypatch, output, algorithm):
    def fail(*args, **kwargs):
        raise AssertionError('Native channel binning must not reconstruct RGB')
    monkeypatch.setattr(reduction, 'demosaicing_CFA_Bayer_bilinear', fail)
    monkeypatch.setattr('exotic.api.demosaicing.reconstruct_bayer', fail)
    image = np.tile([[10., 20.], [30., 40.]], (3, 4))
    mix = pickle.loads(pickle.dumps(reduction.calculate_demosaic_mult(output, algorithm)))
    result = reduction.demosaic_img(image, 'RGGB', output, mix, 1)
    expected = {'green_binned': 50., 'red_binned': 10., 'blue_binned': 40.}[output]
    np.testing.assert_array_equal(result, np.full((3, 4), expected))
