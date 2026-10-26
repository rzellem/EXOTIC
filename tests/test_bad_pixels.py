import json

import numpy as np
import pytest
from astropy.io import fits

from exotic.api.bad_pixels import dark_masks, load_mask, low_pixel_mask, prepare_reference, repair_frame


def write_darks(tmp_path, count, exposures=None):
    rng = np.random.default_rng(72)
    files = []
    for i in range(count):
        image = 400.0 + 100 * i + rng.normal(0, 1, (80, 80))
        image[20, 20] += 100 * (-1) ** i
        image[30, 30] += 500  # hot, but stable after frame-offset removal
        header = fits.Header({'EXPTIME': exposures[i] if exposures else 60.0})
        path = tmp_path / f'dark{i}.fits'
        fits.writeto(path, image, header)
        files.append(path)
    return files


@pytest.mark.parametrize('count', [2, 4, 5, 7])
def test_variability_requires_five_and_ignores_common_offset(tmp_path, count):
    masks, reports = dark_masks(write_darks(tmp_path, count), (80, 80))
    assert reports[0]['variability_evaluated'] == (count >= 5)
    assert masks['dark_variable'][20, 20] == (count >= 5)
    assert not masks['dark_variable'][30, 30]
    assert masks['dark_variable'].sum() <= 2
    assert masks['dark_hot'][30, 30]


def test_mixed_exposures_do_not_fake_five_matching_darks(tmp_path):
    masks, reports = dark_masks(write_darks(tmp_path, 6, [30, 30, 30, 60, 60, 60]), (80, 80))
    assert not masks['dark_variable'].any()
    assert all(not r['variability_evaluated'] for r in reports)


def test_raw_dark_source_survives_master_reuse(tmp_path):
    raw = tmp_path / 'raw'
    raw.mkdir()
    write_darks(raw, 5)
    reference = prepare_reference(
        {'save': str(tmp_path), 'bad_pixel_dark_source': str(raw)}, (80, 80),
        dark_files=[tmp_path / 'MasterDark.fits'],
    )
    assert reference['summary']['dark_groups'][0]['frame_count'] == 5
    assert reference['mask'][20, 20]


def test_combined_dark_products_never_count_toward_temporal_minimum(tmp_path):
    files = write_darks(tmp_path, 5)
    for file in files:
        fits.setval(file, 'NCOMBINE', value=10)
    masks, reports = dark_masks(files, (80, 80))
    assert not masks['dark_variable'].any()
    assert all(not r['variability_evaluated'] for r in reports)


def test_duplicate_paths_do_not_count_as_independent_darks(tmp_path):
    files = write_darks(tmp_path, 1)
    masks, reports = dark_masks(files * 5, (80, 80))
    assert reports[0]['frame_count'] == 1
    assert not reports[0]['variability_evaluated']
    assert not masks['dark_variable'].any()


def test_variability_requires_five_finite_samples_per_pixel(tmp_path):
    files = write_darks(tmp_path, 5)
    with fits.open(files[0], mode='update') as hdus:
        hdus[0].data[20, 20] = np.nan
    masks, reports = dark_masks(files, (80, 80))
    assert reports[0]['variability_evaluated']
    assert not masks['dark_variable'][20, 20]
    assert masks['dark_nonfinite'][20, 20]


def test_input_options_are_read(tmp_path):
    from exotic.inputs import Inputs
    options = {'bad_pixel_map': 'camera-mask.fits', 'detect_bad_pixels_from_darks': False,
               'detect_low_pixels_before_photometry': False}
    file = tmp_path / 'inits.json'
    file.write_text(json.dumps({'user_info': {}, 'planetary_parameters': {}, 'optional_info': options}))
    inputs = Inputs(init_opt='y')
    inputs.comp_params(file, {})
    for key, value in options.items():
        assert inputs.info_dict[key] == value


def test_low_tail_keeps_stars_and_normal_negative_background():
    rng = np.random.default_rng(7)
    frame = rng.normal(-100, 4, (500, 650))
    frame[20, 20] = -1000
    frame[30, 30] = 10000
    mask, threshold = low_pixel_mask(frame)
    assert threshold < -110
    assert mask[20, 20]
    assert not mask[30, 30]
    assert mask.sum() < 200


def test_supplied_mask_shape_and_nonfinite_convention(tmp_path):
    path = tmp_path / 'mask.fits'
    fits.writeto(path, np.array([[0, 1], [np.nan, 0]]))
    np.testing.assert_array_equal(load_mask(path, (2, 2)), [[False, True], [True, False]])
    with pytest.raises(ValueError, match='shape'):
        load_mask(path, (3, 3))


def test_combined_mask_provenance_and_disabled_options(tmp_path):
    darks = write_darks(tmp_path, 4)
    external = tmp_path / 'external.fits'
    mask = np.zeros((80, 80), np.uint8)
    mask[10, 10] = 1
    fits.writeto(external, mask)
    info = {'save': str(tmp_path), 'bad_pixel_map': str(external)}
    reference = prepare_reference(info, mask.shape, dark_files=darks)
    assert reference['mask'][10, 10]
    assert reference['mask'][30, 30]
    assert reference['detect_low']
    summary = json.loads((tmp_path / 'working_artifacts/BadPixelSummary.json').read_text())
    assert summary['dark_variability_min_frames'] == 5
    assert not summary['dark_groups'][0]['variability_evaluated']
    assert prepare_reference({'detect_bad_pixels_from_darks': False,
                              'detect_low_pixels_before_photometry': False}, (80, 80)) is None


def test_repair_masks_all_neighbours_before_infilling_and_preserves_good_pixels():
    original = np.full((20, 20), 100.0)
    original[3:8, 3:8] = 9999.0
    original[0, 0] = np.nan
    mask = original > 1000
    result = repair_frame(original, {'mask': mask})
    np.testing.assert_array_equal(result, np.full_like(result, 100))
    assert original[4, 4] == 9999
    assert np.isnan(original[0, 0])
    with pytest.raises(ValueError, match='no valid pixels'):
        repair_frame(original, {'mask': np.ones_like(mask)})
    with pytest.raises(ValueError, match='shape'):
        repair_frame(original, {'mask': np.zeros((2, 2), bool)})


def test_low_only_reference_repairs_even_without_static_coordinates(tmp_path):
    rng = np.random.default_rng(9)
    frame = rng.normal(400, 3, (300, 300))
    frame[50, 50] = -1000
    reference = prepare_reference({'save': str(tmp_path)}, frame.shape)
    assert not reference['mask'].any()
    from exotic.exotic import repair_bad_pixels_in_frame
    repaired = repair_bad_pixels_in_frame(frame, reference)
    assert 390 < repaired[50, 50] < 410


def test_calibrated_loader_applies_detector_mask(tmp_path):
    from exotic.exotic import load_calibrated_reduction_frame
    path = tmp_path / 'science.fits'
    frame = np.full((30, 30), 450.)
    frame[15, 15] = 9999
    fits.writeto(path, frame, fits.Header({'EXPTIME': 60.0}))
    mask = np.zeros(frame.shape, bool)
    mask[15, 15] = True
    reference = {'mask': mask, 'summary': {}, 'detect_low': False}
    _, repaired = load_calibrated_reduction_frame(
        path, np.full_like(frame, 50), None, None, None, None, None,
        bad_pixel_reference=reference,
    )
    np.testing.assert_allclose(repaired, 400.)
