import numpy as np
import pytest
from astropy.io import fits

import exotic.exotic as exotic_module


def write_frame(path, data, exposure_time=None):
    header = fits.Header()
    if exposure_time is not None:
        header["EXPTIME"] = exposure_time
    fits.writeto(path, np.asarray(data, dtype=float), header=header)
    return str(path)


def test_bias_and_darks_build_exposure_scaled_dark_current(tmp_path):
    master_bias = np.array([[100.0, 101.0], [102.0, 103.0]])
    dark_current = np.array([[0.5, 1.0], [1.5, 2.0]])
    dark_files = [
        write_frame(tmp_path / "dark_10s.fits", master_bias + dark_current * 10.0, 10.0),
        write_frame(tmp_path / "dark_30s.fits", master_bias + dark_current * 30.0, 30.0),
    ]

    result = exotic_module.process_dark_frames(dark_files, master_bias)

    np.testing.assert_allclose(result, dark_current)


def test_darks_without_bias_build_unscaled_master_biasdark(tmp_path):
    dark_1 = np.array([[110.0, 112.0], [114.0, 116.0]])
    dark_2 = np.array([[114.0, 116.0], [118.0, 120.0]])
    dark_files = [
        write_frame(tmp_path / "dark_10s.fits", dark_1, 10.0),
        write_frame(tmp_path / "dark_30s.fits", dark_2, 30.0),
    ]

    result = exotic_module.process_dark_frames(dark_files)

    np.testing.assert_allclose(result, np.median([dark_1, dark_2], axis=0))


def test_flats_are_calibrated_and_normalized_individually_before_stacking(tmp_path):
    master_bias = np.full((2, 2), 10.0)
    dark_current = np.full((2, 2), 2.0)
    corrected_flat_1 = np.array([[100.0, 200.0], [100.0, 200.0]])
    corrected_flat_2 = np.array([[1000.0, 4000.0], [1000.0, 4000.0]])
    flat_files = [
        write_frame(
            tmp_path / "flat_2s.fits",
            corrected_flat_1 + master_bias + dark_current * 2.0,
            2.0,
        ),
        write_frame(
            tmp_path / "flat_5s.fits",
            corrected_flat_2 + master_bias + dark_current * 5.0,
            5.0,
        ),
    ]
    expected = np.median(
        [
            corrected_flat_1 / np.median(corrected_flat_1),
            corrected_flat_2 / np.median(corrected_flat_2),
        ],
        axis=0,
    )

    result = exotic_module.process_flat_frames(flat_files, master_bias, dark_current)

    np.testing.assert_allclose(result, expected)


def test_full_bias_scaled_dark_and_flat_pipeline_recovers_science_counts(tmp_path):
    bias = np.array([[100.0, 101.0], [102.0, 103.0]])
    dark_current = np.array([[0.5, 1.0], [1.5, 2.0]])
    flat_response = np.array([[0.8, 1.2], [0.9, 1.1]])
    true_science = np.array([[200.0, 300.0], [400.0, 500.0]])

    bias_files = [
        write_frame(tmp_path / "bias_1.fits", bias - 1.0),
        write_frame(tmp_path / "bias_2.fits", bias + 1.0),
    ]
    dark_files = [
        write_frame(tmp_path / "dark_10s.fits", bias + dark_current * 10.0, 10.0),
        write_frame(tmp_path / "dark_20s.fits", bias + dark_current * 20.0, 20.0),
    ]
    flat_files = [
        write_frame(
            tmp_path / "flat_2s.fits",
            bias + dark_current * 2.0 + flat_response * 1000.0,
            2.0,
        ),
        write_frame(
            tmp_path / "flat_4s.fits",
            bias + dark_current * 4.0 + flat_response * 2000.0,
            4.0,
        ),
    ]
    science_file = write_frame(
        tmp_path / "science_30s.fits",
        bias + dark_current * 30.0 + true_science * flat_response,
        30.0,
    )

    master_bias = exotic_module.process_bias_frames(bias_files)
    master_dark = exotic_module.process_dark_frames(dark_files, master_bias)
    master_flat = exotic_module.process_flat_frames(flat_files, master_bias, master_dark)
    _, calibrated_science = exotic_module.load_calibrated_reduction_frame(
        science_file,
        master_dark,
        master_bias,
        master_flat,
        None,
        None,
        None,
    )

    np.testing.assert_allclose(master_bias, bias)
    np.testing.assert_allclose(master_dark, dark_current)
    np.testing.assert_allclose(master_flat, flat_response)
    np.testing.assert_allclose(calibrated_science, true_science)


def test_biasdark_only_is_unscaled_for_flats_and_science(tmp_path):
    master_biasdark = np.array([[10.0, 11.0], [12.0, 13.0]])
    flat_response = np.array([[0.8, 1.2], [0.9, 1.1]])
    true_science = np.array([[20.0, 30.0], [40.0, 50.0]])
    flat_files = [
        write_frame(tmp_path / "flat.fits", master_biasdark + flat_response * 1000.0, 2.0),
    ]
    master_flat = exotic_module.process_flat_frames(flat_files, None, master_biasdark)
    raw_science = master_biasdark + true_science * flat_response

    calibrated = exotic_module.apply_cals(
        raw_science,
        master_biasdark,
        None,
        master_flat,
        1,
        exposure_time=999.0,
    )

    np.testing.assert_allclose(master_flat, flat_response)
    np.testing.assert_allclose(calibrated, true_science)


def test_bias_only_corrects_flats_and_science(tmp_path):
    master_bias = np.array([[10.0, 11.0], [12.0, 13.0]])
    flat_response = np.array([[0.8, 1.2], [0.9, 1.1]])
    true_science = np.array([[20.0, 30.0], [40.0, 50.0]])
    flat_files = [
        write_frame(tmp_path / "flat.fits", master_bias + flat_response * 1000.0, 2.0),
    ]
    master_flat = exotic_module.process_flat_frames(flat_files, master_bias, None)
    raw_science = master_bias + true_science * flat_response

    calibrated = exotic_module.apply_cals(
        raw_science,
        None,
        master_bias,
        master_flat,
        1,
    )

    np.testing.assert_allclose(master_flat, flat_response)
    np.testing.assert_allclose(calibrated, true_science)


def test_scaled_dark_requires_positive_exposure_times(tmp_path):
    master_bias = np.full((2, 2), 10.0)
    dark_current = np.full((2, 2), 1.0)
    dark_file = write_frame(tmp_path / "dark_without_exposure.fits", np.full((2, 2), 12.0))
    flat_file = write_frame(tmp_path / "flat_without_exposure.fits", np.full((2, 2), 1000.0))

    with pytest.raises(ValueError, match="positive exposure time.*dark frame"):
        exotic_module.process_dark_frames([dark_file], master_bias)

    with pytest.raises(ValueError, match="positive exposure time.*flat frame"):
        exotic_module.process_flat_frames([flat_file], master_bias, dark_current)

    with pytest.raises(ValueError, match="positive exposure time.*science frame"):
        exotic_module.apply_cals(
            np.full((2, 2), 100.0),
            dark_current,
            master_bias,
            None,
            1,
            exposure_time=0.0,
        )
