import io
import sys
import types
import importlib.util

import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS


def _module_available(name):
    try:
        return importlib.util.find_spec(name) is not None
    except (ModuleNotFoundError, ValueError):
        return False


def _install_stub_module(name, **attrs):
    module = types.ModuleType(name)
    for attr, value in attrs.items():
        setattr(module, attr, value)
    sys.modules[name] = module
    return module


class _DummyDAOStarFinder:
    def __init__(self, *args, **kwargs):
        pass

    def __call__(self, *args, **kwargs):
        return None


if not _module_available("astroalign"):
    _install_stub_module("astroalign", PIXEL_TOL=1)
if not _module_available("astroquery"):
    _install_stub_module("astroquery")
if not _module_available("astroquery.simbad"):
    _install_stub_module("astroquery.simbad", Simbad=object)
if not _module_available("astroquery.gaia"):
    _install_stub_module("astroquery.gaia", Gaia=object)
if not _module_available("barycorrpy"):
    _install_stub_module("barycorrpy")
if not _module_available("barycorrpy.utc_tdb"):
    _install_stub_module("barycorrpy.utc_tdb", JDUTC_to_BJDTDB=lambda *args, **kwargs: None)
if not _module_available("imreg_dft"):
    _install_stub_module("imreg_dft")
if not _module_available("pyvo"):
    _install_stub_module("pyvo")
if not _module_available("photutils"):
    _install_stub_module("photutils")
if not _module_available("photutils.aperture"):
    _install_stub_module("photutils.aperture", CircularAperture=object)
if not _module_available("photutils.detection"):
    _install_stub_module("photutils.detection", DAOStarFinder=_DummyDAOStarFinder)
if not _module_available("colour_demosaicing"):
    _install_stub_module("colour_demosaicing", demosaicing_CFA_Bayer_bilinear=lambda *args, **kwargs: None)
if not _module_available("ldtk"):
    fake_ldtk = _install_stub_module("ldtk")
    fake_ldtk.LDPSet = type("LDPSet", (), {})
    fake_ldtk.ldtk = types.SimpleNamespace(LDPSet=fake_ldtk.LDPSet)
if not _module_available("ldtk.ldmodel"):
    _install_stub_module(
        "ldtk.ldmodel",
        LinearModel=type("LinearModel", (), {}),
        QuadraticModel=type("QuadraticModel", (), {}),
        NonlinearModel=type("NonlinearModel", (), {}),
    )
if not _module_available("lmfit"):
    _install_stub_module("lmfit")

_install_stub_module(
    "exotic.api.elca",
    lc_fitter=lambda *args, **kwargs: None,
    binner=lambda *args, **kwargs: None,
    transit=lambda *args, **kwargs: None,
    get_phase=lambda *args, **kwargs: None,
)

from exotic import exotic as exotic_module


def _gaussian_image(shape=(80, 80), center=(40.0, 35.0), amplitude=5000.0, sigma=2.0, background=100.0):
    y, x = np.indices(shape, dtype=float)
    cx, cy = center
    image = background + amplitude * np.exp(-((x - cx) ** 2 + (y - cy) ** 2) / (2.0 * sigma ** 2))
    return image


def _write_extension_wcs_fits(tmp_path, shape=(100, 120)):
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [shape[1] / 2.0, shape[0] / 2.0]
    wcs.wcs.crval = [210.0, 54.0]
    wcs.wcs.cdelt = np.array([-0.01, 0.01])
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    path = tmp_path / "extension_wcs.fits"
    hdul = fits.HDUList([
        fits.PrimaryHDU(),
        fits.ImageHDU(data=np.zeros(shape, dtype=float), header=wcs.to_header(), name="SCI"),
    ])
    hdul.writeto(path, overwrite=True)
    return path


def test_detect_frame_bad_pixels_flags_isolated_hot_pixel_but_not_broad_star_core():
    image = _gaussian_image(shape=(60, 60), center=(30.0, 30.0), amplitude=1200.0, sigma=1.8, background=100.0)
    image[10, 15] = 8000.0

    mask = exotic_module.detect_frame_bad_pixels(image)

    assert mask[10, 15]
    assert not mask[30, 30]


def test_build_persistent_bad_pixel_map_thresholds_recurrence_and_saves_outputs(tmp_path):
    frames = {}
    for frame_index in range(10):
        frame = np.full((9, 9), 100.0, dtype=float)
        if frame_index < 4:
            frame[2, 3] = 4000.0
        if frame_index < 3:
            frame[6, 5] = 3500.0
        frames[f"frame_{frame_index}.fits"] = frame

    reference = exotic_module.build_persistent_bad_pixel_map(
        list(frames.keys()),
        lambda file_name: frames[file_name],
        save_directory=tmp_path,
    )

    assert reference is not None
    assert reference["required_count"] == 4
    assert reference["mask"][2, 3]
    assert not reference["mask"][6, 5]

    count_image = fits.getdata(tmp_path / "temp" / "BadPixelDetectionCounts.fits")
    mask_image = fits.getdata(tmp_path / "temp" / "BadPixelMask.fits").astype(bool)

    assert count_image[2, 3] == 4
    assert count_image[6, 5] == 3
    assert mask_image[2, 3]
    assert not mask_image[6, 5]


def test_repair_bad_pixels_in_frame_replaces_known_bad_pixel_with_neighbor_median():
    image = np.arange(25, dtype=float).reshape(5, 5)
    image[2, 2] = 9999.0
    reference = {
        "mask": np.zeros((5, 5), dtype=bool),
        "coord_y": np.array([2]),
        "coord_x": np.array([2]),
    }
    reference["mask"][2, 2] = True

    repaired = exotic_module.repair_bad_pixels_in_frame(image, reference)

    assert repaired[2, 2] == pytest.approx(12.0)


def test_fit_centroid_uses_moment_fallback_when_psf_fit_fails(monkeypatch):
    image = _gaussian_image()
    low_flux_warnings = []

    def fail_least_squares(*args, **kwargs):
        raise ValueError("Residuals are not finite in the initial point.")

    monkeypatch.setattr(exotic_module, "least_squares", fail_least_squares)
    monkeypatch.setattr(
        exotic_module.plateStatus,
        "lowFluxAmplitudeWarning",
        lambda star_index, xc, yc: low_flux_warnings.append((star_index, xc, yc)),
    )

    result = exotic_module.fit_centroid(image, [40.0, 35.0], 0)

    assert np.isfinite(result[0])
    assert np.isfinite(result[1])
    assert abs(result[0] - 40.0) < 1.0
    assert abs(result[1] - 35.0) < 1.0
    assert low_flux_warnings == []


def test_fit_centroid_reports_consistent_background_between_fast_and_full_modes():
    image = _gaussian_image(center=(40.3, 35.7), amplitude=5000.0, sigma=2.0, background=123.4)

    fast_result = exotic_module.fit_centroid(image, [40.0, 36.0], 0, fast_mode=True)
    full_result = exotic_module.fit_centroid(image, [40.0, 36.0], 0, fast_mode=False)

    assert np.isfinite(fast_result[6])
    assert np.isfinite(full_result[6])
    assert full_result[6] == pytest.approx(fast_result[6], abs=1e-8)


def test_fit_centroid_or_warn_out_of_frame_skips_centroid_fit(monkeypatch):
    image = np.zeros((40, 50), dtype=float)
    out_of_frame_warnings = []

    monkeypatch.setattr(
        exotic_module,
        "fit_centroid",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("fit_centroid should be skipped")),
    )
    monkeypatch.setattr(
        exotic_module.plateStatus,
        "outOfFrameWarning",
        lambda star_index: out_of_frame_warnings.append(star_index),
    )

    result = exotic_module.fit_centroid_or_warn_out_of_frame(image, [75.0, 20.0], 1)

    assert np.all(np.isnan(result))
    assert out_of_frame_warnings == [1]


def test_skybg_phot_returns_nan_when_annulus_box_is_empty(monkeypatch):
    image = np.zeros((20, 20), dtype=float)
    sky_warnings = []

    monkeypatch.setattr(
        exotic_module,
        "mesh_box",
        lambda *args, **kwargs: (
            np.empty((0, 0), dtype=int),
            np.empty((0, 0), dtype=int),
        ),
    )
    monkeypatch.setattr(
        exotic_module.plateStatus,
        "skyBackgroundWarning",
        lambda star_index, xc, yc: sky_warnings.append((star_index, xc, yc)),
    )

    bgflux, sigmabg, nbg = exotic_module.skybg_phot(image, 0, 30.0, 30.0)

    assert np.isnan(bgflux)
    assert np.isnan(sigmabg)
    assert nbg == 0
    assert sky_warnings == [(0, 30.0, 30.0)]


def test_check_target_pixel_wcs_keeps_input_coords_when_wcs_target_is_off_frame(monkeypatch):
    image = np.zeros((100, 120), dtype=float)
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [60.0, 50.0]
    wcs.wcs.crval = [210.0, 54.0]
    wcs.wcs.cdelt = np.array([-0.01, 0.01])
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    header = wcs.to_header()
    header["NAXIS"] = 2
    header["NAXIS1"] = 120
    header["NAXIS2"] = 100

    ra_list, dec_list = exotic_module.get_ra_dec(header)

    monkeypatch.setattr(
        exotic_module,
        "update_coordinates_with_proper_motion",
        lambda info_dict, obs_time: (212.0, 54.0),
    )
    monkeypatch.setattr(
        exotic_module,
        "get_psf_parameters",
        lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("centroiding should be skipped")),
    )

    x_pixel, y_pixel = exotic_module.check_target_pixel_wcs(
        25.0,
        30.0,
        {"ra": 210.0, "dec": 54.0, "dist": 0.0, "pm_ra": 0.0, "pm_dec": 0.0},
        ra_list,
        dec_list,
        image,
        2461100.5,
        non_interactive_run=True,
        wcs_header=header,
    )

    assert x_pixel == 25.0
    assert y_pixel == 30.0


def test_get_ra_dec_uses_image_shape_when_header_lacks_naxis():
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [60.0, 50.0]
    wcs.wcs.crval = [210.0, 54.0]
    wcs.wcs.cdelt = np.array([-0.01, 0.01])
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    ra_list, dec_list = exotic_module.get_ra_dec(wcs.to_header(), image_shape=(100, 120))

    assert ra_list.shape == (100, 120)
    assert dec_list.shape == (100, 120)


def test_get_first_image_header_skips_empty_primary_hdu(tmp_path):
    wcs_path = _write_extension_wcs_fits(tmp_path)

    header = exotic_module.get_first_image_header(wcs_path)

    assert header["NAXIS1"] == 120
    assert header["NAXIS2"] == 100
    assert header["CTYPE1"] == "RA---TAN"


def test_get_img_scale_uses_first_image_extension_header(tmp_path):
    wcs_path = _write_extension_wcs_fits(tmp_path)

    img_scale_str, img_scale = exotic_module.get_img_scale(fits.Header(), wcs_path, None)

    assert img_scale_str == "Image scale in arcsec/pixel: 36.0"
    assert img_scale == 36.0


def test_should_ignore_header_wcs_defaults_to_false():
    assert exotic_module.should_ignore_header_wcs(None) is False
    assert exotic_module.should_ignore_header_wcs("n") is False
    assert exotic_module.should_ignore_header_wcs("y") is True


def test_get_bad_wcs_threshold_fraction_defaults_to_three_percent():
    assert exotic_module.get_bad_wcs_threshold_fraction(None) == pytest.approx(0.03)
    assert exotic_module.get_bad_wcs_threshold_fraction("") == pytest.approx(0.03)


def test_get_bad_wcs_threshold_fraction_reads_numeric_percent_values():
    assert exotic_module.get_bad_wcs_threshold_fraction(5.5) == pytest.approx(0.055)
    assert exotic_module.get_bad_wcs_threshold_fraction("7.25") == pytest.approx(0.0725)
    assert exotic_module.get_bad_wcs_threshold_fraction("4%") == pytest.approx(0.04)


def test_get_bad_wcs_threshold_fraction_falls_back_for_invalid_values():
    assert exotic_module.get_bad_wcs_threshold_fraction("not-a-number") == pytest.approx(0.03)
    assert exotic_module.get_bad_wcs_threshold_fraction(-1) == pytest.approx(0.03)
    assert exotic_module.get_bad_wcs_threshold_fraction(101) == pytest.approx(0.03)


def test_display_filename_returns_basename_for_unix_and_windows_paths():
    assert (
        exotic_module._display_filename(
            "/content/drive/MyDrive/0.Exoplanets/2.Transits/run/frame_001.fits.fz"
        )
        == "frame_001.fits.fz"
    )
    assert exotic_module._display_filename(r"C:\data\run\frame_002.fits.fz") == "frame_002.fits.fz"


def test_format_plate_solution_reference_uses_basename_only():
    assert (
        exotic_module.format_plate_solution_reference(
            "/mnt/md0/ftp/user_data/psyfitz/DATA_INBOX/Z.good.TOI 2969 b_2026-03-05_ECO1/"
            "NxAst-TOI2969b_rp_2461105d05262731_20260305_1a016_30_eco1.fits.fz"
        )
        == "Here is the filename where we got the WCS from: "
        "NxAst-TOI2969b_rp_2461105d05262731_20260305_1a016_30_eco1.fits.fz"
    )


def test_log_finding_transformation_progress_prints_basename(monkeypatch):
    stdout = io.StringIO()
    debug_messages = []

    monkeypatch.setattr(exotic_module.sys, "stdout", stdout)
    monkeypatch.setattr(exotic_module.log, "debug", lambda message: debug_messages.append(message))

    exotic_module.log_finding_transformation_progress(
        144,
        220,
        "/content/drive/MyDrive/0.Exoplanets/2.Transits/run/frame_145.fits.fz",
        False,
    )

    assert stdout.getvalue() == "Finding transformation 145 of 220 : frame_145.fits.fz\n"
    assert debug_messages == ["Finding transformation 145 of 220 : frame_145.fits.fz\n"]


def test_check_wcs_ignores_header_wcs_when_override_enabled(monkeypatch):
    monkeypatch.setattr(
        exotic_module,
        "search_wcs",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("header WCS should be ignored")),
    )

    wcs_file = exotic_module.check_wcs(
        "frame.fits",
        ".",
        "n",
        ignore_header_wcs=True,
    )

    assert wcs_file is None


def test_check_wcs_keeps_plate_solution_when_override_enabled(monkeypatch):
    monkeypatch.setattr(exotic_module, "get_wcs", lambda *_args, **_kwargs: "solved_wcs.fits")
    monkeypatch.setattr(
        exotic_module,
        "search_wcs",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("header WCS should not be consulted")),
    )

    wcs_file = exotic_module.check_wcs(
        "frame.fits",
        ".",
        "y",
        ignore_header_wcs=True,
    )

    assert wcs_file == "solved_wcs.fits"


def test_check_wcs_prefers_header_wcs_over_plate_solution(monkeypatch):
    monkeypatch.setattr(
        exotic_module,
        "search_wcs",
        lambda *_args, **_kwargs: types.SimpleNamespace(is_celestial=True),
    )
    monkeypatch.setattr(
        exotic_module,
        "get_wcs",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("plate solution should be skipped when header WCS exists")
        ),
    )

    wcs_file = exotic_module.check_wcs(
        "frame.fits",
        ".",
        "y",
    )

    assert wcs_file == "frame.fits"


def test_should_log_plate_solution_path_suppresses_posix_tmp_paths():
    assert exotic_module.should_log_plate_solution_path("/tmp/tmp42fmzrv8/temp/wcs.fits") is False


def test_should_log_plate_solution_path_keeps_non_tmp_paths():
    assert exotic_module.should_log_plate_solution_path("/data/session/temp/wcs.fits") is True
    assert exotic_module.should_log_plate_solution_path("/tmp_backup/temp/wcs.fits") is True


def test_should_use_multiprocess_transform_precompute_respects_header_wcs_override():
    assert exotic_module.should_use_multiprocess_transform_precompute(
        ["frame1.fits", "frame2.fits"],
        requested_processes=2,
        ignore_header_wcs=True,
    ) is True


def test_filter_sparse_missing_wcs_frames_drops_files_below_three_percent(monkeypatch):
    frames = [f"frame_{i}.fits" for i in range(34)]
    missing_frame = frames[7]

    monkeypatch.setattr(exotic_module, "get_first_image_header", lambda file_name: str(file_name))
    monkeypatch.setattr(
        exotic_module,
        "search_wcs_from_header",
        lambda header: types.SimpleNamespace(is_celestial=header != missing_frame),
    )

    filtered, keep_mask, dropped = exotic_module.filter_sparse_missing_wcs_frames(frames)

    assert filtered.tolist() == [frame for frame in frames if frame != missing_frame]
    assert keep_mask.tolist() == [frame != missing_frame for frame in frames]
    assert dropped == [missing_frame]


def test_filter_sparse_missing_wcs_frames_keeps_files_at_three_percent_or_higher(monkeypatch):
    frames = [f"frame_{i}.fits" for i in range(33)]
    missing_frame = frames[5]

    monkeypatch.setattr(exotic_module, "get_first_image_header", lambda file_name: str(file_name))
    monkeypatch.setattr(
        exotic_module,
        "search_wcs_from_header",
        lambda header: types.SimpleNamespace(is_celestial=header != missing_frame),
    )

    filtered, keep_mask, dropped = exotic_module.filter_sparse_missing_wcs_frames(frames)

    assert filtered.tolist() == frames
    assert keep_mask.tolist() == [True] * len(frames)
    assert dropped == []
