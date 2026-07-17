import io
import sys
import types
import importlib.util
import threading
from concurrent.futures import ThreadPoolExecutor

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
    _install_stub_module("photutils.aperture", CircularAperture=object, CircularAnnulus=object)
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

    count_image = fits.getdata(tmp_path / "working_artifacts" / "BadPixelDetectionCounts.fits")
    mask_image = fits.getdata(tmp_path / "working_artifacts" / "BadPixelMask.fits").astype(bool)

    assert count_image[2, 3] == 4
    assert count_image[6, 5] == 3
    assert mask_image[2, 3]
    assert not mask_image[6, 5]


def test_build_persistent_bad_pixel_map_can_scan_with_multiprocessing(tmp_path, monkeypatch):
    paths = []
    for frame_index in range(10):
        frame = np.full((9, 9), 100.0, dtype=float)
        if frame_index < 4:
            frame[2, 3] = 4000.0
        path = tmp_path / f"frame_{frame_index}.fits"
        fits.writeto(path, frame, overwrite=True)
        paths.append(path)

    captured = {}

    class FakeFuture:
        def __init__(self, value):
            self._value = value

        def result(self):
            return self._value

    class FakeExecutor:
        def __init__(self, max_workers, initializer=None, initargs=()):
            captured["max_workers"] = max_workers
            if initializer is not None:
                initializer(*initargs)

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, traceback):
            return False

        def submit(self, fn, task):
            return FakeFuture(fn(task))

    monkeypatch.setattr(exotic_module, "ProcessPoolExecutor", FakeExecutor)
    monkeypatch.setattr(exotic_module, "as_completed", lambda futures: futures)

    reference = exotic_module.build_persistent_bad_pixel_map(
        paths,
        exotic_module.load_image_data,
        save_directory=tmp_path,
        max_processes=2,
    )

    assert captured["max_workers"] == 2
    assert reference is not None
    assert reference["required_count"] == 4
    assert reference["mask"][2, 3]


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


def test_fit_centroid_full_mode_preserves_psf_subpixel_solution():
    rng = np.random.default_rng(7)
    true_center = (40.3, 35.7)
    image = _gaussian_image(
        center=true_center,
        amplitude=120.0,
        sigma=0.8,
        background=1000.0,
    )
    image += rng.normal(0.0, 20.0, size=image.shape)

    full_result = exotic_module.fit_centroid(image, [40.0, 36.0], 0, fast_mode=False)
    psf_result = exotic_module.fit_centroid(
        image,
        [40.0, 36.0],
        0,
        fast_mode=False,
        weightedcenter=False,
    )
    moment_result = exotic_module.fit_centroid(image, [40.0, 36.0], 0, fast_mode=True)

    assert full_result[0] == pytest.approx(psf_result[0], abs=1e-6)
    assert full_result[1] == pytest.approx(psf_result[1], abs=1e-6)

    psf_error = np.hypot(psf_result[0] - true_center[0], psf_result[1] - true_center[1])
    moment_error = np.hypot(moment_result[0] - true_center[0], moment_result[1] - true_center[1])

    assert psf_error < moment_error


def test_fit_psf_photometry_flux_row_preserves_robust_centroid_coordinates():
    image = _gaussian_image(center=(40.3, 35.7), amplitude=180.0, sigma=0.9, background=1000.0)
    centroid_row = np.array([40.1, 35.9, 150.0, 0.8, 0.8, 0.0, 1000.0], dtype=float)

    flux_row = exotic_module.fit_psf_photometry_flux_row(image, centroid_row, 0)

    assert flux_row[0] == pytest.approx(centroid_row[0])
    assert flux_row[1] == pytest.approx(centroid_row[1])
    assert flux_row[2] > 0
    assert 0.5 <= flux_row[3] <= 2.0
    assert 0.5 <= flux_row[4] <= 2.0


def test_fit_centroid_prefers_seed_anchored_solution_in_crowded_field():
    yy, xx = np.mgrid[0:80, 0:80]
    image = np.full((80, 80), 400.0)
    image += 80.0 * np.exp(-((xx - 40.0) ** 2 + (yy - 40.0) ** 2) / (2.0 * 1.0 ** 2))
    image += 220.0 * np.exp(-((xx - 31.5) ** 2 + (yy - 35.0) ** 2) / (2.0 * 1.0 ** 2))

    result = exotic_module.fit_centroid(image, [40.0, 40.0], 0, fast_mode=False)

    assert np.hypot(result[0] - 40.0, result[1] - 40.0) < 1.5
    assert np.hypot(result[0] - 31.5, result[1] - 35.0) > 5.0
    assert result[2] > 0


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

    class _EmptyAnnulusMask:
        data = np.empty((0, 0), dtype=float)

        def cutout(self, *args, **kwargs):
            return None

    class _EmptyCircularAnnulus:
        def __init__(self, *args, **kwargs):
            pass

        def to_mask(self, *args, **kwargs):
            return [_EmptyAnnulusMask()]

    monkeypatch.setattr(exotic_module, "CircularAnnulus", _EmptyCircularAnnulus)
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


@pytest.mark.skipif(not _module_available("photutils.aperture"), reason="requires photutils aperture masks")
def test_skybg_phot_exact_annulus_uses_fractional_pixel_area():
    image = np.ones((80, 80), dtype=float)

    bgflux, sigmabg, nbg = exotic_module.skybg_phot(image, 0, 40.3, 35.7, r=3.0, dr=2.0, fast_mode=False)

    assert bgflux == pytest.approx(1.0, abs=1e-8)
    assert sigmabg == pytest.approx(0.0, abs=1e-8)
    assert nbg == pytest.approx(np.pi * ((3.0 + 2.0) ** 2 - 3.0 ** 2), rel=1e-3)
    assert not np.isclose(nbg, round(nbg), atol=1e-6)


@pytest.mark.skipif(not _module_available("photutils.aperture"), reason="requires photutils aperture masks")
def test_skybg_phot_high_side_clipping_rejects_hot_pixel():
    image = np.full((80, 80), 100.0, dtype=float)
    image[40, 55] = 10000.0

    bgflux, sigmabg, nbg = exotic_module.skybg_phot(image, 0, 40.0, 40.0, r=10.0, dr=10.0, fast_mode=False)

    assert bgflux == pytest.approx(100.0, abs=1e-8)
    assert sigmabg == pytest.approx(0.0, abs=1e-8)
    assert nbg > 250.0


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


def test_get_ra_dec_matches_zero_based_astropy_pixel_coordinates():
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

    for x_pixel, y_pixel in [(0, 0), (59, 49), (119, 99), (23, 71)]:
        expected_ra, expected_dec = wcs.pixel_to_world_values(x_pixel, y_pixel)
        assert ra_list[y_pixel, x_pixel] == pytest.approx(expected_ra, abs=1.0e-12)
        assert dec_list[y_pixel, x_pixel] == pytest.approx(expected_dec, abs=1.0e-12)


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


def test_get_pointing_rejection_sigma_defaults_to_disabled():
    assert exotic_module.get_pointing_rejection_sigma(None) is None
    assert exotic_module.get_pointing_rejection_sigma("") is None


def test_get_pointing_rejection_sigma_reads_positive_numeric_values():
    assert exotic_module.get_pointing_rejection_sigma(3) == pytest.approx(3.0)
    assert exotic_module.get_pointing_rejection_sigma("2.75") == pytest.approx(2.75)


def test_get_pointing_rejection_sigma_disables_for_invalid_text_and_zero():
    assert exotic_module.get_pointing_rejection_sigma("not-a-number") is None
    assert exotic_module.get_pointing_rejection_sigma(-1) is None
    assert exotic_module.get_pointing_rejection_sigma(0) is None
    assert exotic_module.get_pointing_rejection_sigma("off") is None


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


def test_log_alignment_progress_prints_basename(monkeypatch):
    stdout = io.StringIO()
    debug_messages = []

    monkeypatch.setattr(exotic_module.sys, "stdout", stdout)
    monkeypatch.setattr(exotic_module.log, "debug", lambda message: debug_messages.append(message))

    exotic_module.log_alignment_progress(
        144,
        220,
        "/content/drive/MyDrive/0.Exoplanets/2.Transits/run/frame_145.fits.fz",
        False,
    )

    assert stdout.getvalue() == "Aligning frame 145 of 220 : frame_145.fits.fz\n"
    assert debug_messages == ["Aligning frame 145 of 220 : frame_145.fits.fz\n"]


def test_collect_transform_frame_pointings_logs_alignment_progress(monkeypatch):
    progress_messages = []

    monkeypatch.setattr(
        exotic_module,
        "log_info",
        lambda message, warn=False, error=False: progress_messages.append((message, warn, error)),
    )
    monkeypatch.setattr(
        exotic_module,
        "transformation",
        lambda image_data, file_name, report_failure=False, reference_image=None: (
            lambda anchor: np.asarray(anchor, dtype=float)
        ),
    )

    frames = ["frame_0001.fits", "frame_0002.fits", "frame_0003.fits"]
    frame_loader = lambda file_name: np.ones((8, 8), dtype=float)

    positions, usable_mask = exotic_module.collect_transform_frame_pointings(frames, frame_loader=frame_loader)

    assert positions.shape == (3, 2)
    assert usable_mask.tolist() == [True, True, True]
    assert [message for message, _, _ in progress_messages] == [
        "Pointing precheck alignment progress: file 1 of 3 : frame_0001.fits",
        "Pointing precheck alignment progress: file 2 of 3 : frame_0002.fits",
        "Pointing precheck alignment progress: file 3 of 3 : frame_0003.fits",
    ]


def test_collect_transform_frame_pointings_can_return_transform_cache(monkeypatch):
    monkeypatch.setattr(
        exotic_module,
        "log_info",
        lambda *_args, **_kwargs: None,
    )

    expected_transform = exotic_module.SimilarityTransform(scale=1, rotation=0, translation=[2.0, -1.0])
    monkeypatch.setattr(
        exotic_module,
        "transformation",
        lambda image_data, file_name, report_failure=False, reference_image=None: expected_transform,
    )

    frames = ["frame_0001.fits", "frame_0002.fits"]
    frame_loader = lambda file_name: np.ones((8, 8), dtype=float)

    positions, usable_mask, transforms = exotic_module.collect_transform_frame_pointings(
        frames,
        frame_loader=frame_loader,
        return_transforms=True,
    )

    assert usable_mask.tolist() == [True, True]
    assert np.allclose(positions[1], [5.5, 2.5])
    assert set(transforms) == set(frames)
    assert transforms[frames[1]] is expected_transform


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


def test_transformation_pool_initializer_suppresses_inherited_tk_cleanup(monkeypatch):
    calls = []
    reference_image = np.ones((4, 4), dtype=float)
    monkeypatch.setattr(exotic_module, "suppress_inherited_tk_cleanup_in_worker", lambda: calls.append(True))
    monkeypatch.setattr(exotic_module, "load_image_data", lambda _file_name: reference_image)

    exotic_module._TRANSFORM_REFERENCE_IMAGE = None
    exotic_module._TRANSFORM_REFERENCE_CACHE = {"stale": True}

    exotic_module._transformation_pool_initializer("reference.fits")

    assert calls == [True]
    assert exotic_module._TRANSFORM_REFERENCE_IMAGE is reference_image
    assert exotic_module._TRANSFORM_REFERENCE_CACHE is None


def test_apply_parallel_alignment_result_uses_precomputed_fallback_when_wcs_geometry_fails(monkeypatch):
    psf_data = {
        "target": np.zeros((2, 7), dtype=float),
        "comp1": np.zeros((2, 7), dtype=float),
    }
    psf_data["target"][0] = np.array([10.0, 10.0, 100.0, 2.0, 2.0, 0.0, 50.0])
    psf_data["comp1"][0] = np.array([20.0, 10.0, 100.0, 2.0, 2.0, 0.0, 50.0])
    tar_comp_dist = {"comp1": np.array([10, 0], dtype=int)}
    warnings = []

    monkeypatch.setattr(
        exotic_module.plateStatus,
        "lowFluxAmplitudeWarning",
        lambda star_index, xc, yc: warnings.append((star_index, xc, yc)),
    )

    result = {
        "index": 1,
        "file_name": "frame_0002.fits",
        "wcs": {
            "projected_off_frame": False,
            "psf_rows": {
                "target": np.array([10.0, 10.0, 100.0, 2.0, 2.0, 0.0, 50.0]),
                "comp1": np.array([50.0, 50.0, 100.0, 2.0, 2.0, 0.0, 50.0]),
            },
            "warnings": [("low_flux", 1, 50.0, 50.0)],
        },
        "fallback": {
            "psf_rows": {
                "target": np.array([11.0, 10.0, 100.0, 2.0, 2.0, 0.0, 50.0]),
                "comp1": np.array([21.0, 10.0, 90.0, 2.0, 2.0, 0.0, 50.0]),
            },
            "warnings": [("low_flux", 1, 21.0, 10.0)],
        },
    }

    selected = exotic_module.apply_parallel_alignment_result(
        result,
        1,
        psf_data,
        tar_comp_dist,
        ["comp1"],
    )

    assert selected == "fallback"
    assert psf_data["target"][1, 0] == pytest.approx(11.0)
    assert psf_data["comp1"][1, 0] == pytest.approx(21.0)
    assert warnings == [(1, 21.0, 10.0)]


def test_parallel_alignment_task_uses_precomputed_fallback_transform(monkeypatch):
    target_and_comp_pixels = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=float)
    precomputed_transform = exotic_module.SimilarityTransform(
        scale=1,
        rotation=0,
        translation=[5.0, -1.0],
    )

    monkeypatch.setattr(
        exotic_module,
        "_load_alignment_worker_frame",
        lambda _file_name: ({}, np.ones((10, 10), dtype=float)),
    )
    monkeypatch.setattr(
        exotic_module,
        "transformation",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("cached transform should be used")
        ),
    )
    monkeypatch.setattr(
        exotic_module,
        "_fit_alignment_candidate_psfs",
        lambda image_data, predicted_coords, *_args: {
            "coords": np.asarray(predicted_coords, dtype=float),
            "psf_rows": {"target": np.zeros(7, dtype=float)},
            "warnings": [],
        },
    )

    result = exotic_module._parallel_alignment_task((
        1,
        "frame_0002.fits",
        target_and_comp_pixels,
        None,
        True,
        False,
        False,
        True,
        False,
        precomputed_transform,
    ))

    assert np.allclose(result["fallback"]["coords"], [[6.0, 1.0], [8.0, 3.0]])


def test_classify_wcs_fallback_frames_queues_only_missing_or_rejected_candidates():
    def candidate(target_xy, comp_xy):
        coords = np.array([target_xy, comp_xy], dtype=float)
        return {
            "coords": coords,
            "projected_off_frame": False,
            "psf_rows": {
                "target": np.array([*target_xy, 100.0, 2.0, 2.0, 0.0, 50.0]),
                "comp1": np.array([*comp_xy, 90.0, 2.0, 2.0, 0.0, 50.0]),
            },
            "warnings": [],
        }

    results = [
        {"index": 0, "file_name": "frame0.fits", "wcs": candidate((10.0, 10.0), (20.0, 10.0)), "fallback": None},
        {"index": 1, "file_name": "frame1.fits", "wcs": candidate((11.0, 10.0), (50.0, 50.0)), "fallback": None},
        {"index": 2, "file_name": "frame2.fits", "wcs": None, "fallback": None},
        {"index": 3, "file_name": "frame3.fits", "wcs": candidate((13.0, 10.0), (23.0, 10.0)), "fallback": None},
    ]

    missing, rejected = exotic_module.classify_wcs_fallback_frames(
        results,
        np.array([[10.0, 10.0], [20.0, 10.0]]),
    )

    assert missing == [2]
    assert rejected == [1]


def test_build_multiprocess_alignment_results_runs_legacy_batch_only_for_wcs_failures(monkeypatch):
    def candidate(target_xy, comp_xy):
        coords = np.array([target_xy, comp_xy], dtype=float)
        return {
            "coords": coords,
            "projected_off_frame": False,
            "psf_rows": {
                "target": np.array([*target_xy, 100.0, 2.0, 2.0, 0.0, 50.0]),
                "comp1": np.array([*comp_xy, 90.0, 2.0, 2.0, 0.0, 50.0]),
            },
            "warnings": [],
        }

    batches = []

    def fake_run_batch(tasks, *_args, **_kwargs):
        batches.append(tasks)
        if len(batches) == 1:
            return {
                0: {"index": 0, "file_name": "frame0.fits", "wcs": candidate((10.0, 10.0), (20.0, 10.0)), "fallback": None},
                1: {"index": 1, "file_name": "frame1.fits", "wcs": candidate((11.0, 10.0), (21.0, 10.0)), "fallback": None},
                2: {"index": 2, "file_name": "frame2.fits", "wcs": candidate((12.0, 10.0), (50.0, 50.0)), "fallback": None},
                3: {"index": 3, "file_name": "frame3.fits", "wcs": None, "fallback": None},
            }

        return {
            task[0]: {
                "index": task[0],
                "file_name": task[1],
                "wcs": None,
                "fallback": candidate((10.0 + task[0], 10.0), (20.0 + task[0], 10.0)),
            }
            for task in tasks
        }

    monkeypatch.setattr(exotic_module, "_run_multiprocess_alignment_task_batch", fake_run_batch)

    results = exotic_module.build_multiprocess_alignment_results(
        np.array(["frame0.fits", "frame1.fits", "frame2.fits", "frame3.fits"]),
        4,
        np.array([[10.0, 10.0], [20.0, 10.0]]),
        target_and_comp_radec=np.array([[1.0, 2.0], [1.1, 2.1]]),
        compute_fallback_transform=True,
    )

    assert [task[0] for task in batches[0]] == [0, 1, 2, 3]
    assert all(task[7] is False for task in batches[0])
    assert [task[0] for task in batches[1]] == [2, 3]
    assert all(task[4] is True and task[7] is True for task in batches[1])
    assert results[0]["fallback"] is None
    assert results[1]["fallback"] is None
    assert results[2]["fallback"] is not None
    assert results[3]["fallback"] is not None


def test_downsampled_fallback_transformation_restores_full_resolution_translation(monkeypatch):
    calls = []

    def fake_transformation(image_data, _file_name, **kwargs):
        calls.append((image_data.shape, kwargs["reference_image"].shape))
        return exotic_module.SimilarityTransform(
            scale=1.01,
            rotation=0.02,
            translation=[2.0, -3.0],
        )

    monkeypatch.setattr(exotic_module, "transformation", fake_transformation)
    image = np.ones((8, 12), dtype=float)
    reference = np.ones((8, 12), dtype=float)

    result = exotic_module.downsampled_fallback_transformation(
        image,
        "frame.fits",
        reference_image=reference,
        max_dimension=6,
    )

    assert calls == [((4, 6), (4, 6))]
    assert result.scale == pytest.approx(1.01)
    assert result.rotation == pytest.approx(0.02)
    assert np.allclose(result.translation, [4.0, -6.0])


def test_fit_alignment_candidate_psfs_serializes_plate_status_swap(monkeypatch):
    sentinel_status = types.SimpleNamespace(name="original-plate-status")
    started = threading.Event()

    def fake_fit_centroid(_data, pos, starIndex, **_kwargs):
        started.set()
        return np.array([float(starIndex), float(pos[0]), float(pos[1])], dtype=float)

    monkeypatch.setattr(exotic_module, "plateStatus", sentinel_status)
    monkeypatch.setattr(exotic_module, "fit_centroid_or_warn_out_of_frame", fake_fit_centroid)

    lock = exotic_module._PLATE_STATUS_SWAP_LOCK
    lock.acquire()
    executor = ThreadPoolExecutor(max_workers=1)
    future = None
    try:
        future = executor.submit(
            exotic_module._fit_alignment_candidate_psfs,
            np.ones((5, 5), dtype=float),
            np.array([[1.0, 1.0], [2.0, 2.0]], dtype=float),
            False,
            False,
        )
        assert not started.wait(0.2)
    finally:
        lock.release()

    try:
        result = future.result(timeout=2)
    finally:
        executor.shutdown(wait=True)

    assert started.wait(0.2)
    assert result["psf_rows"]["target"].tolist() == [0.0, 1.0, 1.0]
    assert result["psf_rows"]["comp1"].tolist() == [1.0, 2.0, 2.0]
    assert exotic_module.plateStatus is sentinel_status


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


def test_filter_wcs_target_out_of_frame_frames_drops_only_projected_misses(monkeypatch):
    def make_wcs_header(center_ra):
        wcs = WCS(naxis=2)
        wcs.wcs.crpix = [60.0, 50.0]
        wcs.wcs.crval = [center_ra, 54.0]
        wcs.wcs.cdelt = np.array([-0.01, 0.01])
        wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        header = wcs.to_header()
        header["NAXIS"] = 2
        header["NAXIS1"] = 120
        header["NAXIS2"] = 100
        return header

    no_wcs_header = fits.Header()
    no_wcs_header["NAXIS"] = 2
    no_wcs_header["NAXIS1"] = 120
    no_wcs_header["NAXIS2"] = 100

    headers = {
        "target_in_frame.fits": make_wcs_header(210.0),
        "target_off_frame.fits": make_wcs_header(212.0),
        "no_wcs.fits": no_wcs_header,
    }
    messages = []

    monkeypatch.setattr(exotic_module, "get_first_image_header", lambda file_name: headers[file_name])
    monkeypatch.setattr(
        exotic_module,
        "update_coordinates_with_proper_motion",
        lambda info_dict, obs_time: (210.0, 54.0),
    )
    monkeypatch.setattr(
        exotic_module,
        "log_info",
        lambda message, warn=False, error=False: messages.append((message, warn, error)),
    )

    frames = list(headers)
    filtered, keep_mask, dropped = exotic_module.filter_wcs_target_out_of_frame_frames(
        frames,
        {"ra": 210.0, "dec": 54.0},
        obs_times=[2461196.5, 2461196.6, 2461196.7],
    )

    assert filtered.tolist() == ["target_in_frame.fits", "no_wcs.fits"]
    assert keep_mask.tolist() == [True, False, True]
    assert dropped == ["target_off_frame.fits"]
    assert any("Target WCS precheck" in message for message, _, _ in messages)


def test_maybe_reinterpret_decimal_ra_hours_from_wcs_when_only_ra_times_fifteen_matches(monkeypatch):
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [60.0, 50.0]
    wcs.wcs.crval = [16.18494, 74.3313]
    wcs.wcs.cdelt = np.array([-0.01, 0.01])
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    header = wcs.to_header()
    header["NAXIS"] = 2
    header["NAXIS1"] = 120
    header["NAXIS2"] = 100
    messages = []

    monkeypatch.setattr(exotic_module, "get_first_image_header", lambda _file_name: header)
    monkeypatch.setattr(
        exotic_module,
        "log_info",
        lambda message, warn=False, error=False: messages.append((message, warn, error)),
    )

    info = {"ra": 1.078996153, "dec": 74.3313055}

    corrected = exotic_module.maybe_reinterpret_decimal_ra_hours_from_wcs(["frame.fits"], info)

    assert corrected is True
    assert info["ra"] == pytest.approx(16.184942295)
    assert any("interpreted decimal target RA as hours" in message and warn for message, warn, _ in messages)

    already_degrees = {"ra": 16.184942295, "dec": 74.3313055}

    corrected_again = exotic_module.maybe_reinterpret_decimal_ra_hours_from_wcs(
        ["frame.fits"],
        already_degrees,
    )

    assert corrected_again is False
    assert already_degrees["ra"] == pytest.approx(16.184942295)


def test_filter_pointing_outlier_frames_uses_wcs_when_all_frames_have_wcs(monkeypatch):
    frames = [f"frame_{i}.fits" for i in range(6)]
    wcs_positions = np.array(
        [
            [100.0, 100.0],
            [101.0, 100.0],
            [99.0, 100.0],
            [100.0, 101.0],
            [100.0, 99.0],
            [0.0, 0.0],
        ],
        dtype=float,
    )

    monkeypatch.setattr(
        exotic_module,
        "collect_wcs_frame_center_pointings",
        lambda inputfiles: (wcs_positions, np.ones(len(inputfiles), dtype=bool)),
    )
    monkeypatch.setattr(
        exotic_module,
        "collect_transform_frame_pointings",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("transform fallback should not be used")),
    )

    filtered, keep_mask, dropped = exotic_module.filter_pointing_outlier_frames(
        frames,
        pointing_rejection_sigma=3.0,
    )

    assert filtered.tolist() == frames[:-1]
    assert keep_mask.tolist() == [True, True, True, True, True, False]
    assert dropped == [frames[-1]]


def test_filter_pointing_outlier_frames_falls_back_to_transform_when_wcs_is_incomplete(monkeypatch):
    frames = [f"frame_{i}.fits" for i in range(6)]
    transform_positions = np.array(
        [
            [50.0, 50.0],
            [50.5, 49.5],
            [49.5, 50.5],
            [50.0, 51.0],
            [50.0, 49.0],
            [10.0, 10.0],
        ],
        dtype=float,
    )
    transform_calls = []

    monkeypatch.setattr(
        exotic_module,
        "collect_wcs_frame_center_pointings",
        lambda inputfiles: (
            np.full((len(inputfiles), 2), np.nan, dtype=float),
            np.array([True, True, True, True, False, False], dtype=bool),
        ),
    )

    def fake_collect_transform_frame_pointings(inputfiles, frame_loader=None, return_transforms=False, **kwargs):
        transform_calls.append((tuple(inputfiles), frame_loader, return_transforms, kwargs))
        transforms = {
            str(file_name): exotic_module.SimilarityTransform(scale=1, rotation=0, translation=[index, 0])
            for index, file_name in enumerate(inputfiles)
        }
        if return_transforms:
            return transform_positions, np.ones(len(inputfiles), dtype=bool), transforms
        return transform_positions, np.ones(len(inputfiles), dtype=bool)

    monkeypatch.setattr(exotic_module, "collect_transform_frame_pointings", fake_collect_transform_frame_pointings)

    filtered, keep_mask, dropped, cached_transforms = exotic_module.filter_pointing_outlier_frames(
        frames,
        pointing_rejection_sigma=3.0,
        return_alignment_transforms=True,
    )

    assert len(transform_calls) == 1
    assert transform_calls[0][2] is True
    assert filtered.tolist() == frames[:-1]
    assert keep_mask.tolist() == [True, True, True, True, True, False]
    assert dropped == [frames[-1]]
    assert set(cached_transforms) == set(frames[:-1])


def test_reference_frame_rejection_fallback_reports_automatic_removal_and_reprojection(monkeypatch):
    messages = []

    monkeypatch.setattr(
        exotic_module,
        "log_info",
        lambda message, error=False, warn=False: messages.append((message, error, warn)),
    )

    result = exotic_module.reference_frame_rejection_fallback_info(
        "frame_0001.fits",
        ["frame_0001.fits", "frame_0002.fits", "frame_0003.fits"],
        ordered_inputfiles=[
            "frame_0001.fits",
            "frame_0002.fits",
            "frame_0003.fits",
            "frame_0004.fits",
        ],
    )

    assert result["leading_rejected_files"] == ["frame_0001.fits", "frame_0002.fits", "frame_0003.fits"]
    assert result["next_reference_candidate"] == "frame_0004.fits"
    assert any("automatically removing" in message and warn for message, _, warn in messages)
    assert any("target RA/Dec" in message and "nextastro_archive" in message and warn for message, _, warn in messages)
    assert any(
        "Automatically removed leading rejected frame(s)" in message
        and "frame_0001.fits, frame_0002.fits, frame_0003.fits" in message
        and "Continuing from new reference image frame_0004.fits" in message
        and warn
        for message, _, warn in messages
    )


def test_reference_frame_rejection_fallback_only_reports_consecutive_leading_rejections(monkeypatch):
    messages = []

    monkeypatch.setattr(
        exotic_module,
        "log_info",
        lambda message, error=False, warn=False: messages.append((message, error, warn)),
    )

    result = exotic_module.reference_frame_rejection_fallback_info(
        "frame_0001.fits",
        ["frame_0001.fits", "frame_0003.fits"],
        ordered_inputfiles=[
            "frame_0001.fits",
            "frame_0002.fits",
            "frame_0003.fits",
            "frame_0004.fits",
        ],
    )

    assert result["leading_rejected_files"] == ["frame_0001.fits"]
    assert result["next_reference_candidate"] == "frame_0002.fits"
    assert any(
        "Automatically removed leading rejected frame(s)" in message
        and "frame_0001.fits" in message
        and "Continuing from new reference image frame_0002.fits" in message
        and "frame_0003.fits" not in message
        and warn
        for message, _, warn in messages
    )


def test_reference_frame_rejection_fallback_ignores_non_reference_rejections(monkeypatch):
    messages = []

    monkeypatch.setattr(
        exotic_module,
        "log_info",
        lambda message, error=False, warn=False: messages.append((message, error, warn)),
    )

    result = exotic_module.reference_frame_rejection_fallback_info(
        "frame_0001.fits",
        ["frame_0002.fits", "frame_0003.fits"],
    )

    assert result is None
    assert messages == []


def test_reference_fallback_comparison_stars_use_nextastro_archive_image_criteria():
    image = np.zeros((300, 300), dtype=float)

    def add_blob(x_pos, y_pos, value):
        image[y_pos - 1:y_pos + 2, x_pos - 1:x_pos + 2] = value * 0.5
        image[y_pos, x_pos] = value

    add_blob(150, 150, 2000.0)  # target location, excluded by detected-target match
    add_blob(220, 220, 1200.0)
    add_blob(80, 80, 900.0)
    add_blob(180, 180, 1600.0)  # within 50 px of target, excluded
    add_blob(25, 25, 5000.0)    # outside the central 50% frame, excluded

    comp_stars, candidates = exotic_module.select_reference_fallback_comparison_stars(
        image,
        image.shape,
        target_pixel=[150, 150],
        comp_count=2,
    )

    assert comp_stars == [[220.0, 220.0], [80.0, 80.0]]
    assert [candidate["flux"] for candidate in candidates] == sorted(
        [candidate["flux"] for candidate in candidates],
        reverse=True,
    )


def test_automatic_optimal_calibration_selector_filters_flux_and_ranks_color(monkeypatch):
    image = np.zeros((300, 300), dtype=float)

    def add_blob(x_pos, y_pos, value):
        image[y_pos - 1:y_pos + 2, x_pos - 1:x_pos + 2] = value * 0.5
        image[y_pos, x_pos] = value

    add_blob(150, 150, 2000.0)
    add_blob(220, 220, 1700.0)
    add_blob(80, 80, 1800.0)
    add_blob(230, 80, 6000.0)

    ra_wcs = np.tile(np.arange(300, dtype=float), (300, 1))
    dec_wcs = np.tile(np.arange(300, dtype=float)[:, None], (1, 300))
    catalog = {"rows": []}

    def fake_color_match(_catalog, ra, dec, obs_filter, max_separation_arcsec=5.0):
        colors = {
            (150, 150): (12.0, 11.4),
            (220, 220): (13.0, 12.41),
            (80, 80): (13.0, 12.0),
            (230, 80): (10.0, 9.4),
        }
        key = (int(round(float(ra))), int(round(float(dec))))
        if key not in colors:
            return None
        b_mag, v_mag = colors[key]
        return {"catalog_row": {"Bmag": b_mag, "Vmag": v_mag}}

    monkeypatch.setattr(exotic_module, "nextastro_catalog_nearest_color_row", fake_color_match)
    monkeypatch.setattr(
        exotic_module,
        "nextastro_photometry_catalog_match",
        lambda catalog_response, ra, dec, obs_filter: (
            {
                **fake_color_match(catalog_response, ra, dec, obs_filter),
                "mag": 12.0,
                "error": 0.01,
                "mag_band": "V",
            }
            if fake_color_match(catalog_response, ra, dec, obs_filter) is not None
            else None
        ),
    )

    comp_stars, candidates = exotic_module.select_automatic_optimal_calibration_stars(
        image,
        image.shape,
        target_pixel=[150, 150],
        ra_wcs=ra_wcs,
        dec_wcs=dec_wcs,
        obs_filter="V",
        field_catalog=catalog,
        count=2,
        colour_term_metadata={
            "term": 0.2,
            "term_error": 0.01,
            "term_index": "B-V",
        },
    )

    assert comp_stars[0] == [220.0, 220.0]
    assert [candidate["color_delta"] for candidate in candidates] == sorted(
        candidate["color_delta"] for candidate in candidates
    )
    assert candidates[0]["expected_colour_mismatch_mag"] == pytest.approx(
        0.2 * candidates[0]["color_delta"]
    )
    assert candidates[0]["colour_term_uncertainty_mag"] == pytest.approx(
        0.01 * candidates[0]["color_delta"]
    )

    brightest_comp_stars, brightest_candidates = exotic_module.select_automatic_optimal_calibration_stars(
        image,
        image.shape,
        target_pixel=[150, 150],
        ra_wcs=ra_wcs,
        dec_wcs=dec_wcs,
        obs_filter="V",
        field_catalog=catalog,
        count=2,
        brightest_first=True,
        saturation_threshold=5000.0,
    )

    assert brightest_comp_stars[0] == [80.0, 80.0]
    assert [candidate["flux"] for candidate in brightest_candidates] == sorted(
        [candidate["flux"] for candidate in brightest_candidates],
        reverse=True,
    )
    assert all(0.5 <= candidate["brightness_ratio"] <= 2.0 for candidate in candidates)


def test_build_absolute_comp_ensemble_flux_uses_median_normalized_members():
    comp_flux_map = {
        "comp1": np.array([100.0, 102.0, 98.0, 100.0, 101.0, 99.0]),
        "comp2": np.array([200.0, 204.0, 196.0, 200.0, 202.0, 198.0]),
        "comp3": np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]),
    }

    ensemble_flux, member_keys = exotic_module.build_absolute_comp_ensemble_flux(
        comp_flux_map,
        ["comp1", "comp2", "comp3"],
    )

    assert member_keys == ["comp1", "comp2"]
    assert np.nanmedian(ensemble_flux) == pytest.approx(150.0)
    assert ensemble_flux[1] / np.nanmedian(ensemble_flux) == pytest.approx(1.02)
