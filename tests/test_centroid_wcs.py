import io
import sys
import types
import importlib.util

import numpy as np
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


def test_should_ignore_header_wcs_defaults_to_false():
    assert exotic_module.should_ignore_header_wcs(None) is False
    assert exotic_module.should_ignore_header_wcs("n") is False
    assert exotic_module.should_ignore_header_wcs("y") is True


def test_display_filename_returns_basename_for_unix_and_windows_paths():
    assert (
        exotic_module._display_filename(
            "/content/drive/MyDrive/0.Exoplanets/2.Transits/run/frame_001.fits.fz"
        )
        == "frame_001.fits.fz"
    )
    assert exotic_module._display_filename(r"C:\data\run\frame_002.fits.fz") == "frame_002.fits.fz"


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
