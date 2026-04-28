import importlib.util
import sys
import types

import numpy as np
import pytest


def _module_available(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ModuleNotFoundError, ValueError):
        return False


def _set_stub_if_missing(name: str, module: types.ModuleType) -> None:
    if not _module_available(name):
        sys.modules.setdefault(name, module)


fake_barycorrpy = types.ModuleType("barycorrpy")
fake_utc_tdb = types.ModuleType("barycorrpy.utc_tdb")
fake_utc_tdb.JDUTC_to_BJDTDB = lambda *args, **kwargs: None
fake_astroalign = types.ModuleType("astroalign")
fake_astroalign.PIXEL_TOL = 1
fake_astroquery = types.ModuleType("astroquery")
fake_astroquery_simbad = types.ModuleType("astroquery.simbad")
fake_astroquery_simbad.Simbad = type("Simbad", (), {})
fake_astroquery_gaia = types.ModuleType("astroquery.gaia")
fake_astroquery_gaia.Gaia = type("Gaia", (), {})
fake_imreg_dft = types.ModuleType("imreg_dft")
fake_colour_demosaicing = types.ModuleType("colour_demosaicing")
fake_colour_demosaicing.demosaicing_CFA_Bayer_bilinear = lambda *args, **kwargs: None
fake_photutils = types.ModuleType("photutils")
fake_photutils_aperture = types.ModuleType("photutils.aperture")
fake_photutils_aperture.CircularAperture = type("CircularAperture", (), {})
fake_photutils_aperture.CircularAnnulus = type("CircularAnnulus", (), {})
fake_photutils_detection = types.ModuleType("photutils.detection")
fake_photutils_detection.DAOStarFinder = type("DAOStarFinder", (), {})
fake_ldtk = types.ModuleType("ldtk")
fake_ldtk.LDPSet = type("LDPSet", (), {})
fake_ldtk.ldtk = types.SimpleNamespace(LDPSet=fake_ldtk.LDPSet)
fake_ldtk_ldmodel = types.ModuleType("ldtk.ldmodel")
fake_ldtk_ldmodel.LinearModel = type("LinearModel", (), {})
fake_ldtk_ldmodel.QuadraticModel = type("QuadraticModel", (), {})
fake_ldtk_ldmodel.NonlinearModel = type("NonlinearModel", (), {})
fake_lmfit = types.ModuleType("lmfit")
fake_pylightcurve = types.ModuleType("pylightcurve")
fake_pylightcurve_models = types.ModuleType("pylightcurve.models")
fake_pylightcurve_exoplanet = types.ModuleType("pylightcurve.models.exoplanet_lc")
fake_pylightcurve_exoplanet.transit = lambda *args, **kwargs: None
fake_pyvo = types.ModuleType("pyvo")
fake_ultranest = types.ModuleType("ultranest")
fake_ultranest.ReactiveNestedSampler = type("ReactiveNestedSampler", (), {})
fake_elca = types.ModuleType("exotic.api.elca")
fake_elca.lc_fitter = lambda *args, **kwargs: None
fake_elca.binner = lambda *args, **kwargs: None
fake_elca.transit = lambda *args, **kwargs: None
fake_elca.get_phase = lambda *args, **kwargs: None
fake_ld = types.ModuleType("exotic.api.ld")
fake_ld.LimbDarkening = type("LimbDarkening", (), {})
fake_ld.ld_re_punct_p = lambda *args, **kwargs: None

_set_stub_if_missing("astroalign", fake_astroalign)
_set_stub_if_missing("astroquery", fake_astroquery)
_set_stub_if_missing("astroquery.simbad", fake_astroquery_simbad)
_set_stub_if_missing("astroquery.gaia", fake_astroquery_gaia)
_set_stub_if_missing("imreg_dft", fake_imreg_dft)
_set_stub_if_missing("colour_demosaicing", fake_colour_demosaicing)
_set_stub_if_missing("photutils", fake_photutils)
_set_stub_if_missing("photutils.aperture", fake_photutils_aperture)
_set_stub_if_missing("photutils.detection", fake_photutils_detection)
_set_stub_if_missing("ldtk", fake_ldtk)
_set_stub_if_missing("ldtk.ldmodel", fake_ldtk_ldmodel)
_set_stub_if_missing("lmfit", fake_lmfit)
_set_stub_if_missing("pylightcurve", fake_pylightcurve)
_set_stub_if_missing("pylightcurve.models", fake_pylightcurve_models)
_set_stub_if_missing("pylightcurve.models.exoplanet_lc", fake_pylightcurve_exoplanet)
_set_stub_if_missing("pyvo", fake_pyvo)
_set_stub_if_missing("ultranest", fake_ultranest)
_set_stub_if_missing("barycorrpy", fake_barycorrpy)
_set_stub_if_missing("barycorrpy.utc_tdb", fake_utc_tdb)
sys.modules.setdefault("exotic.api.elca", fake_elca)
sys.modules.setdefault("exotic.api.ld", fake_ld)

from exotic.exotic import (  # noqa: E402
    INITIAL_RPRS_BOUND_LOWER_SCALE,
    INITIAL_RPRS_BOUND_UPPER_SCALE,
    RPRS_POSTERIOR_MAX_RETRIES_DEFAULT,
    RPRS_SEARCH_BOUND_MAX,
    RPRS_SEARCH_BOUND_MIN,
    build_single_transit_duration_prior,
    build_initial_rprs_bounds,
    run_nested_lightcurve_fit_with_rprs_posterior_retry,
)


def test_build_initial_rprs_bounds_uses_wider_asymmetric_search_box():
    bounds = build_initial_rprs_bounds(0.1)

    assert bounds == pytest.approx([
        INITIAL_RPRS_BOUND_LOWER_SCALE * 0.1,
        INITIAL_RPRS_BOUND_UPPER_SCALE * 0.1,
    ])


def test_rprs_posterior_retry_walks_bounds_until_retry_cap(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {"calls": []}
    diagnostics_sequence = [
        {"clipped": True, "edge": "upper", "mode": 0.158, "std": 0.006, "bounds": [0.128, 0.188]},
        {"clipped": True, "edge": "upper", "mode": 0.182, "std": 0.005, "bounds": [0.157, 0.207]},
        {"clipped": True, "edge": "upper", "mode": 0.194, "std": 0.004, "bounds": [0.174, 0.214]},
        {"clipped": True, "edge": "upper", "mode": 0.201, "std": 0.003, "bounds": [0.186, 0.216]},
        {"clipped": True, "edge": "upper", "mode": 0.206, "std": 0.003, "bounds": [0.191, 0.221]},
        {"clipped": True, "edge": "upper", "mode": 0.210, "std": 0.003, "bounds": [0.195, 0.225]},
    ]

    def make_fit(diagnostics):
        fit = types.SimpleNamespace(
            parameters={
                "rprs": diagnostics["mode"],
                "tmid": 0.0,
                "inc": 89.0,
                "a2": 0.0,
            }
        )

        def get_parameter_posterior_recenter_diagnostics(key):
            assert key == "rprs"
            return dict(diagnostics)

        fit.get_parameter_posterior_recenter_diagnostics = get_parameter_posterior_recenter_diagnostics
        return fit

    def fake_lc_fitter(
        call_times,
        call_flux,
        call_fluxerr,
        call_airmass,
        call_prior,
        call_bounds,
        jd_times=None,
        mode=None,
        use_impactparameter_rather_than_inclination_to_fit=True,
        duration_prior=None,
    ):
        call_index = len(captured["calls"])
        captured["calls"].append({
            "prior": dict(call_prior),
            "duration_prior": duration_prior,
            "bounds": {
                key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value
                for key, value in call_bounds.items()
            },
        })
        return make_fit(diagnostics_sequence[call_index])

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    times = np.linspace(-0.03, 0.03, 7)
    flux = np.ones(7, dtype=float)
    fluxerr = np.full(7, 0.01, dtype=float)
    airmass = np.ones(7, dtype=float)
    prior = {"tmid": 0.0, "rprs": 0.1, "inc": 89.0, "a2": 0.0}
    bounds = {"rprs": [0.0, 0.125], "tmid": [-0.01, 0.01], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]}

    fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        times,
        flux,
        fluxerr,
        airmass,
        prior,
        bounds,
    )

    assert RPRS_POSTERIOR_MAX_RETRIES_DEFAULT == 5
    assert len(captured["calls"]) == 6
    np.testing.assert_allclose(
        np.asarray([call["bounds"]["rprs"] for call in captured["calls"]], dtype=float),
        np.asarray([
            [0.0, 0.125],
            [0.108, 0.208],
            [0.132, 0.232],
            [0.144, 0.244],
            [0.151, 0.251],
            [0.156, 0.256],
        ], dtype=float),
    )
    assert fit.rprs_posterior_refit_applied is True
    assert fit.rprs_posterior_refit_count == 5
    assert fit.rprs_posterior_refit_edge == "upper"
    assert fit.rprs_posterior_refit_bounds == pytest.approx([0.156, 0.256])
    assert "after 5 retries" in fit.rprs_posterior_refit_note


def test_rprs_posterior_retry_expands_bounds_without_hitting_the_old_0p3_cap(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {"calls": []}
    diagnostics_sequence = [
        {"clipped": True, "edge": "upper", "mode": 0.275, "std": 0.020, "bounds": [0.175, 0.375]},
        {"clipped": False, "edge": None, "mode": 0.278, "std": 0.012, "bounds": [0.175, RPRS_SEARCH_BOUND_MAX]},
    ]

    def make_fit(diagnostics):
        fit = types.SimpleNamespace(
            parameters={
                "rprs": diagnostics["mode"],
                "tmid": 0.0,
                "inc": 89.0,
                "a2": 0.0,
            }
        )

        def get_parameter_posterior_recenter_diagnostics(key):
            assert key == "rprs"
            return dict(diagnostics)

        fit.get_parameter_posterior_recenter_diagnostics = get_parameter_posterior_recenter_diagnostics
        return fit

    def fake_lc_fitter(
        call_times,
        call_flux,
        call_fluxerr,
        call_airmass,
        call_prior,
        call_bounds,
        jd_times=None,
        mode=None,
        use_impactparameter_rather_than_inclination_to_fit=True,
        duration_prior=None,
    ):
        call_index = len(captured["calls"])
        captured["calls"].append({
            "prior": dict(call_prior),
            "duration_prior": duration_prior,
            "bounds": {
                key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value
                for key, value in call_bounds.items()
            },
        })
        return make_fit(diagnostics_sequence[call_index])

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    times = np.linspace(-0.03, 0.03, 7)
    flux = np.ones(7, dtype=float)
    fluxerr = np.full(7, 0.01, dtype=float)
    airmass = np.ones(7, dtype=float)
    prior = {"tmid": 0.0, "rprs": 0.1, "inc": 89.0, "a2": 0.0}
    bounds = {"rprs": [0.0, 0.25], "tmid": [-0.01, 0.01], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]}

    fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        times,
        flux,
        fluxerr,
        airmass,
        prior,
        bounds,
    )

    assert len(captured["calls"]) == 2
    assert captured["calls"][0]["bounds"]["rprs"] == pytest.approx([0.0, 0.25])
    assert captured["calls"][1]["prior"]["rprs"] == pytest.approx(0.275)
    assert captured["calls"][1]["bounds"]["rprs"] == pytest.approx([0.175, 0.375])
    assert fit.rprs_posterior_refit_applied is True
    assert fit.rprs_posterior_refit_count == 1
    assert fit.rprs_posterior_refit_bounds == pytest.approx([0.175, 0.375])


def test_rprs_posterior_retry_can_continue_above_the_old_maximum_exoplanet_range(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {"calls": []}
    diagnostics = {"clipped": True, "edge": "upper", "mode": 0.275, "std": 0.020, "bounds": [0.175, 0.375]}

    def make_fit():
        fit = types.SimpleNamespace(
            parameters={
                "rprs": diagnostics["mode"],
                "tmid": 0.0,
                "inc": 89.0,
                "a2": 0.0,
            }
        )

        def get_parameter_posterior_recenter_diagnostics(key):
            assert key == "rprs"
            return dict(diagnostics)

        fit.get_parameter_posterior_recenter_diagnostics = get_parameter_posterior_recenter_diagnostics
        return fit

    def fake_lc_fitter(
        call_times,
        call_flux,
        call_fluxerr,
        call_airmass,
        call_prior,
        call_bounds,
        jd_times=None,
        mode=None,
        use_impactparameter_rather_than_inclination_to_fit=True,
        duration_prior=None,
    ):
        captured["calls"].append({
            "prior": dict(call_prior),
            "duration_prior": duration_prior,
            "bounds": {
                key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value
                for key, value in call_bounds.items()
            },
        })
        return make_fit()

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    times = np.linspace(-0.03, 0.03, 7)
    flux = np.ones(7, dtype=float)
    fluxerr = np.full(7, 0.01, dtype=float)
    airmass = np.ones(7, dtype=float)
    prior = {"tmid": 0.0, "rprs": 0.35, "inc": 89.0, "a2": 0.0}
    bounds = {"rprs": [RPRS_SEARCH_BOUND_MIN, 0.35], "tmid": [-0.01, 0.01], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]}

    fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        times,
        flux,
        fluxerr,
        airmass,
        prior,
        bounds,
    )

    assert len(captured["calls"]) == 2
    assert captured["calls"][0]["prior"]["rprs"] == pytest.approx(0.35)
    assert captured["calls"][0]["bounds"]["rprs"] == pytest.approx([RPRS_SEARCH_BOUND_MIN, 0.35])
    assert captured["calls"][1]["prior"]["rprs"] == pytest.approx(0.275)
    assert captured["calls"][1]["bounds"]["rprs"] == pytest.approx([0.175, 0.375])
    assert fit.rprs_posterior_refit_applied is True
    assert fit.rprs_posterior_refit_count == 1


def test_ars_posterior_retry_expands_bounds_when_upper_edge_is_truncated(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {"calls": []}
    diagnostics_sequence = [
        {
            "rprs": {"clipped": False, "edge": None, "mode": 0.1, "std": 0.01, "bounds": [0.05, 0.15]},
            "ars": {"clipped": True, "edge": "upper", "mode": 14.45, "std": 0.37, "bounds": [12.60, 16.30]},
        },
        {
            "rprs": {"clipped": False, "edge": None, "mode": 0.1, "std": 0.01, "bounds": [0.05, 0.15]},
            "ars": {"clipped": False, "edge": None, "mode": 14.50, "std": 0.20, "bounds": [12.60, 16.30]},
        },
    ]

    def make_fit(diagnostics):
        fit = types.SimpleNamespace(
            parameters={
                "rprs": diagnostics["rprs"]["mode"],
                "ars": diagnostics["ars"]["mode"],
                "tmid": 0.0,
                "inc": 89.0,
                "a2": 0.0,
            }
        )

        def get_parameter_posterior_recenter_diagnostics(key):
            return dict(diagnostics[key])

        fit.get_parameter_posterior_recenter_diagnostics = get_parameter_posterior_recenter_diagnostics
        return fit

    def fake_lc_fitter(
        call_times,
        call_flux,
        call_fluxerr,
        call_airmass,
        call_prior,
        call_bounds,
        jd_times=None,
        mode=None,
        use_impactparameter_rather_than_inclination_to_fit=True,
        duration_prior=None,
    ):
        call_index = len(captured["calls"])
        captured["calls"].append({
            "prior": dict(call_prior),
            "bounds": {
                key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value
                for key, value in call_bounds.items()
            },
        })
        return make_fit(diagnostics_sequence[call_index])

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        np.linspace(-0.03, 0.03, 7),
        np.ones(7, dtype=float),
        np.full(7, 0.01, dtype=float),
        np.ones(7, dtype=float),
        {"tmid": 0.0, "rprs": 0.1, "ars": 14.0, "inc": 89.0, "a2": 0.0},
        {
            "rprs": [0.0, 0.25],
            "ars": [12.80, 14.80],
            "tmid": [-0.01, 0.01],
            "inc": [84.0, 90.0],
            "a2": [-3.0, 3.0],
        },
    )

    assert len(captured["calls"]) == 2
    assert captured["calls"][0]["bounds"]["ars"] == pytest.approx([12.80, 14.80])
    assert captured["calls"][1]["prior"]["ars"] == pytest.approx(14.45)
    assert captured["calls"][1]["bounds"]["ars"] == pytest.approx([12.60, 16.30])
    assert fit.rprs_posterior_refit_applied is False
    assert fit.ars_posterior_refit_applied is True
    assert fit.ars_posterior_refit_count == 1
    assert fit.ars_posterior_refit_edge == "upper"
    assert fit.ars_posterior_refit_bounds == pytest.approx([12.60, 16.30])


def test_run_nested_lightcurve_fit_passes_duration_prior_when_available(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {}
    diagnostics = {"clipped": False, "edge": None, "mode": 0.1, "std": 0.01, "bounds": [0.05, 0.15]}

    def fake_lc_fitter(
        call_times,
        call_flux,
        call_fluxerr,
        call_airmass,
        call_prior,
        call_bounds,
        jd_times=None,
        mode=None,
        use_impactparameter_rather_than_inclination_to_fit=True,
        duration_prior=None,
    ):
        captured["duration_prior"] = duration_prior
        fit = types.SimpleNamespace(parameters={"rprs": 0.1, "ars": 15.0, "tmid": 0.0, "inc": 89.0, "a2": 0.0})

        def get_parameter_posterior_recenter_diagnostics(key):
            assert key in ("rprs", "ars")
            return dict(diagnostics)

        fit.get_parameter_posterior_recenter_diagnostics = get_parameter_posterior_recenter_diagnostics
        return fit

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    duration_prior = build_single_transit_duration_prior({
        "pPer": 1.0,
        "pPerUnc": 0.001,
        "rprs": 0.1,
        "rprsUnc": 0.01,
        "aRs": 15.0,
        "aRsUnc": 0.1,
        "inc": 89.0,
        "incUnc": 0.1,
        "ecc": 0.0,
        "omega": 0.0,
    })

    fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        np.linspace(-0.03, 0.03, 7),
        np.ones(7, dtype=float),
        np.full(7, 0.01, dtype=float),
        np.ones(7, dtype=float),
        {"tmid": 0.0, "rprs": 0.1, "ars": 15.0, "inc": 89.0, "a2": 0.0},
        {"rprs": [0.0, 0.25], "ars": [14.5, 15.5], "tmid": [-0.01, 0.01], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]},
        duration_prior=duration_prior,
    )

    assert captured["duration_prior"] == duration_prior
    assert fit.duration_prior_applied is True
    assert "expected duration=" in fit.duration_prior_note
