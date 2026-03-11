import importlib.util
import sys
import types
import numpy as np


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

from exotic.exotic import (
    auto_tune_aperture_sigma_grid,
    check_coordinates,
    cheap_lightcurve_prescore,
    comparison_star_stability_summary,
    fit_lightcurve,
    is_comp_star_required,
    is_target_driven_comp_selection_enabled,
    phase_bin_sigma_clip,
    update_coordinates_with_proper_motion,
)


def test_update_coordinates_handles_non_numeric_proper_motion_values():
    info = {
        "ra": 10.0,
        "dec": 20.0,
        "dist": "",
        "pm_ra": "nan-value",
        "pm_dec": None,
    }

    updated_ra, updated_dec = update_coordinates_with_proper_motion(info, 2459945.5)

    assert updated_ra == info["ra"]
    assert updated_dec == info["dec"]


def test_update_coordinates_accepts_numeric_strings():
    info = {
        "ra": 10.0,
        "dec": 20.0,
        "dist": "100",
        "pm_ra": "10.5",
        "pm_dec": "-5.25",
    }

    updated_ra, updated_dec = update_coordinates_with_proper_motion(info, 2459945.5)

    assert isinstance(updated_ra, float)
    assert isinstance(updated_dec, float)


def test_check_coordinates_non_interactive_prefers_wcs_centroid():
    x_pixel, y_pixel = check_coordinates(
        input_x_pixel=5,
        input_y_pixel=5,
        centroid_x=100.25,
        centroid_y=200.75,
        sigma_x=1.0,
        sigma_y=1.0,
        calculated_x_pixel=100,
        calculated_y_pixel=201,
        non_interactive_run=True,
    )

    assert x_pixel == 100.25
    assert y_pixel == 200.75


def test_check_coordinates_non_interactive_uses_wcs_pixel_when_centroid_is_nan():
    x_pixel, y_pixel = check_coordinates(
        input_x_pixel=5,
        input_y_pixel=5,
        centroid_x=float("nan"),
        centroid_y=float("nan"),
        sigma_x=1.0,
        sigma_y=1.0,
        calculated_x_pixel=100,
        calculated_y_pixel=201,
        non_interactive_run=True,
    )

    assert x_pixel == 100
    assert y_pixel == 201


def test_is_comp_star_required_parses_values():
    assert is_comp_star_required(None) is True
    assert is_comp_star_required("y") is True
    assert is_comp_star_required("n") is False


def test_is_target_driven_comp_selection_enabled_parses_values():
    assert is_target_driven_comp_selection_enabled(None) is False
    assert is_target_driven_comp_selection_enabled("y") is True
    assert is_target_driven_comp_selection_enabled("n") is False


def test_auto_tune_aperture_grid_uses_comparison_field_consistency():
    coarse_apertures_sigma = np.array([2.0, 3.0])
    coarse_annuli_sigma = np.array([8.0])
    coarse_aper_data = {
        "target": np.array([[[10.0]], [[11.0]], [[12.0]], [[13.0]], [[14.0]], [[15.0]]]),
        "comp1": np.array([[[5.0], [5.0]], [[5.0], [5.0]], [[5.0], [5.0]], [[5.0], [5.0]], [[5.0], [5.0]], [[5.0], [5.0]]]),
        "comp2": np.array([[[7.5], [7.5]], [[7.5], [7.5]], [[7.5], [12.0]], [[7.5], [7.5]], [[7.5], [7.5]], [[7.5], [7.5]]]),
    }
    subset_airmass = np.arange(1.0, 7.0)

    _, _, best_candidate, _ = auto_tune_aperture_sigma_grid(
        coarse_apertures_sigma,
        coarse_annuli_sigma,
        coarse_aper_data,
        comp_star_count=2,
        subset_airmass=subset_airmass,
        require_comp_star=True,
    )

    assert best_candidate["aper_sigma"] == 2.0
    assert best_candidate["comp_index"] in (0, 1)


def test_comparison_star_stability_summary_penalizes_variable_candidates():
    airmass = np.linspace(1.0, 1.5, 6)
    summary = comparison_star_stability_summary(
        {
            "comp1": np.array([100.0, 101.0, 100.5, 101.5, 100.8, 101.2]),
            "comp2": np.array([80.0, 80.8, 80.4, 81.0, 80.6, 80.9]),
            "comp3": np.array([60.0, 60.4, 84.0, 60.6, 60.5, 60.3]),
        },
        airmass,
    )

    assert np.isfinite(summary["field_score"])
    assert summary["best_comp_index"] in (0, 1)
    assert summary["comp_summaries"][2]["aggregate_score"] > summary["comp_summaries"][0]["aggregate_score"]


def test_cheap_lightcurve_prescore_ignores_relative_flux_above_two():
    tflux = np.array([2.0, 2.0, 2.0, 6.0, 2.0, 2.0])
    cflux = np.full(tflux.shape[0], 2.0)
    airmass = np.linspace(1.0, 1.5, tflux.shape[0])

    score = cheap_lightcurve_prescore(tflux, cflux, airmass)

    assert np.isclose(score, 0.0)


def test_cheap_lightcurve_prescore_keeps_target_only_mode_unfiltered():
    tflux = np.array([10.0, 11.0, 12.0, 13.0, 14.0, 15.0])
    cflux = np.ones(tflux.shape[0])
    airmass = np.linspace(1.0, 1.5, tflux.shape[0])

    score = cheap_lightcurve_prescore(tflux, cflux, airmass)

    assert np.isfinite(score)


def test_phase_bin_sigma_clip_flags_local_phase_outlier():
    phase_centers = np.linspace(-0.045, 0.045, 10)
    phase = np.concatenate([center + np.linspace(-1e-4, 1e-4, 5) for center in phase_centers])
    base_profile = np.array([-0.002, -0.001, 0.0, 0.001, 0.002])
    values = np.concatenate([1.0 + base_profile for _ in phase_centers])
    values[27] = 1.15

    mask = phase_bin_sigma_clip(values, phase, sigma=3, bins=10)

    assert mask.sum() == 1
    assert mask[27]


def test_fit_lightcurve_removes_relative_flux_above_two_before_fit(monkeypatch):
    captured = {}

    def fake_lc_fitter(times, fluxes, flux_unc, airmass, prior, bounds, jd_times=None, mode=None):
        captured["times"] = np.array(times)
        captured["fluxes"] = np.array(fluxes)
        captured["flux_unc"] = np.array(flux_unc)
        captured["airmass"] = np.array(airmass)
        captured["jd_times"] = np.array(jd_times)
        captured["mode"] = mode
        return types.SimpleNamespace()

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr("exotic.exotic.sigma_clip", lambda data, sigma=3, dt=21, po=2: np.zeros(len(data), dtype=bool))

    times = np.linspace(0.0, 0.05, 6)
    tflux = np.array([2.0, 2.0, 2.0, 6.0, 2.0, 2.0])
    cflux = np.full(tflux.shape[0], 2.0)
    airmass = np.linspace(1.0, 1.5, tflux.shape[0])
    jd_times = 2460000.0 + times
    ld = [0.1, 0.1, 0.1, 0.1]
    p_dict = {
        "rprs": 0.1,
        "aRs": 15.0,
        "pPer": 1.0,
        "inc": 89.0,
        "ecc": 0.0,
        "omega": 0.0,
        "midT": 0.02,
        "midTUnc": 0.001,
        "pPerUnc": 0.001,
    }

    myfit, fit_tflux, fit_cflux = fit_lightcurve(times, tflux, cflux, airmass, ld, p_dict, jd_times)

    assert myfit is not None
    assert captured["mode"] == "lm"
    assert len(captured["fluxes"]) == 5
    assert np.all(captured["fluxes"] <= 2.0)
    assert np.allclose(captured["fluxes"], 1.0)
    assert np.allclose(fit_tflux, 2.0)
    assert np.allclose(fit_cflux, 2.0)


def test_fit_lightcurve_refits_after_phase_binned_clip(monkeypatch):
    captured_calls = []

    def fake_lc_fitter(times, fluxes, flux_unc, airmass, prior, bounds, jd_times=None, mode=None):
        call_index = len(captured_calls)
        captured_calls.append({
            "times": np.array(times),
            "fluxes": np.array(fluxes),
            "flux_unc": np.array(flux_unc),
            "airmass": np.array(airmass),
            "jd_times": np.array(jd_times),
            "mode": mode,
        })
        if call_index == 0:
            return types.SimpleNamespace(
                residuals=np.zeros(len(times)),
                phase=np.linspace(-0.05, 0.05, len(times)),
            )
        return types.SimpleNamespace()

    def fake_phase_bin_sigma_clip(values, phase, sigma=3, bins=10, min_points=5, max_iters=3):
        mask = np.zeros(len(values), dtype=bool)
        mask[-1] = True
        return mask

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr("exotic.exotic.sigma_clip", lambda data, sigma=3, dt=21, po=2: np.zeros(len(data), dtype=bool))
    monkeypatch.setattr("exotic.exotic.phase_bin_sigma_clip", fake_phase_bin_sigma_clip)

    times = np.linspace(0.0, 0.08, 8)
    tflux = np.full(times.shape[0], 2.0)
    cflux = np.full(times.shape[0], 2.0)
    airmass = np.linspace(1.0, 1.5, times.shape[0])
    jd_times = 2460000.0 + times
    ld = [0.1, 0.1, 0.1, 0.1]
    p_dict = {
        "rprs": 0.1,
        "aRs": 15.0,
        "pPer": 1.0,
        "inc": 89.0,
        "ecc": 0.0,
        "omega": 0.0,
        "midT": 0.02,
        "midTUnc": 0.001,
        "pPerUnc": 0.001,
    }

    myfit, fit_tflux, fit_cflux = fit_lightcurve(times, tflux, cflux, airmass, ld, p_dict, jd_times)

    assert myfit is not None
    assert len(captured_calls) == 2
    assert captured_calls[0]["mode"] == "lm"
    assert captured_calls[1]["mode"] == "lm"
    assert len(captured_calls[0]["times"]) == 8
    assert len(captured_calls[1]["times"]) == 7
    assert len(fit_tflux) == 7
    assert len(fit_cflux) == 7
