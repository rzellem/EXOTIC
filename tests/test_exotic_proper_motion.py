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
    adaptive_aperture_outlier_mask,
    auto_tune_aperture_sigma_grid,
    check_coordinates,
    cheap_lightcurve_prescore,
    comparison_star_stability_summary,
    fit_lightcurve,
    is_adaptive_aperture_mode_enabled,
    is_comp_star_required,
    is_target_driven_comp_selection_enabled,
    phase_bin_sigma_clip,
    representative_psf_sigma,
    resolve_frame_aperture_radii,
    summarize_adaptive_aperture_usage,
    should_skip_airmass_fit,
    should_use_fast_target_centroid,
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


def test_is_adaptive_aperture_mode_enabled_parses_values():
    assert is_adaptive_aperture_mode_enabled(None) is False
    assert is_adaptive_aperture_mode_enabled("y") is True
    assert is_adaptive_aperture_mode_enabled("n") is False
    assert is_adaptive_aperture_mode_enabled(True) is True


def test_should_use_fast_target_centroid_disables_fast_sigma_path_for_adaptive_runs():
    assert should_use_fast_target_centroid(1, adaptive_apertures=False) is True
    assert should_use_fast_target_centroid(6, adaptive_apertures=False) is False
    assert should_use_fast_target_centroid(1, adaptive_apertures=True) is False


def test_resolve_frame_aperture_radii_scales_sigma_grid():
    apertures, annuli = resolve_frame_aperture_radii(
        np.array([2.0, 3.0]),
        np.array([8.0, 10.0]),
        adaptive_apertures=True,
        frame_sigma=1.5,
        fallback_sigma=1.0,
    )

    assert np.allclose(apertures, np.array([3.0, 4.5]))
    assert np.allclose(annuli, np.array([12.0, 15.0]))


def test_representative_psf_sigma_uses_valid_frames_and_fallback():
    psf_rows = np.array([
        [0.0, 0.0, 1.0, 2.0, 2.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 2.2, 1.8, 0.0, 0.0],
        [0.0, 0.0, 1.0, np.nan, np.nan, 0.0, 0.0],
    ])

    assert np.isclose(representative_psf_sigma(psf_rows, fallback_sigma=1.0), 2.0)
    assert np.isclose(representative_psf_sigma(np.full((0, 7), np.nan), fallback_sigma=1.25), 1.25)


def test_summarize_adaptive_aperture_usage_reports_frame_scaled_stats():
    psf_rows = np.array([
        [0.0, 0.0, 1.0, 2.0, 2.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 3.0, 3.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 4.0, 4.0, 0.0, 0.0],
    ])

    summary = summarize_adaptive_aperture_usage(psf_rows, aperture_scale=2.5, annulus_scale=9.0, fallback_sigma=1.0)

    np.testing.assert_allclose(summary["aperture_series"], np.array([5.0, 7.5, 10.0]))
    np.testing.assert_allclose(summary["annulus_series"], np.array([18.0, 27.0, 36.0]))
    np.testing.assert_allclose(summary["fwhm_series"], np.array([4.71, 7.065, 9.42]))
    assert np.isclose(summary["aperture_median"], 7.5)
    assert np.isclose(summary["aperture_std"], np.std([5.0, 7.5, 10.0]))
    assert np.isclose(summary["aperture_min"], 5.0)
    assert np.isclose(summary["aperture_max"], 10.0)
    assert summary["aperture_sigma"] == 2.5
    assert summary["annulus_sigma"] == 9.0


def test_adaptive_aperture_outlier_mask_rejects_isolated_spike_but_keeps_repeated_lower_mode():
    aperture_series = np.array([10.0, 10.1, 9.4, 10.0, 9.4, 10.1, 10.0, 15.2, 10.1, 9.4, 10.0, 10.1])
    annulus_series = aperture_series * 3.0

    mask = adaptive_aperture_outlier_mask(aperture_series, annulus_series)

    expected = np.zeros_like(aperture_series, dtype=bool)
    expected[7] = True
    np.testing.assert_array_equal(mask, expected)


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


def test_should_skip_airmass_fit_when_airmass_span_is_small():
    airmass = np.array([1.10, 1.12, 1.14, 1.15])

    assert should_skip_airmass_fit(airmass)


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


def test_fit_lightcurve_skips_airmass_term_when_airmass_span_is_small(monkeypatch):
    captured = {}

    def fake_lc_fitter(times, fluxes, flux_unc, airmass, prior, bounds, jd_times=None, mode=None):
        captured["bounds"] = dict(bounds)
        captured["airmass"] = np.array(airmass)
        return types.SimpleNamespace()

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr("exotic.exotic.sigma_clip", lambda data, sigma=3, dt=21, po=2: np.zeros(len(data), dtype=bool))

    times = np.linspace(0.0, 0.05, 6)
    tflux = np.full(times.shape[0], 2.0)
    cflux = np.full(times.shape[0], 2.0)
    airmass = np.array([1.10, 1.11, 1.12, 1.13, 1.14, 1.15])
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

    myfit, _, _ = fit_lightcurve(times, tflux, cflux, airmass, ld, p_dict, jd_times)

    assert myfit is not None
    assert "a2" not in captured["bounds"]
    assert myfit.airmass_fit_skipped is True


def _run_main_until_vertical_flux_bound(monkeypatch, tmp_path, disable_vertical_flux_normalization=Ellipsis):
    import exotic.exotic as exotic_module

    class BoundReached(Exception):
        pass

    prered_file = tmp_path / "prereduced.csv"
    prered_file.write_text(
        "\n".join(
            [
                "2450000.00,1.00,0.01,1.10",
                "2450000.10,1.01,0.01,1.12",
                "2450000.20,0.99,0.01,1.14",
                "2450000.30,1.00,0.01,1.16",
                "2450000.40,1.02,0.01,1.18",
                "2450000.50,1.01,0.01,1.20",
            ]
        )
    )

    user_pdict = {
        "ra": 10.0,
        "dec": 20.0,
        "pName": "Test Planet b",
        "sName": "Test Star",
        "pPer": 1.0,
        "pPerUnc": 0.001,
        "midT": 2450000.25,
        "midTUnc": 0.001,
        "rprs": 0.1,
        "rprsUnc": 0.01,
        "aRs": 15.0,
        "aRsUnc": 0.1,
        "inc": 89.0,
        "incUnc": 0.1,
        "omega": 0.0,
        "ecc": 0.0,
        "teff": 5500.0,
        "teffUncPos": 100.0,
        "teffUncNeg": 100.0,
        "met": 0.0,
        "metUncPos": 0.1,
        "metUncNeg": 0.1,
        "logg": 4.4,
        "loggUncPos": 0.1,
        "loggUncNeg": 0.1,
        "dist": 100.0,
        "pm_ra": 0.0,
        "pm_dec": 0.0,
    }
    exotic_info = {
        "save": tmp_path,
        "prered_file": prered_file,
        "file_time": "BJD_TDB",
        "file_units": "flux",
        "airmass_already_corrected": False,
        "random_seed": 123,
        "date": "2026-03-19",
    }
    if disable_vertical_flux_normalization is not Ellipsis:
        exotic_info["disable_vertical_flux_normalization"] = disable_vertical_flux_normalization

    args = types.SimpleNamespace(
        multiprocess_transformations=None,
        multiprocess_lightcurve_fits=None,
        realtime=None,
        reduce=None,
        prereduced=str(tmp_path / "inits.json"),
        photometry=None,
        override=True,
        nasaexoarch=False,
        non_interactive_run=True,
        use_nextastro_astrometry=False,
        use_nextastro_variability_server=False,
    )

    class FakeInputs:
        def __init__(self, init_opt):
            self.init_opt = init_opt

        def search_init(self, init_path, planet_dict):
            return init_path, dict(user_pdict)

        def prereduced(self, planet):
            return dict(exotic_info), planet or user_pdict["pName"]

    captured = {}

    monkeypatch.setattr(exotic_module, "parse_args", lambda: args)
    monkeypatch.setattr(exotic_module, "Inputs", FakeInputs)
    monkeypatch.setattr(
        exotic_module,
        "get_ld_values",
        lambda *_args, **_kwargs: ([0.1, 0.1, 0.1, 0.1], [0.1], [0.1], [0.1], [0.1]),
    )

    def fake_apply_vertical_flux_normalization_bound(prior, bounds, flux_values, disabled):
        captured["disabled"] = disabled
        raise BoundReached()

    monkeypatch.setattr(
        exotic_module,
        "apply_vertical_flux_normalization_bound",
        fake_apply_vertical_flux_normalization_bound,
    )

    with pytest.raises(BoundReached):
        exotic_module.main()

    return captured["disabled"]


def test_main_prereduced_defaults_vertical_flux_normalization_to_enabled(monkeypatch, tmp_path):
    disabled = _run_main_until_vertical_flux_bound(monkeypatch, tmp_path)

    assert disabled is False


def test_main_prereduced_respects_disable_vertical_flux_normalization_option(monkeypatch, tmp_path):
    disabled = _run_main_until_vertical_flux_bound(monkeypatch, tmp_path, disable_vertical_flux_normalization=True)

    assert disabled is True
