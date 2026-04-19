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
    comparison_candidate_fit_selection_reason,
    comparison_star_coverage_summary,
    comparison_star_stability_summary,
    diagnose_lightcurve_fit_inputs,
    detrend_flux_on_out_of_transit_baseline,
    ensure_lightcurve_fit_failure_reason,
    fit_lightcurve,
    fit_final_lightcurve_with_oot_baseline_detrending,
    fit_lightcurve_to_every_comparison_candidate,
    fit_ranked_comparison_calibration_candidates,
    get_final_fit_baseline_duration_multiplier,
    is_adaptive_aperture_mode_enabled,
    is_comp_star_required,
    is_out_of_transit_baseline_detrending_enabled,
    is_target_driven_comp_selection_enabled,
    log_comparison_calibration_fit_attempt_summaries,
    log_comparison_candidate_fit_summaries,
    log_target_fit_candidate_summaries,
    phase_bin_sigma_clip,
    representative_psf_sigma,
    run_target_driven_photometry_search,
    resolve_frame_aperture_radii,
    summarize_adaptive_aperture_usage,
    should_skip_airmass_fit,
    should_fit_lightcurve_to_every_comparison_candidate,
    should_detect_bad_pixels_before_photometry,
    should_use_aperture_photometry,
    should_use_psf_photometry,
    should_skip_low_comparison_coverage_rejection,
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


def test_should_skip_low_comparison_coverage_rejection_parses_values():
    assert should_skip_low_comparison_coverage_rejection(None) is False
    assert should_skip_low_comparison_coverage_rejection("y") is True
    assert should_skip_low_comparison_coverage_rejection("n") is False


def test_should_fit_lightcurve_to_every_comparison_candidate_parses_values():
    assert should_fit_lightcurve_to_every_comparison_candidate(None) is False
    assert should_fit_lightcurve_to_every_comparison_candidate("y") is True
    assert should_fit_lightcurve_to_every_comparison_candidate("n") is False


def test_should_detect_bad_pixels_before_photometry_parses_values():
    assert should_detect_bad_pixels_before_photometry(None) is True
    assert should_detect_bad_pixels_before_photometry("y") is True
    assert should_detect_bad_pixels_before_photometry("n") is False


def test_is_out_of_transit_baseline_detrending_enabled_parses_values():
    assert is_out_of_transit_baseline_detrending_enabled(None) is True
    assert is_out_of_transit_baseline_detrending_enabled("y") is True
    assert is_out_of_transit_baseline_detrending_enabled("n") is False
    assert is_out_of_transit_baseline_detrending_enabled(True) is True


def test_get_final_fit_baseline_duration_multiplier_parses_values():
    assert get_final_fit_baseline_duration_multiplier(None) == pytest.approx(1.0)
    assert get_final_fit_baseline_duration_multiplier("2.5") == pytest.approx(2.5)
    assert get_final_fit_baseline_duration_multiplier(0) == pytest.approx(0.0)
    assert get_final_fit_baseline_duration_multiplier(-1) == pytest.approx(1.0)


def test_should_use_psf_photometry_parses_values():
    assert should_use_psf_photometry(None) is True
    assert should_use_psf_photometry("y") is True
    assert should_use_psf_photometry("n") is False


def test_should_use_aperture_photometry_parses_values():
    assert should_use_aperture_photometry(None) is True
    assert should_use_aperture_photometry("y") is True
    assert should_use_aperture_photometry("n") is False


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


def test_fit_lightcurve_to_every_comparison_candidate_uses_selected_aperture(monkeypatch):
    calls = []

    class DummyFit:
        def __init__(self, size):
            self.residuals = np.full(size, 0.01)
            self.data = np.ones(size)

    def fake_fit_lightcurve(times, tflux, cflux, airmass, ld, p_dict, jd_times=None, **kwargs):
        calls.append({
            "times": np.asarray(times),
            "tflux": np.asarray(tflux),
            "cflux": np.asarray(cflux),
            "jd_times": np.asarray(jd_times),
            "kwargs": dict(kwargs),
        })
        return DummyFit(len(times)), np.asarray(tflux), np.asarray(cflux)

    monkeypatch.setattr("exotic.exotic.fit_lightcurve", fake_fit_lightcurve)

    times = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    jd_times = np.array([11.0, 12.0, 13.0, 14.0, 15.0, 16.0])
    airmass = np.array([1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
    aper_data = {
        "target": np.array([
            [[1.0], [10.0]],
            [[2.0], [11.0]],
            [[3.0], [12.0]],
            [[4.0], [13.0]],
            [[5.0], [14.0]],
            [[6.0], [15.0]],
        ]),
        "comp1": np.array([
            [[4.0], [20.0]],
            [[5.0], [np.nan]],
            [[6.0], [22.0]],
            [[7.0], [23.0]],
            [[8.0], [24.0]],
            [[9.0], [25.0]],
        ]),
        "comp2": np.array([
            [[7.0], [30.0]],
            [[8.0], [31.0]],
            [[9.0], [32.0]],
            [[10.0], [33.0]],
            [[11.0], [34.0]],
            [[12.0], [35.0]],
        ]),
    }
    photometry_info = {
        "best_fit_lc": object(),
        "comp_star_num": 2,
        "min_aperture": 5.0,
        "min_annulus": 12.0,
        "aperture_index": 1,
        "annulus_index": 0,
    }

    candidate_fits = fit_lightcurve_to_every_comparison_candidate(
        times,
        jd_times,
        airmass,
        ld=np.array([0.1, 0.2, 0.3, 0.4]),
        p_dict={"rprs": 0.1},
        comp_stars=[[100, 200], [300, 400]],
        psf_data={},
        aper_data=aper_data,
        photometry_info=photometry_info,
    )

    assert len(candidate_fits) == 2
    assert candidate_fits[0]["selected"] is False
    assert candidate_fits[1]["selected"] is True
    assert calls[0]["kwargs"]["final_fit_mode"] == "ns"
    assert calls[1]["kwargs"]["final_fit_mode"] == "ns"
    np.testing.assert_array_equal(calls[0]["times"], np.array([1.0, 3.0, 4.0, 5.0, 6.0]))
    np.testing.assert_array_equal(calls[0]["tflux"], np.array([10.0, 12.0, 13.0, 14.0, 15.0]))
    np.testing.assert_array_equal(calls[0]["cflux"], np.array([20.0, 22.0, 23.0, 24.0, 25.0]))
    np.testing.assert_array_equal(calls[0]["jd_times"], np.array([11.0, 13.0, 14.0, 15.0, 16.0]))
    np.testing.assert_array_equal(calls[1]["times"], np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))
    np.testing.assert_array_equal(calls[1]["tflux"], np.array([10.0, 11.0, 12.0, 13.0, 14.0, 15.0]))
    np.testing.assert_array_equal(calls[1]["cflux"], np.array([30.0, 31.0, 32.0, 33.0, 34.0, 35.0]))


def test_fit_lightcurve_to_every_comparison_candidate_records_sparse_candidate_failure(monkeypatch):
    calls = []

    class DummyFit:
        def __init__(self, size):
            self.residuals = np.full(size, 0.01)
            self.data = np.ones(size)

    def fake_fit_lightcurve(times, tflux, cflux, airmass, ld, p_dict, jd_times=None, **kwargs):
        calls.append({
            "times": np.asarray(times),
            "tflux": np.asarray(tflux),
            "cflux": np.asarray(cflux),
            "jd_times": np.asarray(jd_times),
            "kwargs": dict(kwargs),
        })
        return DummyFit(len(times)), np.asarray(tflux), np.asarray(cflux)

    monkeypatch.setattr("exotic.exotic.fit_lightcurve", fake_fit_lightcurve)

    times = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
    jd_times = np.array([11.0, 12.0, 13.0, 14.0, 15.0, 16.0])
    airmass = np.array([1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
    aper_data = {
        "target": np.array([
            [[10.0]],
            [[11.0]],
            [[12.0]],
            [[13.0]],
            [[14.0]],
            [[15.0]],
        ]),
        "comp1": np.array([
            [[20.0]],
            [[np.nan]],
            [[np.nan]],
            [[np.nan]],
            [[np.nan]],
            [[np.nan]],
        ]),
        "comp2": np.array([
            [[30.0]],
            [[31.0]],
            [[32.0]],
            [[33.0]],
            [[34.0]],
            [[35.0]],
        ]),
    }
    photometry_info = {
        "best_fit_lc": object(),
        "comp_star_num": 2,
        "min_aperture": 5.0,
        "min_annulus": 12.0,
        "aperture_index": 0,
        "annulus_index": 0,
    }

    candidate_fits = fit_lightcurve_to_every_comparison_candidate(
        times,
        jd_times,
        airmass,
        ld=np.array([0.1, 0.2, 0.3, 0.4]),
        p_dict={"rprs": 0.1},
        comp_stars=[[100, 200], [300, 400]],
        psf_data={},
        aper_data=aper_data,
        photometry_info=photometry_info,
    )

    assert len(calls) == 1
    assert candidate_fits[0]["fit"] is None
    assert candidate_fits[0]["coverage_rejected"] is True
    assert candidate_fits[0]["fit_diagnostics"]["failed_stage"] == "coverage"
    assert "low-coverage clipping" in candidate_fits[0]["failure_reason"]
    assert candidate_fits[1]["fit"] is not None
    assert candidate_fits[1]["coverage_rejected"] is False
    assert candidate_fits[1]["failure_reason"] is None
    assert calls[0]["kwargs"]["final_fit_mode"] == "ns"


def test_log_comparison_candidate_fit_summaries_includes_reasons(monkeypatch):
    logged = []
    monkeypatch.setattr("exotic.exotic.log_info", lambda message, warn=False, error=False: logged.append(message))

    candidate_fit_summaries = [
        {
            "label": "Comp 1",
            "position": [100, 200],
            "selected": False,
            "fit": None,
            "res_std": np.inf,
            "coverage_count": 1,
            "coverage_total_frame_count": 3,
            "coverage_reference_count": 3.0,
            "coverage_min_required_count": 2,
            "fit_point_count": 0,
            "fit_diagnostics": {"usable_point_count": 0},
            "failure_reason": "comparison candidate rejected after iterative low-coverage clipping (1 < 2 valid frame(s); peer median=3.0).",
        },
        {
            "label": "Comp 2",
            "position": [300, 400],
            "selected": True,
            "fit": object(),
            "res_std": 0.01,
            "coverage_count": 3,
            "coverage_total_frame_count": 3,
            "coverage_reference_count": 3.0,
            "coverage_min_required_count": 2,
            "fit_point_count": 3,
            "fit_diagnostics": {"usable_point_count": 3},
            "failure_reason": None,
            "parameter_summary": "fit_method=ultranest, Tmid=1.0 +/- 0.1",
        },
    ]

    log_comparison_candidate_fit_summaries(
        candidate_fit_summaries,
        {"selection_basis": "comparison_field", "comp_star_num": 2, "min_std": 0.01},
    )

    assert any("Selection basis: comparison-field" in message for message in logged)
    assert any("coverage=1 valid frame(s) out of 3 total; min_required=2; peer_median=3.0" in message for message in logged)
    assert any("Comp 1" in message and "reason=comparison candidate rejected after iterative low-coverage clipping" in message for message in logged)
    assert any("Comp 2 [selected]" in message and "comparison-field calibration ranked this star best" in message for message in logged)
    assert any("parameters: fit_method=ultranest" in message for message in logged)


def test_log_comparison_calibration_fit_attempt_summaries_includes_reasons(monkeypatch):
    logged = []
    monkeypatch.setattr("exotic.exotic.log_info", lambda message, warn=False, error=False: logged.append(message))

    attempts = [
        {
            "label": "Comp 1",
            "position": [100, 200],
            "selected": False,
            "aggregate_score": 0.01,
            "coverage_count": 3,
            "coverage_total_frame_count": 3,
            "coverage_reference_count": 3.0,
            "coverage_min_required_count": 2,
            "fit": None,
            "res_std": np.inf,
            "fit_point_count": 0,
            "fit_diagnostics": {"usable_point_count": 0},
            "failure_reason": "relative-flux filtering left 0 usable point(s); rejected 3/3 frame(s) during target/reference ratio screening (non-finite=0, >2x=3, finite ratio range=3.0000 to 3.0000).",
            "parameter_summary": None,
        },
        {
            "label": "Comp 2",
            "position": [300, 400],
            "selected": True,
            "aggregate_score": 0.02,
            "coverage_count": 3,
            "coverage_total_frame_count": 3,
            "coverage_reference_count": 3.0,
            "coverage_min_required_count": 2,
            "fit": object(),
            "res_std": 0.01,
            "fit_point_count": 3,
            "fit_diagnostics": {"usable_point_count": 3},
            "failure_reason": None,
            "parameter_summary": "fit_method=ultranest, Tmid=1.0 +/- 0.1",
        },
    ]

    log_comparison_calibration_fit_attempt_summaries(attempts, "Aperture photometry (aper=7.05px, annulus=22.73px)")

    assert any("Comparison-star calibration target-fit diagnostics:" in message for message in logged)
    assert any("Photometry method: Aperture photometry (aper=7.05px, annulus=22.73px)" in message for message in logged)
    assert any("Comp 1" in message and "reason=relative-flux filtering left 0 usable point(s)" in message for message in logged)
    assert any("Comp 2 [selected]" in message and "fit_points=3" in message for message in logged)
    assert any("parameters: fit_method=ultranest" in message for message in logged)


def test_log_target_fit_candidate_summaries_includes_methods_and_reasons(monkeypatch):
    logged = []
    monkeypatch.setattr("exotic.exotic.log_info", lambda message, warn=False, error=False: logged.append(message))

    candidate_summaries = [
        {
            "label": "Comp 1",
            "position": [100, 200],
            "selected": False,
            "method_label": "Aperture photometry (aper=7.05px, annulus=22.73px)",
            "prescore": 0.005,
            "fit": None,
            "res_std": np.inf,
            "coverage_count": 3,
            "coverage_total_frame_count": 3,
            "coverage_reference_count": 3.0,
            "coverage_min_required_count": 2,
            "fit_point_count": 0,
            "fit_diagnostics": {"usable_point_count": 0},
            "failure_reason": "relative-flux filtering left 0 usable point(s); rejected 3/3 frame(s) during target/reference ratio screening (non-finite=0, >2x=3, finite ratio range=3.0000 to 3.0000).",
            "parameter_summary": None,
        },
    ]

    log_target_fit_candidate_summaries(candidate_summaries)

    assert any("Target-fit candidate diagnostics:" in message for message in logged)
    assert any(
        "Comp 1" in message
        and "with Aperture photometry (aper=7.05px, annulus=22.73px)" in message
        and "reason=relative-flux filtering left 0 usable point(s)" in message
        for message in logged
    )


def test_comparison_candidate_fit_selection_reason_describes_comparison_field_retry():
    reason = comparison_candidate_fit_selection_reason(
        {
            "selected": True,
            "failure_reason": None,
            "res_std": 0.01,
        },
        {
            "selection_basis": "comparison_field_retry",
            "comp_star_num": 2,
            "min_std": 0.01,
        },
    )

    assert "fell back to this star" in reason


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


def test_comparison_star_coverage_summary_rejects_sparse_candidates():
    coverage = comparison_star_coverage_summary(
        {
            "comp1": np.array([100.0, 101.0, 100.5, 101.5, 100.8, 101.2]),
            "comp2": np.array([80.0, 80.8, 80.4, 81.0, 80.6, 80.9]),
            "comp3": np.array([60.0, np.nan, np.nan, np.nan, np.nan, 60.3]),
        }
    )

    assert not coverage["comp1"]["coverage_rejected"]
    assert not coverage["comp2"]["coverage_rejected"]
    assert coverage["comp3"]["coverage_rejected"]
    assert coverage["comp3"]["coverage_count"] == 2
    assert coverage["comp3"]["coverage_total_frame_count"] == 6


def test_comparison_star_coverage_summary_iteratively_rejects_low_count_tail():
    coverage = comparison_star_coverage_summary(
        {
            "comp1": np.array([10.0] * 10),
            "comp2": np.array([11.0] * 10),
            "comp3": np.array([12.0] * 10),
            "comp4": np.array([13.0] * 7 + [np.nan] * 3),
            "comp5": np.array([14.0] * 6 + [np.nan] * 4),
            "comp6": np.array([15.0] + [np.nan] * 9),
        }
    )

    assert not coverage["comp1"]["coverage_rejected"]
    assert not coverage["comp2"]["coverage_rejected"]
    assert not coverage["comp3"]["coverage_rejected"]
    assert coverage["comp4"]["coverage_rejected"]
    assert coverage["comp5"]["coverage_rejected"]
    assert coverage["comp6"]["coverage_rejected"]
    assert coverage["comp1"]["coverage_total_frame_count"] == 10
    assert coverage["comp1"]["coverage_reference_count"] == pytest.approx(10.0)
    assert coverage["comp1"]["coverage_min_required_count"] == 8


def test_comparison_star_stability_summary_rejects_low_coverage_candidates():
    airmass = np.linspace(1.0, 1.5, 6)
    summary = comparison_star_stability_summary(
        {
            "comp1": np.array([100.0, 101.0, 100.5, 101.5, 100.8, 101.2]),
            "comp2": np.array([80.0, 80.8, 80.4, 81.0, 80.6, 80.9]),
            "comp3": np.array([60.0, np.nan, np.nan, np.nan, np.nan, 60.3]),
        },
        airmass,
    )

    assert np.isfinite(summary["field_score"])
    assert summary["best_comp_index"] in (0, 1)
    assert summary["comp_summaries"][2]["coverage_rejected"]
    assert np.isinf(summary["comp_summaries"][2]["aggregate_score"])


def test_cheap_lightcurve_prescore_ignores_large_ratios_when_requested():
    tflux = np.array([2.0, 2.0, 2.0, 6.0, 2.0, 2.0])
    cflux = np.full(tflux.shape[0], 2.0)
    airmass = np.linspace(1.0, 1.5, tflux.shape[0])

    score = cheap_lightcurve_prescore(tflux, cflux, airmass, enforce_relative_flux_max=True)

    assert np.isclose(score, 0.0)


def test_cheap_lightcurve_prescore_allows_large_raw_target_reference_ratios():
    tflux = np.full(6, 30.0)
    cflux = np.full(6, 10.0)
    airmass = np.linspace(1.0, 1.5, 6)

    score = cheap_lightcurve_prescore(tflux, cflux, airmass, enforce_relative_flux_max=False)

    assert np.isfinite(score)


def test_cheap_lightcurve_prescore_keeps_target_only_mode_unfiltered():
    tflux = np.array([10.0, 11.0, 12.0, 13.0, 14.0, 15.0])
    cflux = np.ones(tflux.shape[0])
    airmass = np.linspace(1.0, 1.5, tflux.shape[0])

    score = cheap_lightcurve_prescore(tflux, cflux, airmass)

    assert np.isfinite(score)


def test_should_skip_airmass_fit_when_airmass_span_is_small():
    airmass = np.array([1.10, 1.12, 1.14, 1.15])

    assert should_skip_airmass_fit(airmass)


def test_detrend_flux_on_out_of_transit_baseline_removes_linear_slope():
    times = np.array([-2.0, -1.0, -0.25, 0.0, 0.25, 1.0, 2.0])
    baseline = 1.0 + 0.02 * times
    transit_profile = np.array([1.0, 1.0, 1.0, 0.99, 1.0, 1.0, 1.0])
    flux = baseline * transit_profile
    fluxerr = np.full_like(times, 0.01)
    fit = types.SimpleNamespace(
        transit=transit_profile,
        parameters={"tmid": 0.0},
    )

    result = detrend_flux_on_out_of_transit_baseline(times, flux, fluxerr, fit)

    assert result["applied"] is True
    assert np.allclose(result["flux"][[0, 1, 2, 4, 5, 6]], 1.0, atol=1e-8)
    assert result["flux"][3] == pytest.approx(0.99, abs=1e-8)
    assert result["slope"] == pytest.approx(0.02, abs=1e-8)


def test_fit_final_lightcurve_with_oot_baseline_detrending_refits_with_flattened_flux(monkeypatch):
    import exotic.exotic as exotic_module

    times = np.array([-2.0, -1.0, -0.25, 0.0, 0.25, 1.0, 2.0])
    flux = (1.0 + 0.02 * times) * np.array([1.0, 1.0, 1.0, 0.99, 1.0, 1.0, 1.0])
    fluxerr = np.full_like(times, 0.01)
    airmass = np.ones_like(times)
    prior = {"rprs": 0.1, "tmid": 0.0, "inc": 89.0, "a2": 0.0}
    bounds = {"rprs": [0.0, 0.2], "tmid": [-0.1, 0.1], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]}

    captured = {"calls": []}

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
    ):
        captured["calls"].append(np.array(call_flux, dtype=float))
        return types.SimpleNamespace(
            transit=np.array([1.0, 1.0, 1.0, 0.99, 1.0, 1.0, 1.0]),
            parameters={"tmid": 0.0, "rprs": 0.1, "inc": 89.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a2": 0.01},
            data=np.array(call_flux, dtype=float),
            residuals=np.zeros_like(call_flux, dtype=float),
        )

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    fit, refit_flux, refit_unc = fit_final_lightcurve_with_oot_baseline_detrending(
        times,
        flux,
        fluxerr,
        airmass,
        prior,
        bounds,
        skip_airmass_fit=False,
        disable_vertical_flux_normalization=False,
        detrend_on_outoftransit_baseline=True,
    )

    assert len(captured["calls"]) == 2
    assert np.allclose(captured["calls"][0], flux)
    assert np.allclose(captured["calls"][1][[0, 1, 2, 4, 5, 6]], 1.0, atol=1e-8)
    assert refit_flux[3] == pytest.approx(0.99, abs=1e-8)
    assert np.allclose(refit_unc[[0, 1, 2, 4, 5, 6]], 0.01 / (1.0 + 0.02 * times[[0, 1, 2, 4, 5, 6]]))
    assert fit.oot_baseline_detrending_applied is True
    assert fit.oot_baseline_pre_points == 3
    assert fit.oot_baseline_post_points == 3


def test_phase_bin_sigma_clip_flags_local_phase_outlier():
    phase_centers = np.linspace(-0.045, 0.045, 10)
    phase = np.concatenate([center + np.linspace(-1e-4, 1e-4, 5) for center in phase_centers])
    base_profile = np.array([-0.002, -0.001, 0.0, 0.001, 0.002])
    values = np.concatenate([1.0 + base_profile for _ in phase_centers])
    values[27] = 1.15

    mask = phase_bin_sigma_clip(values, phase, sigma=3, bins=10)

    assert mask.sum() == 1
    assert mask[27]


def test_fit_final_lightcurve_preserves_explicit_plot_time_range(monkeypatch):
    import exotic.exotic as exotic_module

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
    ):
        return types.SimpleNamespace(
            parameters={"tmid": 0.0, "rprs": 0.1, "inc": 89.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a2": 0.01},
            data=np.array(call_flux, dtype=float),
            residuals=np.zeros_like(call_flux, dtype=float),
        )

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    times = np.linspace(0.0, 0.05, 6)
    flux = np.ones(6, dtype=float)
    fluxerr = np.full(6, 0.01, dtype=float)
    airmass = np.linspace(1.0, 1.5, 6)
    prior = {"tmid": 0.0, "rprs": 0.1, "inc": 89.0, "a2": 0.0}
    bounds = {"rprs": [0.05, 0.15], "tmid": [-0.01, 0.01], "inc": [84.0, 90.0]}
    plot_time_range = (-0.12, 0.18)

    fit, _, _ = fit_final_lightcurve_with_oot_baseline_detrending(
        times,
        flux,
        fluxerr,
        airmass,
        prior,
        bounds,
        detrend_on_outoftransit_baseline=False,
        plot_time_range=plot_time_range,
    )

    assert fit.plot_time_range == pytest.approx(plot_time_range)


def test_fit_final_lightcurve_retries_nested_fit_when_rprs_posterior_is_clipped(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {"calls": []}

    def make_fit(call_flux, diagnostics):
        fit = types.SimpleNamespace(
            parameters={"tmid": 0.0, "rprs": 0.152, "inc": 89.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.002, "inc": 0.1, "a2": 0.01},
            data=np.array(call_flux, dtype=float),
            residuals=np.zeros_like(call_flux, dtype=float),
        )
        fit.get_parameter_posterior_recenter_diagnostics = lambda key: diagnostics if key == "rprs" else None
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
    ):
        captured["calls"].append({
            "prior": dict(call_prior),
            "bounds": {key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value for key, value in call_bounds.items()},
        })
        if len(captured["calls"]) == 1:
            return make_fit(
                call_flux,
                {
                    "clipped": True,
                    "edge": "upper",
                    "mode": 0.158,
                    "std": 0.006,
                    "bounds": [0.128, 0.188],
                    "reason": "posterior peaks against the upper search bound.",
                },
            )
        return make_fit(
            call_flux,
            {
                "clipped": False,
                "edge": None,
                "mode": 0.159,
                "std": 0.005,
                "bounds": [0.128, 0.188],
                "reason": "posterior support is comfortably inside the sampled bounds.",
            },
        )

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    times = np.linspace(-0.03, 0.03, 7)
    flux = np.ones(7, dtype=float)
    fluxerr = np.full(7, 0.01, dtype=float)
    airmass = np.ones(7, dtype=float)
    prior = {"tmid": 0.0, "rprs": 0.1, "inc": 89.0, "a2": 0.0}
    bounds = {"rprs": [0.0, 0.125], "tmid": [-0.01, 0.01], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]}

    fit, _, _ = fit_final_lightcurve_with_oot_baseline_detrending(
        times,
        flux,
        fluxerr,
        airmass,
        prior,
        bounds,
        detrend_on_outoftransit_baseline=False,
    )

    assert len(captured["calls"]) == 2
    assert captured["calls"][0]["bounds"]["rprs"] == pytest.approx([0.0, 0.125])
    assert captured["calls"][1]["prior"]["rprs"] == pytest.approx(0.158)
    assert captured["calls"][1]["bounds"]["rprs"] == pytest.approx([0.128, 0.188])
    assert fit.rprs_posterior_refit_applied is True
    assert fit.rprs_posterior_refit_count == 1
    assert fit.rprs_posterior_refit_edge == "upper"
    assert fit.rprs_posterior_refit_bounds == pytest.approx([0.128, 0.188])


def test_fit_final_lightcurve_prefit_refinement_trims_baseline_and_recenters_tmid(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {"calls": []}

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
    ):
        captured["calls"].append({
            "times": np.array(call_times, dtype=float),
            "bounds": {
                key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value
                for key, value in call_bounds.items()
            },
        })
        return types.SimpleNamespace(
            duration_expected=2.0,
            duration_measured=2.0,
            parameters={"tmid": 0.0, "rprs": 0.1, "inc": 89.0, "a2": 0.0, "per": 10.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a2": 0.01},
            data=np.array(call_flux, dtype=float),
            residuals=np.zeros_like(call_flux, dtype=float),
        )

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    times = np.array([-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0], dtype=float)
    flux = np.ones(times.shape[0], dtype=float)
    fluxerr = np.full(times.shape[0], 0.01, dtype=float)
    airmass = np.ones(times.shape[0], dtype=float)
    prior = {"tmid": 0.0, "rprs": 0.1, "inc": 89.0, "a2": 0.0, "per": 10.0}
    bounds = {"rprs": [0.0, 0.2], "tmid": [-2.0, 2.0], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]}

    fit, trimmed_flux, trimmed_unc = fit_final_lightcurve_with_oot_baseline_detrending(
        times,
        flux,
        fluxerr,
        airmass,
        prior,
        bounds,
        detrend_on_outoftransit_baseline=False,
        baseline_duration_multiplier=0.5,
    )

    assert len(captured["calls"]) == 2
    assert captured["calls"][0]["times"] == pytest.approx(times)
    assert captured["calls"][1]["times"] == pytest.approx(np.array([-2.0, -1.0, 0.0, 1.0, 2.0]))
    assert captured["calls"][1]["bounds"]["tmid"] == pytest.approx([-1.0, 1.0])
    assert trimmed_flux == pytest.approx(np.ones(5))
    assert trimmed_unc == pytest.approx(np.full(5, 0.01))
    assert fit.prefit_refinement_applied is True
    assert fit.prefit_refinement_trimmed_pre_points == 1
    assert fit.prefit_refinement_trimmed_post_points == 1
    assert fit.prefit_refinement_tmid_bounds == pytest.approx([-1.0, 1.0])


def test_fit_lightcurve_keeps_large_raw_target_reference_ratios(monkeypatch):
    captured = {}

    def fake_lc_fitter(
        times,
        fluxes,
        flux_unc,
        airmass,
        prior,
        bounds,
        jd_times=None,
        mode=None,
        use_impactparameter_rather_than_inclination_to_fit=True,
    ):
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
    assert len(captured["fluxes"]) == 6
    assert np.allclose(captured["fluxes"], np.array([1.0, 1.0, 1.0, 3.0, 1.0, 1.0]))
    assert np.allclose(fit_tflux, tflux)
    assert np.allclose(fit_cflux, 2.0)


def test_fit_lightcurve_preserves_explicit_plot_time_range(monkeypatch):
    captured = {}

    def fake_lc_fitter(
        times,
        fluxes,
        flux_unc,
        airmass,
        prior,
        bounds,
        jd_times=None,
        mode=None,
        use_impactparameter_rather_than_inclination_to_fit=True,
    ):
        fit = types.SimpleNamespace()
        captured["fit"] = fit
        return fit

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr("exotic.exotic.sigma_clip", lambda data, sigma=3, dt=21, po=2: np.zeros(len(data), dtype=bool))

    times = np.linspace(0.0, 0.05, 6)
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
    plot_time_range = (-0.12, 0.18)

    myfit, _, _ = fit_lightcurve(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times,
        plot_time_range=plot_time_range,
    )

    assert myfit is captured["fit"]
    assert myfit.plot_time_range == pytest.approx(plot_time_range)


def test_fit_lightcurve_centers_vertical_flux_bound_on_raw_flux_ratio(monkeypatch):
    captured = {}

    def fake_lc_fitter(
        times,
        fluxes,
        flux_unc,
        airmass,
        prior,
        bounds,
        jd_times=None,
        mode=None,
        use_impactparameter_rather_than_inclination_to_fit=True,
    ):
        fit = types.SimpleNamespace()
        captured["fit"] = fit
        captured["prior"] = dict(prior)
        captured["bounds"] = {
            key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value
            for key, value in bounds.items()
        }
        return fit

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr("exotic.exotic.sigma_clip", lambda data, sigma=3, dt=21, po=2: np.zeros(len(data), dtype=bool))

    times = np.linspace(0.0, 0.05, 6)
    tflux = np.full(times.shape[0], 100.0)
    cflux = np.full(times.shape[0], 2000.0)
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

    myfit, _, _ = fit_lightcurve(times, tflux, cflux, airmass, ld, p_dict, jd_times)

    assert myfit is captured["fit"]
    assert captured["prior"]["a0"] == pytest.approx(0.05)
    assert captured["prior"]["a1"] == pytest.approx(0.05)
    assert captured["bounds"]["a0"] == pytest.approx([0.0375, 0.0625])


def test_fit_lightcurve_rejects_undersampled_series(monkeypatch):
    called = {"count": 0}

    def fake_lc_fitter(*args, **kwargs):
        called["count"] += 1
        return types.SimpleNamespace()

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr("exotic.exotic.sigma_clip", lambda data, sigma=3, dt=21, po=2: np.zeros(len(data), dtype=bool))

    times = np.linspace(0.0, 0.03, 4)
    tflux = np.full(times.shape[0], 2.0)
    cflux = np.full(times.shape[0], 2.0)
    airmass = np.linspace(1.0, 1.3, times.shape[0])
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

    assert myfit is None
    assert fit_tflux is None
    assert fit_cflux is None
    assert called["count"] == 0


def test_fit_lightcurve_refits_after_phase_binned_clip(monkeypatch):
    captured_calls = []

    def fake_lc_fitter(
        times,
        fluxes,
        flux_unc,
        airmass,
        prior,
        bounds,
        jd_times=None,
        mode=None,
        use_impactparameter_rather_than_inclination_to_fit=True,
    ):
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


def test_fit_lightcurve_runs_nested_fit_when_requested(monkeypatch):
    captured_modes = []

    def fake_lc_fitter(
        times,
        fluxes,
        flux_unc,
        airmass,
        prior,
        bounds,
        jd_times=None,
        mode=None,
        use_impactparameter_rather_than_inclination_to_fit=True,
    ):
        captured_modes.append(mode)
        return types.SimpleNamespace()

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr("exotic.exotic.sigma_clip", lambda data, sigma=3, dt=21, po=2: np.zeros(len(data), dtype=bool))

    times = np.linspace(0.0, 0.05, 6)
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

    myfit, _, _ = fit_lightcurve(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times,
        final_fit_mode="ns",
    )

    assert myfit is not None
    assert captured_modes == ["lm", "ns"]


def test_run_target_driven_photometry_search_selects_best_method_across_psf_and_aperture(monkeypatch):
    evaluated = []

    class DummyFit:
        def __init__(self, residual_level):
            self.residuals = np.full(6, residual_level)
            self.data = np.ones(6)

    def fake_evaluate(task):
        _, tflux, cflux, _, _, _, _, _, _, _ = task
        evaluated.append(np.asarray(cflux))
        cflux = np.asarray(cflux)
        tflux = np.asarray(tflux)
        if np.allclose(cflux, 20.0):
            return {"myfit": DummyFit(0.02), "res_std": 0.02}, tflux, cflux
        if np.allclose(cflux, 40.0):
            return {"myfit": DummyFit(0.01), "res_std": 0.01}, tflux, cflux
        raise AssertionError("Unexpected candidate flux passed to evaluator.")

    monkeypatch.setattr("exotic.exotic.evaluate_lightcurve_candidate", fake_evaluate)

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.5, 6)
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
    psf_target_amp = 20.0 / (2.0 * np.pi)
    psf_comp_amp = 20.0 / (2.0 * np.pi)
    psf_data = {
        "target": np.column_stack([
            np.zeros(6),
            np.zeros(6),
            np.full(6, psf_target_amp),
            np.ones(6),
            np.ones(6),
        ]),
        "comp1": np.column_stack([
            np.ones(6),
            np.ones(6),
            np.full(6, psf_comp_amp),
            np.ones(6),
            np.ones(6),
        ]),
    }
    aper_data = {
        "target": np.full((6, 1, 1), 40.0),
        "comp1": np.full((6, 1, 1), 40.0),
    }

    result = run_target_driven_photometry_search(
        times,
        jd_times,
        airmass,
        ld,
        p_dict,
        comp_stars=[[100.0, 200.0]],
        psf_data=psf_data,
        aper_data=aper_data,
        apers=np.array([5.0]),
        annuli=np.array([12.0]),
        sigma=1.0,
        require_comp_star=True,
        use_psf_photometry=True,
        use_aperture_photometry=True,
    )

    assert len(evaluated) == 2
    assert {tuple(np.unique(values)) for values in evaluated} == {(20.0,), (40.0,)}
    assert result["best_candidate"]["method"] == "aperture"
    assert result["best_candidate"]["comp_index"] == 0
    assert result["min_std"] == pytest.approx(0.01)


def test_fit_ranked_comparison_calibration_candidates_retries_next_best_candidate(monkeypatch):
    class DummyFit:
        def __init__(self):
            self.residuals = np.full(6, 0.01)
            self.data = np.ones(6)

    def fake_fit_lightcurve(times, tflux, cflux, airmass, ld, p_dict, jd_times, **kwargs):
        cflux = np.asarray(cflux, dtype=float)
        if np.allclose(cflux, 0.0):
            return None, None, None
        return DummyFit(), np.asarray(tflux, dtype=float), cflux

    monkeypatch.setattr("exotic.exotic.fit_lightcurve", fake_fit_lightcurve)

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.5, 6)
    comparison_calibration = {
        "method": "aperture",
        "a": 0,
        "an": 0,
        "aper": 5.0,
        "annulus": 12.0,
        "best_comp_index": 0,
        "comp_summaries": [
            {
                "comp_index": 0,
                "key": "comp1",
                "aggregate_score": 0.01,
                "coverage_rejected": False,
            },
            {
                "comp_index": 1,
                "key": "comp2",
                "aggregate_score": 0.02,
                "coverage_rejected": False,
            },
        ],
    }
    aper_data = {
        "target": np.full((6, 1, 1), 10.0),
        "comp1": np.zeros((6, 1, 1)),
        "comp2": np.full((6, 1, 1), 5.0),
    }

    result = fit_ranked_comparison_calibration_candidates(
        times,
        jd_times,
        airmass,
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={},
        comparison_calibration=comparison_calibration,
        psf_data={},
        aper_data=aper_data,
        target_psf_flux=np.ones(6),
    )

    assert [attempt["comp_index"] for attempt in result["attempts"]] == [0, 1]
    assert result["selected_result"]["comp_index"] == 1
    assert "relative-flux filtering left 0 usable point(s)" in result["attempts"][0]["fit_diagnostics"]["failure_reason"]
    assert "non-finite=6" in result["attempts"][0]["fit_diagnostics"]["failure_reason"]
    assert result["attempts"][1]["fit"] is not None


def test_diagnose_lightcurve_fit_inputs_reports_relative_flux_breakdown():
    diagnostics = diagnose_lightcurve_fit_inputs(
        np.linspace(0.0, 0.05, 6),
        np.full(6, 30.0),
        np.full(6, 10.0),
        np.linspace(1.0, 1.5, 6),
    )

    assert diagnostics["failed_stage"] == "relative_flux_filter"
    assert "relative-flux filtering left 0 usable point(s)" in diagnostics["failure_reason"]
    assert "non-finite=0" in diagnostics["failure_reason"]
    assert ">2x=6" in diagnostics["failure_reason"]
    assert "finite ratio range=3.0000 to 3.0000" in diagnostics["failure_reason"]


def test_run_target_driven_photometry_search_returns_failed_candidate_summaries(monkeypatch):
    monkeypatch.setattr("exotic.exotic.fit_lightcurve", lambda *args, **kwargs: (None, None, None))

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.5, 6)
    aper_data = {
        "target": np.full((6, 1, 1), 30.0),
        "comp1": np.full((6, 1, 1), 10.0),
        "comp2": np.full((6, 1, 1), 12.0),
    }

    result = run_target_driven_photometry_search(
        times,
        jd_times,
        airmass,
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={},
        comp_stars=[[100.0, 200.0], [300.0, 400.0]],
        psf_data={},
        aper_data=aper_data,
        apers=np.array([7.05]),
        annuli=np.array([22.73]),
        sigma=1.0,
        require_comp_star=True,
        use_psf_photometry=False,
        use_aperture_photometry=True,
        multiprocess_lightcurve_fits=0,
    )

    assert result["best_candidate"] is None
    assert len(result["candidate_summaries"]) == 2
    assert all(
        summary["failure_reason"] is not None
        for summary in result["candidate_summaries"]
    )
    assert all(
        ">2x=" not in summary["failure_reason"]
        for summary in result["candidate_summaries"]
    )
    assert result["candidate_summaries"][0]["method_label"] == "Aperture photometry (aper=7.05px, annulus=22.73px)"


def test_fit_lightcurve_to_every_comparison_candidate_forwards_full_plot_time_range(monkeypatch):
    captured_plot_ranges = []

    class DummyFit:
        def __init__(self, plot_time_range):
            self.plot_time_range = plot_time_range
            self.parameters = {"tmid": 0.0, "rprs": 0.1, "inc": 89.0, "a0": 1.0, "a2": 0.0}
            self.errors = {"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.01, "a2": 0.01}
            self.residuals = np.full(6, 0.01, dtype=float)
            self.data = np.ones(6, dtype=float)

    def fake_fit_lightcurve(times, tflux, cflux, airmass, ld, p_dict, jd_times=None, **kwargs):
        plot_time_range = kwargs.get("plot_time_range")
        captured_plot_ranges.append(plot_time_range)
        return DummyFit(plot_time_range), np.asarray(tflux, dtype=float), np.asarray(cflux, dtype=float)

    monkeypatch.setattr("exotic.exotic.fit_lightcurve", fake_fit_lightcurve)
    monkeypatch.setattr(
        "exotic.exotic.diagnose_lightcurve_fit_inputs",
        lambda *args, **kwargs: {
            "input_point_count": 6,
            "has_reference_flux": True,
            "relative_flux_point_count": 6,
            "sigma_clip_point_count": 6,
            "usable_point_count": 6,
            "failed_stage": None,
            "failure_reason": None,
        },
    )

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.5, 6)
    psf_series = np.ones((6, 7), dtype=float)
    psf_data = {
        "target": psf_series.copy(),
        "comp1": psf_series.copy(),
    }
    photometry_info = {
        "best_fit_lc": object(),
        "comp_star_num": 1,
        "min_aperture": 0,
    }
    plot_time_range = (-0.12, 0.18)

    summaries = fit_lightcurve_to_every_comparison_candidate(
        times,
        jd_times,
        airmass,
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={},
        comp_stars=[[100.0, 200.0]],
        psf_data=psf_data,
        aper_data=None,
        photometry_info=photometry_info,
        plot_time_range=plot_time_range,
    )

    assert captured_plot_ranges == [plot_time_range]
    assert summaries[0]["fit"].plot_time_range == pytest.approx(plot_time_range)


def test_ensure_lightcurve_fit_failure_reason_preserves_existing_diagnostic_reason():
    diagnostics = {
        "failed_stage": "minimum_points",
        "failure_reason": "only 4 usable point(s) remained after filtering; need at least 5 for a lightcurve fit.",
    }

    result = ensure_lightcurve_fit_failure_reason(
        diagnostics,
        fit_result=None,
        failed_stage="lightcurve_fit",
        failure_reason="the lightcurve fitter did not converge to a usable solution.",
    )

    assert result["failed_stage"] == "minimum_points"
    assert result["failure_reason"] == diagnostics["failure_reason"]


def test_ensure_lightcurve_fit_failure_reason_adds_generic_reason_when_missing():
    diagnostics = {
        "failed_stage": None,
        "failure_reason": None,
    }

    result = ensure_lightcurve_fit_failure_reason(
        diagnostics,
        fit_result=None,
        failed_stage="lightcurve_fit",
        failure_reason="the lightcurve fitter did not converge to a usable solution.",
    )

    assert result["failed_stage"] == "lightcurve_fit"
    assert result["failure_reason"] == "the lightcurve fitter did not converge to a usable solution."


def test_fit_lightcurve_can_disable_impact_parameter_parameterization(monkeypatch):
    captured = {"flags": []}

    def fake_lc_fitter(
        times,
        fluxes,
        flux_unc,
        airmass,
        prior,
        bounds,
        jd_times=None,
        mode=None,
        use_impactparameter_rather_than_inclination_to_fit=True,
    ):
        captured["flags"].append(use_impactparameter_rather_than_inclination_to_fit)
        return types.SimpleNamespace()

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr("exotic.exotic.sigma_clip", lambda data, sigma=3, dt=21, po=2: np.zeros(len(data), dtype=bool))

    times = np.linspace(0.0, 0.05, 6)
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

    fit_lightcurve(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times,
        use_impactparameter_rather_than_inclination_to_fit=False,
    )

    assert captured["flags"] == [False]


def test_fit_lightcurve_skips_airmass_term_when_airmass_span_is_small(monkeypatch):
    captured = {}

    def fake_lc_fitter(
        times,
        fluxes,
        flux_unc,
        airmass,
        prior,
        bounds,
        jd_times=None,
        mode=None,
        use_impactparameter_rather_than_inclination_to_fit=True,
    ):
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
