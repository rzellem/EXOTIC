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
    SPARSE_POSTERIOR_LIVE_POINT_RETRY_FACTOR_DEFAULT,
    build_fast_ultranest_lightcurve_series,
    build_expected_transit_coverage_assessment,
    evaluate_sparse_posterior_sample_support,
    extend_sparse_posterior_live_points_if_needed,
    fit_final_lightcurve_with_oot_baseline_detrending,
    finalize_comparison_candidate_full_reduction,
    refit_selected_fast_comparison_on_full_lightcurve,
    configure_rprs_search_bound_max,
    should_run_fast_ultranest_before_final_run,
    build_single_transit_duration_prior,
    build_initial_rprs_bounds,
    run_nested_lightcurve_fit_with_rprs_posterior_retry,
)


def test_build_initial_rprs_bounds_allows_zero_depth_search_box():
    bounds = build_initial_rprs_bounds(0.1)

    assert bounds == pytest.approx([
        RPRS_SEARCH_BOUND_MIN,
        INITIAL_RPRS_BOUND_UPPER_SCALE * 0.1,
    ])
    assert INITIAL_RPRS_BOUND_LOWER_SCALE == pytest.approx(0.0)


def test_build_initial_rprs_bounds_clamps_to_configured_search_ceiling(monkeypatch):
    import exotic.exotic as exotic_module

    monkeypatch.setattr(exotic_module, "RPRS_SEARCH_BOUND_MAX", 0.5)

    assert build_initial_rprs_bounds(0.2) == pytest.approx([RPRS_SEARCH_BOUND_MIN, 0.5])
    assert build_initial_rprs_bounds(0.7) == pytest.approx([RPRS_SEARCH_BOUND_MIN, 0.5])


def test_configure_rprs_search_bound_max_updates_retry_ceiling(monkeypatch):
    import exotic.exotic as exotic_module

    monkeypatch.setattr(exotic_module, "RPRS_SEARCH_BOUND_MAX", 0.5)

    assert configure_rprs_search_bound_max("0.4") == pytest.approx(0.4)
    assert exotic_module.RPRS_SEARCH_BOUND_MAX == pytest.approx(0.4)


def test_rprs_posterior_retry_clamps_to_configured_search_ceiling(monkeypatch):
    import exotic.exotic as exotic_module

    monkeypatch.setattr(exotic_module, "RPRS_SEARCH_BOUND_MAX", 0.5)
    captured = {"calls": []}
    diagnostics_sequence = [
        {"clipped": True, "edge": "upper", "mode": 0.49, "std": 0.10, "bounds": [0.29, 0.89]},
        {"clipped": True, "edge": "upper", "mode": 0.49, "std": 0.08, "bounds": [0.38, 0.78]},
    ]

    def make_fit(diagnostics):
        fit = types.SimpleNamespace(
            parameters={"tmid": 0.0, "rprs": diagnostics["mode"], "inc": 89.0, "a2": 0.0}
        )
        fit.get_parameter_posterior_recenter_diagnostics = (
            lambda key: dict(diagnostics) if key == "rprs" else None
        )
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
        return make_fit(diagnostics_sequence[min(call_index, len(diagnostics_sequence) - 1)])

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    times = np.linspace(-0.03, 0.03, 7)
    flux = np.ones(7, dtype=float)
    fluxerr = np.full(7, 0.01, dtype=float)
    airmass = np.ones(7, dtype=float)
    prior = {"tmid": 0.0, "rprs": 0.4, "inc": 89.0, "a2": 0.0}
    bounds = {"rprs": [0.0, 0.4], "tmid": [-0.01, 0.01], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]}

    fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        times,
        flux,
        fluxerr,
        airmass,
        prior,
        bounds,
    )

    assert len(captured["calls"]) == 2
    assert captured["calls"][1]["bounds"]["rprs"][1] == pytest.approx(0.5)
    assert fit.rprs_posterior_refit_bounds[1] == pytest.approx(0.5)


def test_fast_ultranest_option_defaults_enabled_and_parses_false_values():
    assert should_run_fast_ultranest_before_final_run(None) is True
    assert should_run_fast_ultranest_before_final_run("n") is False
    assert should_run_fast_ultranest_before_final_run(False) is False


def test_fast_ultranest_binning_reduces_large_light_curve_to_twenty_points():
    times = np.linspace(0.0, 1.0, 80)
    flux = 1.0 + 0.01 * np.sin(np.linspace(0.0, 2.0 * np.pi, 80))
    unc = np.full(80, 0.01)
    airmass = np.linspace(1.0, 1.5, 80)

    result = build_fast_ultranest_lightcurve_series(times, flux, unc, airmass)

    assert result["applied"] is True
    assert result["original_point_count"] == 80
    assert result["binned_point_count"] <= 20
    assert result["time"].shape == result["flux"].shape == result["unc"].shape == result["airmass"].shape


def test_fast_ultranest_binning_skips_short_light_curve():
    times = np.linspace(0.0, 1.0, 60)
    result = build_fast_ultranest_lightcurve_series(
        times,
        np.ones(60),
        np.full(60, 0.01),
        np.linspace(1.0, 1.2, 60),
    )

    assert result["applied"] is False
    assert result["binned_point_count"] == 60


def test_expected_transit_coverage_assessment_flags_ingress_only_as_very_low():
    prior = {"tmid": 10.0, "per": 2.0, "rprs": 0.1, "ars": 12.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0}
    duration_prior = {"applied": True, "expected_duration": 0.1}
    times = np.linspace(9.90, 9.955, 12)

    assessment = build_expected_transit_coverage_assessment(
        times,
        prior,
        flux_values=np.ones(times.shape[0]),
        flux_errors=np.full(times.shape[0], 0.001),
        duration_prior=duration_prior,
    )

    assert assessment["valid"] is True
    assert assessment["observed_segment"] == "pre-ingress baseline plus ingress"
    assert assessment["transit_fraction_observed"] == pytest.approx(0.05)
    assert assessment["success_label"] == "very low"
    assert assessment["expected_successful"] is False


def test_final_fit_logs_partial_coverage_before_first_ultranest_call(monkeypatch):
    import exotic.exotic as exotic_module

    events = []

    def fake_log_info(message, warn=False, error=False):
        events.append(("log", str(message), warn))
        return True

    def fake_run_nested(times, flux_values, flux_errors, airmass, prior, bounds, **kwargs):
        events.append(("run_nested", "", False))
        local_times = np.asarray(times, dtype=float)
        model = np.ones(local_times.shape[0], dtype=float)
        model[-1:] -= 0.01
        fit = types.SimpleNamespace(
            time=local_times,
            data=np.asarray(flux_values, dtype=float),
            dataerr=np.asarray(flux_errors, dtype=float),
            model=model,
            residuals=np.zeros(local_times.shape[0], dtype=float),
            airmass=np.asarray(airmass, dtype=float),
            parameters={"rprs": prior["rprs"], "tmid": prior["tmid"], "inc": prior["inc"], "a2": 0.0, "per": prior["per"]},
            errors={"rprs": 0.01, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
            bounds=dict(bounds),
            duration_expected=0.1,
            duration_measured=0.1,
        )
        fit.get_parameter_posterior_recenter_diagnostics = (
            lambda key: {"clipped": False, "edge": None, "mode": fit.parameters.get(key, np.nan), "std": 0.01}
        )
        return fit

    monkeypatch.setattr(exotic_module, "log_info", fake_log_info)
    monkeypatch.setattr(exotic_module, "run_nested_lightcurve_fit_with_rprs_posterior_retry", fake_run_nested)
    monkeypatch.setattr(exotic_module, "apply_plot_time_range", lambda fit, plot_time_range: fit)
    monkeypatch.setattr(
        exotic_module,
        "build_final_fit_prefit_refinement_plan",
        lambda times, flux_values, flux_errors, airmass, prior, bounds, fit, **kwargs: {
            "applied": False,
            "note": "not needed",
            "times": np.asarray(times, dtype=float),
            "flux": np.asarray(flux_values, dtype=float),
            "unc": np.asarray(flux_errors, dtype=float),
            "airmass": np.asarray(airmass, dtype=float),
            "jd_times": None,
            "prior": dict(prior),
            "bounds": dict(bounds),
            "duration": 0.1,
            "original_point_count": len(times),
            "refined_point_count": len(times),
            "trimmed_pre_points": 0,
            "trimmed_post_points": 0,
            "original_tmid_bounds": bounds["tmid"],
            "refined_tmid_bounds": bounds["tmid"],
        },
    )

    times = np.linspace(9.90, 9.955, 12)
    fit, _, _ = fit_final_lightcurve_with_oot_baseline_detrending(
        times,
        np.ones(times.shape[0], dtype=float),
        np.full(times.shape[0], 0.001, dtype=float),
        np.linspace(1.0, 1.1, times.shape[0]),
        {"rprs": 0.1, "tmid": 10.0, "inc": 89.0, "a2": 0.0, "per": 2.0, "ars": 12.0, "ecc": 0.0, "omega": 0.0},
        {"rprs": [0.0, 0.5], "tmid": [9.95, 10.05], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]},
        detrend_on_outoftransit_baseline=False,
        duration_prior={"applied": True, "expected_duration": 0.1},
    )

    coverage_index = next(i for i, event in enumerate(events) if "Pre-UltraNest transit coverage assessment" in event[1])
    nested_index = next(i for i, event in enumerate(events) if event[0] == "run_nested")
    assert coverage_index < nested_index
    assert any("Estimated fit success: VERY LOW" in event[1] and event[2] for event in events)
    assert fit.pre_ultranest_transit_coverage_status == "very low"


def test_finalize_comparison_candidate_runs_pre_final_ultranest_on_binned_series(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {}

    def fake_fit_final(
        times,
        flux_values,
        flux_errors,
        airmass,
        prior,
        bounds,
        jd_times=None,
        **kwargs,
    ):
        captured["point_count"] = len(times)
        captured["bounds"] = dict(bounds)
        captured["fix_baseline_terms_for_final"] = kwargs.get("fix_baseline_terms_for_final")
        fit = types.SimpleNamespace(
            time=np.asarray(times, dtype=float),
            data=np.asarray(flux_values, dtype=float),
            dataerr=np.asarray(flux_errors, dtype=float),
            airmass=np.asarray(airmass, dtype=float),
            parameters={
                **dict(prior),
                "tmid": 0.5,
                "rprs": 0.1,
                "ars": 10.0,
                "inc": 89.0,
                "a0": 1.0,
                "a1": 1.0,
                "a2": 0.02,
            },
            errors={"tmid": 0.001, "rprs": 0.001, "ars": 0.1, "inc": 0.1, "a0": 0.01, "a2": 0.01},
            transit=np.ones(len(times), dtype=float),
            residuals=np.zeros(len(times), dtype=float),
            duration_measured=0.04,
            duration_expected=0.04,
        )
        return fit, np.asarray(flux_values, dtype=float), np.asarray(flux_errors, dtype=float)

    monkeypatch.setattr(exotic_module, "fit_final_lightcurve_with_oot_baseline_detrending", fake_fit_final)

    times = np.linspace(0.0, 1.0, 80)
    target_flux = 100.0 * (1.0 + 0.002 * np.sin(np.linspace(0.0, 2.0 * np.pi, 80)))
    result = finalize_comparison_candidate_full_reduction(
        times,
        target_flux,
        np.full(80, 100.0),
        np.linspace(1.0, 1.3, 80),
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={
            "pName": "Test b",
            "midT": 0.5,
            "midTUnc": 0.001,
            "pPer": 1.0,
            "pPerUnc": 0.001,
            "rprs": 0.1,
            "aRs": 10.0,
            "aRsUnc": 0.1,
            "inc": 89.0,
            "ecc": 0.0,
            "omega": 0.0,
        },
        jd_times=2460000.0 + times,
        run_fast_ultranest_before_final_run=True,
    )

    assert result["applied"] is True
    assert captured["point_count"] <= 20
    assert captured["fix_baseline_terms_for_final"] is False
    assert "a0" in captured["bounds"]
    assert "a2" in captured["bounds"]
    assert len(result["good_times"]) == result["fast_ultranest_binning"]["original_point_count"]
    assert len(result["good_times"]) > 60
    assert result["fast_ultranest_binning"]["applied"] is True
    assert result["fit"].fast_ultranest_binning_applied is True


def test_selected_fast_candidate_final_refit_uses_full_series_and_fixed_baseline(monkeypatch):
    import exotic.exotic as exotic_module

    monkeypatch.setenv("EXOTIC_ULTRANEST_MIN_NUM_LIVE_POINTS", "200")
    monkeypatch.setenv("EXOTIC_SPARSE_POSTERIOR_LIVE_POINT_RETRY", "1")
    captured = {}

    def fake_run_nested(
        times,
        flux_values,
        flux_errors,
        airmass,
        prior,
        bounds,
        jd_times=None,
        **kwargs,
    ):
        captured["point_count"] = len(times)
        captured["prior"] = dict(prior)
        captured["bounds"] = dict(bounds)
        captured["fixed_parameter_errors"] = dict(kwargs.get("fixed_parameter_errors", {}))
        captured["fixed_flux_baseline"] = kwargs.get("fixed_flux_baseline")
        captured["ultranest_min_num_live_points"] = kwargs.get("ultranest_min_num_live_points")
        captured["max_rprs_retries"] = kwargs.get("max_rprs_retries")
        captured["max_ars_retries"] = kwargs.get("max_ars_retries")
        captured["max_impact_parameter_retries"] = kwargs.get("max_impact_parameter_retries")
        fit = types.SimpleNamespace(
            time=np.asarray(times, dtype=float),
            data=np.asarray(flux_values, dtype=float),
            dataerr=np.asarray(flux_errors, dtype=float),
            airmass=np.asarray(airmass, dtype=float),
            parameters=dict(prior),
            errors=dict(kwargs.get("fixed_parameter_errors", {})),
            residuals=np.zeros(len(times), dtype=float),
            transit=np.ones(len(times), dtype=float),
            duration_measured=0.04,
            duration_expected=0.04,
            transit_qc={"status": "pass", "summary": "ok"},
            transit_qc_status="pass",
        )
        fit.get_parameter_posterior_samples = lambda key: np.linspace(0.0, 1.0, 1500)
        return fit

    monkeypatch.setattr(exotic_module, "run_nested_lightcurve_fit_with_rprs_posterior_retry", fake_run_nested)

    previous_fit = types.SimpleNamespace(
        fast_ultranest_binning_applied=True,
        parameters={
            "rprs": 0.1,
            "ars": 10.0,
            "per": 1.0,
            "tmid": 0.5,
            "inc": 89.0,
            "u0": 0.1,
            "u1": 0.1,
            "u2": 0.1,
            "u3": 0.1,
            "ecc": 0.0,
            "omega": 0.0,
            "a0": 1.03,
            "a1": 1.03,
            "a2": 0.12,
        },
        errors={"a0": 0.02, "a1": 0.02, "a2": 0.03, "rprs": 0.001, "tmid": 0.001, "ars": 0.1},
        bounds={
            "rprs": [0.05, 0.15],
            "tmid": [0.49, 0.51],
            "ars": [9.0, 11.0],
            "inc": [85.0, 90.0],
            "a0": [0.95, 1.05],
            "a2": [-3.0, 3.0],
        },
    )
    times = np.linspace(0.0, 1.0, 80)
    selected_result = {
        "fit": previous_fit,
        "good_times": times,
        "good_flux": np.ones(80),
        "good_unc": np.full(80, 0.01),
        "good_airmass": np.linspace(1.0, 1.3, 80),
        "good_jd_times": 2460000.0 + times,
    }

    returned = refit_selected_fast_comparison_on_full_lightcurve(
        selected_result,
        {
            "midT": 0.5,
            "midTUnc": 0.001,
            "pPer": 1.0,
            "rprs": 0.1,
            "aRs": 10.0,
            "inc": 89.0,
            "ecc": 0.0,
            "omega": 0.0,
        },
        detrend_on_outoftransit_baseline=False,
    )

    assert returned is not None
    assert captured["point_count"] == 80
    assert captured["fixed_flux_baseline"] is True
    assert captured["ultranest_min_num_live_points"] == 1200
    assert captured["max_rprs_retries"] == 0
    assert captured["max_ars_retries"] == 0
    assert captured["max_impact_parameter_retries"] == 0
    assert captured["prior"]["a0"] == pytest.approx(1.03)
    assert captured["prior"]["a2"] == pytest.approx(0.12)
    assert captured["fixed_parameter_errors"]["a0"] == pytest.approx(0.02)
    assert captured["fixed_parameter_errors"]["a2"] == pytest.approx(0.03)
    assert "a0" not in captured["bounds"]
    assert "a2" not in captured["bounds"]


def test_selected_fast_candidate_final_refit_resets_fixed_baseline_after_linear_detrend(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {}

    def fake_run_nested(
        times,
        flux_values,
        flux_errors,
        airmass,
        prior,
        bounds,
        jd_times=None,
        **kwargs,
    ):
        captured["flux"] = np.asarray(flux_values, dtype=float)
        captured["prior"] = dict(prior)
        captured["bounds"] = dict(bounds)
        captured["fixed_parameter_errors"] = dict(kwargs.get("fixed_parameter_errors", {}))
        captured["fixed_flux_baseline"] = kwargs.get("fixed_flux_baseline")
        fit = types.SimpleNamespace(
            time=np.asarray(times, dtype=float),
            data=np.asarray(flux_values, dtype=float),
            dataerr=np.asarray(flux_errors, dtype=float),
            airmass=np.asarray(airmass, dtype=float),
            parameters=dict(prior),
            errors=dict(kwargs.get("fixed_parameter_errors", {})),
            residuals=np.zeros(len(times), dtype=float),
            transit=np.ones(len(times), dtype=float),
            duration_measured=0.2,
            duration_expected=0.2,
            transit_qc={"status": "pass", "summary": "ok"},
            transit_qc_status="pass",
        )
        return fit

    monkeypatch.setattr(exotic_module, "run_nested_lightcurve_fit_with_rprs_posterior_retry", fake_run_nested)
    monkeypatch.setattr(exotic_module, "selected_final_live_point_target", lambda *args, **kwargs: (200, None))

    times = np.array([-2.0, -1.0, -0.25, 0.0, 0.25, 1.0, 2.0])
    transit_profile = np.array([1.0, 1.0, 1.0, 0.99, 1.0, 1.0, 1.0])
    baseline = 1.03 + 0.02 * times
    previous_fit = types.SimpleNamespace(
        fast_ultranest_binning_applied=True,
        transit=transit_profile,
        parameters={
            "rprs": 0.1,
            "ars": 10.0,
            "per": 1.0,
            "tmid": 0.0,
            "inc": 89.0,
            "u0": 0.1,
            "u1": 0.1,
            "u2": 0.1,
            "u3": 0.1,
            "ecc": 0.0,
            "omega": 0.0,
            "a0": 1.03,
            "a1": 1.03,
            "a2": 0.12,
        },
        errors={"a0": 0.02, "a1": 0.02, "a2": 0.03, "rprs": 0.001, "tmid": 0.001, "ars": 0.1},
        bounds={
            "rprs": [0.05, 0.15],
            "tmid": [-0.1, 0.1],
            "ars": [9.0, 11.0],
            "inc": [85.0, 90.0],
            "a0": [0.95, 1.05],
            "a2": [-3.0, 3.0],
        },
    )
    selected_result = {
        "fit": previous_fit,
        "good_times": times,
        "good_flux": baseline * transit_profile,
        "good_unc": np.full(times.shape, 0.01),
        "good_airmass": np.linspace(1.0, 1.3, times.size),
        "good_jd_times": 2460000.0 + times,
        "fast_fit_bounds": previous_fit.bounds,
    }

    returned, fit_flux, _ = refit_selected_fast_comparison_on_full_lightcurve(
        selected_result,
        {
            "midT": 0.0,
            "midTUnc": 0.001,
            "pPer": 1.0,
            "rprs": 0.1,
            "aRs": 10.0,
            "inc": 89.0,
            "ecc": 0.0,
            "omega": 0.0,
        },
        detrend_on_outoftransit_baseline=True,
    )

    assert returned is not None
    assert captured["fixed_flux_baseline"] is True
    assert captured["prior"]["a0"] == pytest.approx(1.0)
    assert captured["prior"]["a1"] == pytest.approx(1.0)
    assert captured["prior"]["a2"] == pytest.approx(0.0)
    assert captured["fixed_parameter_errors"]["a0"] == pytest.approx(0.02)
    assert captured["fixed_parameter_errors"]["a2"] == pytest.approx(0.03)
    assert np.allclose(captured["flux"][[0, 1, 2, 4, 5, 6]], 1.0, atol=1e-8)
    assert captured["flux"][3] == pytest.approx(0.99, abs=1e-8)
    assert np.allclose(fit_flux, captured["flux"])
    assert returned.oot_baseline_parameter_fit_applied is False
    assert "instead of reusing" in returned.oot_baseline_parameter_fit_note


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


def test_rprs_posterior_retry_expands_lower_edge_down_to_zero(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {"calls": []}
    diagnostics_sequence = [
        {"clipped": True, "edge": "lower", "mode": 0.030, "std": 0.006, "bounds": [0.000, 0.100]},
        {"clipped": False, "edge": None, "mode": 0.031, "std": 0.005, "bounds": [0.000, 0.120]},
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

    fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        np.linspace(-0.03, 0.03, 7),
        np.ones(7, dtype=float),
        np.full(7, 0.01, dtype=float),
        np.ones(7, dtype=float),
        {"tmid": 0.0, "rprs": 0.1, "inc": 89.0, "a2": 0.0},
        {"rprs": [0.025, 0.300], "tmid": [-0.01, 0.01], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]},
    )

    assert len(captured["calls"]) == 2
    assert captured["calls"][0]["bounds"]["rprs"] == pytest.approx([0.025, 0.300])
    assert captured["calls"][1]["prior"]["rprs"] == pytest.approx(0.030)
    assert captured["calls"][1]["bounds"]["rprs"] == pytest.approx([0.000, 0.120])
    assert fit.rprs_posterior_refit_applied is True
    assert fit.rprs_posterior_refit_count == 1
    assert fit.rprs_posterior_refit_edge == "lower"
    assert fit.rprs_posterior_refit_bounds == pytest.approx([0.000, 0.120])


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


def test_partial_coverage_suppresses_open_geometry_posterior_retries(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {"calls": []}
    diagnostics = {
        "rprs": {"clipped": False, "edge": None, "mode": 0.1, "std": 0.01, "bounds": [0.05, 0.15]},
        "ars": {"clipped": True, "edge": "lower", "mode": 5.0, "std": 3.0, "bounds": [0.000001, 20.0]},
        "b": {"clipped": True, "edge": "upper", "mode": 1.6, "std": 0.3, "bounds": [0.5, 2.5]},
    }

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
        captured["calls"].append({"prior": dict(call_prior), "bounds": dict(call_bounds)})
        fit = types.SimpleNamespace(
            sampled_keys=["rprs", "ars", "b", "tmid"],
            sample_bounds={"rprs": [0.0, 0.25], "ars": [0.000001, 20.0], "b": [0.0, 2.5]},
            parameters={"rprs": 0.1, "ars": 5.0, "tmid": 0.0, "inc": 80.0, "a2": 0.0},
        )
        fit.get_parameter_posterior_recenter_diagnostics = lambda key: dict(diagnostics[key])
        return fit

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        np.linspace(-0.03, 0.03, 7),
        np.ones(7, dtype=float),
        np.full(7, 0.01, dtype=float),
        np.ones(7, dtype=float),
        {"tmid": 0.0, "rprs": 0.1, "ars": 10.0, "inc": 89.0, "a2": 0.0},
        {
            "rprs": [0.0, 0.25],
            "ars": [5.0, 15.0],
            "tmid": [-0.01, 0.01],
            "inc": [70.0, 90.0],
            "a2": [-3.0, 3.0],
        },
        pre_ultranest_coverage_assessment={
            "valid": True,
            "success_label": "low",
            "expected_successful": False,
            "pre_ingress_points": 0,
            "post_egress_points": 8,
        },
    )

    assert len(captured["calls"]) == 1
    assert fit.ars_posterior_refit_applied is False
    assert "one-sided/LOW" in fit.ars_posterior_refit_note
    assert fit.b_posterior_refit_applied is False
    assert "one-sided/LOW" in fit.b_posterior_refit_note


def test_impact_parameter_posterior_retry_expands_inclination_bounds(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {"calls": []}
    diagnostics_sequence = [
        {
            "rprs": {"clipped": False, "edge": None, "mode": 0.1, "std": 0.01, "bounds": [0.05, 0.15]},
            "b": {"clipped": True, "edge": "upper", "mode": 1.6, "std": 0.18, "bounds": [0.7, 2.5]},
        },
        {
            "rprs": {"clipped": False, "edge": None, "mode": 0.1, "std": 0.01, "bounds": [0.05, 0.15]},
            "b": {"clipped": False, "edge": None, "mode": 1.6, "std": 0.12, "bounds": [0.7, 2.5]},
        },
    ]

    def make_fit(diagnostics):
        fit = types.SimpleNamespace(
            sampled_keys=["rprs", "b", "tmid"],
            sample_bounds={"rprs": [0.0, 0.25], "b": [0.0, 2.5], "tmid": [-0.01, 0.01]},
            parameters={
                "rprs": diagnostics["rprs"]["mode"],
                "ars": 10.0,
                "tmid": 0.0,
                "inc": 80.5,
                "a2": 0.0,
            },
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
            "use_impactparameter": use_impactparameter_rather_than_inclination_to_fit,
        })
        return make_fit(diagnostics_sequence[call_index])

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        np.linspace(-0.03, 0.03, 7),
        np.ones(7, dtype=float),
        np.full(7, 0.01, dtype=float),
        np.ones(7, dtype=float),
        {"tmid": 0.0, "rprs": 0.1, "ars": 10.0, "inc": 85.0, "a2": 0.0},
        {"rprs": [0.0, 0.25], "tmid": [-0.01, 0.01], "inc": [80.0, 90.0], "a2": [-3.0, 3.0]},
    )

    assert len(captured["calls"]) == 2
    assert captured["calls"][0]["bounds"]["inc"] == pytest.approx([80.0, 90.0])
    assert captured["calls"][1]["bounds"]["inc"][0] == pytest.approx(np.degrees(np.arccos(0.25)))
    assert captured["calls"][1]["bounds"]["inc"][1] == pytest.approx(90.0)
    assert captured["calls"][1]["prior"]["inc"] == pytest.approx(80.5)
    assert fit.b_posterior_refit_applied is True
    assert fit.b_posterior_refit_count == 1
    assert fit.b_posterior_refit_edge == "upper"
    assert fit.b_posterior_refit_bounds == pytest.approx([np.degrees(np.arccos(0.25)), 90.0])


def test_impact_parameter_retry_is_skipped_when_b_is_sampled_directly(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {"calls": []}
    diagnostics = {
        "rprs": {"clipped": False, "edge": None, "mode": 0.1, "std": 0.01, "bounds": [0.05, 0.15]},
        "b": {"clipped": True, "edge": "upper", "mode": 1.1, "std": 0.04, "bounds": [0.0, 1.12]},
    }

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
        captured["calls"].append({"prior": dict(call_prior), "bounds": dict(call_bounds)})
        fit = types.SimpleNamespace(
            sampled_keys=["rprs", "b", "tmid"],
            sample_bounds={"rprs": [0.0, 0.25], "b": [0.0, 1.12], "tmid": [-0.01, 0.01]},
            impact_parameter_sampled_directly=True,
            parameters={"rprs": 0.1, "ars": 10.0, "tmid": 0.0, "inc": 83.5, "a2": 0.0},
        )
        fit.get_parameter_posterior_recenter_diagnostics = lambda key: dict(diagnostics[key])
        return fit

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    fit = run_nested_lightcurve_fit_with_rprs_posterior_retry(
        np.linspace(-0.03, 0.03, 7),
        np.ones(7, dtype=float),
        np.full(7, 0.01, dtype=float),
        np.ones(7, dtype=float),
        {"tmid": 0.0, "rprs": 0.1, "ars": 10.0, "inc": 85.0, "a2": 0.0},
        {"rprs": [0.0, 0.25], "tmid": [-0.01, 0.01], "inc": [80.0, 90.0], "a2": [-3.0, 3.0]},
    )

    assert len(captured["calls"]) == 1
    assert fit.b_posterior_refit_applied is False
    assert fit.b_posterior_refit_count == 0


def test_impact_parameter_posterior_retry_expands_toward_face_on_boundary(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {"calls": []}
    diagnostics_sequence = [
        {
            "rprs": {"clipped": False, "edge": None, "mode": 0.1, "std": 0.01, "bounds": [0.05, 0.15]},
            "b": {"clipped": True, "edge": "lower", "mode": 0.55, "std": 0.10, "bounds": [0.0, 1.0]},
        },
        {
            "rprs": {"clipped": False, "edge": None, "mode": 0.1, "std": 0.01, "bounds": [0.05, 0.15]},
            "b": {"clipped": False, "edge": None, "mode": 0.55, "std": 0.08, "bounds": [0.0, 1.0]},
        },
    ]

    def make_fit(diagnostics):
        fit = types.SimpleNamespace(
            sampled_keys=["rprs", "b", "tmid"],
            sample_bounds={"rprs": [0.0, 0.25], "b": [0.0, 1.0], "tmid": [-0.01, 0.01]},
            parameters={"rprs": 0.1, "ars": 10.0, "tmid": 0.0, "inc": 86.0, "a2": 0.0},
        )
        fit.get_parameter_posterior_recenter_diagnostics = lambda key: dict(diagnostics[key])
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
        {"tmid": 0.0, "rprs": 0.1, "ars": 10.0, "inc": 84.0, "a2": 0.0},
        {"rprs": [0.0, 0.25], "tmid": [-0.01, 0.01], "inc": [80.0, 87.0], "a2": [-3.0, 3.0]},
    )

    assert len(captured["calls"]) == 2
    assert captured["calls"][1]["bounds"]["inc"] == pytest.approx([80.0, 90.0])
    assert fit.b_posterior_refit_applied is True
    assert fit.b_posterior_refit_edge == "lower"


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


def test_sparse_posterior_metric_flags_under_sampled_key_parameters():
    fit = types.SimpleNamespace()
    fit.get_parameter_posterior_samples = lambda key: np.linspace(0.0, 1.0, 100)

    diagnostics = evaluate_sparse_posterior_sample_support(
        fit,
        base_live_points=200,
    )

    assert diagnostics["sparse"] is True
    assert diagnostics["minimum_effective_samples"] == 1000
    assert diagnostics["parameters"]["rprs"]["effective_sample_count"] == pytest.approx(100)
    assert "rprs" in diagnostics["reason"]


def test_sparse_posterior_extension_continues_existing_ultranest_sampler(monkeypatch):
    monkeypatch.setenv("EXOTIC_ULTRANEST_MIN_NUM_LIVE_POINTS", "200")

    class SparseFit:
        def __init__(self):
            self.samples = {
                "rprs": np.linspace(0.09, 0.11, 100),
                "tmid": np.linspace(-0.001, 0.001, 100),
                "ars": np.linspace(9.5, 10.5, 100),
            }
            self.max_ncalls = 1000
            self.extension_calls = []
            self.cleared = False

        def get_parameter_posterior_samples(self, key):
            return self.samples[key]

        def extend_ultranest_fit(self, min_num_live_points=None, max_ncalls=None):
            self.extension_calls.append({
                "min_num_live_points": min_num_live_points,
                "max_ncalls": max_ncalls,
            })
            self.samples = {
                "rprs": np.linspace(0.09, 0.11, 1000),
                "tmid": np.linspace(-0.001, 0.001, 1000),
                "ars": np.linspace(9.5, 10.5, 1000),
            }
            return True

        def clear_ultranest_resume_state(self):
            self.cleared = True

    fit = SparseFit()
    returned = extend_sparse_posterior_live_points_if_needed(
        fit,
        enabled=True,
        extension_factor=SPARSE_POSTERIOR_LIVE_POINT_RETRY_FACTOR_DEFAULT,
    )

    assert returned is fit
    assert fit.extension_calls == [{
        "min_num_live_points": 1200,
        "max_ncalls": 6000,
    }]
    assert fit.sparse_posterior_live_point_extension_applied is True
    assert "200->1200" in fit.sparse_posterior_live_point_extension_note
    assert fit.cleared is True


def test_run_nested_lightcurve_fit_can_retain_sampler_for_final_extension(monkeypatch):
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
        keep_ultranest_sampler=False,
    ):
        captured["keep_ultranest_sampler"] = keep_ultranest_sampler
        fit = types.SimpleNamespace(parameters={"rprs": 0.1, "ars": 15.0, "tmid": 0.0, "inc": 89.0, "a2": 0.0})
        fit.get_parameter_posterior_recenter_diagnostics = lambda key: dict(diagnostics)
        return fit

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    run_nested_lightcurve_fit_with_rprs_posterior_retry(
        np.linspace(-0.03, 0.03, 7),
        np.ones(7, dtype=float),
        np.full(7, 0.01, dtype=float),
        np.ones(7, dtype=float),
        {"tmid": 0.0, "rprs": 0.1, "ars": 15.0, "inc": 89.0, "a2": 0.0},
        {"rprs": [0.0, 0.25], "ars": [14.5, 15.5], "tmid": [-0.01, 0.01], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]},
        keep_ultranest_sampler=True,
    )

    assert captured["keep_ultranest_sampler"] is True
