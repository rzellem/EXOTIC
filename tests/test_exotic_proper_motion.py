import importlib
import importlib.util
import json
import sys
import types
from pathlib import Path
import numpy as np
import pytest
from astropy.wcs import WCS


def _module_available(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except Exception:
        return False


def _set_stub_if_missing(name: str, module: types.ModuleType) -> None:
    if not _module_available(name):
        sys.modules.setdefault(name, module)


fake_barycorrpy = types.ModuleType("barycorrpy")
fake_utc_tdb = types.ModuleType("barycorrpy.utc_tdb")
fake_utc_tdb.JDUTC_to_BJDTDB = lambda *args, **kwargs: None
fake_barycorrpy.utc_tdb = fake_utc_tdb
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
try:
    importlib.import_module("barycorrpy.utc_tdb")
except Exception:
    sys.modules["barycorrpy"] = fake_barycorrpy
    sys.modules["barycorrpy.utc_tdb"] = fake_utc_tdb
sys.modules.setdefault("exotic.api.elca", fake_elca)
sys.modules.setdefault("exotic.api.ld", fake_ld)

from exotic.exotic import (
    APERTURE_MAX_FWHM_MULTIPLIER,
    APERTURE_MIN_FWHM_MULTIPLIER,
    APERTURE_SIGMA_MAX,
    APERTURE_SIGMA_MIN,
    GAUSSIAN_SIGMA_TO_FWHM,
    adaptive_aperture_outlier_mask,
    annotate_transit_qc_expected_values,
    aperture_contains_overexposed_pixel,
    auto_tune_aperture_sigma_grid,
    build_aperture_correction_profile,
    build_initial_ars_bounds,
    build_single_transit_duration_prior,
    build_target_fit_candidate_jobs,
    build_time_rejection_diagnostic,
    check_coordinates,
    cheap_lightcurve_prescore,
    centroid_offset_matches_reference,
    choose_centroid_seed_position,
    compute_star_aperture_grid,
    compute_transit_qc_ktmf,
    detect_aperture_correction_star_candidates,
    transit_qc_residual_scatter_score,
    transit_qc_residual_flatness_summary,
    transit_qc_sampling_summary,
    transit_qc_tmid_gaussianity_summary,
    apply_comparison_star_suitability_outlier_rejection,
    comparison_calibration_selection_reason,
    comparison_candidate_triangle_plot_output_path,
    comparison_candidate_fit_selection_reason,
    comparison_star_coverage_summary,
    comparison_star_stability_summary,
    compute_photometry_noise_budget,
    configure_windows_multiprocessing_main_spec,
    deduplicate_comparison_star_coords,
    diagnose_lightcurve_fit_inputs,
    detrend_flux_on_out_of_transit_baseline,
    alignment_candidate_quality_score,
    aperture_estimation_comparison_stars,
    aperture_frame_sigma_from_psf_data,
    apply_raw_target_photometry_selection,
    build_tracked_comparison_pool,
    collapse_aperture_data_to_selected_grid_cell,
    ensure_lightcurve_fit_failure_reason,
    evaluate_lightcurve_candidate,
    evaluate_transit_detection_qc,
    finalize_comparison_candidate_full_reduction,
    fit_lightcurve,
    fit_final_lightcurve_with_oot_baseline_detrending,
    fit_lightcurve_to_every_comparison_candidate,
    fit_ranked_comparison_calibration_candidates,
    get_final_fit_baseline_duration_multiplier,
    get_multiprocess_bad_pixel_precheck_processes,
    estimate_ephemeris_tmid_and_bounds,
    estimate_tmid_and_bounds_with_eebls,
    initialize_aperture_data_store,
    is_adaptive_aperture_mode_enabled,
    is_comp_star_required,
    is_out_of_transit_baseline_detrending_enabled,
    is_target_driven_comp_selection_enabled,
    limited_ensemble_comparison_keys,
    log_comparison_calibration_fit_attempt_summaries,
    log_comparison_candidate_fit_summaries,
    log_target_fit_candidate_summaries,
    noise_budget_config_from_info,
    normalize_flux_series_to_approximate_unity,
    parse_overexposure_threshold_fraction,
    parse_saturation_value,
    saturation_value_from_header,
    phase_bin_sigma_clip,
    parse_deviation_from_expected_transit_in_qc_sigma,
    parse_maximum_number_of_ensemble_comparisons_for_stellar_variability,
    parse_maximum_number_of_ensemble_comparisons_for_transit,
    prepare_final_fit_lightcurve_series,
    prepare_lightcurve_fit_input_series,
    project_comparison_radec_to_pixels,
    psf_frame_quality_components,
    psf_frame_quality_mask,
    psf_solution_quality_score,
    target_psf_shape_quality_components,
    target_psf_shape_quality_mask,
    target_comp_flux_scatter,
    fitted_lightcurve_scatter_on_dataset,
    populate_aperture_data_for_frame,
    rank_comparison_candidate_preflight_plans,
    refit_selected_fast_comparison_on_full_lightcurve,
    representative_psf_sigma,
    ranked_comparison_calibration_summaries,
    resolve_sky_annulus_geometry,
    run_target_driven_photometry_search,
    resolve_frame_aperture_radii,
    robust_flux_floor_mask,
    robust_target_reference_flux_mask,
    save_final_triangle_plot,
    save_selected_photometry_debug_series,
    select_comparison_calibrated_photometry,
    select_alignment_candidate,
    select_preferred_comparison_attempt,
    should_keep_header_wcs_alignment,
    should_prefer_pixel_values_over_wcs_for_target,
    sigma_clip,
    summarize_adaptive_aperture_usage,
    summarize_prior_transit_coverage,
    should_skip_airmass_fit,
    should_require_apparent_magnitudes,
    should_use_exactly_the_comps_provided,
    should_use_eebls_to_initialize_tmid_and_bounds,
    should_use_ensemble_photometry_for_stellar_variability,
    should_photometer_fortuitous_variables,
    should_reject_overexposed_stars,
    should_fit_lightcurve_to_every_comparison_candidate,
    should_detect_bad_pixels_before_photometry,
    should_use_aperture_photometry,
    should_use_aperture_corrections_and_full_image_fwhm,
    should_exit_at_first_qc_pass_solution,
    should_pick_comparison_by_eebls_snr,
    should_stop_after_promising_partial_comparison_attempt,
    should_use_psf_photometry,
    should_skip_low_comparison_coverage_rejection,
    should_use_fast_target_centroid,
    should_use_deviation_from_expected_transit_in_qc,
    update_coordinates_with_proper_motion,
    zoomed_final_triangle_plot_output_path,
)


def test_save_final_triangle_plot_regenerates_even_when_selected_candidate_artifact_exists(tmp_path):
    class DummyFigure:
        def savefig(self, path):
            Path(path).write_bytes(b"regenerated-final")

    class DummyFit:
        def __init__(self):
            self.called = False

        def plot_triangle(self):
            self.called = True
            return DummyFigure()

    planet_name = "TOI-1728 b"
    observation_date = "2024-12-14"
    source_dir = tmp_path / "comp6"
    final_dir = tmp_path / "final"
    source_temp = source_dir / "working_artifacts"
    source_temp.mkdir(parents=True)
    source_plot = source_temp / "Triangle_TOI-1728b_2024-12-14.png"
    source_plot.write_bytes(b"stale-selected-comp-6")

    fit = DummyFit()
    output_path = save_final_triangle_plot(
        fit,
        final_dir,
        planet_name,
        observation_date,
        source_dir=source_dir,
    )

    assert output_path == final_dir / "Diagnostics" / "FinalTriangle_TOI-1728b_2024-12-14.png"
    assert output_path.read_bytes() == b"regenerated-final"
    assert (final_dir / "Diagnostics" / "Triangle_TOI-1728b_2024-12-14.png").read_bytes() == b"regenerated-final"
    assert fit.called is True


def test_comparison_candidate_triangle_plot_uses_candidate_specific_name(tmp_path):
    output_path = comparison_candidate_triangle_plot_output_path(
        tmp_path / "comp7",
        "WASP-80 b",
        "2025-06-22",
        6,
    )

    assert output_path.name == "Comp7_Triangle_WASP-80b_2025-06-22.png"
    assert output_path.parent == tmp_path / "comp7" / "working_artifacts"


def test_comparison_candidate_triangle_plot_uses_date_only_from_timestamp(tmp_path):
    output_path = comparison_candidate_triangle_plot_output_path(
        tmp_path / "comp7",
        "XO-1/b",
        "2026-05-06T19:51:13.964-0700",
        6,
    )

    assert output_path.name == "Comp7_Triangle_XO-1-b_2026-05-06.png"


def test_zoomed_final_triangle_plot_uses_named_artifact(tmp_path):
    output_path = zoomed_final_triangle_plot_output_path(
        tmp_path / "final",
        "XO-1/b",
        "2026-05-06T19:51:13.964-0700",
    )

    assert output_path == tmp_path / "final" / "Diagnostics" / "ZoomedTrianglePlot_XO-1-b_2026-05-06.png"


def test_save_final_triangle_plot_creates_zoomed_companion_when_supported(tmp_path):
    class DummyFigure:
        def __init__(self, content):
            self.content = content

        def savefig(self, path):
            Path(path).write_bytes(self.content)

    class DummyFit:
        def __init__(self):
            self.calls = []

        def plot_triangle(self, plot_title=None, zoom_sigma=None):
            self.calls.append({"plot_title": plot_title, "zoom_sigma": zoom_sigma})
            content = b"zoomed" if zoom_sigma == 5.0 else b"full"
            return DummyFigure(content)

    fit = DummyFit()
    output_path = save_final_triangle_plot(
        fit,
        tmp_path / "final",
        "TOI-1728 b",
        "2024-12-14",
        source_dir=tmp_path / "comp4",
    )

    zoomed_output_path = zoomed_final_triangle_plot_output_path(
        tmp_path / "final",
        "TOI-1728 b",
        "2024-12-14",
    )
    assert output_path.read_bytes() == b"full"
    assert zoomed_output_path.read_bytes() == b"zoomed"
    assert fit.calls == [
        {"plot_title": "Final selected fit (comparison candidate #4)", "zoom_sigma": None},
        {"plot_title": "Final selected fit (comparison candidate #4) (5-sigma zoom)", "zoom_sigma": 5.0},
    ]


def test_save_final_triangle_plot_regenerates_when_selected_artifact_missing(tmp_path):
    class DummyFigure:
        def savefig(self, path):
            Path(path).write_bytes(b"regenerated")

    class DummyFit:
        def __init__(self):
            self.called = False

        def plot_triangle(self):
            self.called = True
            return DummyFigure()

    fit = DummyFit()
    output_path = save_final_triangle_plot(
        fit,
        tmp_path / "final",
        "TOI-1728 b",
        "2024-12-14",
        source_dir=tmp_path / "missing-comp",
    )

    assert fit.called is True
    assert output_path.read_bytes() == b"regenerated"
    assert output_path.parent == tmp_path / "final" / "Diagnostics"
    assert output_path.name == "FinalTriangle_TOI-1728b_2024-12-14.png"
    assert (
        tmp_path
        / "final"
        / "Diagnostics"
        / "Triangle_TOI-1728b_2024-12-14.png"
    ).read_bytes() == b"regenerated"


def test_save_final_triangle_plot_labels_selected_candidate_when_supported(tmp_path):
    class DummyFigure:
        def savefig(self, path):
            Path(path).write_bytes(b"regenerated")

    class DummyFit:
        def __init__(self):
            self.plot_title = None

        def plot_triangle(self, plot_title=None):
            self.plot_title = plot_title
            return DummyFigure()

    fit = DummyFit()
    output_path = save_final_triangle_plot(
        fit,
        tmp_path / "final",
        "TOI-1728 b",
        "2024-12-14",
        source_dir=tmp_path / "comp4",
    )

    assert output_path.read_bytes() == b"regenerated"
    assert fit.plot_title == "Final selected fit (comparison candidate #4)"


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


def test_prepare_final_fit_lightcurve_series_uses_two_sided_modeled_oot():
    times = np.linspace(-0.04, 0.04, 9)
    detrended = np.array([1.0, 1.0, 1.0, 0.99, 0.98, 0.99, 1.0, 1.0, 1.0], dtype=float)
    fit = types.SimpleNamespace(
        time=times,
        data=detrended.copy(),
        dataerr=np.full(times.shape, 0.01, dtype=float),
        detrended=detrended.copy(),
        detrendederr=np.full(times.shape, 0.01, dtype=float),
        airmass_model=np.ones(times.shape, dtype=float),
        transit=np.array([1.0, 1.0, 1.0, 0.99, 0.98, 0.99, 1.0, 1.0, 1.0], dtype=float),
        parameters={"tmid": 0.0},
    )

    prepared = prepare_final_fit_lightcurve_series(fit)

    assert prepared["applied"] is True
    assert prepared["used_two_sided_oot"] is True
    assert "modeled out-of-transit" in prepared["note"]
    assert prepared["flux"] == pytest.approx(detrended)
    assert np.all(prepared["unc"] > 0)


def test_prepare_final_fit_lightcurve_series_avoids_one_sided_oot_raw_flux_bias():
    times = np.linspace(-0.05, 0.05, 11)
    detrended = 1.0 - 0.02 * np.exp(-0.5 * (times / 0.012) ** 2)
    airmass_model = 1.0 + 2.0 * times
    raw_flux = detrended * airmass_model
    transit_model = np.where(times < 0.0, 1.0, 0.985)
    fit = types.SimpleNamespace(
        time=times,
        data=raw_flux,
        dataerr=np.full(times.shape, 0.01, dtype=float),
        detrended=detrended.copy(),
        detrendederr=np.full(times.shape, 0.01, dtype=float),
        airmass_model=airmass_model,
        transit=transit_model,
        parameters={"tmid": 0.0},
    )

    prepared = prepare_final_fit_lightcurve_series(fit)
    old_oot_mask = transit_model == 1.0
    legacy_flux = raw_flux / np.nanmedian(raw_flux[old_oot_mask])
    expected_flux = detrended / np.nanmedian(detrended)

    assert prepared["applied"] is True
    assert prepared["used_two_sided_oot"] is False
    assert "only bracketed one side of transit" in prepared["note"]
    assert prepared["flux"] == pytest.approx(expected_flux)
    assert prepared["flux"][-1] != pytest.approx(legacy_flux[-1], abs=1e-3)
    assert np.all(prepared["unc"] > 0)


def test_save_selected_photometry_debug_series_writes_stage_masks(tmp_path):
    fit = types.SimpleNamespace(
        selected_photometry_debug={
            "times": np.array([1.0, 2.0, 3.0], dtype=float),
            "target_flux": np.array([10.0, 11.0, 12.0], dtype=float),
            "comp_flux": np.array([5.0, 5.0, 6.0], dtype=float),
            "raw_ratio": np.array([2.0, 2.2, 2.0], dtype=float),
            "initial_sigma_keep_mask": np.array([True, False, True], dtype=bool),
            "prefit_raw_ratio_keep_mask": np.array([True, True, True], dtype=bool),
            "phase_clip_keep_mask_on_sigma_filtered": np.array([True, False], dtype=bool),
        }
    )

    output_path = save_selected_photometry_debug_series(tmp_path, "Qatar-10 b", "20260420", fit)

    assert output_path is not None
    assert output_path.exists()

    rows = np.loadtxt(output_path, delimiter=",", skiprows=1)
    assert rows.shape == (3, 10)
    assert np.isnan(rows[:, 4:7]).all()
    assert rows[:, 7].astype(int).tolist() == [1, 0, 1]
    assert rows[:, 8].astype(int).tolist() == [1, 1, 1]
    assert rows[:, 9].astype(int).tolist() == [1, 0, 0]


def test_finalize_comparison_candidate_phase_clips_before_nested_fit(monkeypatch):
    captured = {}

    def fake_lc_fitter(times, flux, unc, airmass, prior, bounds, jd_times=None, mode=None, **kwargs):
        assert mode == "lm"
        return types.SimpleNamespace(
            residuals=np.linspace(-0.01, 0.01, len(times)),
            phase=np.linspace(-0.5, 0.5, len(times)),
        )

    def fake_phase_clip(residuals, phase, sigma=3, bins=10):
        mask = np.zeros(len(residuals), dtype=bool)
        mask[3] = True
        return mask

    def fake_final_fit(
        times,
        flux,
        unc,
        airmass,
        prior,
        bounds,
        jd_times=None,
        **kwargs,
    ):
        captured["times"] = np.asarray(times, dtype=float).copy()
        captured["jd_times"] = np.asarray(jd_times, dtype=float).copy()
        fit = types.SimpleNamespace(
            time=np.asarray(times, dtype=float),
            airmass=np.asarray(airmass, dtype=float),
            data=np.asarray(flux, dtype=float),
            dataerr=np.asarray(unc, dtype=float),
            detrended=np.asarray(flux, dtype=float),
            detrendederr=np.asarray(unc, dtype=float),
            airmass_model=np.ones(len(times), dtype=float),
            transit=np.ones(len(times), dtype=float),
            phase=np.linspace(-0.5, 0.5, len(times)),
            residuals=np.zeros(len(times), dtype=float),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a1": 1.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a1": 0.01, "a2": 0.01},
        )
        return fit, np.asarray(flux, dtype=float), np.asarray(unc, dtype=float)

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr("exotic.exotic.phase_bin_sigma_clip", fake_phase_clip)
    monkeypatch.setattr("exotic.exotic.fit_final_lightcurve_with_oot_baseline_detrending", fake_final_fit)
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

    times = np.linspace(0.0, 0.09, 10)
    result = finalize_comparison_candidate_full_reduction(
        times,
        np.full(10, 100.0, dtype=float),
        np.full(10, 100.0, dtype=float),
        np.linspace(1.0, 1.2, 10),
        [0.1, 0.1, 0.1, 0.1],
        {
            "midT": 0.045,
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
    )

    assert result["applied"] is True
    assert captured["times"].tolist() == pytest.approx(np.delete(times, 3).tolist())
    assert result["source_indices"].tolist() == [0, 1, 2, 4, 5, 6, 7, 8, 9]
    assert any(
        diagnostic["stage"] == "Final-fit phase residual clip"
        and diagnostic["dropped_point_count"] == 1
        for diagnostic in result["fit"].frame_filter_diagnostics
    )
    assert result["fit"].selected_photometry_debug[
        "phase_clip_keep_mask_on_sigma_filtered"
    ].tolist() == [True, True, True, False, True, True, True, True, True, True]


def test_finalize_comparison_candidate_can_disable_phase_residual_clip(monkeypatch):
    captured = {"phase_clip_called": False}

    def fake_lc_fitter(times, flux, unc, airmass, prior, bounds, jd_times=None, mode=None, **kwargs):
        return types.SimpleNamespace(
            residuals=np.linspace(-0.01, 0.01, len(times)),
            phase=np.linspace(-0.5, 0.5, len(times)),
        )

    def fake_phase_clip(residuals, phase, sigma=3, bins=10):
        captured["phase_clip_called"] = True
        mask = np.zeros(len(residuals), dtype=bool)
        mask[3] = True
        return mask

    def fake_final_fit(times, flux, unc, airmass, prior, bounds, jd_times=None, **kwargs):
        captured["times"] = np.asarray(times, dtype=float).copy()
        fit = types.SimpleNamespace(
            time=np.asarray(times, dtype=float),
            airmass=np.asarray(airmass, dtype=float),
            data=np.asarray(flux, dtype=float),
            dataerr=np.asarray(unc, dtype=float),
            detrended=np.asarray(flux, dtype=float),
            detrendederr=np.asarray(unc, dtype=float),
            airmass_model=np.ones(len(times), dtype=float),
            transit=np.ones(len(times), dtype=float),
            phase=np.linspace(-0.5, 0.5, len(times)),
            residuals=np.zeros(len(times), dtype=float),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a1": 1.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a1": 0.01, "a2": 0.01},
        )
        return fit, np.asarray(flux, dtype=float), np.asarray(unc, dtype=float)

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr("exotic.exotic.phase_bin_sigma_clip", fake_phase_clip)
    monkeypatch.setattr("exotic.exotic.fit_final_lightcurve_with_oot_baseline_detrending", fake_final_fit)
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

    times = np.linspace(0.0, 0.09, 10)
    result = finalize_comparison_candidate_full_reduction(
        times,
        np.full(10, 100.0, dtype=float),
        np.full(10, 100.0, dtype=float),
        np.linspace(1.0, 1.2, 10),
        [0.1, 0.1, 0.1, 0.1],
        {
            "midT": 0.045,
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
        run_final_fit_phase_residual_clip=False,
    )

    assert result["applied"] is True
    assert captured["phase_clip_called"] is False
    assert captured["times"].tolist() == pytest.approx(times.tolist())
    assert result["source_indices"].tolist() == list(range(10))
    assert not any(
        diagnostic["stage"] == "Final-fit phase residual clip"
        and diagnostic["dropped_point_count"] > 0
        for diagnostic in result["fit"].frame_filter_diagnostics
    )
    assert result["fit"].selected_photometry_debug[
        "phase_clip_keep_mask_on_sigma_filtered"
    ].tolist() == [True] * 10


def test_finalize_comparison_candidate_keeps_flux_aligned_after_final_fit_subsets_times(monkeypatch):
    fit_calls = []
    initial_subset_mask = np.ones(38, dtype=bool)
    initial_subset_mask[[1, 3, 5, 7, 9, 11]] = False

    def fake_lc_fitter(times, flux, unc, airmass, prior, bounds, jd_times=None, mode=None, **kwargs):
        return types.SimpleNamespace(
            residuals=np.zeros(len(times), dtype=float),
            phase=np.linspace(-0.5, 0.5, len(times)),
        )

    def fake_final_fit(times, flux, unc, airmass, prior, bounds, jd_times=None, **kwargs):
        times = np.asarray(times, dtype=float)
        flux = np.asarray(flux, dtype=float)
        unc = np.asarray(unc, dtype=float)
        airmass = np.asarray(airmass, dtype=float)
        fit_calls.append({
            "times": times.copy(),
            "flux": flux.copy(),
            "unc": unc.copy(),
            "airmass": airmass.copy(),
        })

        if len(fit_calls) == 1:
            keep_mask = initial_subset_mask
            residuals = np.zeros(np.count_nonzero(keep_mask), dtype=float)
            residuals[10] = 1.0
        else:
            keep_mask = np.ones(len(times), dtype=bool)
            residuals = np.zeros(len(times), dtype=float)

        retained_times = times[keep_mask]
        retained_flux = flux[keep_mask]
        retained_unc = unc[keep_mask]
        retained_airmass = airmass[keep_mask]
        fit = types.SimpleNamespace(
            time=retained_times,
            airmass=retained_airmass,
            data=retained_flux,
            dataerr=retained_unc,
            detrended=retained_flux,
            detrendederr=retained_unc,
            airmass_model=np.ones(len(retained_times), dtype=float),
            transit=np.ones(len(retained_times), dtype=float),
            phase=np.linspace(-0.5, 0.5, len(retained_times)),
            residuals=residuals,
            parameters={
                "tmid": 0.5,
                "rprs": 0.1,
                "ars": 10.0,
                "inc": 89.0,
                "a0": 1.0,
                "a1": 1.0,
                "a2": 0.0,
            },
            errors={
                "tmid": 0.001,
                "rprs": 0.001,
                "ars": 0.1,
                "inc": 0.1,
                "a0": 0.01,
                "a1": 0.01,
                "a2": 0.01,
            },
        )
        return fit, retained_flux, retained_unc

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr("exotic.exotic.fit_final_lightcurve_with_oot_baseline_detrending", fake_final_fit)
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

    times = np.linspace(0.0, 1.0, 38)
    result = finalize_comparison_candidate_full_reduction(
        times,
        np.full(38, 100.0, dtype=float),
        np.full(38, 100.0, dtype=float),
        np.linspace(1.0, 1.2, 38),
        [0.1, 0.1, 0.1, 0.1],
        {
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
        run_fast_ultranest_before_final_run=False,
        run_final_fit_phase_residual_clip=False,
        run_final_residual_rejection=True,
        use_eebls_to_initialize_tmid_and_bounds=False,
    )

    initially_retained_indices = np.flatnonzero(initial_subset_mask)
    expected_source_indices = np.delete(initially_retained_indices, 10)
    assert result["applied"] is True
    assert [len(call["times"]) for call in fit_calls] == [38, 31]
    assert all(call["times"].shape == call["flux"].shape == call["unc"].shape for call in fit_calls)
    assert result["source_indices"].tolist() == expected_source_indices.tolist()
    assert result["good_times"].tolist() == pytest.approx(times[expected_source_indices].tolist())
    assert all(
        len(result[key]) == 31
        for key in (
            "good_times",
            "good_flux",
            "good_unc",
            "good_airmass",
            "good_jd_times",
            "good_target_flux",
            "good_comp_flux",
            "good_target_flux_error",
            "good_comp_flux_error",
            "source_indices",
        )
    )


def test_detrend_flux_on_out_of_transit_baseline_falls_back_to_prior_ephemeris():
    times = np.linspace(-0.08, 0.08, 17)
    baseline = 1.0 + 0.25 * times
    transit = np.ones_like(times)
    transit[(times >= 0.03) & (times <= 0.05)] = 0.99
    flux = baseline.copy()
    unc = np.full_like(times, 0.01)
    fit = types.SimpleNamespace(
        transit=transit,
        parameters={"tmid": 0.10},
    )
    prior = {
        "tmid": 0.0,
        "per": 1.0,
        "rprs": 0.1,
        "ars": 12.0,
        "inc": 88.0,
        "ecc": 0.0,
        "omega": 0.0,
    }

    modeled = detrend_flux_on_out_of_transit_baseline(times, flux, unc, fit, min_side_points=2)
    fallback = detrend_flux_on_out_of_transit_baseline(times, flux, unc, fit, prior=prior, min_side_points=2)
    prior_coverage = summarize_prior_transit_coverage(times, prior, flux_values=flux, flux_errors=unc)

    assert modeled["applied"] is False
    assert prior_coverage["valid"] is True
    assert prior_coverage["has_two_sided_oot"] is True
    assert fallback["applied"] is True
    assert fallback["used_prior_ephemeris"] is True
    assert "ephemeris-centered transit window" in fallback["note"]


def test_detrend_flux_on_out_of_transit_baseline_requires_default_side_coverage():
    times = np.linspace(-0.08, 0.08, 17)
    baseline = 1.0 + 0.25 * times
    transit_profile = np.ones_like(times)
    transit_profile[(times >= -0.01) & (times <= 0.01)] = 0.99
    flux = baseline * transit_profile
    unc = np.full_like(times, 0.01)
    fit = types.SimpleNamespace(
        transit=transit_profile,
        parameters={"tmid": 0.0},
    )

    result = detrend_flux_on_out_of_transit_baseline(times, flux, unc, fit)

    assert result["applied"] is False
    assert result["pre_points"] <= 12
    assert result["post_points"] <= 12
    assert "need more than 12 on each side" in result["note"]


def test_deduplicate_comparison_star_coords_preserves_nearby_user_stars():
    preserved_coords, duplicate_messages = deduplicate_comparison_star_coords(
        [
            [1826.0, 1499.0],
            [1827.0, 1511.0],
            [1828.0, 1487.0],
            [842.0, 1810.0],
        ],
        min_separation_pixels=15.0,
    )

    assert preserved_coords == [
        [1826.0, 1499.0],
        [1827.0, 1511.0],
        [1828.0, 1487.0],
        [842.0, 1810.0],
    ]
    assert duplicate_messages == []


def test_robust_flux_floor_mask_rejects_tiny_positive_outliers():
    flux = np.full(30, 100.0)
    flux[[5, 17]] = [1.0, 0.5]

    mask = robust_flux_floor_mask(flux)

    assert mask.sum() == 28
    assert not mask[5]
    assert not mask[17]


def test_robust_target_reference_flux_mask_requires_both_series_to_be_plausible():
    target_flux = np.full(30, 100.0)
    reference_flux = np.full(30, 120.0)
    target_flux[7] = 1.0
    reference_flux[13] = 1.0

    mask = robust_target_reference_flux_mask(target_flux, reference_flux)

    assert mask.sum() == 28
    assert not mask[7]
    assert not mask[13]


def test_psf_frame_quality_mask_rejects_high_seeing_and_low_amplitude_outliers():
    frame_count = 30
    phase = np.linspace(0.0, 2.0 * np.pi, frame_count)
    psf_rows = np.zeros((frame_count, 7), dtype=float)
    psf_rows[:, 0] = 10.0
    psf_rows[:, 1] = 20.0
    psf_rows[:, 2] = 200.0 * (1.0 + 0.02 * np.sin(phase))
    psf_rows[:, 3] = 1.1 * (1.0 + 0.01 * np.cos(phase))
    psf_rows[:, 4] = 1.0 * (1.0 + 0.01 * np.sin(phase))
    psf_rows[5, 2] = 45.0
    psf_rows[12, 3:5] = 6.5

    components = psf_frame_quality_components(psf_rows)
    mask = psf_frame_quality_mask(psf_rows)

    assert mask.sum() == frame_count - 2
    assert not mask[5]
    assert not mask[12]
    assert components["amplitude_outlier_mask"][5]
    assert components["seeing_outlier_mask"][12]


def test_target_psf_shape_quality_rejects_broad_target_but_preserves_amplitude_dips():
    frame_count = 30
    phase = np.linspace(0.0, 2.0 * np.pi, frame_count)
    target_rows = np.zeros((frame_count, 7), dtype=float)
    target_rows[:, 0] = 10.0
    target_rows[:, 1] = 20.0
    target_rows[:, 2] = 200.0 * (1.0 + 0.02 * np.sin(phase))
    target_rows[:, 3] = 1.0
    target_rows[:, 4] = 1.0
    target_rows[5, 2] = 45.0
    target_rows[12, 3:5] = 6.5

    reference_rows = target_rows.copy()
    reference_rows[:, 2] = 250.0
    reference_rows[:, 3:5] = 1.0

    components = target_psf_shape_quality_components(target_rows, reference_rows)
    mask = target_psf_shape_quality_mask(target_rows, reference_rows)

    assert mask.sum() == frame_count - 1
    assert mask[5]
    assert not mask[12]
    assert not components["invalid_mask"][5]
    assert components["seeing_outlier_mask"][12]
    assert components["reference_width_outlier_mask"][12]


def test_build_target_fit_candidate_jobs_masks_pairwise_psf_failures_but_preserves_target_dips():
    frame_count = 30

    def build_psf_rows(amplitudes):
        psf_rows = np.zeros((frame_count, 7), dtype=float)
        psf_rows[:, 0] = 10.0
        psf_rows[:, 1] = 20.0
        psf_rows[:, 2] = amplitudes
        psf_rows[:, 3] = 1.0
        psf_rows[:, 4] = 1.0
        return psf_rows

    target_amplitudes = np.full(frame_count, 100.0)
    comp_amplitudes = np.full(frame_count, 120.0)
    target_amplitudes[7] = 80.0
    comp_amplitudes[13] = 1.0

    psf_data = {
        "target": build_psf_rows(target_amplitudes),
        "comp1": build_psf_rows(comp_amplitudes),
    }
    psf_data["target"][19, 3:5] = 6.5

    candidate_jobs = build_target_fit_candidate_jobs(
        psf_data,
        aper_data=None,
        apers=None,
        annuli=None,
        airmass=np.linspace(1.0, 1.3, frame_count),
        comp_stars=[[1827.0, 1511.0]],
        sigma=3.0,
        require_comp_star=True,
        skip_low_comparison_coverage_rejection=False,
        use_psf_photometry=True,
        use_aperture_photometry=False,
    )

    assert len(candidate_jobs) == 1
    assert candidate_jobs[0]["method"] == "psf"
    assert candidate_jobs[0]["mask"].sum() == 28
    assert candidate_jobs[0]["mask"][7]
    assert not candidate_jobs[0]["mask"][13]
    assert not candidate_jobs[0]["mask"][19]
    assert candidate_jobs[0]["coverage_count"] == 29


def test_build_target_fit_candidate_jobs_uses_psf_flux_rows_for_psf_quality():
    frame_count = 30

    def build_psf_rows(amplitudes):
        psf_rows = np.zeros((frame_count, 7), dtype=float)
        psf_rows[:, 0] = 10.0
        psf_rows[:, 1] = 20.0
        psf_rows[:, 2] = amplitudes
        psf_rows[:, 3] = 1.0
        psf_rows[:, 4] = 1.0
        return psf_rows

    psf_data = {
        "target": build_psf_rows(np.full(frame_count, 100.0)),
        "comp1": build_psf_rows(np.full(frame_count, 120.0)),
    }
    psf_data["target"][19, 3:5] = 6.5

    psf_flux_data = {
        "target": build_psf_rows(np.full(frame_count, 100.0)),
        "comp1": build_psf_rows(np.full(frame_count, 120.0)),
    }
    psf_flux_data["target"][7, 2] = 80.0
    psf_flux_data["target"][21, 3:5] = 6.5
    psf_flux_data["comp1"][13, 2] = 1.0

    candidate_jobs = build_target_fit_candidate_jobs(
        psf_data,
        aper_data=None,
        apers=None,
        annuli=None,
        airmass=np.linspace(1.0, 1.3, frame_count),
        comp_stars=[[1827.0, 1511.0]],
        sigma=3.0,
        require_comp_star=True,
        skip_low_comparison_coverage_rejection=False,
        use_psf_photometry=True,
        use_aperture_photometry=False,
        psf_flux_data=psf_flux_data,
    )

    assert len(candidate_jobs) == 1
    assert candidate_jobs[0]["method"] == "psf"
    assert candidate_jobs[0]["mask"].sum() == 28
    assert candidate_jobs[0]["mask"][7]
    assert candidate_jobs[0]["mask"][19]
    assert not candidate_jobs[0]["mask"][13]
    assert not candidate_jobs[0]["mask"][21]
    assert candidate_jobs[0]["coverage_count"] == 29


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_legacy_psf_photometry_flux_row_uses_weighted_centroid_override():
    import exotic.exotic as exotic_module

    y_grid, x_grid = np.mgrid[0:31, 0:31]
    data = exotic_module.gaussian_psf(
        x_grid,
        y_grid,
        15.25,
        14.65,
        200.0,
        1.8,
        2.2,
        0.05,
        30.0,
    )
    seed_row = np.array([15.0, 15.0, 100.0, 1.0, 1.0, 0.0, 30.0], dtype=float)

    row = exotic_module.fit_legacy_psf_photometry_flux_row(data, seed_row, 0, box=8)

    xv, yv = exotic_module.mesh_box(seed_row[:2], 8, maxx=data.shape[1], maxy=data.shape[0])
    subarray = data[yv, xv]
    expected_wx = np.sum(xv[0] * subarray.sum(0)) / subarray.sum(0).sum()
    expected_wy = np.sum(yv[:, 0] * subarray.sum(1)) / subarray.sum(1).sum()

    assert row[0] == pytest.approx(expected_wx)
    assert row[1] == pytest.approx(expected_wy)
    assert row[2] > 0
    assert row[3] > 0
    assert row[4] > 0


def test_load_psf_flux_seed_tracks_accepts_legacy_selected_comp_file(tmp_path):
    import exotic.exotic as exotic_module

    run_dir = tmp_path / "old_run"
    temp_dir = run_dir / "temp"
    temp_dir.mkdir(parents=True)
    target_rows = np.tile(np.array([[10.0, 20.0, 100.0, 1.0, 1.1, 0.0, 30.0]]), (3, 1))
    comp_rows = np.tile(np.array([[30.0, 40.0, 150.0, 1.2, 1.3, 0.0, 31.0]]), (3, 1))
    np.savetxt(temp_dir / "psf_data_target.txt", target_rows)
    np.savetxt(temp_dir / "psf_data_comp.txt", comp_rows)

    seed_tracks = exotic_module.load_psf_flux_seed_tracks(str(run_dir), 3, ["comp1"])

    assert set(seed_tracks) == {"target", "comp1"}
    assert seed_tracks["target"].shape == (3, 7)
    assert seed_tracks["comp1"].shape == (3, 7)
    assert seed_tracks["target"][0, 0] == pytest.approx(10.0)
    assert seed_tracks["comp1"][0, 0] == pytest.approx(30.0)


def test_centroid_offset_matches_reference_uses_float_geometry_tolerance():
    target = np.array([2383.27, 867.04, 6.3, 9.0, 0.7, 0.0, 223.0])
    comp = np.array([1821.90, 549.21, 131.0, 3.0, 4.9, 0.0, 226.0])

    assert centroid_offset_matches_reference(comp, target, 562.88, 315.42)


def test_should_keep_header_wcs_alignment_prefers_geometry_over_flux_swings():
    decision = should_keep_header_wcs_alignment(
        projected_off_frame=False,
        frame_index=5,
        target_psf_row=np.array([2383.27, 867.04, 6.3, 9.0, 0.7, 0.0, 223.0]),
        previous_target_psf_row=np.array([2382.88, 454.14, 154.0, 2.8, 2.7, 0.0, 223.0]),
        comp_psf_rows={
            "comp1": np.array([1821.90, 549.21, 131.0, 3.0, 4.9, 0.0, 226.0]),
        },
        previous_comp_psf_rows={
            "comp1": np.array([1822.57, 135.19, 3645.8, 2.9, 5.8, 0.0, 226.0]),
        },
        expected_offsets={
            "comp1": np.array([562.88, 315.42]),
        },
    )

    assert decision["use_wcs_alignment"] is True
    assert decision["reason"] == "geometry_match"
    assert not decision["target_flux_change_ok"]
    assert not decision["comp_flux_change_ok"]
    assert decision["geometry_match_count"] == 1


def test_should_keep_header_wcs_alignment_ignores_one_bad_comp_when_majority_match():
    decision = should_keep_header_wcs_alignment(
        projected_off_frame=False,
        frame_index=8,
        target_psf_row=np.array([2387.50, 868.60, 50.0, 3.0, 3.0, 0.0, 223.0]),
        previous_target_psf_row=np.array([2382.88, 454.14, 154.0, 2.8, 2.7, 0.0, 223.0]),
        comp_psf_rows={
            "comp1": np.array([1827.0, 558.1, 200.0, 3.0, 5.0, 0.0, 226.0]),
            "comp2": np.array([1300.0, 900.0, 300.0, 3.0, 5.0, 0.0, 226.0]),
        },
        previous_comp_psf_rows={
            "comp1": np.array([1822.6, 135.2, 3645.8, 2.9, 5.8, 0.0, 226.0]),
            "comp2": np.array([1600.0, 700.0, 280.0, 3.0, 5.0, 0.0, 226.0]),
        },
        expected_offsets={
            "comp1": np.array([560.7, 310.2]),
            "comp2": np.array([1187.5, 31.4]),
        },
    )

    assert decision["use_wcs_alignment"] is True
    assert decision["geometry_test_count"] == 2
    assert decision["geometry_match_count"] == 1


def test_should_keep_header_wcs_alignment_rejects_geometry_mismatch():
    decision = should_keep_header_wcs_alignment(
        projected_off_frame=False,
        frame_index=8,
        target_psf_row=np.array([1576.0, 1317.0, 1.1, 20.0, 1.9, 0.0, 223.0]),
        previous_target_psf_row=np.array([2382.88, 454.14, 154.0, 2.8, 2.7, 0.0, 223.0]),
        comp_psf_rows={
            "comp1": np.array([1826.7, 1510.9, 1.9, 15.4, 1.5, 0.0, 226.0]),
        },
        previous_comp_psf_rows={
            "comp1": np.array([1822.57, 135.19, 3645.8, 2.9, 5.8, 0.0, 226.0]),
        },
        expected_offsets={
            "comp1": np.array([562.88, 315.42]),
        },
    )

    assert decision["use_wcs_alignment"] is False
    assert decision["reason"] == "geometry_mismatch"
    assert decision["geometry_match_count"] == 0


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


def test_check_coordinates_non_interactive_keeps_input_when_wcs_psf_is_implausible():
    x_pixel, y_pixel = check_coordinates(
        input_x_pixel=246,
        input_y_pixel=271,
        centroid_x=238.5,
        centroid_y=266.1,
        sigma_x=4.6,
        sigma_y=0.7,
        calculated_x_pixel=245,
        calculated_y_pixel=270,
        non_interactive_run=True,
        wcs_psf_quality_score=np.inf,
        input_psf_quality_score=0.1,
    )

    assert x_pixel == 246
    assert y_pixel == 271


def test_check_coordinates_non_interactive_keeps_plausible_input_when_wcs_finds_other_source():
    x_pixel, y_pixel = check_coordinates(
        input_x_pixel=246,
        input_y_pixel=271,
        centroid_x=236.6,
        centroid_y=265.2,
        sigma_x=1.2,
        sigma_y=0.8,
        calculated_x_pixel=236,
        calculated_y_pixel=265,
        non_interactive_run=True,
        wcs_psf_quality_score=0.05,
        input_psf_quality_score=0.25,
    )

    assert x_pixel == 246
    assert y_pixel == 271


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


def test_check_coordinates_can_prefer_input_pixels_over_wcs_conflict():
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
        prefer_pixel_values_over_wcs_for_target="y",
    )

    assert x_pixel == 5
    assert y_pixel == 5


def test_should_prefer_pixel_values_over_wcs_for_target_parses_values():
    assert should_prefer_pixel_values_over_wcs_for_target(None) is False
    assert should_prefer_pixel_values_over_wcs_for_target("n") is False
    assert should_prefer_pixel_values_over_wcs_for_target("y") is True
    assert should_prefer_pixel_values_over_wcs_for_target(True) is True


def test_psf_solution_quality_score_rejects_offset_or_elongated_solutions():
    good = np.array([246.2, 270.8, 140.0, 1.1, 0.9, 0.0, 40.0])
    offset = np.array([238.5, 266.1, 140.0, 1.1, 0.9, 0.0, 40.0])
    elongated = np.array([246.2, 270.8, 140.0, 4.6, 0.7, 0.0, 40.0])

    assert np.isfinite(psf_solution_quality_score(good, seed_pos=[246.0, 271.0]))
    assert not np.isfinite(psf_solution_quality_score(offset, seed_pos=[246.0, 271.0]))
    assert not np.isfinite(psf_solution_quality_score(elongated, seed_pos=[246.0, 271.0]))


def test_alignment_candidate_selection_rejects_broad_offset_wcs_target_solution():
    psf_data = {
        "target": np.zeros((2, 7), dtype=float),
        "comp1": np.zeros((2, 7), dtype=float),
    }
    psf_data["target"][0] = [245.8, 270.7, 110.0, 1.1, 0.8, 0.0, 40.0]
    psf_data["comp1"][0] = [360.3, 443.3, 174.0, 1.1, 1.0, 0.0, 40.0]
    tar_comp_dist = {"comp1": np.array([115.0, 173.0])}

    wcs_candidate = {
        "coords": np.array([[244.0, 268.7], [360.4, 443.2]], dtype=float),
        "projected_off_frame": False,
        "psf_rows": {
            "target": np.array([244.0, 268.7, 210.0, 7.3, 5.6, 0.0, 40.0]),
            "comp1": np.array([360.4, 443.2, 174.0, 1.1, 1.0, 0.0, 40.0]),
        },
        "warnings": [],
    }
    fallback_candidate = {
        "coords": np.array([[246.0, 270.8], [360.2, 443.3]], dtype=float),
        "psf_rows": {
            "target": np.array([245.9, 270.8, 111.0, 1.1, 0.8, 0.0, 40.0]),
            "comp1": np.array([360.2, 443.3, 173.0, 1.1, 1.0, 0.0, 40.0]),
        },
        "warnings": [],
    }

    assert not np.isfinite(
        alignment_candidate_quality_score(
            wcs_candidate,
            comp_keys=["comp1"],
            previous_target_psf_row=psf_data["target"][0],
            previous_comp_psf_rows={"comp1": psf_data["comp1"][0]},
            expected_offsets=tar_comp_dist,
        )
    )

    selected_source, selected_candidate, diagnostics = select_alignment_candidate(
        {"wcs": wcs_candidate, "fallback": fallback_candidate, "file_name": "frame.fits"},
        frame_index=1,
        psf_data=psf_data,
        tar_comp_dist=tar_comp_dist,
        comp_keys=["comp1"],
    )

    assert selected_source == "fallback"
    assert selected_candidate is fallback_candidate
    assert not np.isfinite(diagnostics["wcs_score"])
    assert np.isfinite(diagnostics["fallback_score"])


def test_is_comp_star_required_parses_values():
    assert is_comp_star_required(None) is True
    assert is_comp_star_required("y") is True
    assert is_comp_star_required("n") is False


def test_mixed_exposure_times_keep_target_only_allowed_with_scaling_warning(monkeypatch):
    import exotic.exotic as exotic_module

    messages = []
    monkeypatch.setattr(
        exotic_module,
        "log_info",
        lambda message, **kwargs: messages.append((message, kwargs)),
    )

    assert exotic_module.exposure_time_spread_fraction([60.0, 60.3, 60.5]) < 0.01
    assert not exotic_module.exposure_variation_requires_comp_star([60.0, 60.3, 60.5])
    assert exotic_module.resolve_require_comp_star_for_exposure_times("n", [60.0, 60.3, 60.5]) is False

    assert exotic_module.exposure_time_spread_fraction([60.0, 61.0]) > 0.01
    assert exotic_module.exposure_variation_requires_comp_star([60.0, 61.0])
    assert exotic_module.resolve_require_comp_star_for_exposure_times("n", [60.0, 61.0]) is False
    assert exotic_module.resolve_require_comp_star_for_exposure_times("y", [60.0, 61.0]) is True
    assert any("scale source counts to a common exposure time" in message for message, _ in messages)


def test_img_time_bjd_tdb_prefers_direct_mid_exposure_bjd(monkeypatch):
    import exotic.exotic as exotic_module

    header = exotic_module.fits.Header()
    header["BJD_TDB"] = 2461152.1287422837
    header["DATE-AVG"] = "2026-04-22T15:05:23.333333"
    header["DATE-UTC"] = "2026-04-22T15:05:08.333333"
    header["EXPTIME"] = 30.0

    monkeypatch.setattr(
        exotic_module,
        "convert_jd_to_bjd",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("conversion should not run")),
    )

    assert exotic_module.img_time_bjd_tdb(header, {}, {}) == pytest.approx(2461152.1287422837)


def test_img_time_bjd_tdb_uses_nina_date_avg_before_start_time(monkeypatch):
    import exotic.exotic as exotic_module

    header = exotic_module.fits.Header()
    header["DATE-AVG"] = "2026-04-22T15:05:23.333333"
    header["DATE-UTC"] = "2026-04-22T15:05:08.333333"
    header["EXPTIME"] = 10.0

    converted_inputs = []

    def fake_convert_jd_to_bjd(values, _p_dict, _info_dict):
        converted_inputs.extend(values)
        return np.asarray(values, dtype=float) + 0.25

    monkeypatch.setattr(exotic_module, "convert_jd_to_bjd", fake_convert_jd_to_bjd)

    expected_midpoint_jd = exotic_module.Time("2026-04-22T15:05:23.333333", scale="utc").jd
    expected_start_plus_exposure_jd = (
        exotic_module.Time("2026-04-22T15:05:08.333333", scale="utc").jd
        + 5.0 / 86400.0
    )

    assert exotic_module.img_time_jd(header) == pytest.approx(expected_midpoint_jd)
    assert exotic_module.img_time_bjd_tdb(header, {}, {}) == pytest.approx(expected_midpoint_jd + 0.25)
    assert converted_inputs == pytest.approx([expected_midpoint_jd])
    assert converted_inputs[0] != pytest.approx(expected_start_plus_exposure_jd, abs=1e-9)


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


def test_stellar_variability_ensemble_config_defaults_on_and_supports_opt_out():
    assert should_use_ensemble_photometry_for_stellar_variability(None) is True
    assert should_use_ensemble_photometry_for_stellar_variability("y") is True
    assert should_use_ensemble_photometry_for_stellar_variability("n") is False
    assert should_use_ensemble_photometry_for_stellar_variability(False) is False


def test_apparent_and_exact_comparison_config_defaults_and_values():
    assert should_require_apparent_magnitudes(None) is True
    assert should_require_apparent_magnitudes("n") is False
    assert should_use_exactly_the_comps_provided(None) is False
    assert should_use_exactly_the_comps_provided("y") is True


def test_project_comparison_radec_to_reference_pixels_and_reject_out_of_frame():
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [50.0, 50.0]
    wcs.wcs.cdelt = np.array([-0.001, 0.001])
    wcs.wcs.crval = [31.04125, 46.68972]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]

    projected = project_comparison_radec_to_pixels(
        [[31.04125, 46.68972]],
        wcs.to_header(),
        (100, 100),
    )

    assert projected[0] == pytest.approx([49.0, 49.0])
    with pytest.raises(ValueError, match="projects outside the reference image"):
        project_comparison_radec_to_pixels(
            [[32.04125, 46.68972]],
            wcs.to_header(),
            (100, 100),
        )


def test_independent_ensemble_comparison_limits_default_to_five_and_have_no_upper_cap():
    assert parse_maximum_number_of_ensemble_comparisons_for_transit(None) == 5
    assert parse_maximum_number_of_ensemble_comparisons_for_transit(1) == 5
    assert parse_maximum_number_of_ensemble_comparisons_for_transit("250") == 250
    assert parse_maximum_number_of_ensemble_comparisons_for_stellar_variability(None) == 5
    assert parse_maximum_number_of_ensemble_comparisons_for_stellar_variability(2) == 2
    assert parse_maximum_number_of_ensemble_comparisons_for_stellar_variability("125") == 125


def test_limited_ensemble_comparison_keys_uses_configured_maximum():
    ranked_summaries = [
        {'key': f'comp{index}'}
        for index in range(1, 13)
    ]

    assert limited_ensemble_comparison_keys(ranked_summaries, None) == [
        f'comp{index}' for index in range(1, 6)
    ]
    assert limited_ensemble_comparison_keys(ranked_summaries, 12) == [
        f'comp{index}' for index in range(1, 13)
    ]


def test_transit_ensemble_fit_uses_configured_maximum(monkeypatch):
    frame_count = 6
    quality_mask = np.ones(frame_count, dtype=bool)
    ranked_summaries = [
        {
            'key': f'comp{index}',
            'comp_index': index - 1,
            'aggregate_score': index / 1000.0,
            'coverage_rejected': False,
            'suitability_outlier_rejected': False,
            'psf_quality_keep_mask': quality_mask,
        }
        for index in range(1, 13)
    ]
    comparison_calibration = {
        'method': 'aperture',
        'a': 0,
        'an': 0,
        'aper': 5.0,
        'annulus': 12.0,
        'comp_summaries': ranked_summaries,
    }
    aper_data = {
        'target': np.full((frame_count, 1, 1), 1000.0),
        **{
            f'comp{index}': np.full((frame_count, 1, 1), 100.0 + index)
            for index in range(1, 13)
        },
    }
    captured = {}

    def capture_active_keys(comp_flux_map, active_keys, validity_mask_func):
        captured['active_keys'] = list(active_keys)
        raise RuntimeError('captured configured transit ensemble')

    monkeypatch.setattr(
        'exotic.exotic.build_absolute_comp_ensemble_flux',
        capture_active_keys,
    )

    with pytest.raises(RuntimeError, match='captured configured transit ensemble'):
        fit_ranked_comparison_calibration_candidates(
            np.linspace(0.0, 0.05, frame_count),
            np.linspace(2460000.0, 2460000.05, frame_count),
            np.linspace(1.0, 1.2, frame_count),
            ld=[0.1, 0.1, 0.1, 0.1],
            p_dict={},
            comparison_calibration=comparison_calibration,
            psf_data={},
            aper_data=aper_data,
            target_psf_flux=np.ones(frame_count),
            use_ensemble_photometry_rather_than_single_comp=True,
            maximum_number_of_ensemble_comparisons_for_transit=9,
        )

    assert captured['active_keys'] == [
        f'comp{index}' for index in range(1, 10)
    ]


def test_exact_transit_ensemble_uses_every_supplied_comparison(monkeypatch):
    frame_count = 6
    quality_mask = np.ones(frame_count, dtype=bool)
    comparison_calibration = {
        'method': 'aperture',
        'a': 0,
        'an': 0,
        'aper': 5.0,
        'annulus': 12.0,
        'comp_summaries': [
            {
                'key': f'comp{index}',
                'comp_index': index - 1,
                'aggregate_score': index / 1000.0,
                'coverage_rejected': False,
                'suitability_outlier_rejected': False,
                'psf_quality_keep_mask': quality_mask,
            }
            for index in range(1, 6)
        ],
    }
    aper_data = {
        'target': np.full((frame_count, 1, 1), 1000.0),
        **{
            f'comp{index}': np.full((frame_count, 1, 1), 100.0 + index)
            for index in range(1, 6)
        },
    }
    captured = {}

    def capture_active_keys(comp_flux_map, active_keys, validity_mask_func):
        captured['active_keys'] = list(active_keys)
        raise RuntimeError('captured exact transit ensemble')

    monkeypatch.setattr(
        'exotic.exotic.build_absolute_comp_ensemble_flux',
        capture_active_keys,
    )

    with pytest.raises(RuntimeError, match='captured exact transit ensemble'):
        fit_ranked_comparison_calibration_candidates(
            np.linspace(0.0, 0.05, frame_count),
            np.linspace(2460000.0, 2460000.05, frame_count),
            np.linspace(1.0, 1.2, frame_count),
            ld=[0.1, 0.1, 0.1, 0.1],
            p_dict={},
            comparison_calibration=comparison_calibration,
            psf_data={},
            aper_data=aper_data,
            target_psf_flux=np.ones(frame_count),
            use_ensemble_photometry_rather_than_single_comp=True,
            maximum_number_of_ensemble_comparisons_for_transit=1,
            use_exactly_the_comps_provided=True,
        )

    assert captured['active_keys'] == [f'comp{index}' for index in range(1, 6)]


def test_exact_comparison_mode_does_not_expand_tracked_pool():
    supplied = [[10.0, 20.0], [30.0, 40.0]]
    automatic = [[50.0, 60.0], [70.0, 80.0]]

    tracked, messages = build_tracked_comparison_pool(
        supplied,
        automatic,
        use_exactly_the_comps_provided=True,
    )

    assert tracked == supplied
    assert messages == []


def test_fortuitous_variable_photometry_config_defaults_on_and_supports_opt_out():
    assert should_photometer_fortuitous_variables(None) is True
    assert should_photometer_fortuitous_variables("y") is True
    assert should_photometer_fortuitous_variables("n") is False
    assert should_photometer_fortuitous_variables(False) is False


def test_should_detect_bad_pixels_before_photometry_parses_values():
    assert should_detect_bad_pixels_before_photometry(None) is False
    assert should_detect_bad_pixels_before_photometry("y") is True
    assert should_detect_bad_pixels_before_photometry("n") is False


def test_get_multiprocess_bad_pixel_precheck_processes_parses_values():
    assert get_multiprocess_bad_pixel_precheck_processes(None) is None
    assert get_multiprocess_bad_pixel_precheck_processes("n") is None
    assert get_multiprocess_bad_pixel_precheck_processes("0") is None
    assert get_multiprocess_bad_pixel_precheck_processes("y") >= 1
    assert get_multiprocess_bad_pixel_precheck_processes("3") == 3
    assert get_multiprocess_bad_pixel_precheck_processes(2) == 2


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


def test_should_use_aperture_corrections_and_full_image_fwhm_parses_values():
    assert should_use_aperture_corrections_and_full_image_fwhm(None) is False
    assert should_use_aperture_corrections_and_full_image_fwhm("y") is True
    assert should_use_aperture_corrections_and_full_image_fwhm("n") is False
    assert should_use_aperture_corrections_and_full_image_fwhm(True) is True


def test_overexposure_rejection_config_parsers_default_and_override():
    assert should_reject_overexposed_stars(None) is True
    assert should_reject_overexposed_stars("y") is True
    assert should_reject_overexposed_stars("n") is False
    assert should_reject_overexposed_stars(False) is False

    assert parse_saturation_value(None) == pytest.approx(65535.0)
    assert parse_saturation_value("") == pytest.approx(65535.0)
    assert parse_saturation_value("42000") == pytest.approx(42000.0)
    assert parse_saturation_value(-1) == pytest.approx(65535.0)
    assert parse_saturation_value("not-a-number") == pytest.approx(65535.0)

    assert parse_overexposure_threshold_fraction(None) == pytest.approx(0.9)
    assert parse_overexposure_threshold_fraction("0.75") == pytest.approx(0.75)
    assert parse_overexposure_threshold_fraction(1.0) == pytest.approx(1.0)
    assert parse_overexposure_threshold_fraction(0) == pytest.approx(0.9)
    assert parse_overexposure_threshold_fraction(1.5) == pytest.approx(0.9)


def test_saturation_value_from_header_uses_cecilia_microobservatory_value():
    assert saturation_value_from_header({"TELESCOP": "Cecilia "}) == pytest.approx(4096.0)
    assert saturation_value_from_header({
        "TELESCOP": "Cecilia ",
        "SATURATE": 65535.0,
    }) == pytest.approx(4096.0)
    assert saturation_value_from_header({"SATURATE": 76500.0}) == pytest.approx(76500.0)
    assert saturation_value_from_header({}) is None


def test_should_use_eebls_to_initialize_tmid_and_bounds_parses_values():
    assert should_use_eebls_to_initialize_tmid_and_bounds(None) is True
    assert should_use_eebls_to_initialize_tmid_and_bounds("y") is True
    assert should_use_eebls_to_initialize_tmid_and_bounds("n") is False
    assert should_use_eebls_to_initialize_tmid_and_bounds(True) is True


def test_should_pick_comparison_by_eebls_snr_parses_values():
    assert should_pick_comparison_by_eebls_snr(None) is True
    assert should_pick_comparison_by_eebls_snr("y") is True
    assert should_pick_comparison_by_eebls_snr("n") is False
    assert should_pick_comparison_by_eebls_snr(True) is True


def test_should_use_deviation_from_expected_transit_in_qc_parses_values():
    assert should_use_deviation_from_expected_transit_in_qc(None) is True
    assert should_use_deviation_from_expected_transit_in_qc("y") is True
    assert should_use_deviation_from_expected_transit_in_qc("n") is False
    assert should_use_deviation_from_expected_transit_in_qc(True) is True


def test_parse_deviation_from_expected_transit_in_qc_sigma_parses_values():
    assert parse_deviation_from_expected_transit_in_qc_sigma(None) == pytest.approx(5.0)
    assert parse_deviation_from_expected_transit_in_qc_sigma("7.5") == pytest.approx(7.5)
    assert parse_deviation_from_expected_transit_in_qc_sigma(3) == pytest.approx(3.0)
    assert parse_deviation_from_expected_transit_in_qc_sigma(-1) == pytest.approx(5.0)


def test_should_exit_at_first_qc_pass_solution_parses_values():
    assert should_exit_at_first_qc_pass_solution(None) is True
    assert should_exit_at_first_qc_pass_solution("y") is True
    assert should_exit_at_first_qc_pass_solution("n") is False
    assert should_exit_at_first_qc_pass_solution(True) is True


def test_validate_ultranest_mpi_runtime_rejects_whole_program_mpi(monkeypatch):
    import exotic.exotic as exotic_module

    monkeypatch.setattr(
        exotic_module,
        "get_mpi_status",
        lambda: {"available": True, "size": 72, "rank": 0, "source": "mpi4py", "error": None},
    )

    with pytest.raises(RuntimeError, match="duplicates the full reduction"):
        exotic_module.validate_ultranest_mpi_runtime()


def test_configure_windows_multiprocessing_main_spec_retargets_console_launcher(monkeypatch):
    import exotic.exotic as exotic_module

    fake_main = types.SimpleNamespace(
        __spec__=types.SimpleNamespace(name="exotic"),
        __file__=r"C:\Python312\Scripts\exotic.exe",
        __package__="",
    )
    spawn_executables = []

    monkeypatch.setattr(exotic_module.sys, "platform", "win32")
    monkeypatch.setattr(exotic_module.sys, "_base_executable", r"C:\Python312\python.exe", raising=False)
    monkeypatch.setattr(exotic_module.sys, "executable", r"C:\Python312\Scripts\exotic.exe")
    monkeypatch.setattr(exotic_module.sys, "frozen", True, raising=False)
    monkeypatch.setattr(exotic_module.multiprocessing, "set_executable", spawn_executables.append)
    monkeypatch.setitem(sys.modules, "__main__", fake_main)

    assert configure_windows_multiprocessing_main_spec() is True
    assert fake_main.__spec__ is None
    assert fake_main.__file__ is None
    assert fake_main.__package__ is None
    assert spawn_executables == [r"C:\Python312\python.exe"]
    assert exotic_module.sys.frozen is False


def test_configure_windows_multiprocessing_main_spec_skips_non_windows(monkeypatch):
    import exotic.exotic as exotic_module

    fake_main = types.SimpleNamespace(__spec__=types.SimpleNamespace(name="exotic"), __file__="exotic.exe")

    monkeypatch.setattr(exotic_module.sys, "platform", "linux")
    monkeypatch.setitem(sys.modules, "__main__", fake_main)

    assert configure_windows_multiprocessing_main_spec() is False
    assert fake_main.__spec__.name == "exotic"


def test_configure_windows_multiprocessing_main_spec_preserves_regular_script(monkeypatch):
    import exotic.exotic as exotic_module

    fake_spec = types.SimpleNamespace(name="run_exotic")
    fake_main = types.SimpleNamespace(
        __spec__=fake_spec,
        __file__=r"C:\work\run_exotic.py",
        __package__="",
    )

    monkeypatch.setattr(exotic_module.sys, "platform", "win32")
    monkeypatch.setattr(exotic_module.sys, "_base_executable", r"C:\Python312\python.exe", raising=False)
    monkeypatch.setattr(exotic_module.sys, "executable", r"C:\Python312\python.exe")
    monkeypatch.setattr(exotic_module.multiprocessing, "set_executable", lambda _path: None)
    monkeypatch.setitem(sys.modules, "__main__", fake_main)

    assert configure_windows_multiprocessing_main_spec() is True
    assert fake_main.__spec__ is fake_spec
    assert fake_main.__file__ == r"C:\work\run_exotic.py"


def test_windows_python_spawn_executable_falls_back_to_exec_prefix(monkeypatch):
    import exotic.exotic as exotic_module

    monkeypatch.setattr(exotic_module.sys, "_base_executable", r"C:\Python312\Scripts\exotic.exe", raising=False)
    monkeypatch.setattr(exotic_module.sys, "executable", r"C:\Python312\Scripts\exotic.exe")
    monkeypatch.setattr(exotic_module.sys, "exec_prefix", r"C:\Python312")
    monkeypatch.setattr(exotic_module.sys, "base_exec_prefix", r"C:\Python312", raising=False)

    assert exotic_module._windows_python_spawn_executable() == r"C:\Python312\python.exe"


def test_process_pool_executor_uses_threads_on_windows(monkeypatch):
    import exotic.exotic as exotic_module

    captured = {}

    class FakeThreadPoolExecutor:
        def __init__(self, *args, **kwargs):
            captured["args"] = args
            captured["kwargs"] = kwargs

    monkeypatch.setattr(exotic_module.sys, "platform", "win32")
    monkeypatch.setattr(exotic_module, "ThreadPoolExecutor", FakeThreadPoolExecutor)

    executor = exotic_module.ProcessPoolExecutor(max_workers=3, initializer=lambda: None)

    assert isinstance(executor, FakeThreadPoolExecutor)
    assert captured["kwargs"]["max_workers"] == 3
    assert "initializer" in captured["kwargs"]


def test_build_time_rejection_diagnostic_groups_contiguous_ranges():
    times = np.array([1.0, 1.1, 1.2, 1.5, 1.6, 2.0], dtype=float)
    keep_mask = np.array([True, False, False, True, False, True], dtype=bool)

    diagnostic = build_time_rejection_diagnostic("Example stage", times, keep_mask)

    assert diagnostic["stage"] == "Example stage"
    assert diagnostic["input_point_count"] == 6
    assert diagnostic["kept_point_count"] == 3
    assert diagnostic["dropped_point_count"] == 3
    assert diagnostic["dropped_ranges"] == [
        {"start": pytest.approx(1.1), "end": pytest.approx(1.2), "count": 2},
        {"start": pytest.approx(1.6), "end": pytest.approx(1.6), "count": 1},
    ]


def test_estimate_tmid_and_bounds_with_eebls_identifies_box_like_transit():
    times = np.linspace(0.0, 0.2, 240)
    tmid = 0.101
    duration = 0.028
    flux = np.ones(times.shape[0], dtype=float)
    in_transit = np.abs(times - tmid) <= duration / 2.0
    flux[in_transit] -= 0.018
    flux += 0.0015 * (times - np.nanmean(times))
    flux_errors = np.full(times.shape[0], 0.002, dtype=float)
    prior = {
        "tmid": 0.08,
        "per": 1.0,
        "rprs": np.sqrt(0.018),
        "ars": 12.0,
        "inc": 88.5,
        "ecc": 0.0,
        "omega": 0.0,
    }

    summary = estimate_tmid_and_bounds_with_eebls(
        times,
        flux,
        flux_errors,
        prior,
        [0.04, 0.12],
    )

    assert summary["applied"] is True
    assert summary["method"] == "eebls"
    assert summary["tmid"] == pytest.approx(tmid, abs=0.01)
    assert summary["bounds"][0] < summary["tmid"] < summary["bounds"][1]
    assert summary["depth"] > 0
    assert summary["depth_snr"] > 0


def test_estimate_tmid_and_bounds_with_eebls_keeps_depth_snr_for_one_sided_event():
    times = np.linspace(0.0, 0.11, 160)
    tmid = 0.101
    duration = 0.028
    flux = np.ones(times.shape[0], dtype=float)
    in_transit = np.abs(times - tmid) <= duration / 2.0
    flux[in_transit] -= 0.018
    flux_errors = np.full(times.shape[0], 0.002, dtype=float)
    prior = {
        "tmid": 0.08,
        "per": 1.0,
        "rprs": np.sqrt(0.018),
        "ars": 12.0,
        "inc": 88.5,
        "ecc": 0.0,
        "omega": 0.0,
    }

    summary = estimate_tmid_and_bounds_with_eebls(
        times,
        flux,
        flux_errors,
        prior,
        [0.04, 0.12],
    )

    assert summary["method"] == "eebls"
    assert summary["applied"] is False
    assert summary["depth"] > 0
    assert summary["depth_snr"] > 0
    assert "keeping the EEBLS depth SNR only" in summary["note"]


def test_estimate_ephemeris_tmid_and_bounds_caps_bracketed_runs_to_duration_scale():
    prior = {
        "tmid": 2458247.90746,
        "per": 1.645321,
        "rprs": 0.1265,
        "ars": 4.9,
        "inc": 85.87,
        "ecc": 0.0,
        "omega": 0.0,
    }
    times = np.linspace(2461151.8012, 2461151.9966, 220)
    expected_duration = 0.049

    summary = estimate_ephemeris_tmid_and_bounds(
        times,
        prior["tmid"],
        prior["per"],
        midt_unc=0.00036,
        per_unc=1.0e-5,
        expected_duration=expected_duration,
        sigma_multiplier=35.0,
    )

    assert summary["duration_capped"] is True
    assert summary["observed_window_capped"] is True
    assert summary["observations_bracket_expected_transit"] is True
    assert summary["tmid"] == pytest.approx(2461151.899025, abs=1e-6)
    assert summary["half_width"] > 0
    assert summary["bounds"][0] == pytest.approx(times.min() + 0.5 * expected_duration)
    assert summary["bounds"][1] == pytest.approx(times.max() - 0.5 * expected_duration)


def test_estimate_ephemeris_tmid_and_bounds_keeps_wider_bounds_for_one_sided_runs():
    prior = {
        "tmid": 2458247.90746,
        "per": 1.645321,
    }
    times = np.linspace(2461151.92, 2461152.02, 120)

    summary = estimate_ephemeris_tmid_and_bounds(
        times,
        prior["tmid"],
        prior["per"],
        midt_unc=0.00036,
        per_unc=1.0e-5,
        expected_duration=0.049,
        sigma_multiplier=35.0,
    )

    assert summary["duration_capped"] is False
    assert summary["observations_bracket_expected_transit"] is False
    assert summary["half_width"] == pytest.approx(0.25 * prior["per"])


def test_estimate_ephemeris_tmid_and_bounds_selects_nearest_epoch_for_ingress_only_runs():
    # Regression for issue #1387: an ingress-only partial transit whose true mid
    # falls minutes AFTER the last surviving frame. floor(phases).max() snapped to
    # the previous cycle and reported Tmid one full period early; the epoch nearest
    # the data is the correct one. Geometry taken from the 2026-07-27 TOI-1516 b
    # MicroObservatory night that surfaced the bug (two independent reductions
    # reported 2461246.93, one period before the actual night of the frames).
    prior_tmid = 2458765.325
    period = 2.056014
    times = np.linspace(2461248.8938, 2461248.9847, 34)
    expected_mid = prior_tmid + 1208 * period  # 2461248.9899, ~7 min after times.max()

    summary = estimate_ephemeris_tmid_and_bounds(
        times,
        prior_tmid,
        period,
        midt_unc=0.00023,
        per_unc=2.1e-6,
        expected_duration=0.1177,
        sigma_multiplier=25.0,
    )

    assert summary["cycle_index"] == pytest.approx(1208.0)
    assert summary["tmid"] == pytest.approx(expected_mid, abs=1e-6)
    # The search bounds must be able to reach the true mid.
    assert summary["bounds"][0] <= expected_mid <= summary["bounds"][1]
    assert summary["observations_bracket_expected_transit"] is False


def test_is_adaptive_aperture_mode_enabled_parses_values():
    assert is_adaptive_aperture_mode_enabled(None) is False
    assert is_adaptive_aperture_mode_enabled("y") is True
    assert is_adaptive_aperture_mode_enabled("n") is False
    assert is_adaptive_aperture_mode_enabled(True) is True


def test_aperture_sigma_bounds_match_physical_fwhm_limits():
    assert APERTURE_SIGMA_MIN == pytest.approx(
        APERTURE_MIN_FWHM_MULTIPLIER * GAUSSIAN_SIGMA_TO_FWHM
    )
    assert APERTURE_SIGMA_MAX == pytest.approx(
        APERTURE_MAX_FWHM_MULTIPLIER * GAUSSIAN_SIGMA_TO_FWHM
    )


def test_aperture_correction_profile_recovers_gaussian_curve_of_growth():
    sigma = 2.0
    fwhm = GAUSSIAN_SIGMA_TO_FWHM * sigma
    y, x = np.mgrid[0:120, 0:120]
    image = np.full((120, 120), 10.0, dtype=float)
    positions = np.array([
        [25.0, 25.0],
        [25.0, 70.0],
        [70.0, 25.0],
        [70.0, 70.0],
        [95.0, 95.0],
    ])
    for xc, yc in positions:
        image += 1200.0 * np.exp(-((x - xc) ** 2 + (y - yc) ** 2) / (2.0 * sigma ** 2))

    field_star_psfs = np.column_stack([
        positions[:, 0],
        positions[:, 1],
        np.full(positions.shape[0], 1200.0),
        np.full(positions.shape[0], sigma),
        np.full(positions.shape[0], sigma),
        np.zeros(positions.shape[0]),
        np.full(positions.shape[0], 10.0),
    ])
    radii = np.array([
        APERTURE_MIN_FWHM_MULTIPLIER * fwhm,
        fwhm,
        APERTURE_MAX_FWHM_MULTIPLIER * fwhm,
    ])

    profile = build_aperture_correction_profile(
        image,
        radii,
        fwhm_hint=fwhm,
        field_star_psfs=field_star_psfs,
    )

    assert profile["applied"] is True
    assert profile["star_count"] == positions.shape[0]
    assert profile["image_fwhm"] == pytest.approx(fwhm)
    assert profile["correction_factors"][0] == pytest.approx(2.0, rel=0.15)
    assert profile["correction_factors"][1] == pytest.approx(1.066, rel=0.08)
    assert profile["correction_factors"][2] == pytest.approx(1.0, abs=0.02)


def test_detect_aperture_correction_star_candidates_finds_numpy_local_peaks():
    sigma = 1.8
    fwhm = GAUSSIAN_SIGMA_TO_FWHM * sigma
    y, x = np.mgrid[0:140, 0:140]
    image = np.full((140, 140), 10.0, dtype=float)
    positions = np.array([
        [30.0, 35.0],
        [95.0, 42.0],
        [58.0, 108.0],
    ])
    amplitudes = np.array([1000.0, 850.0, 700.0])
    for (xc, yc), amplitude in zip(positions, amplitudes):
        image += amplitude * np.exp(-((x - xc) ** 2 + (y - yc) ** 2) / (2.0 * sigma ** 2))

    candidates = detect_aperture_correction_star_candidates(image, fwhm_hint=fwhm)

    assert candidates.shape[0] >= positions.shape[0]
    for xc, yc in positions:
        nearest = np.min(np.hypot(candidates[:, 0] - xc, candidates[:, 1] - yc))
        assert nearest < 1.5


def test_compute_star_aperture_grid_applies_aperture_correction_factors():
    y, x = np.mgrid[0:41, 0:41]
    image = 100.0 * np.exp(-((x - 20.0) ** 2 + (y - 20.0) ** 2) / (2.0 * 2.0 ** 2))
    apertures = np.array([2.5, 4.0])
    annuli = np.array([0.0])

    raw_flux, _ = compute_star_aperture_grid(
        image,
        0,
        20.0,
        20.0,
        apertures,
        annuli,
    )
    corrected_flux, _ = compute_star_aperture_grid(
        image,
        0,
        20.0,
        20.0,
        apertures,
        annuli,
        aperture_correction_factors=np.array([2.0, 1.25]),
    )

    np.testing.assert_allclose(corrected_flux[:, 0], raw_flux[:, 0] * np.array([2.0, 1.25]))


def test_compute_star_aperture_grid_ignores_invalid_aperture_geometry():
    image = np.ones((20, 20), dtype=float)

    flux, bg = compute_star_aperture_grid(
        image,
        0,
        10.0,
        10.0,
        np.array([np.nan]),
        np.array([0.0]),
    )

    assert flux.shape == (1, 1)
    assert bg.shape == (1, 1)
    assert np.isnan(flux[0, 0])
    assert np.isnan(bg[0, 0])


def test_aperture_contains_overexposed_pixel_checks_aperture_only(monkeypatch):
    import exotic.exotic as exotic_module

    class FakeMask:
        def __init__(self, xc, yc, radius):
            self.x0 = int(np.floor(xc - radius))
            self.x1 = int(np.ceil(xc + radius)) + 1
            self.y0 = int(np.floor(yc - radius))
            self.y1 = int(np.ceil(yc + radius)) + 1
            y, x = np.mgrid[self.y0:self.y1, self.x0:self.x1]
            self.data = (((x - xc) ** 2 + (y - yc) ** 2) <= radius ** 2).astype(float)

        def cutout(self, data):
            return np.asarray(data)[self.y0:self.y1, self.x0:self.x1]

    class FakeCircularAperture:
        def __init__(self, positions, r):
            self.xc, self.yc = positions[0]
            self.r = r

        def to_mask(self, method="exact"):
            return [FakeMask(self.xc, self.yc, self.r)]

    monkeypatch.setattr(exotic_module, "CircularAperture", FakeCircularAperture)

    data = np.zeros((20, 20), dtype=float)
    data[10, 10] = 90.0
    data[2, 2] = 100.0

    assert aperture_contains_overexposed_pixel(data, 10.0, 10.0, 2.5, 80.0) is True
    assert aperture_contains_overexposed_pixel(data, 10.0, 10.0, 2.5, 95.0) is False
    assert aperture_contains_overexposed_pixel(data, 10.0, 10.0, 2.5, 90.0) is False


def test_populate_aperture_data_skips_field_star_corrections_when_disabled(monkeypatch):
    import exotic.exotic as exotic_module

    def fail_field_star_estimate(*_args, **_kwargs):
        raise AssertionError("field-star FWHM estimation should be opt-in")

    monkeypatch.setattr(exotic_module, "estimate_isolated_field_star_psfs", fail_field_star_estimate)
    y, x = np.mgrid[0:41, 0:41]
    image = 100.0 * np.exp(-((x - 20.0) ** 2 + (y - 20.0) ** 2) / (2.0 * 2.0 ** 2))
    psf_data = {
        "target": np.array([[20.0, 20.0, 100.0, 2.0, 2.0, 0.0, 0.0]]),
    }
    aper_data = initialize_aperture_data_store(
        frame_count=1,
        aperture_count=1,
        annulus_count=1,
        comp_star_count=0,
    )

    profile = populate_aperture_data_for_frame(
        image,
        0,
        psf_data,
        0,
        aper_data,
        np.array([4.0]),
        np.array([0.0]),
        fast_aperture_mask=False,
        use_aperture_corrections_and_full_image_fwhm=False,
    )

    assert profile["applied"] is False
    assert np.isfinite(aper_data["target"][0, 0, 0])


def test_stellar_variability_aperture_estimation_uses_first_five_vetted_comparisons():
    science_comp_stars = [[float(index), float(index + 100)] for index in range(8)]

    assert aperture_estimation_comparison_stars(
        science_comp_stars,
        stellar_variability_only=False,
    ) == science_comp_stars
    assert aperture_estimation_comparison_stars(
        science_comp_stars,
        stellar_variability_only=True,
    ) == science_comp_stars[:5]


def test_stellar_variability_aperture_grid_excludes_variable_target_and_uses_comp_seeing(monkeypatch):
    import exotic.exotic as exotic_module

    measured = []

    def fake_compute_star_aperture_grid(
        _data,
        star_index,
        _xc,
        _yc,
        apertures,
        annuli,
        **_kwargs,
    ):
        aperture_values = np.asarray(apertures, dtype=float).reshape(-1)
        annulus_values = np.asarray(annuli, dtype=float).reshape(-1)
        measured.append((star_index, aperture_values.copy(), annulus_values.copy()))
        shape = (len(aperture_values), len(annulus_values))
        flux = np.full(shape, 100.0 + star_index, dtype=float)
        background = np.full(shape, 10.0 + star_index, dtype=float)
        noise = {
            component: np.ones(shape, dtype=float)
            for component in exotic_module.NOISE_BUDGET_COMPONENT_KEYS
        }
        return flux, background, noise

    monkeypatch.setattr(exotic_module, "compute_star_aperture_grid", fake_compute_star_aperture_grid)
    psf_data = {
        # The VSX science target deliberately has very different seeing. It must not
        # determine the stellar-variability aperture grid.
        "target": np.array([[10.0, 10.0, 100.0, 9.0, 9.0, 0.0, 0.0]]),
        "comp1": np.array([[12.0, 10.0, 90.0, 2.0, 2.0, 0.0, 0.0]]),
        "comp2": np.array([[14.0, 10.0, 80.0, 4.0, 4.0, 0.0, 0.0]]),
    }
    assert aperture_frame_sigma_from_psf_data(
        psf_data,
        0,
        comparison_indices=[0, 1],
    ) == pytest.approx(3.0)

    full_grid = initialize_aperture_data_store(1, 1, 1, 2)
    populate_aperture_data_for_frame(
        np.zeros((25, 25), dtype=float),
        0,
        psf_data,
        2,
        full_grid,
        np.array([2.0]),
        np.array([8.0]),
        fast_aperture_mask=False,
        adaptive_apertures=True,
        comp_indices=[0, 1],
        include_target=False,
        frame_sigma_comp_indices=[0, 1],
    )

    assert [star_index for star_index, _apers, _annuli in measured] == [1, 2]
    assert all(apers[0] == pytest.approx(6.0) for _index, apers, _annuli in measured)
    assert np.all(np.isnan(full_grid["target"]))

    frozen = collapse_aperture_data_to_selected_grid_cell(full_grid, 0, 0)
    measured.clear()
    populate_aperture_data_for_frame(
        np.zeros((25, 25), dtype=float),
        0,
        psf_data,
        2,
        frozen,
        np.array([2.0]),
        np.array([8.0]),
        fast_aperture_mask=False,
        adaptive_apertures=True,
        comp_indices=[],
        include_target=True,
        frame_sigma_comp_indices=[0, 1],
    )

    assert len(measured) == 1
    assert measured[0][0] == 0
    assert measured[0][1][0] == pytest.approx(6.0)
    assert frozen["target"][0, 0, 0] == pytest.approx(100.0)


def test_frozen_aperture_path_only_grids_estimators_then_backfills_additional_stars(monkeypatch):
    import exotic.exotic as exotic_module

    measured_star_indices = []

    def fake_compute_star_aperture_grid(
        _data,
        star_index,
        _xc,
        _yc,
        apertures,
        annuli,
        **_kwargs,
    ):
        measured_star_indices.append(star_index)
        shape = (len(np.asarray(apertures).reshape(-1)), len(np.asarray(annuli).reshape(-1)))
        flux = np.full(shape, 100.0 + star_index, dtype=float)
        background = np.full(shape, 10.0 + star_index, dtype=float)
        noise = {
            component: np.full(shape, 1.0 + star_index, dtype=float)
            for component in exotic_module.NOISE_BUDGET_COMPONENT_KEYS
        }
        return flux, background, noise

    monkeypatch.setattr(exotic_module, "compute_star_aperture_grid", fake_compute_star_aperture_grid)
    psf_data = {
        "target": np.array([[10.0, 10.0, 100.0, 2.0, 2.0, 0.0, 0.0]]),
        "comp1": np.array([[12.0, 10.0, 90.0, 2.0, 2.0, 0.0, 0.0]]),
        "comp2": np.array([[14.0, 10.0, 80.0, 2.0, 2.0, 0.0, 0.0]]),
        "comp3": np.array([[16.0, 10.0, 70.0, 2.0, 2.0, 0.0, 0.0]]),
        "comp4": np.array([[18.0, 10.0, 60.0, 2.0, 2.0, 0.0, 0.0]]),
    }
    full_grid = initialize_aperture_data_store(1, 2, 2, 4)
    populate_aperture_data_for_frame(
        np.zeros((25, 25), dtype=float),
        0,
        psf_data,
        4,
        full_grid,
        np.array([3.0, 4.0]),
        np.array([8.0, 10.0]),
        fast_aperture_mask=False,
        comp_indices=[0, 1],
    )

    assert measured_star_indices == [0, 1, 2]
    assert np.all(np.isfinite(full_grid["target"]))
    assert np.all(np.isfinite(full_grid["comp1"]))
    assert np.all(np.isfinite(full_grid["comp2"]))
    assert np.all(np.isnan(full_grid["comp3"]))
    assert np.all(np.isnan(full_grid["comp4"]))

    frozen = collapse_aperture_data_to_selected_grid_cell(full_grid, 1, 0)
    measured_star_indices.clear()
    populate_aperture_data_for_frame(
        np.zeros((25, 25), dtype=float),
        0,
        psf_data,
        4,
        frozen,
        np.array([4.0]),
        np.array([8.0]),
        fast_aperture_mask=False,
        comp_indices=[2, 3],
        include_target=False,
    )

    assert measured_star_indices == [3, 4]
    assert frozen["target"].shape == (1, 1, 1)
    assert frozen["target"][0, 0, 0] == pytest.approx(100.0)
    assert frozen["comp1"][0, 0, 0] == pytest.approx(101.0)
    assert frozen["comp2"][0, 0, 0] == pytest.approx(102.0)
    assert frozen["comp3"][0, 0, 0] == pytest.approx(103.0)
    assert frozen["comp4"][0, 0, 0] == pytest.approx(104.0)


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


def test_resolve_frame_aperture_radii_accepts_scalar_values():
    apertures, annuli = resolve_frame_aperture_radii(
        2.0,
        8.0,
        adaptive_apertures=True,
        frame_sigma=1.5,
        fallback_sigma=1.0,
    )

    assert apertures.shape == (1,)
    assert annuli.shape == (1,)
    assert apertures[0] == pytest.approx(3.0)
    assert annuli[0] == pytest.approx(12.0)


def test_resolve_sky_annulus_geometry_enforces_fwhm_floor_and_min_sky_pixels():
    geometry = resolve_sky_annulus_geometry(aperture_radius=1.5, annulus_width=2.0, psf_sigma=1.0)

    assert geometry["inner_radius"] == pytest.approx(3.0 * 2.355)
    assert geometry["effective_sky_pixels"] == pytest.approx(250.0, abs=1e-9)
    assert geometry["annulus_width"] > 2.0


def test_choose_centroid_seed_position_prefers_previous_fit_for_small_predicted_jumps():
    seed = choose_centroid_seed_position([100.8, 200.2], previous_psf_row=np.array([100.2, 199.9, 1, 1, 1, 0, 0]))
    np.testing.assert_allclose(seed, np.array([100.2, 199.9]))


def test_choose_centroid_seed_position_falls_back_to_predicted_for_large_jump_or_invalid_previous():
    seed_far = choose_centroid_seed_position([110.0, 210.0], previous_psf_row=np.array([100.0, 200.0, 1, 1, 1, 0, 0]))
    seed_nan = choose_centroid_seed_position([110.0, 210.0], previous_psf_row=np.array([np.nan, 200.0, 1, 1, 1, 0, 0]))

    np.testing.assert_allclose(seed_far, np.array([110.0, 210.0]))
    np.testing.assert_allclose(seed_nan, np.array([110.0, 210.0]))


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
            "eebls_snr": np.nan,
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
            "eebls_snr": 6.5,
            "transit_delta_bic": 18.4,
            "residual_scatter": 0.0042,
            "ktmf_metric": 4.35,
            "ktmf_contributions": [
                {
                    "label": "EEBLS Depth SNR",
                    "available": True,
                    "points": 1.25,
                    "max_points": 1.40,
                    "score": 0.89,
                    "detail": "6.50",
                }
            ],
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
        {
            "selection_basis": "comparison_field",
            "selection_metric": "ktmf",
            "comp_star_num": 2,
            "comparison_ktmf_metric": 4.35,
            "comparison_eebls_snr": 6.5,
            "comparison_transit_delta_bic": 18.4,
        },
    )

    assert any("Selection basis: comparison-field" in message for message in logged)
    assert any("Selection metric: KTMF" in message for message in logged)
    assert any("coverage=1 valid frame(s) out of 3 total; min_required=2; peer_median=3.0" in message for message in logged)
    assert any("Comp 1" in message and "reason=comparison candidate rejected after iterative low-coverage clipping" in message for message in logged)
    assert any("Comp 2 [selected]" in message and "ktmf=4.35/5.00" in message and "comparison-field calibration ranked this star best" in message for message in logged)
    assert any("KTMF contribution: EEBLS Depth SNR +1.25/1.40" in message for message in logged)
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
            "eebls_snr": np.nan,
            "fit_point_count": 0,
            "fit_diagnostics": {"usable_point_count": 0},
            "failure_reason": "relative-flux filtering left 0 usable point(s); rejected 3/3 frame(s) during invalid target/reference ratio screening (non-finite=0, non-positive=3, finite ratio range=-1.0000 to -1.0000).",
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
            "eebls_snr": 5.2,
            "transit_delta_bic": 18.4,
            "residual_scatter": 0.0035,
            "ktmf_metric": 4.60,
            "ktmf_contributions": [
                {
                    "label": "Residual Scatter Around Full Model Fit",
                    "available": True,
                    "points": 0.63,
                    "max_points": 0.80,
                    "score": 0.79,
                    "detail": "0.3500%",
                }
            ],
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
    assert any("Comp 2 [selected]" in message and "ktmf=4.60/5.00" in message and "fit_points=3" in message for message in logged)
    assert any("KTMF contribution: Residual Scatter Around Full Model Fit +0.63/0.80" in message for message in logged)
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
            "residual_scatter": np.inf,
            "eebls_snr": np.nan,
            "ktmf_metric": 0.0,
            "ktmf_contributions": [
                {
                    "label": "Deviation From Expected Value",
                    "available": False,
                    "points": 0.0,
                    "max_points": 0.0,
                    "score": np.nan,
                    "detail": "expected-value deviation disabled or unavailable",
                }
            ],
            "coverage_count": 3,
            "coverage_total_frame_count": 3,
            "coverage_reference_count": 3.0,
            "coverage_min_required_count": 2,
            "fit_point_count": 0,
            "fit_diagnostics": {"usable_point_count": 0},
            "failure_reason": "relative-flux filtering left 0 usable point(s); rejected 3/3 frame(s) during invalid target/reference ratio screening (non-finite=0, non-positive=3, finite ratio range=-1.0000 to -1.0000).",
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
    assert any("KTMF contribution: Deviation From Expected Value +0.00/0.00 (unavailable;" in message for message in logged)


def test_compute_transit_qc_ktmf_uses_rebalanced_component_weights():
    summary = {
        "delta_bic": 10.0,
        "delta_chi2": 50.0,
        "deviation_from_expected_value": 0.6,
        "tmid_deviation_sigma": 1.0,
        "rprs_deviation_sigma": 2.0,
        "residual_scatter": 0.005,
        "transit_depth_for_residual_scatter": 0.01,
        "residual_scatter_to_depth_ratio": 0.5,
        "residual_flatness_score": 0.5,
        "residual_flatness_detail": "curve=0.50",
        "tmid_gaussianity_score": 0.8,
        "tmid_gaussianity_score_uncertainty": 0.04,
        "tmid_gaussianity_detail": "strongly Gaussian-like",
        "rprs_sigma": 6.0,
        "duration_ratio": 1.0,
        "eebls_depth_snr": 8.0,
    }

    ktmf_metric, contributions = compute_transit_qc_ktmf(summary)
    contributions_by_label = {contribution["label"]: contribution for contribution in contributions}

    assert "Model Evidence" not in contributions_by_label
    assert "Delta BIC" not in contributions_by_label
    assert "Delta chi2" not in contributions_by_label
    scale = 5.0 / (2.0 + 0.7 + 1.0 + 1.0 + 0.75 + 1.3)
    assert contributions_by_label["Deviation From Expected Value"]["max_points"] == pytest.approx(2.0 * scale)
    assert contributions_by_label["Residual Scatter Around Full Model Fit"]["max_points"] == pytest.approx(0.7 * scale)
    assert contributions_by_label["Residual Flatness"]["max_points"] == pytest.approx(1.0 * scale)
    assert contributions_by_label["Tmid Posterior Gaussianity"]["max_points"] == pytest.approx(1.0 * scale)
    assert contributions_by_label["Tmid Posterior Gaussianity"]["score_uncertainty"] == pytest.approx(0.04)
    assert "Rp/R* Significance" not in contributions_by_label
    assert contributions_by_label["Duration Consistency"]["max_points"] == pytest.approx(0.75 * scale)
    assert contributions_by_label["EEBLS Depth SNR"]["max_points"] == pytest.approx(1.3 * scale)
    assert "Rp/R* sigma=2.00" in contributions_by_label["Deviation From Expected Value"]["detail"]
    assert "Tmid" not in contributions_by_label["Deviation From Expected Value"]["detail"]
    assert contributions_by_label["Residual Flatness"]["score"] == pytest.approx(0.5)

    expected_ktmf = scale * (
        2.0 * 0.6
        + 0.7 * 1.0
        + 1.0 * 0.5
        + 1.0 * 0.8
        + 0.75 * 1.0
        + 1.3 * (1.0 - np.exp(-2.0))
    )
    assert ktmf_metric == pytest.approx(expected_ktmf)


class _TmidPosteriorFit:
    def __init__(self, values, weights=None, sampled_keys=("tmid",)):
        self.values = np.asarray(values, dtype=float)
        self.weights = None if weights is None else np.asarray(weights, dtype=float)
        self.sampled_keys = list(sampled_keys)
        self.bounds = {"tmid": [float(np.min(self.values)), float(np.max(self.values))]}

    def _get_triangle_plot_samples(self):
        return self.values[:, np.newaxis], np.zeros(self.values.size), self.weights


def test_tmid_posterior_gaussianity_distinguishes_gaussian_flat_skewed_and_multimodal_shapes():
    rng = np.random.default_rng(20260715)
    center = 2460835.82621
    gaussian_values = center + rng.normal(0.0, 0.0015, 5000)
    flat_values = np.linspace(center - 0.006, center + 0.006, 5000)
    skewed_values = center + 0.002 * (rng.lognormal(-1.0, 0.5, 5000) - 0.42)
    multimodal_values = center + np.concatenate([
        rng.normal(-0.003, 0.0005, 2500),
        rng.normal(0.003, 0.0005, 2500),
    ])

    gaussian = transit_qc_tmid_gaussianity_summary(_TmidPosteriorFit(gaussian_values))
    flat = transit_qc_tmid_gaussianity_summary(_TmidPosteriorFit(flat_values))
    skewed = transit_qc_tmid_gaussianity_summary(_TmidPosteriorFit(skewed_values))
    multimodal = transit_qc_tmid_gaussianity_summary(_TmidPosteriorFit(multimodal_values))

    assert gaussian["available"] is True
    assert gaussian["score"] > 0.90
    assert np.isfinite(gaussian["score_uncertainty"])
    assert "strongly Gaussian-like" in gaussian["detail"]
    assert flat["score"] < 0.05
    assert skewed["score"] < 0.60
    assert multimodal["score"] < 0.20


def test_tmid_posterior_gaussianity_uses_ultranest_sample_weights():
    rng = np.random.default_rng(717)
    center = 2460835.82621
    flat_values = np.linspace(center - 0.01, center + 0.01, 6000)
    gaussian_values = center + rng.normal(0.0, 0.001, 2500)
    values = np.concatenate([flat_values, gaussian_values])
    weights = np.concatenate([
        np.full(flat_values.size, 1e-8),
        np.ones(gaussian_values.size),
    ])

    summary = transit_qc_tmid_gaussianity_summary(_TmidPosteriorFit(values, weights=weights))

    assert summary["available"] is True
    assert summary["effective_sample_count"] == pytest.approx(2500.0, rel=1e-4)
    assert summary["score"] > 0.85


def test_tmid_posterior_gaussianity_is_unavailable_when_tmid_was_fixed():
    summary = transit_qc_tmid_gaussianity_summary(
        _TmidPosteriorFit(np.linspace(0.0, 1.0, 500), sampled_keys=())
    )

    assert summary["available"] is False
    assert not np.isfinite(summary["score"])
    assert "fixed rather than sampled" in summary["detail"]


def test_transit_qc_residual_scatter_score_full_credit_floor_and_zero_ceiling():
    transit_depth = 0.02
    assert transit_qc_residual_scatter_score(0.0, transit_depth) == pytest.approx(1.0)
    assert transit_qc_residual_scatter_score(0.01, transit_depth) == pytest.approx(1.0)
    assert transit_qc_residual_scatter_score(0.08, transit_depth) == pytest.approx(0.0)
    assert transit_qc_residual_scatter_score(0.09, transit_depth) == pytest.approx(0.0)
    assert not np.isfinite(transit_qc_residual_scatter_score(0.005))

    mid_score = transit_qc_residual_scatter_score(0.03, transit_depth)
    assert 0.0 < mid_score < 1.0
    assert mid_score < transit_qc_residual_scatter_score(0.02, transit_depth)


def test_evaluate_transit_detection_qc_uses_transit_component_depth_for_residual_scatter():
    transit_component = np.ones(21, dtype=float)
    transit_component[8:13] = 0.984
    baseline_trend = np.linspace(0.0, 0.04, transit_component.size)
    full_model = transit_component + baseline_trend
    data = full_model + np.array(
        [
            0.0002, -0.0001, 0.0001, -0.0002, 0.0000, 0.0001, -0.0001,
            0.0002, -0.0002, 0.0001, -0.0001, 0.0002, -0.0002, 0.0001,
            0.0000, -0.0001, 0.0002, -0.0001, 0.0001, 0.0000, -0.0001,
        ],
        dtype=float,
    )
    fit = types.SimpleNamespace(
        data=data,
        dataerr=np.full(data.shape[0], 0.0015, dtype=float),
        model=full_model,
        transit=transit_component,
        airmass=np.ones(data.shape[0], dtype=float),
        airmass_fit_skipped=True,
        parameters={"rprs": 0.10, "tmid": 0.5, "inc": 89.0, "a2": 0.0},
        errors={"rprs": 0.01, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
        bounds={"rprs": [0.0, 1.0], "tmid": [0.4, 0.6], "inc": [80.0, 90.0]},
        duration_expected=5.0,
        duration_measured=5.0,
    )

    summary = evaluate_transit_detection_qc(fit)

    assert summary["computed"] is True
    assert summary["transit_depth_for_residual_scatter"] == pytest.approx(0.016, abs=5e-4)
    assert summary["residual_scatter_to_depth_ratio"] == pytest.approx(
        summary["residual_scatter"] / summary["transit_depth_for_residual_scatter"]
    )


def test_transit_qc_residual_flatness_summary_penalizes_residual_structure():
    phase = np.linspace(-0.05, 0.05, 80)
    alternating_noise = 0.001 * np.where(np.arange(phase.size) % 2 == 0, -1.0, 1.0)

    flat_summary = transit_qc_residual_flatness_summary(alternating_noise, phase)
    trend_summary = transit_qc_residual_flatness_summary(
        alternating_noise + 0.004 * np.linspace(-1.0, 1.0, phase.size),
        phase,
    )
    curve_summary = transit_qc_residual_flatness_summary(
        alternating_noise + 0.004 * np.sin(2.0 * np.pi * np.linspace(0.0, 1.0, phase.size)),
        phase,
    )
    smooth_bowl_summary = transit_qc_residual_flatness_summary(
        alternating_noise + 0.003 * np.maximum(0.0, 1.0 - (phase / 0.02) ** 2),
        phase,
    )
    heteroscedastic_summary = transit_qc_residual_flatness_summary(
        alternating_noise * np.r_[np.ones(40), np.full(40, 4.0)],
        phase,
    )
    one_sided_summary = transit_qc_residual_flatness_summary(
        alternating_noise - 0.003,
        phase,
    )

    assert flat_summary["available"] is True
    assert flat_summary["score"] > 0.9
    assert trend_summary["score"] < flat_summary["score"]
    assert trend_summary["score"] < 0.5
    assert curve_summary["score"] < flat_summary["score"]
    assert curve_summary["score"] < 0.6
    assert smooth_bowl_summary["score"] < flat_summary["score"]
    assert smooth_bowl_summary["score"] < 0.5
    assert smooth_bowl_summary["dominant"] == "curvature/sinusoid"
    assert heteroscedastic_summary["score"] < flat_summary["score"]
    assert heteroscedastic_summary["score"] < 0.8
    assert one_sided_summary["score"] < flat_summary["score"]
    assert one_sided_summary["score"] < 0.4
    assert one_sided_summary["dominant"] == "zero bias"
    assert one_sided_summary["sign_imbalance"] > 0.8


def test_transit_qc_residual_flatness_summary_tolerates_one_quiet_patch():
    phase = np.linspace(-0.05, 0.05, 84)
    residuals = 0.001 * np.sin(np.arange(phase.size) * 2.3999632)
    residuals[-10:] *= 0.12

    summary = transit_qc_residual_flatness_summary(residuals, phase)

    assert summary["available"] is True
    assert summary["score"] > 0.5
    assert summary["scatter_ratio"] < 3.0


def test_transit_qc_sampling_summary_scores_ingress_egress_and_baseline_counts():
    fit = types.SimpleNamespace(
        time=np.array([
            -0.090, -0.075, -0.060,
            -0.050, -0.045, -0.040, -0.035, -0.030,
            -0.020, -0.010, 0.000, 0.010, 0.020,
            0.030, 0.035, 0.040, 0.045, 0.050,
            0.060, 0.075, 0.090,
        ]),
        parameters={"tmid": 0.0, "rprs": 0.1},
        duration_expected=0.1,
    )

    summary = transit_qc_sampling_summary(fit)

    assert summary["available"] is True
    assert summary["ingress_count"] == 5
    assert summary["egress_count"] == 5
    assert summary["in_transit_count"] == 15
    assert summary["pre_baseline_count"] == 3
    assert summary["post_baseline_count"] == 3
    assert 0.0 < summary["score"] < 1.0
    assert "ingress=5, egress=5" in summary["detail"]


def test_compute_transit_qc_ktmf_omits_prior_assumed_rprs_component():
    summary = {
        "delta_bic": 10.0,
        "delta_chi2": 50.0,
        "deviation_from_expected_value": 1.0,
        "rprs_deviation_sigma": 0.0,
        "rprs_prior_assumed": True,
        "rprs_prior_assumed_note": "Rp/R* was fixed to the input prior.",
        "residual_scatter": 0.005,
        "transit_depth_for_residual_scatter": 0.01,
        "residual_scatter_to_depth_ratio": 0.5,
        "point_count": 86,
        "duration_ratio": 1.0,
        "eebls_depth_snr": 8.0,
    }

    ktmf_metric, contributions = compute_transit_qc_ktmf(summary)
    contributions_by_label = {contribution["label"]: contribution for contribution in contributions}

    omitted = contributions_by_label["Deviation From Expected Value"]
    assert omitted["available"] is False
    assert omitted["points"] == pytest.approx(0.0)
    assert omitted["max_points"] == pytest.approx(0.0)
    assert "fixed to the input prior" in omitted["detail"]

    assert "Model Evidence" not in contributions_by_label
    scale = 5.0 / (0.7 + 0.75 + 1.3)
    assert contributions_by_label["Residual Scatter Around Full Model Fit"]["max_points"] == pytest.approx(0.7 * scale)
    assert contributions_by_label["Duration Consistency"]["max_points"] == pytest.approx(0.75 * scale)
    assert contributions_by_label["EEBLS Depth SNR"]["max_points"] == pytest.approx(1.3 * scale)

    expected_ktmf = scale * (
        0.7 * 1.0
        + 0.75 * 1.0
        + 1.3 * (1.0 - np.exp(-2.0))
    )
    assert ktmf_metric == pytest.approx(expected_ktmf)


def test_compute_transit_qc_ktmf_adds_tmid_gaussianity_for_prior_assumed_geometry():
    summary = {
        "geometry_prior_assumed": True,
        "geometry_prior_assumed_note": "Transit geometry was fixed to priors.",
        "deviation_from_expected_value": 1.0,
        "residual_scatter": 0.005,
        "transit_depth_for_residual_scatter": 0.01,
        "residual_scatter_to_depth_ratio": 0.5,
        "point_count": 86,
        "duration_ratio": 1.0,
        "sampling_score": 1.0,
        "sampling_detail": "ingress=4, egress=4",
        "eebls_depth_snr": 8.0,
        "tmid_gaussianity_score": 0.75,
        "tmid_gaussianity_score_uncertainty": 0.05,
        "tmid_gaussianity_detail": "broadly Gaussian-like",
    }

    ktmf_metric, contributions = compute_transit_qc_ktmf(summary)
    contributions_by_label = {contribution["label"]: contribution for contribution in contributions}

    assert contributions_by_label["Deviation From Expected Value"]["available"] is False
    assert contributions_by_label["Duration Consistency"]["available"] is False
    assert contributions_by_label["Sampling / Cadence"]["available"] is False
    assert "fixed to priors" in contributions_by_label["Duration Consistency"]["detail"]

    scale = 5.0 / (0.7 + 1.0 + 1.3)
    assert contributions_by_label["Residual Scatter Around Full Model Fit"]["max_points"] == pytest.approx(0.7 * scale)
    assert contributions_by_label["Tmid Posterior Gaussianity"]["max_points"] == pytest.approx(1.0 * scale)
    assert contributions_by_label["Tmid Posterior Gaussianity"]["score"] == pytest.approx(0.75)
    assert contributions_by_label["EEBLS Depth SNR"]["max_points"] == pytest.approx(1.3 * scale)

    expected_ktmf = scale * (
        0.7 * 1.0
        + 1.0 * 0.75
        + 1.3 * (1.0 - np.exp(-2.0))
    )
    assert ktmf_metric == pytest.approx(expected_ktmf)


def test_comparison_candidate_fit_selection_reason_describes_comparison_field_retry():
    reason = comparison_candidate_fit_selection_reason(
        {
            "selected": True,
            "failure_reason": None,
            "transit_delta_bic": 18.4,
        },
        {
            "selection_basis": "comparison_field_retry",
            "comp_star_num": 2,
            "comparison_transit_delta_bic": 18.4,
        },
    )

    assert "fell back to this star" in reason


def test_comparison_calibration_selection_reason_reports_suitability_outlier_rejection():
    reason = comparison_calibration_selection_reason(
        {
            "selected": False,
            "coverage_rejected": False,
            "aggregate_score": 0.139668,
            "suitability_outlier_rejected": True,
            "suitability_high_threshold": 0.0398471675,
        },
        best_comp_score=0.021654,
    )

    assert "rejected by high-side sigma clipping" in reason
    assert "13.9668%" in reason


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


def test_exact_comparison_mode_bypasses_star_and_frame_vetting():
    airmass = np.linspace(1.0, 1.5, 8)
    summary = comparison_star_stability_summary(
        {
            "comp1": np.array([100.0, 101.0, 100.0, 101.0, 100.0, 101.0, 100.0, 101.0]),
            "comp2": np.array([80.0, 80.0, 80.0, 160.0, 80.0, 80.0, 80.0, 80.0]),
            "comp3": np.array([60.0, 60.0, 60.0, 60.0, 60.0, 60.0, np.nan, np.nan]),
        },
        airmass,
        bypass_vetting=True,
    )

    assert [row["key"] for row in summary["comp_summaries"]] == [
        "comp1",
        "comp2",
        "comp3",
    ]
    assert not any(row["coverage_rejected"] for row in summary["comp_summaries"])
    assert not any(row["suitability_outlier_rejected"] for row in summary["comp_summaries"])
    assert np.all(summary["field_image_keep_mask"])


def test_apply_comparison_star_suitability_outlier_rejection_rejects_high_tail():
    comp_summaries = [
        {"label": "Comp 1", "aggregate_score": 0.139668, "coverage_rejected": False},
        {"label": "Comp 2", "aggregate_score": 0.051809, "coverage_rejected": False},
        {"label": "Comp 3", "aggregate_score": 0.024802, "coverage_rejected": False},
        {"label": "Comp 4", "aggregate_score": 0.027680, "coverage_rejected": False},
        {"label": "Comp 5", "aggregate_score": 0.036997, "coverage_rejected": False},
        {"label": "Comp 6", "aggregate_score": 0.037970, "coverage_rejected": False},
        {"label": "Comp 7", "aggregate_score": 0.023458, "coverage_rejected": False},
        {"label": "Comp 8", "aggregate_score": 0.021654, "coverage_rejected": False},
        {"label": "Comp 9", "aggregate_score": 0.025084, "coverage_rejected": False},
        {"label": "Comp 10", "aggregate_score": 0.024918, "coverage_rejected": False},
    ]

    result = apply_comparison_star_suitability_outlier_rejection(comp_summaries)

    assert result["rejected_indices"] == [0, 1]
    assert comp_summaries[0]["suitability_outlier_rejected"] is True
    assert comp_summaries[1]["suitability_outlier_rejected"] is True
    assert comp_summaries[4]["suitability_outlier_rejected"] is False
    assert comp_summaries[5]["suitability_outlier_rejected"] is False
    assert 0.037970 < result["high_threshold"] < 0.051809


def test_comparison_star_stability_summary_iterates_after_suitability_outlier_rejection(monkeypatch):
    monkeypatch.setattr(
        "exotic.exotic.normalize_flux_series",
        lambda flux_values, validity_mask_func=None: np.asarray(flux_values, dtype=float),
    )
    monkeypatch.setattr(
        "exotic.exotic.normalized_ratio_series",
        lambda flux_a, flux_b: np.array([1.0], dtype=float),
    )

    def fake_build_normalized_comp_ensemble(normalized_flux_map, exclude_key):
        comp_index = float(exclude_key.replace("comp", ""))
        return np.array([-float(len(normalized_flux_map)), comp_index], dtype=float)

    def fake_prescore(tflux, cflux, airmass, enforce_relative_flux_max=False):
        comp_index = int(np.rint(np.asarray(tflux, dtype=float).flat[0]))
        reference = np.asarray(cflux, dtype=float).reshape(-1)
        if reference.size == 0:
            return np.inf
        if np.allclose(reference, 1.0):
            return 0.0
        if reference[0] < 0:
            active_count = int(np.rint(abs(reference[0])))
            ensemble_scores = {
                6: {1: 10.0, 2: 2.5, 3: 1.0, 4: 1.1, 5: 1.2, 6: 1.4},
                5: {2: 4.0, 3: 1.0, 4: 1.1, 5: 1.2, 6: 1.4},
                4: {3: 1.0, 4: 1.1, 5: 1.2, 6: 1.4},
            }
            return ensemble_scores.get(active_count, {}).get(comp_index, 1.0)
        return 0.5

    monkeypatch.setattr(
        "exotic.exotic.build_normalized_comp_ensemble",
        fake_build_normalized_comp_ensemble,
    )
    monkeypatch.setattr("exotic.exotic.cheap_lightcurve_prescore", fake_prescore)

    summary = comparison_star_stability_summary(
        {
            "comp1": np.array([1.0], dtype=float),
            "comp2": np.array([2.0], dtype=float),
            "comp3": np.array([3.0], dtype=float),
            "comp4": np.array([4.0], dtype=float),
            "comp5": np.array([5.0], dtype=float),
            "comp6": np.array([6.0], dtype=float),
        },
        np.array([1.0], dtype=float),
    )

    assert summary["suitability_outlier_rejected_count"] == 2
    assert summary["comp_summaries"][0]["suitability_outlier_rejected"] is True
    assert summary["comp_summaries"][1]["suitability_outlier_rejected"] is True
    assert summary["comp_summaries"][2]["suitability_outlier_rejected"] is False
    assert summary["best_comp_index"] == 2


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


def test_comparison_star_coverage_summary_keeps_nearly_complete_candidates():
    frame_count = 146
    coverage = comparison_star_coverage_summary(
        {
            **{
                f"comp{comp_index + 1}": np.ones(frame_count, dtype=float)
                for comp_index in range(8)
            },
            "comp9": np.concatenate([np.ones(145, dtype=float), [np.nan]]),
            "comp10": np.concatenate([np.ones(142, dtype=float), np.full(4, np.nan)]),
        }
    )

    assert coverage["comp9"]["coverage_count"] == 145
    assert coverage["comp10"]["coverage_count"] == 142
    assert coverage["comp9"]["coverage_min_required_count"] == 117
    assert coverage["comp10"]["coverage_min_required_count"] == 117
    assert coverage["comp9"]["coverage_rejected"] is False
    assert coverage["comp10"]["coverage_rejected"] is False


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


def test_comparison_star_stability_summary_rejects_noisy_nearly_complete_candidates_as_outliers():
    frame_count = 146
    airmass = np.linspace(1.0, 1.5, frame_count)
    phase = np.linspace(0.0, 4.0 * np.pi, frame_count)
    stable_flux = 100.0 * (1.0 + 0.001 * np.sin(phase))
    comp_flux_map = {
        f"comp{comp_index + 1}": stable_flux * (1.0 + 0.0001 * comp_index)
        for comp_index in range(8)
    }
    noisy_flux = 100.0 * (1.0 + 0.35 * np.sin(np.linspace(0.0, 14.0 * np.pi, frame_count)))
    noisy_flux[-1] = np.nan
    choppy_flux = 100.0 * (1.0 + 0.25 * np.sign(np.sin(np.linspace(0.0, 20.0 * np.pi, frame_count))))
    choppy_flux[-4:] = np.nan
    comp_flux_map["comp9"] = noisy_flux
    comp_flux_map["comp10"] = choppy_flux

    summary = comparison_star_stability_summary(comp_flux_map, airmass)
    comp9_summary = summary["comp_summaries"][8]
    comp10_summary = summary["comp_summaries"][9]

    assert comp9_summary["coverage_count"] == 145
    assert comp10_summary["coverage_count"] == 142
    assert comp9_summary["coverage_rejected"] is False
    assert comp10_summary["coverage_rejected"] is False
    assert comp9_summary["suitability_outlier_rejected"] is True
    assert comp10_summary["suitability_outlier_rejected"] is True
    reason = comparison_calibration_selection_reason(
        comp9_summary,
        summary["best_comp_score"],
    )
    assert "high-side sigma clipping" in reason
    assert "low coverage" not in reason


def test_comparison_star_stability_summary_rejects_shared_bad_frame():
    airmass = np.linspace(1.0, 1.5, 6)
    summary = comparison_star_stability_summary(
        {
            "comp1": np.array([100.0, 100.8, 99.6, 100.4, 100.1, 140.0]),
            "comp2": np.array([80.0, 79.5, 80.6, 80.2, 79.8, 40.0]),
            "comp3": np.array([120.0, 121.0, 119.2, 120.5, 119.7, 100.0]),
        },
        airmass,
    )

    np.testing.assert_array_equal(
        summary["field_image_keep_mask"],
        np.array([True, True, True, True, True, False], dtype=bool),
    )
    assert summary["image_outlier_rejected_count"] == 1
    assert summary["image_outlier_required_valid_pairs"] == 2
    assert summary["image_outlier_available_pairs"] == 3
    assert summary["image_outlier_valid_pair_counts"][-1] == 3
    assert summary["image_outlier_outlier_pair_counts"][-1] == 3


def test_comparison_star_stability_summary_flags_candidate_specific_bad_frame():
    airmass = np.linspace(1.0, 1.5, 12)
    comp1 = np.full(12, 100.0, dtype=float)
    comp2 = np.full(12, 80.0, dtype=float)
    comp3 = np.full(12, 120.0, dtype=float)
    comp4 = np.full(12, 90.0, dtype=float)
    comp1[7] = 60.0

    summary = comparison_star_stability_summary(
        {
            "comp1": comp1,
            "comp2": comp2,
            "comp3": comp3,
            "comp4": comp4,
        },
        airmass,
    )

    assert summary["field_image_keep_mask"].all()
    comp1_summary = summary["comp_summaries"][0]
    comp2_summary = summary["comp_summaries"][1]
    assert comp1_summary["ensemble_frame_rejected_indices"] == [7]
    assert comp1_summary["ensemble_frame_rejected_count"] == 1
    assert comp1_summary["ensemble_frame_valid_pair_counts"][7] == 3
    assert comp1_summary["ensemble_frame_outlier_pair_counts"][7] == 3
    assert comp2_summary["ensemble_frame_rejected_count"] == 0


def test_comparison_star_stability_summary_clips_candidate_psf_spikes_before_suitability_rejection():
    frame_count = 89
    airmass = np.linspace(1.25, 1.06, frame_count)
    phase = np.linspace(0.0, 4.0 * np.pi, frame_count)
    comp_flux_map = {
        f"comp{index + 1}": 100.0 * (1.0 + 0.002 * np.sin(phase + index))
        for index in range(9)
    }
    spike_indices = np.array([2, 5, 11, 12, 13, 39, 43, 45, 56], dtype=int)
    comp_flux_map["comp6"] = comp_flux_map["comp6"].copy()
    comp_flux_map["comp6"][spike_indices] *= 0.35

    summary = comparison_star_stability_summary(comp_flux_map, airmass)
    comp6_summary = summary["comp_summaries"][5]

    assert comp6_summary["coverage_count"] == frame_count
    assert comp6_summary["suitability_outlier_rejected"] is False
    assert comp6_summary["aggregate_score"] < 0.01
    assert comp6_summary["ensemble_frame_rejected_count"] == len(spike_indices)
    assert comp6_summary["ensemble_frame_rejected_indices"] == spike_indices.tolist()


def test_select_comparison_calibrated_photometry_masks_psf_quality_before_aperture_ensemble():
    frame_count = 30
    airmass = np.linspace(1.2, 1.0, frame_count)

    def build_psf_rows():
        rows = np.zeros((frame_count, 7), dtype=float)
        rows[:, 0] = 10.0
        rows[:, 1] = 20.0
        rows[:, 2] = 200.0
        rows[:, 3] = 1.0
        rows[:, 4] = 1.0
        return rows

    psf_data = {
        "target": build_psf_rows(),
        "comp1": build_psf_rows(),
        "comp2": build_psf_rows(),
        "comp3": build_psf_rows(),
    }
    psf_data["comp1"][7, 2] = 40.0
    psf_data["comp1"][7, 3:5] = 6.0

    aper_data = {
        "target": np.full((frame_count, 1, 1), 1000.0),
        "target_bg": np.full((frame_count, 1, 1), 10.0),
    }
    for key in ("comp1", "comp2", "comp3"):
        aper_data[key] = np.full((frame_count, 1, 1), 100.0)
        aper_data[f"{key}_bg"] = np.full((frame_count, 1, 1), 10.0)
    aper_data["comp1"][7, 0, 0] = 1.0

    calibration = select_comparison_calibrated_photometry(
        psf_data,
        aper_data,
        apers=np.array([2.5]),
        annuli=np.array([10.0]),
        airmass=airmass,
        comp_stars=[[10.0, 20.0], [30.0, 40.0], [50.0, 60.0]],
        sigma=1.0,
        use_psf_photometry=False,
        use_aperture_photometry=True,
    )
    comp1_summary = calibration["comp_summaries"][0]

    assert comp1_summary["psf_quality_rejected_count"] == 1
    assert comp1_summary["coverage_count"] == frame_count - 1
    assert comp1_summary["ensemble_frame_rejected_count"] == 0
    assert np.isnan(comp1_summary["ensemble_ratio_series"][7])


def test_exact_comparison_calibration_does_not_reject_an_infinite_stability_score():
    frame_count = 6
    psf_rows = np.ones((frame_count, 7), dtype=float)
    psf_rows[:, 3:5] = 1.0
    calibration = select_comparison_calibrated_photometry(
        {
            "target": psf_rows.copy(),
            "comp1": psf_rows.copy(),
        },
        {
            "target": np.full((frame_count, 1, 1), 1000.0),
            "comp1": np.full((frame_count, 1, 1), np.nan),
        },
        apers=np.array([2.5]),
        annuli=np.array([10.0]),
        airmass=np.linspace(1.0, 1.5, frame_count),
        comp_stars=[[10.0, 20.0]],
        sigma=1.0,
        use_psf_photometry=False,
        use_aperture_photometry=True,
        use_exactly_the_comps_provided=True,
    )

    assert calibration is not None
    assert calibration["best_comp_index"] == 0
    assert calibration["comp_summaries"][0]["coverage_rejected"] is False


def test_select_comparison_calibrated_photometry_masks_overexposed_comp_measurements():
    frame_count = 24
    airmass = np.linspace(1.2, 1.0, frame_count)

    def build_psf_rows():
        rows = np.zeros((frame_count, 7), dtype=float)
        rows[:, 0] = 10.0
        rows[:, 1] = 20.0
        rows[:, 2] = 200.0
        rows[:, 3] = 1.0
        rows[:, 4] = 1.0
        return rows

    psf_data = {
        "target": build_psf_rows(),
        "comp1": build_psf_rows(),
        "comp2": build_psf_rows(),
    }
    aper_data = {
        "target": np.full((frame_count, 1, 1), 1000.0),
        "target_bg": np.full((frame_count, 1, 1), 10.0),
    }
    for key in ("comp1", "comp2"):
        aper_data[key] = np.full((frame_count, 1, 1), 100.0)
        aper_data[f"{key}_bg"] = np.full((frame_count, 1, 1), 10.0)

    comp_overexposed_masks = {
        "comp1": np.zeros(frame_count, dtype=bool),
        "comp2": np.zeros(frame_count, dtype=bool),
    }
    comp_overexposed_masks["comp2"][5] = True
    aper_data["comp2"][5, 0, 0] = np.nan

    calibration = select_comparison_calibrated_photometry(
        psf_data,
        aper_data,
        apers=np.array([2.5]),
        annuli=np.array([10.0]),
        airmass=airmass,
        comp_stars=[[10.0, 20.0], [30.0, 40.0]],
        sigma=1.0,
        use_psf_photometry=False,
        use_aperture_photometry=True,
        comp_overexposed_masks=comp_overexposed_masks,
    )
    comp2_summary = calibration["comp_summaries"][1]

    assert comp2_summary["overexposure_rejected_count"] == 1
    assert comp2_summary["coverage_count"] == frame_count - 1
    assert np.isnan(comp2_summary["ensemble_ratio_series"][5])
    assert calibration["best_comp_index"] == 0
    assert calibration["comp_summaries"][0]["coverage_count"] == frame_count


def test_cheap_lightcurve_prescore_treats_large_ratio_flag_as_noop():
    tflux = np.array([2.0, 2.0, 2.0, 6.0, 2.0, 2.0])
    cflux = np.full(tflux.shape[0], 2.0)
    airmass = np.linspace(1.0, 1.5, tflux.shape[0])

    score_with_flag = cheap_lightcurve_prescore(tflux, cflux, airmass, enforce_relative_flux_max=True)
    score_without_flag = cheap_lightcurve_prescore(tflux, cflux, airmass, enforce_relative_flux_max=False)

    assert np.isfinite(score_with_flag)
    assert np.isclose(score_with_flag, score_without_flag)


def test_normalize_flux_series_to_approximate_unity_scales_by_robust_baseline():
    flux = np.array([3.0, 3.3, 2.7, 3.0, 30.0], dtype=float)
    unc = np.full(flux.shape[0], 0.3, dtype=float)

    normalized_flux, normalized_unc, baseline = normalize_flux_series_to_approximate_unity(flux, unc)

    assert baseline == pytest.approx(3.0)
    assert np.nanmedian(normalized_flux[:4]) == pytest.approx(1.0)
    assert np.nanmedian(normalized_unc[:4]) == pytest.approx(0.1)


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

    result = detrend_flux_on_out_of_transit_baseline(times, flux, fluxerr, fit, min_side_points=2)

    assert result["applied"] is True
    assert np.allclose(result["flux"][[0, 1, 2, 4, 5, 6]], 1.0, atol=1e-8)
    assert result["flux"][3] == pytest.approx(0.99, abs=1e-8)
    assert result["slope"] == pytest.approx(0.02, abs=1e-8)
    assert result["reference_time_bjd_tdb"] == pytest.approx(0.0)


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
        oot_baseline_min_points_per_side=2,
    )

    assert len(captured["calls"]) == 2
    assert np.allclose(captured["calls"][0], flux)
    assert np.allclose(captured["calls"][1][[0, 1, 2, 4, 5, 6]], 1.0, atol=1e-8)
    assert refit_flux[3] == pytest.approx(0.99, abs=1e-8)
    assert np.allclose(refit_unc[[0, 1, 2, 4, 5, 6]], 0.01 / (1.0 + 0.02 * times[[0, 1, 2, 4, 5, 6]]))
    assert fit.oot_baseline_detrending_applied is True
    assert fit.oot_baseline_reference_time_bjd_tdb == pytest.approx(0.0)
    assert fit.oot_baseline_pre_points == 3
    assert fit.oot_baseline_post_points == 3


def test_fit_final_lightcurve_linear_detrend_does_not_reapply_fixed_airmass_baseline(monkeypatch):
    import exotic.exotic as exotic_module

    times = np.array([-2.0, -1.0, -0.25, 0.0, 0.25, 1.0, 2.0])
    transit_profile = np.array([1.0, 1.0, 1.0, 0.99, 1.0, 1.0, 1.0])
    flux = (1.03 + 0.02 * times) * transit_profile
    fluxerr = np.full_like(times, 0.01)
    airmass = np.linspace(1.0, 1.3, times.size)
    prior = {"rprs": 0.1, "tmid": 0.0, "inc": 89.0, "a0": 1.03, "a1": 1.03, "a2": 0.2}
    bounds = {
        "rprs": [0.0, 0.2],
        "tmid": [-0.1, 0.1],
        "inc": [84.0, 90.0],
        "a0": [0.95, 1.05],
        "a2": [-3.0, 3.0],
    }
    captured = {"calls": []}

    def fake_run_nested(
        call_times,
        call_flux,
        call_fluxerr,
        call_airmass,
        call_prior,
        call_bounds,
        **kwargs,
    ):
        captured["calls"].append({
            "flux": np.asarray(call_flux, dtype=float),
            "prior": dict(call_prior),
            "bounds": dict(call_bounds),
            "fixed_flux_baseline": kwargs.get("fixed_flux_baseline"),
            "fixed_parameter_errors": dict(kwargs.get("fixed_parameter_errors", {})),
        })
        return types.SimpleNamespace(
            time=np.asarray(call_times, dtype=float),
            data=np.asarray(call_flux, dtype=float),
            dataerr=np.asarray(call_fluxerr, dtype=float),
            airmass=np.asarray(call_airmass, dtype=float),
            transit=transit_profile.copy(),
            parameters={
                "tmid": 0.0,
                "rprs": 0.1,
                "inc": 89.0,
                "a0": call_prior.get("a0", 1.0),
                "a2": call_prior.get("a2", 0.2),
            },
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.001, "a2": 0.01},
            residuals=np.zeros_like(call_flux, dtype=float),
            duration_expected=0.5,
            duration_measured=0.5,
        )

    monkeypatch.setattr(exotic_module, "run_nested_lightcurve_fit_with_rprs_posterior_retry", fake_run_nested)
    monkeypatch.setattr(
        exotic_module,
        "build_final_fit_prefit_refinement_plan",
        lambda call_times, call_flux, call_fluxerr, call_airmass, call_prior, call_bounds, fit, **kwargs: {
            "applied": False,
            "note": "test no prefit refinement",
            "times": np.asarray(call_times, dtype=float),
            "flux": np.asarray(call_flux, dtype=float),
            "unc": np.asarray(call_fluxerr, dtype=float),
            "airmass": np.asarray(call_airmass, dtype=float),
            "jd_times": None,
            "prior": dict(call_prior),
            "bounds": dict(call_bounds),
            "duration": 0.5,
            "original_point_count": len(call_times),
            "refined_point_count": len(call_times),
            "trimmed_pre_points": 0,
            "trimmed_post_points": 0,
            "original_tmid_bounds": call_bounds["tmid"],
            "refined_tmid_bounds": call_bounds["tmid"],
        },
    )
    monkeypatch.setattr(exotic_module, "annotate_transit_detection_qc", lambda fit: None)

    fit, refit_flux, _ = fit_final_lightcurve_with_oot_baseline_detrending(
        times,
        flux,
        fluxerr,
        airmass,
        prior,
        bounds,
        detrend_on_outoftransit_baseline=True,
        oot_baseline_min_points_per_side=2,
        extend_sparse_posterior_live_points=False,
    )

    assert len(captured["calls"]) == 2
    final_call = captured["calls"][1]
    assert final_call["fixed_flux_baseline"] is True
    assert final_call["prior"]["a0"] == pytest.approx(1.0)
    assert final_call["prior"]["a2"] == pytest.approx(0.0)
    assert final_call["fixed_parameter_errors"]["a0"] > 0
    assert final_call["fixed_parameter_errors"]["a1"] == pytest.approx(
        final_call["fixed_parameter_errors"]["a0"],
    )
    assert final_call["fixed_parameter_errors"]["a2"] == pytest.approx(0.0)
    assert "a0" not in final_call["bounds"]
    assert "a2" not in final_call["bounds"]
    assert np.allclose(final_call["flux"][[0, 1, 2, 4, 5, 6]], 1.0, atol=1e-8)
    assert final_call["flux"][3] == pytest.approx(0.99, abs=1e-8)
    assert np.allclose(refit_flux, final_call["flux"])
    assert fit.oot_baseline_parameter_fit_applied is False
    assert "already flattened" in fit.oot_baseline_parameter_fit_note
    assert np.isfinite(fit.oot_baseline_parameter_fit_a0)
    assert fit.oot_baseline_parameter_fit_a0_error > 0
    assert np.isfinite(fit.oot_baseline_parameter_fit_a2)
    assert fit.oot_baseline_parameter_fit_a2_error > 0
    assert fit.pre_detrending_baseline_source.startswith("out-of-transit airmass/baseline")
    assert fit.pre_detrending_baseline_scale_parameter == "a0"
    assert fit.pre_detrending_baseline_scale_value == pytest.approx(
        fit.oot_baseline_parameter_fit_a0
    )
    assert fit.pre_detrending_baseline_scale_error == pytest.approx(
        fit.oot_baseline_parameter_fit_a0_error
    )
    assert fit.pre_detrending_baseline_a2_value == pytest.approx(
        fit.oot_baseline_parameter_fit_a2
    )
    assert fit.pre_detrending_baseline_a2_error == pytest.approx(
        fit.oot_baseline_parameter_fit_a2_error
    )


def test_fit_final_lightcurve_uses_oot_baseline_parameter_refit_when_linear_detrend_skips(monkeypatch):
    import exotic.exotic as exotic_module

    times = np.array([-0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03])
    transit_profile = np.array([0.99, 0.99, 0.99, 0.99, 1.0, 1.0, 1.0])
    airmass = np.linspace(1.0, 1.6, times.size)
    flux = np.exp(0.2 * (airmass - np.mean(airmass))) * transit_profile
    fluxerr = np.full_like(times, 0.01)
    prior = {"rprs": 0.1, "tmid": -0.01, "inc": 89.0, "a2": 0.0}
    bounds = {
        "rprs": [0.0, 0.2],
        "tmid": [-0.03, 0.01],
        "inc": [84.0, 90.0],
        "a0": [0.95, 1.05],
        "a2": [-3.0, 3.0],
    }
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
        baseline_fit_mask=None,
        fixed_parameter_errors=None,
    ):
        captured["calls"].append({
            "bounds": dict(call_bounds),
            "baseline_fit_mask": None if baseline_fit_mask is None else np.asarray(baseline_fit_mask, dtype=bool),
            "fixed_parameter_errors": dict(fixed_parameter_errors or {}),
            "prior": dict(call_prior),
        })
        return types.SimpleNamespace(
            transit=transit_profile.copy(),
            parameters={
                "tmid": -0.01,
                "rprs": 0.1,
                "inc": 89.0,
                "a2": call_prior.get("a2", 0.0),
                "a0": call_prior.get("a0", 1.0),
                "a1": call_prior.get("a0", 1.0),
            },
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a2": 0.01, "a0": 0.001, "a1": 0.001},
            data=np.array(call_flux, dtype=float),
            residuals=np.zeros_like(call_flux, dtype=float),
        )

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    fit, _, _ = fit_final_lightcurve_with_oot_baseline_detrending(
        times,
        flux,
        fluxerr,
        airmass,
        prior,
        bounds,
        detrend_on_outoftransit_baseline=True,
        oot_baseline_min_points_per_side=0,
    )

    assert len(captured["calls"]) == 2
    assert captured["calls"][0]["baseline_fit_mask"] is None
    assert captured["calls"][1]["baseline_fit_mask"].tolist() == [False, False, False, False, True, True, True]
    assert "a0" not in captured["calls"][1]["bounds"]
    assert "a2" not in captured["calls"][1]["bounds"]
    assert captured["calls"][1]["fixed_parameter_errors"]["a0"] > 0
    assert captured["calls"][1]["fixed_parameter_errors"]["a1"] == pytest.approx(
        captured["calls"][1]["fixed_parameter_errors"]["a0"],
    )
    assert "a2" in captured["calls"][1]["fixed_parameter_errors"]
    assert fit.oot_baseline_parameter_fit_applied is True
    assert fit.oot_baseline_detrending_applied is False


def test_phase_bin_sigma_clip_flags_local_phase_outlier():
    phase_centers = np.linspace(-0.045, 0.045, 10)
    phase = np.concatenate([center + np.linspace(-1e-4, 1e-4, 5) for center in phase_centers])
    base_profile = np.array([-0.002, -0.001, 0.0, 0.001, 0.002])
    values = np.concatenate([1.0 + base_profile for _ in phase_centers])
    values[27] = 1.15

    mask = phase_bin_sigma_clip(values, phase, sigma=3, bins=10)

    assert mask.sum() == 1
    assert mask[27]


def test_sigma_clip_preserves_an_exactly_constant_series():
    values = np.full(6, 3.0)
    times = np.linspace(0.0, 0.05, values.size)

    mask = sigma_clip(values, sigma=3, dt=5, times=times)

    assert not mask.any()


def test_sigma_clip_does_not_consume_process_global_random_state():
    values = 1.0 + 0.001 * np.sin(np.linspace(0.0, 4.0 * np.pi, 21))
    values[10] += 0.2
    state_before = np.random.get_state()

    sigma_clip(values, sigma=3, dt=11)

    state_after = np.random.get_state()
    assert state_after[0] == state_before[0]
    assert np.array_equal(state_after[1], state_before[1])
    assert state_after[2:] == state_before[2:]


def test_sigma_clip_respects_large_time_gaps_between_segments():
    times_pre = 2461151.80 + np.arange(50, dtype=float) * 0.00075
    times_post = 2461151.98 + np.arange(8, dtype=float) * 0.00075
    times = np.concatenate([times_pre, times_post])

    rng = np.random.default_rng(42)
    values_pre = (
        0.0235
        + 0.0006 * np.sin(np.linspace(0, 8 * np.pi, times_pre.size))
        + 0.0004 * np.linspace(0, 1, times_pre.size)
        + rng.normal(0, 5e-5, times_pre.size)
    )
    values_post = np.full(times_post.size, np.nanmedian(values_pre[-20:]) - 8e-4)
    values_post += rng.normal(0, 5e-5, times_post.size)
    values = np.concatenate([values_pre, values_post])
    values[54] += 0.8

    old_mask = sigma_clip(values, sigma=3, dt=37, times=None)
    gap_aware_mask = sigma_clip(values, sigma=3, dt=37, times=times)

    assert old_mask[54:58].all()
    assert not gap_aware_mask[54:58].any()
    assert gap_aware_mask.sum() < old_mask.sum()


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

    monkeypatch.setattr(exotic_module, "RPRS_RANGE_RESTRICTION_ENABLED", False)
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
    assert captured["calls"][1]["bounds"]["rprs"] == pytest.approx([0.0, 0.208])
    assert fit.rprs_posterior_refit_applied is True
    assert fit.rprs_posterior_refit_count == 1
    assert fit.rprs_posterior_refit_edge == "upper"
    assert fit.rprs_posterior_refit_bounds == pytest.approx([0.0, 0.208])


def test_fit_final_lightcurve_carries_retry_bounds_into_oot_baseline_refit(monkeypatch):
    import exotic.exotic as exotic_module

    monkeypatch.setattr(exotic_module, "RPRS_RANGE_RESTRICTION_ENABLED", False)
    times = np.array([-2.0, -1.0, -0.25, 0.0, 0.25, 1.0, 2.0])
    flux = (1.0 + 0.02 * times) * np.array([1.0, 1.0, 1.0, 0.99, 1.0, 1.0, 1.0])
    fluxerr = np.full_like(times, 0.01)
    airmass = np.ones_like(times)
    prior = {"rprs": 0.1, "tmid": 0.0, "inc": 89.0, "a2": 0.0}
    bounds = {"rprs": [0.0, 0.125], "tmid": [-0.1, 0.1], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]}
    transit_model = np.array([1.0, 1.0, 1.0, 0.99, 1.0, 1.0, 1.0])
    diagnostics = [
        {
            "clipped": True,
            "edge": "upper",
            "mode": 0.158,
            "std": 0.006,
            "bounds": [0.128, 0.188],
            "reason": "posterior peaks against the upper search bound.",
        },
        {
            "clipped": False,
            "edge": None,
            "mode": 0.159,
            "std": 0.005,
            "bounds": [0.108, 0.208],
            "reason": "posterior support is comfortably inside the sampled bounds.",
        },
        {
            "clipped": False,
            "edge": None,
            "mode": 0.160,
            "std": 0.005,
            "bounds": [0.108, 0.208],
            "reason": "posterior support is comfortably inside the sampled bounds.",
        },
    ]
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
        duration_prior=None,
    ):
        call_index = len(captured["calls"])
        call_diagnostics = diagnostics[min(call_index, len(diagnostics) - 1)]
        captured["calls"].append({
            "flux": np.array(call_flux, dtype=float),
            "prior": dict(call_prior),
            "bounds": {
                key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value
                for key, value in call_bounds.items()
            },
        })
        fit = types.SimpleNamespace(
            transit=transit_model,
            parameters={"tmid": 0.0, "rprs": call_diagnostics["mode"], "inc": 89.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a2": 0.01},
            data=np.array(call_flux, dtype=float),
            residuals=np.zeros_like(call_flux, dtype=float),
        )
        fit.get_parameter_posterior_recenter_diagnostics = (
            lambda key: dict(call_diagnostics) if key == "rprs" else None
        )
        return fit

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    fit, _, _ = fit_final_lightcurve_with_oot_baseline_detrending(
        times,
        flux,
        fluxerr,
        airmass,
        prior,
        bounds,
        detrend_on_outoftransit_baseline=True,
        oot_baseline_min_points_per_side=2,
    )

    assert len(captured["calls"]) == 3
    assert captured["calls"][0]["bounds"]["rprs"] == pytest.approx([0.0, 0.125])
    assert captured["calls"][1]["bounds"]["rprs"] == pytest.approx([0.0, 0.208])
    assert captured["calls"][2]["bounds"]["rprs"] == pytest.approx([0.0, 0.208])
    assert np.allclose(captured["calls"][2]["flux"][[0, 1, 2, 4, 5, 6]], 1.0, atol=1e-8)
    assert fit.oot_baseline_detrending_applied is True


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


def test_fit_final_lightcurve_prefit_refinement_skips_one_sided_transit_solution(monkeypatch):
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
        call_times = np.array(call_times, dtype=float)
        captured["calls"].append({
            "times": call_times,
            "bounds": {
                key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value
                for key, value in call_bounds.items()
            },
        })
        transit = np.ones_like(call_times, dtype=float)
        transit[call_times >= 0.0] = 0.95
        return types.SimpleNamespace(
            duration_expected=2.0,
            duration_measured=2.0,
            parameters={"tmid": 2.0, "rprs": 0.1, "inc": 89.0, "a2": 0.0, "per": 10.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a2": 0.01},
            data=np.array(call_flux, dtype=float),
            residuals=np.zeros_like(call_flux, dtype=float),
            time=call_times,
            transit=transit,
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

    assert len(captured["calls"]) == 1
    assert captured["calls"][0]["times"] == pytest.approx(times)
    assert trimmed_flux == pytest.approx(flux)
    assert trimmed_unc == pytest.approx(fluxerr)
    assert fit.prefit_refinement_applied is False
    assert "one side of the modeled transit" in fit.prefit_refinement_note


def test_fit_final_lightcurve_prefit_refinement_does_not_expand_tmid_past_original_bounds(monkeypatch):
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
        call_times = np.array(call_times, dtype=float)
        captured["calls"].append({
            "times": call_times,
            "bounds": {
                key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value
                for key, value in call_bounds.items()
            },
        })
        transit = np.where(np.abs(call_times - 0.4) <= 0.25, 0.98, 1.0)
        return types.SimpleNamespace(
            duration_expected=0.5,
            duration_measured=0.5,
            parameters={"tmid": 0.4, "rprs": 0.1, "inc": 89.0, "a2": 0.0, "per": 10.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a2": 0.01},
            data=np.array(call_flux, dtype=float),
            residuals=np.zeros_like(call_flux, dtype=float),
            time=call_times,
            transit=transit,
        )

    monkeypatch.setattr(exotic_module, "lc_fitter", fake_lc_fitter)

    times = np.array([-0.8, -0.4, 0.0, 0.4, 0.8], dtype=float)
    flux = np.ones(times.shape[0], dtype=float)
    fluxerr = np.full(times.shape[0], 0.01, dtype=float)
    airmass = np.ones(times.shape[0], dtype=float)
    prior = {"tmid": 0.0, "rprs": 0.1, "inc": 89.0, "a2": 0.0, "per": 10.0}
    bounds = {"rprs": [0.0, 0.2], "tmid": [-0.5, 0.5], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]}

    fit, _, _ = fit_final_lightcurve_with_oot_baseline_detrending(
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
    assert captured["calls"][1]["bounds"]["tmid"] == pytest.approx([0.15, 0.5])
    assert fit.prefit_refinement_tmid_bounds == pytest.approx([0.15, 0.5])


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
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

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
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

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


def test_fit_lightcurve_centers_vertical_flux_bound_on_normalized_flux(monkeypatch):
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
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

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
    assert captured["prior"]["a0"] == pytest.approx(1.0)
    assert captured["prior"]["a1"] == pytest.approx(1.0)
    assert captured["bounds"]["a0"] == pytest.approx([0.95, 1.05])


def test_fit_lightcurve_rejects_undersampled_series(monkeypatch):
    called = {"count": 0}

    def fake_lc_fitter(*args, **kwargs):
        called["count"] += 1
        return types.SimpleNamespace()

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

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
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )
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
    captured_duration_priors = []

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
        duration_prior=None,
    ):
        captured_modes.append(mode)
        captured_duration_priors.append(duration_prior)
        return types.SimpleNamespace()

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

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
    assert captured_duration_priors[0] is None
    assert captured_duration_priors[1] is not None
    assert captured_duration_priors[1]["applied"] is True
    assert captured_duration_priors[1]["expected_duration"] > 0
    assert myfit.pre_ultranest_transit_coverage_valid is True
    assert myfit.pre_ultranest_transit_coverage_expected_successful is True


def test_fit_lightcurve_attaches_frame_filter_diagnostics(monkeypatch):
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
        return types.SimpleNamespace(
            time=np.asarray(times, dtype=float),
            data=np.asarray(fluxes, dtype=float),
            dataerr=np.asarray(flux_unc, dtype=float),
            detrended=np.asarray(fluxes, dtype=float),
            detrendederr=np.asarray(flux_unc, dtype=float),
            airmass=np.asarray(airmass, dtype=float),
            airmass_model=np.ones(len(times), dtype=float),
            residuals=np.zeros(len(times), dtype=float),
            phase=np.linspace(-0.1, 0.1, len(times)),
            transit=np.ones(len(times), dtype=float),
            model=np.asarray(fluxes, dtype=float),
            wf=np.ones(len(times), dtype=float),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a1": 1.0, "a2": 0.0, "per": 1.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a1": 0.01, "a2": 0.01},
        )

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.array([False] * (len(data) - 1) + [True], dtype=bool),
    )
    monkeypatch.setattr(
        "exotic.exotic.phase_bin_sigma_clip",
        lambda values, phase, sigma=3, bins=10, min_points=5, max_iters=3: np.zeros(len(values), dtype=bool),
    )

    times = np.arange(12, dtype=float)
    tflux = np.full(times.shape[0], 100.0, dtype=float)
    cflux = np.full(times.shape[0], 50.0, dtype=float)
    cflux[1] = 0.0
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
        "midT": 0.5,
        "midTUnc": 0.001,
        "pPerUnc": 0.001,
    }

    myfit, _, _ = fit_lightcurve(times, tflux, cflux, airmass, ld, p_dict, jd_times)

    assert myfit is not None
    diagnostics = myfit.frame_filter_diagnostics
    assert [diagnostic["stage"] for diagnostic in diagnostics] == [
        "Target/reference ratio filter",
        "Initial sigma clip",
        "Pre-fit raw-ratio outlier clip",
        "Finite/positive photometry filter",
    ]
    assert diagnostics[0]["dropped_point_count"] == 1
    assert diagnostics[0]["first_dropped_time"] == pytest.approx(1.0)
    assert diagnostics[1]["dropped_point_count"] == 1
    assert diagnostics[1]["first_dropped_time"] == pytest.approx(11.0)
    assert diagnostics[2]["dropped_point_count"] == 0
    assert diagnostics[3]["dropped_point_count"] == 0


def test_evaluate_transit_detection_qc_prefers_transit_model():
    transit_model = np.ones(21, dtype=float)
    transit_model[8:13] = 0.99
    data = transit_model + np.array(
        [
            0.0002, -0.0001, 0.0001, -0.0002, 0.0000, 0.0001, -0.0001,
            0.0002, -0.0002, 0.0001, -0.0001, 0.0002, -0.0002, 0.0001,
            0.0000, -0.0001, 0.0002, -0.0001, 0.0001, 0.0000, -0.0001,
        ],
        dtype=float,
    )
    fit = types.SimpleNamespace(
        data=data,
        dataerr=np.full(data.shape[0], 0.0015, dtype=float),
        model=transit_model,
        airmass=np.ones(data.shape[0], dtype=float),
        airmass_fit_skipped=True,
        parameters={"rprs": 0.10, "tmid": 0.5, "inc": 89.0, "a2": 0.0},
        errors={"rprs": 0.01, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
        bounds={"rprs": [0.0, 1.0], "tmid": [0.4, 0.6], "inc": [80.0, 90.0]},
        duration_expected=5.0,
        duration_measured=5.0,
    )

    summary = evaluate_transit_detection_qc(fit)

    assert summary["computed"] is True
    assert summary["preferred_model"] == "transit"
    assert summary["status"] == "pass"
    assert summary["delta_bic"] > 10.0
    assert summary["delta_chi2"] > 0.0


def test_evaluate_transit_detection_qc_passes_strong_model_with_low_rprs_precision():
    transit_model = np.ones(21, dtype=float)
    transit_model[8:13] = 0.99
    data = transit_model + np.array(
        [
            0.0002, -0.0001, 0.0001, -0.0002, 0.0000, 0.0001, -0.0001,
            0.0002, -0.0002, 0.0001, -0.0001, 0.0002, -0.0002, 0.0001,
            0.0000, -0.0001, 0.0002, -0.0001, 0.0001, 0.0000, -0.0001,
        ],
        dtype=float,
    )
    fit = types.SimpleNamespace(
        data=data,
        dataerr=np.full(data.shape[0], 0.0015, dtype=float),
        model=transit_model,
        airmass=np.ones(data.shape[0], dtype=float),
        airmass_fit_skipped=True,
        parameters={"rprs": 0.10, "tmid": 0.5, "inc": 89.0, "a2": 0.0},
        errors={"rprs": 0.20, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
        bounds={"rprs": [0.0, 1.0], "tmid": [0.4, 0.6], "inc": [80.0, 90.0]},
        duration_expected=5.0,
        duration_measured=5.0,
    )

    summary = evaluate_transit_detection_qc(fit)

    assert summary["computed"] is True
    assert summary["status"] == "pass"
    assert summary["rprs_sigma"] == pytest.approx(0.5)
    assert summary["ktmf_metric"] >= 3.5
    assert "not used as a transit-detection veto" in " ".join(summary["notes"])


def test_evaluate_transit_detection_qc_uses_ktmf_marginal_band_despite_weak_bic(monkeypatch):
    monkeypatch.setattr(
        "exotic.exotic.compute_transit_qc_ktmf",
        lambda summary: (3.34, []),
    )
    transit_model = np.ones(21, dtype=float)
    transit_model[8:13] = 0.99
    data = np.ones(21, dtype=float)
    fit = types.SimpleNamespace(
        data=data,
        dataerr=np.full(data.shape[0], 0.02, dtype=float),
        model=transit_model,
        airmass=np.ones(data.shape[0], dtype=float),
        airmass_fit_skipped=True,
        parameters={"rprs": 0.10, "tmid": 0.5, "inc": 89.0, "a2": 0.0},
        errors={"rprs": 0.02, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
        bounds={"rprs": [0.0, 1.0], "tmid": [0.4, 0.6], "inc": [80.0, 90.0]},
        duration_expected=5.0,
        duration_measured=5.0,
    )

    summary = evaluate_transit_detection_qc(fit)

    assert summary["computed"] is True
    assert summary["delta_bic"] < 6.0
    assert summary["ktmf_metric"] == pytest.approx(3.34)
    assert summary["status"] == "marginal"
    assert "KTMF indicates a marginal transit fit" in summary["summary"]


def test_evaluate_transit_detection_qc_uses_midpoint_anchored_duration_for_partial():
    times = np.linspace(0.0, 3.0, 13)
    transit_model = np.ones(times.shape[0], dtype=float)
    transit_model[times <= 2.25] = 0.99
    data = transit_model + np.array(
        [
            0.0002, -0.0001, 0.0001, -0.0002, 0.0000, 0.0001, -0.0001,
            0.0002, -0.0002, 0.0001, -0.0001, 0.0002, -0.0002,
        ],
        dtype=float,
    )
    fit = types.SimpleNamespace(
        time=times,
        data=data,
        dataerr=np.full(data.shape[0], 0.0015, dtype=float),
        transit=transit_model,
        model=transit_model,
        airmass=np.ones(data.shape[0], dtype=float),
        airmass_fit_skipped=True,
        parameters={"rprs": 0.10, "tmid": 0.0, "inc": 89.0, "a2": 0.0},
        errors={"rprs": 0.01, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
        bounds={"rprs": [0.0, 1.0], "tmid": [-0.1, 0.1], "inc": [80.0, 90.0]},
        duration_expected=5.0,
        duration_measured=2.5,
        pre_ultranest_transit_coverage={
            "valid": True,
            "covers_ingress": False,
            "covers_mid_transit": True,
            "covers_egress": True,
            "observed_segment": "mid-transit to egress",
            "expected_tmid": 0.0,
        },
    )

    summary = evaluate_transit_detection_qc(fit)
    contributions_by_label = {
        contribution["label"]: contribution
        for contribution in summary["ktmf_contributions"]
    }

    assert summary["duration_measured_for_qc"] == pytest.approx(4.75)
    assert summary["duration_ratio"] == pytest.approx(0.95)
    assert contributions_by_label["Duration Consistency"]["available"] is True
    assert "midpoint-anchored partial estimate" in contributions_by_label["Duration Consistency"]["detail"]


def test_evaluate_transit_detection_qc_skips_duration_for_edge_only_partial():
    times = np.linspace(-3.0, -0.25, 12)
    transit_model = np.ones(times.shape[0], dtype=float)
    transit_model[times >= -2.5] = 0.99
    data = transit_model + np.array(
        [
            0.0002, -0.0001, 0.0001, -0.0002, 0.0000, 0.0001,
            -0.0001, 0.0002, -0.0002, 0.0001, -0.0001, 0.0002,
        ],
        dtype=float,
    )
    fit = types.SimpleNamespace(
        time=times,
        data=data,
        dataerr=np.full(data.shape[0], 0.0015, dtype=float),
        transit=transit_model,
        model=transit_model,
        airmass=np.ones(data.shape[0], dtype=float),
        airmass_fit_skipped=True,
        parameters={"rprs": 0.10, "tmid": 0.0, "inc": 89.0, "a2": 0.0},
        errors={"rprs": 0.01, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
        bounds={"rprs": [0.0, 1.0], "tmid": [-0.1, 0.1], "inc": [80.0, 90.0]},
        duration_expected=5.0,
        duration_measured=2.75,
        pre_ultranest_transit_coverage={
            "valid": True,
            "covers_ingress": True,
            "covers_mid_transit": False,
            "covers_egress": False,
            "observed_segment": "ingress-only partial",
            "expected_tmid": 0.0,
        },
    )

    summary = evaluate_transit_detection_qc(fit)
    contributions_by_label = {
        contribution["label"]: contribution
        for contribution in summary["ktmf_contributions"]
    }

    assert np.isnan(summary["duration_ratio"])
    assert summary["duration_consistency_applicable"] is False
    assert contributions_by_label["Duration Consistency"]["available"] is False
    assert "only partially observed" in contributions_by_label["Duration Consistency"]["detail"]


def test_evaluate_transit_detection_qc_fails_when_flat_model_is_better():
    transit_model = np.ones(21, dtype=float)
    transit_model[8:13] = 0.99
    data = np.ones(21, dtype=float) + np.array(
        [
            0.0002, -0.0001, 0.0001, -0.0002, 0.0000, 0.0001, -0.0001,
            0.0002, -0.0002, 0.0001, -0.0001, 0.0002, -0.0002, 0.0001,
            0.0000, -0.0001, 0.0002, -0.0001, 0.0001, 0.0000, -0.0001,
        ],
        dtype=float,
    )
    fit = types.SimpleNamespace(
        data=data,
        dataerr=np.full(data.shape[0], 0.0015, dtype=float),
        model=transit_model,
        airmass=np.ones(data.shape[0], dtype=float),
        airmass_fit_skipped=True,
        parameters={"rprs": 0.10, "tmid": 0.5, "inc": 89.0, "a2": 0.0},
        errors={"rprs": 0.01, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
        bounds={"rprs": [0.0, 1.0], "tmid": [0.4, 0.6], "inc": [80.0, 90.0]},
        duration_expected=5.0,
        duration_measured=5.0,
    )

    summary = evaluate_transit_detection_qc(fit)

    assert summary["computed"] is True
    assert summary["status"] == "fail"
    assert summary["preferred_model"] == "flat"
    assert summary["delta_chi2"] < 0.0


def test_evaluate_transit_detection_qc_marks_large_expected_value_deviation_fail_via_ktmf():
    transit_model = np.ones(21, dtype=float)
    transit_model[8:13] = 0.99
    data = transit_model + np.array(
        [
            0.0002, -0.0001, 0.0001, -0.0002, 0.0000, 0.0001, -0.0001,
            0.0002, -0.0002, 0.0001, -0.0001, 0.0002, -0.0002, 0.0001,
            0.0000, -0.0001, 0.0002, -0.0001, 0.0001, 0.0000, -0.0001,
        ],
        dtype=float,
    )
    fit = types.SimpleNamespace(
        data=data,
        dataerr=np.full(data.shape[0], 0.0015, dtype=float),
        model=transit_model,
        airmass=np.ones(data.shape[0], dtype=float),
        airmass_fit_skipped=True,
        parameters={"rprs": 0.18, "tmid": 0.506, "inc": 89.0, "a2": 0.0},
        errors={"rprs": 0.01, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
        bounds={"rprs": [0.0, 1.0], "tmid": [0.4, 0.6], "inc": [80.0, 90.0]},
        duration_expected=5.0,
        duration_measured=5.0,
        transit_qc_expected_tmid=0.5,
        transit_qc_expected_tmid_unc=0.001,
        transit_qc_expected_rprs=0.10,
        transit_qc_expected_rprs_unc=0.01,
        transit_qc_use_deviation_from_expected_transit_in_qc=True,
        transit_qc_deviation_sigma_threshold=5.0,
    )

    summary = evaluate_transit_detection_qc(fit)

    assert summary["computed"] is True
    assert summary["status"] == "fail"
    expected_comparison_unc = np.sqrt(0.01 ** 2 + 0.01 ** 2 + (0.05 * 0.10) ** 2)
    assert summary["rprs_deviation_unc"] == pytest.approx(expected_comparison_unc)
    assert summary["rprs_deviation_systematic_floor"] == pytest.approx(0.05 * 0.10)
    assert summary["rprs_deviation_sigma"] == pytest.approx(abs(0.18 - 0.10) / expected_comparison_unc)
    assert summary["deviation_from_expected_value"] == pytest.approx(0.0)
    assert summary["ktmf_metric"] <= 5.0
    assert summary["ktmf_metric"] < 3.0
    assert np.isnan(summary["tmid_deviation_sigma"])
    assert np.isnan(summary["tmid_deviation_minutes"])
    assert summary["rprs_deviation_fit_unc"] == pytest.approx(0.01)


def test_annotate_transit_qc_expected_values_prefers_propagated_epoch_tmid():
    fit = types.SimpleNamespace(
        initial_tmid_search_tmid=2460658.8654321,
        initial_tmid_search_uncertainty=0.0025,
    )

    annotate_transit_qc_expected_values(
        fit,
        {
            "midT": 2455867.402743,
            "midTUnc": 4.9e-05,
            "rprs": 0.1488,
            "rprsUnc": 0.00055,
        },
    )

    assert fit.transit_qc_expected_tmid == pytest.approx(2460658.8654321)
    assert fit.transit_qc_expected_tmid_unc == pytest.approx(0.0025)
    assert fit.transit_qc_expected_rprs == pytest.approx(0.1488)
    assert fit.transit_qc_expected_rprs_unc == pytest.approx(0.00055)


def test_annotate_transit_qc_expected_values_coerces_scalar_like_inputs():
    fit = types.SimpleNamespace(
        initial_tmid_search_tmid=np.array(["2460658.8654321"]),
        initial_tmid_search_uncertainty="0.0025",
    )

    annotate_transit_qc_expected_values(
        fit,
        {
            "midT": "2455867.402743",
            "midTUnc": ["4.9e-05"],
            "rprs": "0.1488",
            "rprsUnc": np.array(["0.00055"]),
            "use_deviation_from_expected_transit_in_qc": "n",
            "deviation_from_expected_transit_in_qc_sigma": "7.5",
        },
    )

    assert fit.transit_qc_expected_tmid == pytest.approx(2460658.8654321)
    assert fit.transit_qc_expected_tmid_unc == pytest.approx(0.0025)
    assert fit.transit_qc_expected_rprs == pytest.approx(0.1488)
    assert fit.transit_qc_expected_rprs_unc == pytest.approx(0.00055)
    assert fit.transit_qc_use_deviation_from_expected_transit_in_qc is False
    assert fit.transit_qc_deviation_sigma_threshold == pytest.approx(7.5)


def test_evaluate_transit_detection_qc_does_not_calculate_tmid_expected_value_deviation():
    times = np.linspace(0.0, 1.0, 21)
    transit_model = np.ones(times.shape[0], dtype=float)
    transit_model[9:12] -= 0.02
    data = transit_model + np.array(
        [
            0.0001, -0.0001, 0.0002, -0.0002, 0.0000, 0.0001, -0.0001,
            0.0002, -0.0002, 0.0001, 0.0000, -0.0001, 0.0002, -0.0002,
            0.0001, 0.0000, -0.0001, 0.0001, -0.0001, 0.0000, 0.0001,
        ],
        dtype=float,
    )
    fit = types.SimpleNamespace(
        time=times,
        data=data,
        dataerr=np.full(data.shape[0], 0.0015, dtype=float),
        model=transit_model,
        airmass=np.ones(data.shape[0], dtype=float),
        airmass_fit_skipped=True,
        parameters={"rprs": 0.10, "tmid": 0.528, "inc": 89.0, "a2": 0.0, "per": 2.0},
        errors={"rprs": 0.01, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
        bounds={"rprs": [0.0, 1.0], "tmid": [0.4, 0.6], "inc": [80.0, 90.0]},
        prior={"per": 2.0, "tmid": 0.5},
        duration_expected=5.0,
        duration_measured=5.0,
        transit_qc_expected_tmid=0.5,
        transit_qc_expected_tmid_unc=0.001,
        transit_qc_expected_rprs=0.10,
        transit_qc_expected_rprs_unc=0.01,
        transit_qc_use_deviation_from_expected_transit_in_qc=True,
        transit_qc_deviation_sigma_threshold=5.0,
    )

    summary = evaluate_transit_detection_qc(fit)

    assert summary["status"] == "pass"
    assert np.isnan(summary["tmid_deviation_minutes"])
    assert np.isnan(summary["tmid_deviation_sigma"])
    assert summary["rprs_deviation_sigma"] == pytest.approx(0.0)
    assert summary["deviation_from_expected_value"] == pytest.approx(1.0)
    assert not any("Expected-value Tmid" in note for note in summary["notes"])
    assert "QC rejected the fit because" not in summary["summary"]
    assert "Tmid of the fit is 40.32 minutes away from the ephemeris Tmid" not in summary["summary"]
    assert "not supported strongly enough against a flat/null model" not in summary["summary"]


def test_expected_value_rprs_deviation_uses_combined_uncertainty_with_systematic_floor():
    transit_model = np.ones(21, dtype=float)
    transit_model[9:12] -= 0.0287
    data = transit_model.copy()
    fit = types.SimpleNamespace(
        data=data,
        dataerr=np.full(data.shape[0], 0.0015, dtype=float),
        model=transit_model,
        airmass=np.ones(data.shape[0], dtype=float),
        airmass_fit_skipped=True,
        parameters={"rprs": 0.1694, "tmid": 0.5, "inc": 89.0, "a2": 0.0},
        errors={"rprs": 0.0046, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
        bounds={"rprs": [0.0, 0.5], "tmid": [0.4, 0.6], "inc": [80.0, 90.0]},
        duration_expected=5.0,
        duration_measured=5.0,
        transit_qc_expected_tmid=0.5,
        transit_qc_expected_tmid_unc=0.001,
        transit_qc_expected_rprs=0.1589,
        transit_qc_expected_rprs_unc=0.0001,
        transit_qc_use_deviation_from_expected_transit_in_qc=True,
        transit_qc_deviation_sigma_threshold=5.0,
    )

    summary = evaluate_transit_detection_qc(fit)

    expected_comparison_unc = np.sqrt(0.0046 ** 2 + 0.0001 ** 2 + (0.05 * 0.1589) ** 2)
    assert summary["rprs_deviation_fit_unc"] == pytest.approx(0.0046)
    assert summary["rprs_deviation_expected_unc"] == pytest.approx(0.0001)
    assert summary["rprs_deviation_systematic_floor"] == pytest.approx(0.05 * 0.1589)
    assert summary["rprs_deviation_unc"] == pytest.approx(expected_comparison_unc)
    assert summary["rprs_deviation_sigma"] == pytest.approx(abs(0.1694 - 0.1589) / expected_comparison_unc)
    assert summary["deviation_from_expected_value"] == pytest.approx(
        1.0 - summary["rprs_deviation_sigma"] / 5.0
    )
    assert summary["deviation_from_expected_value"] > 0.0


def test_expected_value_rprs_deviation_uses_model_data_fit_uncertainty():
    transit_model = np.ones(31, dtype=float)
    transit_model[12:19] -= 0.0287
    residual_pattern = np.array(
        [
            0.0, 0.006, -0.005, 0.004, -0.006, 0.005, -0.004, 0.006,
            -0.005, 0.004, -0.006, 0.005, -0.004, 0.006, -0.005, 0.004,
            -0.006, 0.005, -0.004, 0.006, -0.005, 0.004, -0.006, 0.005,
            -0.004, 0.006, -0.005, 0.004, -0.006, 0.005, 0.0,
        ],
        dtype=float,
    )
    fit = types.SimpleNamespace(
        time=np.linspace(0.0, 1.0, transit_model.size),
        data=transit_model + residual_pattern,
        dataerr=np.full(transit_model.size, 0.003, dtype=float),
        model=transit_model,
        transit=transit_model,
        airmass=np.ones(transit_model.size, dtype=float),
        airmass_fit_skipped=True,
        parameters={"rprs": 0.1694, "tmid": 0.5, "inc": 89.0, "a2": 0.0},
        errors={"rprs": 0.0046, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
        bounds={"rprs": [0.0, 0.5], "tmid": [0.4, 0.6], "inc": [80.0, 90.0]},
        duration_expected=5.0,
        duration_measured=5.0,
        transit_qc_expected_tmid=0.5,
        transit_qc_expected_tmid_unc=0.001,
        transit_qc_expected_rprs=0.1589,
        transit_qc_expected_rprs_unc=0.0001,
        transit_qc_use_deviation_from_expected_transit_in_qc=True,
        transit_qc_deviation_sigma_threshold=5.0,
    )

    summary = evaluate_transit_detection_qc(fit)

    expected_fit_unc = np.sqrt(
        summary["rprs_deviation_model_fit_unc"] ** 2
        + summary["rprs_deviation_data_fit_unc"] ** 2
    )
    expected_comparison_unc = np.sqrt(
        expected_fit_unc ** 2
        + summary["rprs_deviation_expected_unc"] ** 2
        + summary["rprs_deviation_systematic_floor"] ** 2
    )
    assert summary["rprs_deviation_model_fit_unc"] == pytest.approx(0.0046)
    assert summary["rprs_deviation_data_fit_unc"] > 0.0
    assert summary["rprs_deviation_fit_unc"] == pytest.approx(expected_fit_unc)
    assert summary["rprs_deviation_unc"] == pytest.approx(expected_comparison_unc)
    contribution = next(
        item for item in summary["ktmf_contributions"]
        if item["label"] == "Deviation From Expected Value"
    )
    assert "model uncertainty=0.004600" in contribution["detail"]
    assert "data/red-noise uncertainty=" in contribution["detail"]


def test_selected_full_resolution_refit_keeps_expected_value_context(monkeypatch):
    import exotic.exotic as exotic_module

    times = np.linspace(0.0, 1.0, 21)
    transit_model = np.ones(times.shape[0], dtype=float)
    transit_model[9:12] -= 0.0287
    errors = np.full(times.shape[0], 0.0015, dtype=float)
    airmass = np.ones(times.shape[0], dtype=float)

    previous_fit = types.SimpleNamespace(
        fast_ultranest_binning_applied=True,
        parameters={
            "rprs": 0.1694,
            "tmid": 0.5,
            "ars": 5.0,
            "inc": 89.0,
            "per": 1.0,
            "u0": 0.1,
            "u1": 0.1,
            "u2": 0.1,
            "u3": 0.1,
            "ecc": 0.0,
            "omega": 0.0,
            "a0": 1.0,
            "a1": 1.0,
            "a2": 0.0,
        },
        errors={"rprs": 0.0046, "tmid": 0.001, "ars": 0.1, "inc": 0.1, "a0": 0.01, "a2": 0.01},
        bounds={"rprs": [0.0, 0.5], "tmid": [0.4, 0.6], "ars": [1.0, 10.0], "inc": [80.0, 90.0]},
    )

    def fake_run_nested(*args, **kwargs):
        return types.SimpleNamespace(
            time=times,
            data=transit_model.copy(),
            dataerr=errors.copy(),
            model=transit_model.copy(),
            airmass=airmass.copy(),
            prior={"per": 1.0, "tmid": 0.5},
            parameters={
                "rprs": 0.1694,
                "tmid": 0.5,
                "ars": 5.0,
                "inc": 89.0,
                "per": 1.0,
                "a0": 1.0,
                "a1": 1.0,
                "a2": 0.0,
            },
            errors={"rprs": 0.0046, "tmid": 0.001, "ars": 0.1, "inc": 0.1, "a0": 0.01, "a2": 0.01},
            bounds={"rprs": [0.0, 0.5], "tmid": [0.4, 0.6], "ars": [1.0, 10.0], "inc": [80.0, 90.0]},
            airmass_fit_skipped=True,
            eebls_diagnostic_depth_snr=50.0,
            duration_expected=0.1,
            duration_measured=0.1,
        )

    monkeypatch.setattr(exotic_module, "run_nested_lightcurve_fit_with_rprs_posterior_retry", fake_run_nested)
    monkeypatch.setattr(exotic_module, "build_expected_transit_coverage_assessment", lambda *args, **kwargs: {})
    monkeypatch.setattr(exotic_module, "log_expected_transit_coverage_assessment", lambda *args, **kwargs: None)
    monkeypatch.setattr(exotic_module, "annotate_pre_ultranest_transit_coverage", lambda *args, **kwargs: None)
    monkeypatch.setattr(exotic_module, "selected_final_live_point_target", lambda *args, **kwargs: (200, None))

    p_dict = {
        "rprs": 0.1589,
        "rprsUnc": 0.0001,
        "midT": 0.5,
        "midTUnc": 0.001,
        "pPer": 1.0,
        "pPerUnc": 0.0,
        "aRs": 5.0,
        "aRsUnc": 0.1,
        "inc": 89.0,
        "ecc": 0.0,
        "omega": 0.0,
        "use_deviation_from_expected_transit_in_qc": True,
        "deviation_from_expected_transit_in_qc_sigma": 5.0,
    }
    selected_result = {
        "fit": previous_fit,
        "good_times": times,
        "good_flux": transit_model.copy(),
        "good_unc": errors.copy(),
        "good_airmass": airmass.copy(),
        "good_jd_times": times.copy(),
        "fast_fit_bounds": previous_fit.bounds,
    }

    refit, _, _ = refit_selected_fast_comparison_on_full_lightcurve(
        selected_result,
        p_dict,
        detrend_on_outoftransit_baseline=False,
        duration_prior={"duration": 0.1},
    )

    expected_comparison_unc = np.sqrt(0.0046 ** 2 + 0.0001 ** 2 + (0.05 * 0.1589) ** 2)
    assert refit.transit_qc_rprs_deviation_fit_unc == pytest.approx(0.0046)
    assert refit.transit_qc_rprs_deviation_expected_unc == pytest.approx(0.0001)
    assert refit.transit_qc_rprs_deviation_systematic_floor == pytest.approx(0.05 * 0.1589)
    assert refit.transit_qc_rprs_deviation_unc == pytest.approx(expected_comparison_unc)
    assert refit.transit_qc_rprs_deviation_sigma == pytest.approx(abs(0.1694 - 0.1589) / expected_comparison_unc)
    contribution = next(
        item for item in refit.transit_qc_ktmf_contributions
        if item["label"] == "Deviation From Expected Value"
    )
    assert contribution["score"] > 0.0
    assert contribution["max_points"] > 0.0
    assert "fit uncertainty=0.004600" in contribution["detail"]
    assert "expected uncertainty=0.000100" in contribution["detail"]
    assert "comparison uncertainty=" in contribution["detail"]
    assert "systematic floor=" in contribution["detail"]
    assert "Tmid" not in contribution["detail"]


def test_evaluate_transit_detection_qc_computes_missing_eebls_depth_snr(monkeypatch):
    def fake_eebls(times, flux_values, flux_errors, prior, fallback_bounds):
        return {
            "method": "eebls",
            "applied": True,
            "tmid": 0.5,
            "bounds": [0.45, 0.55],
            "duration": 0.1,
            "depth": 0.01,
            "depth_snr": 7.25,
            "note": "test eebls diagnostic",
        }

    monkeypatch.setattr("exotic.exotic.estimate_tmid_and_bounds_with_eebls", fake_eebls)

    times = np.linspace(0.0, 1.0, 21)
    transit_model = np.ones(times.shape[0], dtype=float)
    transit_model[9:12] -= 0.02
    data = transit_model.copy()
    fit = types.SimpleNamespace(
        time=times,
        data=data,
        dataerr=np.full(data.shape[0], 0.0015, dtype=float),
        model=transit_model,
        airmass=np.ones(data.shape[0], dtype=float),
        airmass_fit_skipped=True,
        parameters={"rprs": 0.10, "tmid": 0.5, "inc": 89.0, "a2": 0.0, "per": 2.0},
        errors={"rprs": 0.01, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
        bounds={"rprs": [0.0, 1.0], "tmid": [0.4, 0.6], "inc": [80.0, 90.0]},
        prior={"per": 2.0, "tmid": 0.5},
        duration_expected=5.0,
        duration_measured=5.0,
        transit_qc_expected_tmid=0.5,
        transit_qc_expected_tmid_unc=0.01,
        transit_qc_expected_rprs=0.10,
        transit_qc_expected_rprs_unc=0.05,
        transit_qc_use_deviation_from_expected_transit_in_qc=False,
        transit_qc_deviation_sigma_threshold=5.0,
    )

    summary = evaluate_transit_detection_qc(fit)

    assert summary["eebls_depth_snr"] == pytest.approx(7.25)
    assert fit.eebls_diagnostic_depth_snr == pytest.approx(7.25)


def test_fit_final_lightcurve_with_oot_baseline_detrending_preserves_expected_tmid_context(monkeypatch):
    run_count = {"value": 0, "duration_priors": []}

    def fake_run_nested(times, flux_values, flux_errors, airmass, prior, bounds, **kwargs):
        run_count["value"] += 1
        run_count["duration_priors"].append(kwargs.get("duration_prior"))
        local_times = np.asarray(times, dtype=float)
        model = np.ones(local_times.shape[0], dtype=float)
        model[1:-1] -= 0.01
        return types.SimpleNamespace(
            time=local_times,
            data=model.copy(),
            dataerr=np.full(local_times.shape[0], 0.001, dtype=float),
            model=model.copy(),
            residuals=np.zeros(local_times.shape[0], dtype=float),
            airmass=np.asarray(airmass, dtype=float),
            prior=dict(prior),
            parameters={"rprs": 0.1, "tmid": prior["tmid"], "inc": 89.0, "a2": 0.0, "per": prior["per"]},
            errors={"rprs": 0.01, "tmid": 0.001, "inc": 0.1, "a2": 0.01},
            bounds=dict(bounds),
            duration_expected=0.1,
            duration_measured=0.1,
        )

    monkeypatch.setattr(
        "exotic.exotic.run_nested_lightcurve_fit_with_rprs_posterior_retry",
        fake_run_nested,
    )
    monkeypatch.setattr("exotic.exotic.apply_plot_time_range", lambda fit, plot_time_range: fit)
    monkeypatch.setattr("exotic.exotic.apply_vertical_flux_normalization_bound", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        "exotic.exotic.build_final_fit_prefit_refinement_plan",
        lambda times, flux_values, flux_errors, airmass, prior, bounds, fit, **kwargs: {
            "applied": True,
            "note": "test prefit refinement",
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

    times = np.linspace(2460000.45, 2460000.55, 8)
    fit, _, _ = fit_final_lightcurve_with_oot_baseline_detrending(
        times,
        np.ones(times.shape[0], dtype=float),
        np.full(times.shape[0], 0.001, dtype=float),
        np.linspace(1.0, 1.1, times.shape[0]),
        {"rprs": 0.1, "tmid": 2460000.5, "inc": 89.0, "a2": 0.0, "per": 2.0},
        {"rprs": [0.0, 1.0], "tmid": [2460000.45, 2460000.55], "inc": [84.0, 90.0], "a2": [-3.0, 3.0]},
        detrend_on_outoftransit_baseline=False,
        expected_planet_dict={
            "midT": 2455000.0,
            "midTUnc": 0.0001,
            "pPer": 2.0,
            "pPerUnc": 0.001,
            "rprs": 0.1,
            "rprsUnc": 0.01,
            "aRs": 15.0,
            "aRsUnc": 0.1,
            "inc": 89.0,
            "incUnc": 0.1,
            "ecc": 0.0,
            "omega": 0.0,
        },
        expected_tmid_search_summary={
            "method": "ephemeris",
            "applied": True,
            "tmid": 2460000.5,
            "uncertainty": 0.002,
            "bounds": [2460000.45, 2460000.55],
            "duration": 0.1,
            "depth": np.nan,
            "depth_snr": np.nan,
            "note": "test propagated tmid",
        },
        eebls_search_summary={
            "method": "eebls",
            "applied": True,
            "tmid": 2460000.5,
            "bounds": [2460000.47, 2460000.53],
            "duration": 0.1,
            "depth": 0.01,
            "depth_snr": 6.5,
            "note": "test eebls",
        },
    )

    assert run_count["value"] == 2
    assert all(prior is not None and prior.get("applied") for prior in run_count["duration_priors"])
    assert fit.initial_tmid_search_tmid == pytest.approx(2460000.5)
    assert fit.transit_qc_expected_tmid == pytest.approx(2460000.5)
    assert fit.transit_qc_expected_tmid_unc == pytest.approx(0.002)
    assert fit.eebls_diagnostic_depth_snr == pytest.approx(6.5)


def test_fit_ranked_comparison_calibration_candidates_selects_highest_ktmf_success(monkeypatch):
    def fake_diagnostics(*args, **kwargs):
        return {"usable_point_count": 6}

    def fake_finalize(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times=None,
        **kwargs,
    ):
        comp_marker = int(np.nanmedian(cflux))
        residual_scale_map = {50: 0.05, 40: 0.02, 30: 0.03}
        delta_bic_map = {50: 8.0, 40: 18.0, 30: 12.0}
        ktmf_map = {50: 2.40, 40: 4.70, 30: 3.90}
        residual_scale = residual_scale_map[comp_marker]
        residuals = residual_scale * np.array([-1.0, 1.0, -1.0, 1.0, -1.0, 1.0], dtype=float)
        fit = types.SimpleNamespace(
            residuals=residuals,
            data=np.ones_like(residuals),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a0": 1.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.01, "a2": 0.01},
            transit_qc_delta_bic=delta_bic_map[comp_marker],
            transit_qc_ktmf_metric=ktmf_map[comp_marker],
        )
        return {
            "applied": True,
            "fit": fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": np.asarray(cflux, dtype=float),
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.diagnose_lightcurve_fit_inputs", fake_diagnostics)
    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.2, 6)
    ld = [0.1, 0.1, 0.1, 0.1]
    p_dict = {"midT": 0.5, "pPer": 1.0, "rprs": 0.1, "aRs": 10.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0}
    aper_data = {
        "target": np.full((6, 1, 1), 100.0, dtype=float),
        "comp1": np.full((6, 1, 1), 50.0, dtype=float),
        "comp2": np.full((6, 1, 1), 40.0, dtype=float),
        "comp3": np.full((6, 1, 1), 30.0, dtype=float),
    }
    comparison_calibration = {
        "method": "aperture",
        "a": 0,
        "an": 0,
        "comp_summaries": [
            {"label": "Comp 1", "position": (10.0, 10.0), "aggregate_score": 0.01, "coverage_count": 6, "coverage_total_frame_count": 6, "coverage_reference_count": 6.0, "coverage_min_required_count": 5, "coverage_rejected": False, "comp_index": 0},
            {"label": "Comp 2", "position": (20.0, 20.0), "aggregate_score": 0.02, "coverage_count": 6, "coverage_total_frame_count": 6, "coverage_reference_count": 6.0, "coverage_min_required_count": 5, "coverage_rejected": False, "comp_index": 1},
            {"label": "Comp 3", "position": (30.0, 30.0), "aggregate_score": 0.03, "coverage_count": 6, "coverage_total_frame_count": 6, "coverage_reference_count": 6.0, "coverage_min_required_count": 5, "coverage_rejected": False, "comp_index": 2},
        ],
    }

    result = fit_ranked_comparison_calibration_candidates(
        times,
        jd_times,
        airmass,
        ld,
        p_dict,
        comparison_calibration,
        psf_data={},
        aper_data=aper_data,
        target_psf_flux=np.full(6, 100.0, dtype=float),
    )

    assert len(result["attempts"]) == 3
    assert result["selection_metric"] == "ktmf_combined_quality"
    assert result["selected_result"]["comp_index"] == 1
    assert result["selected_result"]["rank"] == 1
    assert result["selected_result"]["selected"] is True
    assert result["selected_result"]["ktmf_metric"] == pytest.approx(4.70)
    assert "highest KTMF/projected-scatter" in result["selected_result"]["selection_reason"]
    assert result["attempts"][0]["selection_reason"].startswith(
        "not selected: full-resolution UltraNest model residual scatter"
    )


def test_select_preferred_comparison_attempt_rejects_noisy_high_ktmf_before_ranking():
    attempts = [
        {
            "label": "Comp 1",
            "rank": 0,
            "ktmf_metric": 4.8,
            "residual_scatter": 0.040,
            "eebls_snr": 3.0,
            "transit_delta_bic": 10.0,
        },
        {
            "label": "Comp 2",
            "rank": 1,
            "ktmf_metric": 3.6,
            "residual_scatter": 0.010,
            "eebls_snr": 2.8,
            "transit_delta_bic": 8.0,
        },
        {
            "label": "Comp 3",
            "rank": 2,
            "ktmf_metric": 3.8,
            "residual_scatter": 0.014,
            "eebls_snr": 2.6,
            "transit_delta_bic": 7.0,
        },
    ]

    selected, metric = select_preferred_comparison_attempt(attempts)

    assert metric == "ktmf_combined_quality"
    assert selected["label"] == "Comp 2"
    assert attempts[0]["scatter_gate_passed"] is False
    assert attempts[0]["scatter_gate_threshold"] == pytest.approx(0.015)
    assert selected["scatter_adjusted_ktmf_metric"] == pytest.approx(3.6)
    assert attempts[2]["scatter_adjusted_ktmf_metric"] == pytest.approx(3.8 * 0.010 / 0.014)


def test_select_preferred_comparison_attempt_uses_ktmf_and_projected_selection_scatter_only():
    attempts = [
        {
            "label": "Comp 1",
            "rank": 0,
            "ktmf_metric": 2.59,
            "selection_scatter": 0.022659,
            "target_comp_scatter": 0.005693,
            "aggregate_score": 0.010,
            "eebls_snr": 3.27,
            "transit_delta_bic": 2.71,
        },
        {
            "label": "Comp 2",
            "rank": 1,
            "ktmf_metric": 3.23,
            "selection_scatter": 0.027437,
            "target_comp_scatter": 0.003000,
            "aggregate_score": 0.001,
            "eebls_snr": 3.53,
            "transit_delta_bic": 2.28,
        },
        {
            "label": "Comp 9",
            "rank": 2,
            "ktmf_metric": 3.25,
            "selection_scatter": 0.023904,
            "target_comp_scatter": 0.030258,
            "aggregate_score": 0.100,
            "eebls_snr": 1.61,
            "transit_delta_bic": -8.36,
        },
    ]

    selected, metric = select_preferred_comparison_attempt(attempts)

    assert metric == "ktmf_combined_quality"
    assert selected["label"] == "Comp 9"
    assert attempts[0]["combined_quality_ktmf_metric"] == pytest.approx(2.59 / 2.2659)
    assert attempts[1]["combined_quality_ktmf_metric"] == pytest.approx(3.23 / 2.7437)
    assert attempts[2]["combined_quality_ktmf_metric"] == pytest.approx(3.25 / 2.3904)
    assert selected["combined_quality_ktmf_metric"] > attempts[1]["combined_quality_ktmf_metric"]
    assert selected["combined_quality_ktmf_metric"] > attempts[0]["combined_quality_ktmf_metric"]


def test_target_comp_flux_scatter_measures_normalized_target_reference_ratio():
    comp_flux = np.full(8, 100.0, dtype=float)
    ratio = np.array([1.00, 1.01, 0.99, 1.02, 0.98, 1.00, 1.01, 0.99], dtype=float)
    target_flux = comp_flux * ratio

    scatter = target_comp_flux_scatter(target_flux, comp_flux, min_points=5)

    assert scatter == pytest.approx(0.014826, rel=1.0e-3)


def test_fitted_lightcurve_scatter_on_dataset_projects_fit_to_full_flux(monkeypatch):
    def fake_transit(times, parameters):
        return np.ones_like(np.asarray(times, dtype=float))

    monkeypatch.setattr("exotic.exotic.transit", fake_transit)
    fit = types.SimpleNamespace(
        parameters={"a0": 1.0, "a2": 0.0},
        airmass_reference=1.0,
    )
    times = np.arange(8, dtype=float)
    flux_values = np.array([1.0, 1.01, 0.99, 1.02, 0.98, 1.0, 1.01, 0.99], dtype=float)
    airmass = np.ones_like(times)

    scatter = fitted_lightcurve_scatter_on_dataset(fit, times, flux_values, airmass)

    assert scatter == pytest.approx(np.std(flux_values - 1.0) / np.median(flux_values))


def test_fit_ranked_comparison_calibration_candidates_extends_only_selected_final_fit(monkeypatch):
    monkeypatch.setenv("EXOTIC_ULTRANEST_MIN_NUM_LIVE_POINTS", "200")
    monkeypatch.setenv("EXOTIC_SPARSE_POSTERIOR_LIVE_POINT_RETRY", "1")

    def fake_diagnostics(*args, **kwargs):
        return {"usable_point_count": 6}

    created_fits = {}

    class RetainedSamplerFit:
        def __init__(self, comp_marker, ktmf_metric):
            self.comp_marker = comp_marker
            self.extension_calls = []
            self.cleared = False
            self.max_ncalls = 1000
            self.time = np.linspace(0.0, 0.05, 6)
            self.data = np.ones(6, dtype=float)
            self.residuals = np.full(6, 0.01, dtype=float)
            self.parameters = {
                "tmid": 0.5,
                "rprs": 0.1,
                "inc": 89.0,
                "ars": 10.0,
                "a0": 1.0,
                "a2": 0.0,
            }
            self.errors = {
                "tmid": 0.001,
                "rprs": 0.001,
                "inc": 0.1,
                "ars": 0.1,
                "a0": 0.01,
                "a2": 0.01,
            }
            self.bounds = {
                "rprs": [0.08, 0.12],
                "tmid": [0.49, 0.51],
                "ars": [9.0, 11.0],
                "inc": [85.0, 90.0],
                "a2": [-3.0, 3.0],
            }
            self.transit_qc_delta_bic = 12.0 + comp_marker / 100.0
            self.transit_qc_ktmf_metric = ktmf_metric

        def get_parameter_posterior_samples(self, key):
            ranges = {
                "rprs": (0.09, 0.11),
                "tmid": (0.499, 0.501),
                "ars": (9.5, 10.5),
            }
            low, high = ranges[key]
            return np.linspace(low, high, 1500)

        def extend_ultranest_fit(self, min_num_live_points=None, max_ncalls=None):
            self.extension_calls.append({
                "min_num_live_points": min_num_live_points,
                "max_ncalls": max_ncalls,
                "bounds": self.bounds.copy(),
            })
            return True

        def clear_ultranest_resume_state(self):
            self.cleared = True

    def fake_finalize(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times=None,
        **kwargs,
    ):
        comp_marker = int(np.nanmedian(cflux))
        ktmf_map = {50: 2.40, 40: 4.70, 30: 3.90}
        fit = RetainedSamplerFit(comp_marker, ktmf_map[comp_marker])
        created_fits[comp_marker] = fit
        return {
            "applied": True,
            "fit": fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": np.asarray(cflux, dtype=float),
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.diagnose_lightcurve_fit_inputs", fake_diagnostics)
    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.2, 6)
    aper_data = {
        "target": np.full((6, 1, 1), 100.0, dtype=float),
        "comp1": np.full((6, 1, 1), 50.0, dtype=float),
        "comp2": np.full((6, 1, 1), 40.0, dtype=float),
        "comp3": np.full((6, 1, 1), 30.0, dtype=float),
    }
    comparison_calibration = {
        "method": "aperture",
        "a": 0,
        "an": 0,
        "comp_summaries": [
            {"label": "Comp 1", "aggregate_score": 0.01, "coverage_rejected": False, "comp_index": 0},
            {"label": "Comp 2", "aggregate_score": 0.02, "coverage_rejected": False, "comp_index": 1},
            {"label": "Comp 3", "aggregate_score": 0.03, "coverage_rejected": False, "comp_index": 2},
        ],
    }

    result = fit_ranked_comparison_calibration_candidates(
        times,
        jd_times,
        airmass,
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={"midT": 0.5, "pPer": 1.0, "rprs": 0.1, "aRs": 10.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0},
        comparison_calibration=comparison_calibration,
        psf_data={},
        aper_data=aper_data,
        target_psf_flux=np.full(6, 100.0, dtype=float),
    )

    assert result["selected_result"]["comp_index"] == 1
    assert created_fits[40].extension_calls == [{
        "min_num_live_points": 1200,
        "max_ncalls": 6000,
        "bounds": created_fits[40].bounds,
    }]
    assert created_fits[50].extension_calls == []
    assert created_fits[30].extension_calls == []
    assert created_fits[40].cleared is True
    assert created_fits[50].cleared is True
    assert created_fits[30].cleared is True
    assert "selected comparison-star final" in created_fits[40].sparse_posterior_live_point_extension_note


def test_fit_ranked_comparison_calibration_candidates_stops_at_first_qc_pass_by_default(monkeypatch):
    def fake_diagnostics(*args, **kwargs):
        return {"usable_point_count": 6}

    call_markers = []

    def fake_finalize(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times=None,
        **kwargs,
    ):
        comp_marker = int(np.nanmedian(cflux))
        call_markers.append(comp_marker)
        status_map = {50: "marginal", 40: "pass", 30: "pass"}
        ktmf_map = {50: 4.90, 40: 3.20, 30: 5.00}
        status = status_map[comp_marker]
        ktmf_metric = ktmf_map[comp_marker]
        fit = types.SimpleNamespace(
            residuals=np.full(6, 0.01, dtype=float),
            data=np.ones(6, dtype=float),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a0": 1.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.01, "a2": 0.01},
            transit_qc={
                "status": status,
                "summary": "ok",
                "delta_bic": 12.0 + comp_marker / 100.0,
                "ktmf_metric": ktmf_metric,
                "ktmf_contributions": [],
            },
            transit_qc_status=status,
            transit_qc_ktmf_metric=ktmf_metric,
            transit_qc_delta_bic=12.0 + comp_marker / 100.0,
        )
        return {
            "applied": True,
            "fit": fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": np.asarray(cflux, dtype=float),
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.diagnose_lightcurve_fit_inputs", fake_diagnostics)
    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.2, 6)
    aper_data = {
        "target": np.full((6, 1, 1), 100.0, dtype=float),
        "comp1": np.full((6, 1, 1), 50.0, dtype=float),
        "comp2": np.full((6, 1, 1), 40.0, dtype=float),
        "comp3": np.full((6, 1, 1), 30.0, dtype=float),
    }
    comparison_calibration = {
        "method": "aperture",
        "a": 0,
        "an": 0,
        "comp_summaries": [
            {"label": "Comp 1", "aggregate_score": 0.01, "coverage_rejected": False, "comp_index": 0},
            {"label": "Comp 2", "aggregate_score": 0.02, "coverage_rejected": False, "comp_index": 1},
            {"label": "Comp 3", "aggregate_score": 0.03, "coverage_rejected": False, "comp_index": 2},
        ],
    }

    result = fit_ranked_comparison_calibration_candidates(
        times,
        jd_times,
        airmass,
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={"midT": 0.5, "pPer": 1.0, "rprs": 0.1, "aRs": 10.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0},
        comparison_calibration=comparison_calibration,
        psf_data={},
        aper_data=aper_data,
        target_psf_flux=np.full(6, 100.0, dtype=float),
    )

    assert call_markers == [50, 40]
    assert len(result["attempts"]) == 2
    assert result["stopped_after_first_qc_pass"] is True
    assert result["selection_metric"] == "first_qc_pass"
    assert result["selected_result"]["comp_index"] == 1
    assert result["selected_result"]["search_stopped_after_qc_pass"] is True
    assert "first completed comparison-star candidate" in result["selected_result"]["selection_reason"]


def test_fit_ranked_comparison_calibration_candidates_stops_at_promising_partial_marginal(monkeypatch):
    def fake_diagnostics(*args, **kwargs):
        return {"usable_point_count": 6}

    monkeypatch.setattr(
        "exotic.exotic.build_comparison_candidate_preflight",
        lambda *args, **kwargs: {
            "prepared_series": None,
            "coverage_priority": 2,
            "scout": {"score": np.nan},
        },
    )

    call_markers = []

    def fake_finalize(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times=None,
        **kwargs,
    ):
        comp_marker = int(np.nanmedian(cflux))
        call_markers.append(comp_marker)
        status = "marginal" if comp_marker == 50 else "pass"
        fit = types.SimpleNamespace(
            residuals=np.full(6, 0.01, dtype=float),
            data=np.ones(6, dtype=float),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a0": 1.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.01, "a2": 0.01},
            transit_qc={"status": status, "summary": "ok"},
            transit_qc_status=status,
            transit_qc_ktmf_metric=3.5,
            transit_qc_delta_bic=15.1,
        )
        return {
            "applied": True,
            "fit": fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": np.asarray(cflux, dtype=float),
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.diagnose_lightcurve_fit_inputs", fake_diagnostics)
    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.2, 6)
    aper_data = {
        "target": np.full((6, 1, 1), 100.0, dtype=float),
        "comp1": np.full((6, 1, 1), 50.0, dtype=float),
        "comp2": np.full((6, 1, 1), 40.0, dtype=float),
    }
    comparison_calibration = {
        "method": "aperture",
        "a": 0,
        "an": 0,
        "comp_summaries": [
            {"label": "Comp 1", "aggregate_score": 0.01, "coverage_rejected": False, "comp_index": 0},
            {"label": "Comp 2", "aggregate_score": 0.02, "coverage_rejected": False, "comp_index": 1},
        ],
    }

    result = fit_ranked_comparison_calibration_candidates(
        times,
        jd_times,
        airmass,
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={"midT": 0.5, "pPer": 1.0, "rprs": 0.1, "aRs": 10.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0},
        comparison_calibration=comparison_calibration,
        psf_data={},
        aper_data=aper_data,
        target_psf_flux=np.full(6, 100.0, dtype=float),
    )

    assert call_markers == [50]
    assert result["selection_metric"] == "promising_partial"
    assert result["selected_result"]["search_stopped_after_promising_partial"] is True
    assert result["stopped_after_promising_partial"] is True


def test_fit_ranked_comparison_calibration_candidates_can_evaluate_all_qc_passes_when_exit_disabled(monkeypatch):
    def fake_diagnostics(*args, **kwargs):
        return {"usable_point_count": 6}

    call_markers = []

    def fake_finalize(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times=None,
        **kwargs,
    ):
        comp_marker = int(np.nanmedian(cflux))
        call_markers.append(comp_marker)
        ktmf_metric = {50: 3.10, 40: 4.00, 30: 4.80}[comp_marker]
        fit = types.SimpleNamespace(
            residuals=np.full(6, 0.01, dtype=float),
            data=np.ones(6, dtype=float),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a0": 1.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.01, "a2": 0.01},
            transit_qc={
                "status": "pass",
                "summary": "ok",
                "delta_bic": 12.0 + comp_marker / 100.0,
                "ktmf_metric": ktmf_metric,
                "ktmf_contributions": [],
            },
            transit_qc_status="pass",
            transit_qc_ktmf_metric=ktmf_metric,
            transit_qc_delta_bic=12.0 + comp_marker / 100.0,
        )
        return {
            "applied": True,
            "fit": fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": np.asarray(cflux, dtype=float),
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.diagnose_lightcurve_fit_inputs", fake_diagnostics)
    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.2, 6)
    aper_data = {
        "target": np.full((6, 1, 1), 100.0, dtype=float),
        "comp1": np.full((6, 1, 1), 50.0, dtype=float),
        "comp2": np.full((6, 1, 1), 40.0, dtype=float),
        "comp3": np.full((6, 1, 1), 30.0, dtype=float),
    }
    comparison_calibration = {
        "method": "aperture",
        "a": 0,
        "an": 0,
        "comp_summaries": [
            {"label": "Comp 1", "aggregate_score": 0.01, "coverage_rejected": False, "comp_index": 0},
            {"label": "Comp 2", "aggregate_score": 0.02, "coverage_rejected": False, "comp_index": 1},
            {"label": "Comp 3", "aggregate_score": 0.03, "coverage_rejected": False, "comp_index": 2},
        ],
    }

    result = fit_ranked_comparison_calibration_candidates(
        times,
        jd_times,
        airmass,
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={"midT": 0.5, "pPer": 1.0, "rprs": 0.1, "aRs": 10.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0},
        comparison_calibration=comparison_calibration,
        psf_data={},
        aper_data=aper_data,
        target_psf_flux=np.full(6, 100.0, dtype=float),
        exit_at_first_qc_pass_solution=False,
    )

    assert call_markers == [50, 40, 30]
    assert result["stopped_after_first_qc_pass"] is False
    assert result["selection_metric"] == "ktmf_combined_quality"
    assert result["selected_result"]["comp_index"] == 2
    assert result["selected_result"]["ktmf_metric"] == pytest.approx(4.80)


def test_ranked_comparison_calibration_summaries_skip_suitability_outliers():
    ranked = ranked_comparison_calibration_summaries(
        {
            "comp_summaries": [
                {"comp_index": 0, "aggregate_score": 0.139668, "coverage_rejected": False, "suitability_outlier_rejected": True},
                {"comp_index": 1, "aggregate_score": 0.051809, "coverage_rejected": False, "suitability_outlier_rejected": True},
                {"comp_index": 2, "aggregate_score": 0.024802, "coverage_rejected": False, "suitability_outlier_rejected": False},
                {"comp_index": 3, "aggregate_score": 0.027680, "coverage_rejected": False, "suitability_outlier_rejected": False},
            ]
        }
    )

    assert [summary["comp_index"] for summary in ranked] == [2, 3]


def test_comparison_preflight_ranking_prioritizes_full_coverage_then_scout_score():
    plans = [
        {
            "field_rank": 0,
            "summary": {"comp_index": 7, "aggregate_score": 0.002175, "label": "Comp 8"},
            "preflight": {"coverage_priority": 2, "scout": {"score": 0.42}},
        },
        {
            "field_rank": 1,
            "summary": {"comp_index": 0, "aggregate_score": 0.002331, "label": "Comp 1"},
            "preflight": {"coverage_priority": 2, "scout": {"score": 0.91}},
        },
        {
            "field_rank": 4,
            "summary": {"comp_index": 2, "aggregate_score": 0.002804, "label": "Comp 3"},
            "preflight": {"coverage_priority": 0, "scout": {"score": 0.25}},
        },
    ]

    ranked = rank_comparison_candidate_preflight_plans(plans)

    assert [plan["summary"]["comp_index"] for plan in ranked] == [2, 0, 7]


def test_promising_partial_comparison_attempt_can_stop_candidate_search():
    attempt = {
        "fit": object(),
        "full_reduction_applied": True,
        "rejected_by_transit_qc": False,
        "transit_qc_status": "marginal",
        "preflight_coverage_priority": 2,
        "ktmf_metric": 3.50,
        "transit_delta_bic": 15.09,
    }

    assert should_stop_after_promising_partial_comparison_attempt(attempt) is True

    attempt["preflight_coverage_priority"] = 4
    assert should_stop_after_promising_partial_comparison_attempt(attempt) is False


def test_fit_ranked_comparison_calibration_candidates_applies_field_image_clip(monkeypatch):
    observed_lengths = []

    def fake_diagnostics(times, *args, **kwargs):
        observed_lengths.append(len(times))
        return {"usable_point_count": len(times)}

    def fake_finalize(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times=None,
        **kwargs,
    ):
        observed_lengths.append(len(times))
        fit = types.SimpleNamespace(
            residuals=np.full(len(times), 0.01, dtype=float),
            data=np.ones(len(times), dtype=float),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a0": 1.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.01, "a2": 0.01},
            transit_qc={"status": "pass", "summary": "ok", "ktmf_metric": 4.2},
            transit_qc_status="pass",
            transit_qc_summary="ok",
            transit_qc_ktmf_metric=4.2,
            transit_qc_delta_bic=16.0,
            frame_filter_diagnostics=[{"stage": "Comparison-field image clip", "dropped_point_count": 2}],
        )
        return {
            "applied": True,
            "fit": fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": np.asarray(cflux, dtype=float),
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.diagnose_lightcurve_fit_inputs", fake_diagnostics)
    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.2, 6)
    comparison_calibration = {
        "method": "aperture",
        "method_label": "Aperture photometry (aper=5.00px, annulus=12.00px)",
        "a": 0,
        "an": 0,
        "aper": 5.0,
        "annulus": 12.0,
        "field_image_keep_mask": np.array([True, False, True, True, False, True], dtype=bool),
        "image_outlier_sigma": 4.25,
        "image_outlier_required_valid_pairs": 2,
        "comp_summaries": [
            {"label": "Comp 1", "position": (10.0, 10.0), "aggregate_score": 0.01, "coverage_count": 6, "coverage_total_frame_count": 6, "coverage_reference_count": 6.0, "coverage_min_required_count": 5, "coverage_rejected": False, "comp_index": 0},
        ],
    }
    aper_data = {
        "target": np.full((6, 1, 1), 100.0, dtype=float),
        "comp1": np.full((6, 1, 1), 50.0, dtype=float),
    }

    result = fit_ranked_comparison_calibration_candidates(
        times,
        jd_times,
        airmass,
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={"midT": 0.5, "pPer": 1.0, "rprs": 0.1, "aRs": 10.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0},
        comparison_calibration=comparison_calibration,
        psf_data={},
        aper_data=aper_data,
        target_psf_flux=np.full(6, 100.0, dtype=float),
    )

    assert observed_lengths == [4, 4]
    diagnostic = result["attempts"][0]["fit"].frame_filter_diagnostics[0]
    assert diagnostic["stage"] == "Comparison-field image clip"
    assert diagnostic["dropped_point_count"] == 2
    assert result["attempts"][0]["fit_point_count"] == 4


def test_fit_ranked_comparison_calibration_candidates_masks_target_psf_shape(monkeypatch):
    observed_lengths = []

    def fake_diagnostics(times, *args, **kwargs):
        observed_lengths.append(("diagnostics", len(times)))
        return {"usable_point_count": len(times), "failure_reason": None}

    def fake_preflight(*args, **kwargs):
        return {"coverage_priority": 1, "prepared_series": None}

    def fake_finalize(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times=None,
        **kwargs,
    ):
        observed_lengths.append(("finalize", len(times)))
        fit = types.SimpleNamespace(
            residuals=np.full(len(times), 0.01, dtype=float),
            data=np.ones(len(times), dtype=float),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a0": 1.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.01, "a2": 0.01},
            transit_qc={"status": "pass", "summary": "ok", "ktmf_metric": 4.2},
            transit_qc_status="pass",
            transit_qc_summary="ok",
            transit_qc_ktmf_metric=4.2,
            transit_qc_delta_bic=16.0,
            frame_filter_diagnostics=[],
        )
        return {
            "applied": True,
            "fit": fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": np.asarray(cflux, dtype=float),
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.diagnose_lightcurve_fit_inputs", fake_diagnostics)
    monkeypatch.setattr("exotic.exotic.build_comparison_candidate_preflight", fake_preflight)
    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    frame_count = 30
    times = np.linspace(0.0, 0.2, frame_count)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.2, frame_count)

    psf_rows = np.zeros((frame_count, 7), dtype=float)
    psf_rows[:, 0] = 10.0
    psf_rows[:, 1] = 20.0
    psf_rows[:, 2] = 100.0
    psf_rows[:, 3] = 1.0
    psf_rows[:, 4] = 1.0
    psf_data = {
        "target": psf_rows.copy(),
        "comp1": psf_rows.copy(),
    }
    psf_data["comp1"][:, 2] = 120.0
    psf_data["target"][12, 3:5] = 6.5
    psf_flux_data = {
        "target": psf_data["target"].copy(),
        "comp1": psf_data["comp1"].copy(),
    }
    psf_flux_data["comp1"][:, 2] = 240.0

    target_psf_flux = 2 * np.pi * psf_data["target"][:, 2] * psf_data["target"][:, 3] * psf_data["target"][:, 4]
    comparison_calibration = {
        "method": "psf",
        "method_label": "PSF photometry",
        "a": None,
        "an": None,
        "aper": 0.0,
        "annulus": 15.0,
        "comp_summaries": [
            {
                "label": "Comp 1",
                "position": (10.0, 10.0),
                "aggregate_score": 0.01,
                "coverage_count": frame_count,
                "coverage_total_frame_count": frame_count,
                "coverage_reference_count": float(frame_count),
                "coverage_min_required_count": 5,
                "coverage_rejected": False,
                "comp_index": 0,
                "key": "comp1",
            },
        ],
    }

    result = fit_ranked_comparison_calibration_candidates(
        times,
        jd_times,
        airmass,
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={"midT": 0.5, "pPer": 1.0, "rprs": 0.1, "aRs": 10.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0},
        comparison_calibration=comparison_calibration,
        psf_data=psf_data,
        aper_data=None,
        target_psf_flux=target_psf_flux,
        psf_flux_data=psf_flux_data,
    )

    assert observed_lengths == [("diagnostics", frame_count - 1), ("finalize", frame_count - 1)]
    assert result["selected_result"]["fit_point_count"] == frame_count - 1
    assert np.nanmax(result["selected_result"]["tflux_fit"]) < 1000.0
    assert np.nanmedian(result["selected_result"]["cflux_fit"]) == pytest.approx(2.0 * np.pi * 240.0)


def test_fit_ranked_comparison_calibration_candidates_applies_candidate_intercomparison_clip(monkeypatch):
    observed_lengths = []

    def fake_diagnostics(times, *args, **kwargs):
        observed_lengths.append(("diagnostics", len(times)))
        return {"usable_point_count": len(times), "failure_reason": None}

    def fake_preflight(*args, **kwargs):
        return {"coverage_priority": 1, "prepared_series": None}

    def fake_finalize(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times=None,
        **kwargs,
    ):
        observed_lengths.append(("finalize", len(times)))
        fit = types.SimpleNamespace(
            residuals=np.full(len(times), 0.01, dtype=float),
            data=np.ones(len(times), dtype=float),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a0": 1.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.01, "a2": 0.01},
            transit_qc={"status": "pass", "summary": "ok", "ktmf_metric": 4.2},
            transit_qc_status="pass",
            transit_qc_summary="ok",
            transit_qc_ktmf_metric=4.2,
            transit_qc_delta_bic=16.0,
            frame_filter_diagnostics=[],
        )
        return {
            "applied": True,
            "fit": fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": np.asarray(cflux, dtype=float),
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.diagnose_lightcurve_fit_inputs", fake_diagnostics)
    monkeypatch.setattr("exotic.exotic.build_comparison_candidate_preflight", fake_preflight)
    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.2, 6)
    comparison_calibration = {
        "method": "aperture",
        "method_label": "Aperture photometry (aper=5.00px, annulus=12.00px)",
        "a": 0,
        "an": 0,
        "aper": 5.0,
        "annulus": 12.0,
        "field_image_keep_mask": np.ones(6, dtype=bool),
        "comp_summaries": [
            {
                "label": "Comp 1",
                "position": (10.0, 10.0),
                "aggregate_score": 0.01,
                "coverage_count": 6,
                "coverage_total_frame_count": 6,
                "coverage_reference_count": 6.0,
                "coverage_min_required_count": 5,
                "coverage_rejected": False,
                "suitability_outlier_rejected": False,
                "comp_index": 0,
                "ensemble_frame_keep_mask": np.array([True, True, False, True, True, True], dtype=bool),
                "ensemble_frame_required_valid_pairs": 2,
                "ensemble_frame_sigma": 4.25,
            },
        ],
    }
    aper_data = {
        "target": np.full((6, 1, 1), 100.0, dtype=float),
        "comp1": np.full((6, 1, 1), 50.0, dtype=float),
    }

    result = fit_ranked_comparison_calibration_candidates(
        times,
        jd_times,
        airmass,
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={"midT": 0.5, "pPer": 1.0, "rprs": 0.1, "aRs": 10.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0},
        comparison_calibration=comparison_calibration,
        psf_data={},
        aper_data=aper_data,
        target_psf_flux=np.full(6, 100.0, dtype=float),
    )

    assert observed_lengths == [("diagnostics", 5), ("finalize", 5)]
    diagnostic = result["attempts"][0]["fit"].frame_filter_diagnostics[0]
    assert diagnostic["stage"] == "Comparison-candidate intercomparison clip"
    assert diagnostic["dropped_point_count"] == 1
    assert result["attempts"][0]["fit_point_count"] == 5


def test_fit_ranked_comparison_calibration_candidates_saves_outputs_for_completed_candidates(
    monkeypatch, tmp_path
):
    def fake_diagnostics(*args, **kwargs):
        return {"usable_point_count": 6}

    call_markers = []

    def fake_finalize(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times=None,
        **kwargs,
    ):
        comp_marker = int(np.nanmedian(cflux))
        call_markers.append(comp_marker)
        fit = types.SimpleNamespace(
            residuals=np.full(6, 0.01, dtype=float),
            data=np.ones(6, dtype=float),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a0": 1.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.01, "a2": 0.01},
            transit_qc_ktmf_metric={50: 3.2, 40: 4.4}[comp_marker],
            transit_qc_delta_bic=12.0 + comp_marker / 100.0,
        )
        return {
            "applied": True,
            "fit": fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": np.asarray(cflux, dtype=float),
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.diagnose_lightcurve_fit_inputs", fake_diagnostics)
    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    saved_dirs = []

    def fake_save(save_dir, provisional_fit, final_fit, p_dict, observation_date, comp_index, **kwargs):
        candidate_dir = Path(save_dir) / f"comp{comp_index + 1}"
        candidate_dir.mkdir(parents=True, exist_ok=True)
        saved_dirs.append(candidate_dir)
        return candidate_dir

    monkeypatch.setattr(
        "exotic.exotic.save_comparison_candidate_full_reduction_outputs",
        fake_save,
    )

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.2, 6)
    aper_data = {
        "target": np.full((6, 1, 1), 100.0, dtype=float),
        "comp1": np.full((6, 1, 1), 50.0, dtype=float),
        "comp2": np.full((6, 1, 1), 40.0, dtype=float),
    }
    comparison_calibration = {
        "method": "aperture",
        "a": 0,
        "an": 0,
        "comp_summaries": [
            {"label": "Comp 1", "aggregate_score": 0.01, "coverage_rejected": False, "comp_index": 0},
            {"label": "Comp 2", "aggregate_score": 0.02, "coverage_rejected": False, "comp_index": 1},
        ],
    }

    result = fit_ranked_comparison_calibration_candidates(
        times,
        jd_times,
        airmass,
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={"midT": 0.5, "pPer": 1.0, "rprs": 0.1, "aRs": 10.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0},
        comparison_calibration=comparison_calibration,
        psf_data={},
        aper_data=aper_data,
        target_psf_flux=np.full(6, 100.0, dtype=float),
        save_dir=tmp_path,
        planet_name="HAT-P-32 b",
        observation_date="2026-04-28",
    )

    assert call_markers == [50, 40]
    assert len(result["attempts"]) == 2
    assert result["selection_metric"] == "ktmf_combined_quality"
    assert result["selected_result"]["comp_index"] == 1
    assert [attempt["final_output_dir"] for attempt in result["attempts"]] == [
        str(tmp_path / "comp1"),
        str(tmp_path / "comp2"),
    ]
    assert saved_dirs == [tmp_path / "comp1", tmp_path / "comp2"]


def test_fit_ranked_comparison_calibration_candidates_can_prefer_highest_eebls_snr(monkeypatch):
    def fake_diagnostics(*args, **kwargs):
        return {"usable_point_count": 6}

    def fake_finalize(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times=None,
        **kwargs,
    ):
        comp_marker = int(np.nanmedian(cflux))
        if comp_marker == 50:
            residual_level = 0.01
            eebls_snr = 4.0
        else:
            residual_level = 0.012
            eebls_snr = 7.5

        residuals = residual_level * np.array([-1.0, 1.0, -1.0, 1.0, -1.0, 1.0], dtype=float)
        fit = types.SimpleNamespace(
            residuals=residuals,
            data=np.ones_like(residuals),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a0": 1.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.01, "a2": 0.01},
            eebls_diagnostic_depth_snr=eebls_snr,
            transit_qc_delta_bic=(10.0 if comp_marker == 50 else 12.0),
        )
        return {
            "applied": True,
            "fit": fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": np.asarray(cflux, dtype=float),
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.diagnose_lightcurve_fit_inputs", fake_diagnostics)
    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.2, 6)
    ld = [0.1, 0.1, 0.1, 0.1]
    p_dict = {"midT": 0.5, "pPer": 1.0, "rprs": 0.1, "aRs": 10.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0}
    aper_data = {
        "target": np.full((6, 1, 1), 100.0, dtype=float),
        "comp1": np.full((6, 1, 1), 50.0, dtype=float),
        "comp2": np.full((6, 1, 1), 40.0, dtype=float),
    }
    comparison_calibration = {
        "method": "aperture",
        "a": 0,
        "an": 0,
        "comp_summaries": [
            {"label": "Comp 1", "position": (10.0, 10.0), "aggregate_score": 0.01, "coverage_count": 6, "coverage_total_frame_count": 6, "coverage_reference_count": 6.0, "coverage_min_required_count": 5, "coverage_rejected": False, "comp_index": 0},
            {"label": "Comp 2", "position": (20.0, 20.0), "aggregate_score": 0.02, "coverage_count": 6, "coverage_total_frame_count": 6, "coverage_reference_count": 6.0, "coverage_min_required_count": 5, "coverage_rejected": False, "comp_index": 1},
        ],
    }

    result = fit_ranked_comparison_calibration_candidates(
        times,
        jd_times,
        airmass,
        ld,
        p_dict,
        comparison_calibration,
        psf_data={},
        aper_data=aper_data,
        target_psf_flux=np.full(6, 100.0, dtype=float),
        pick_comparison_by_eebls_snr=True,
    )

    assert result["selection_metric"] == "eebls_snr"
    assert result["selected_result"]["comp_index"] == 1
    assert result["selected_result"]["eebls_snr"] == pytest.approx(7.5)
    assert "highest selection-pass EEBLS SNR" in result["selected_result"]["selection_reason"]
    assert result["attempts"][0]["selection_reason"].startswith("not selected: selection-pass EEBLS SNR")


def test_fit_ranked_comparison_calibration_candidates_logs_per_comp_run_reporting(monkeypatch):
    logged = []

    def fake_diagnostics(*args, **kwargs):
        return {"usable_point_count": 6}

    final_fit = types.SimpleNamespace(
        residuals=np.full(6, 0.01, dtype=float),
        data=np.ones(6, dtype=float),
        parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a0": 1.0, "a2": 0.0},
        errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.01, "a2": 0.01},
        ns_type="ultranest",
        transit_qc={
            "status": "pass",
            "summary": "Transit model strongly preferred over flat/null model.",
            "delta_bic": 18.4,
            "residual_scatter": 0.0035,
            "ktmf_metric": 4.6,
            "ktmf_contributions": [],
        },
        transit_qc_status="pass",
        transit_qc_summary="Transit model strongly preferred over flat/null model.",
        transit_qc_delta_bic=18.4,
        transit_qc_residual_scatter=0.0035,
        transit_qc_ktmf_metric=4.6,
        transit_qc_ktmf_contributions=[],
        rprs_posterior_refit_applied=True,
        rprs_posterior_refit_count=1,
        rprs_posterior_refit_note="Applied 1 automatic Rp/R* posterior range refit(s).",
        prefit_refinement_applied=True,
        prefit_refinement_note="Applied a focused final-fit prefit refinement window.",
        oot_baseline_detrending_applied=False,
        oot_baseline_detrending_note="Skipped; need out-of-transit coverage on both sides of transit to fit a linear baseline.",
    )

    def fake_finalize(times, tflux, cflux, airmass, ld, p_dict, jd_times=None, **kwargs):
        return {
            "applied": True,
            "fit": final_fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": np.asarray(cflux, dtype=float),
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "completed the full comparison-candidate reduction.",
        }

    monkeypatch.setattr("exotic.exotic.log_info", lambda message, warn=False, error=False: logged.append(message))
    monkeypatch.setattr("exotic.exotic.diagnose_lightcurve_fit_inputs", fake_diagnostics)
    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.2, 6)
    aper_data = {
        "target": np.full((6, 1, 1), 100.0, dtype=float),
        "comp1": np.full((6, 1, 1), 50.0, dtype=float),
    }
    comparison_calibration = {
        "method": "aperture",
        "method_label": "Aperture photometry (aper=5.00px, annulus=12.00px)",
        "a": 0,
        "an": 0,
        "comp_summaries": [
            {
                "label": "Comp 1",
                "position": (10.0, 10.0),
                "aggregate_score": 0.01,
                "coverage_count": 6,
                "coverage_total_frame_count": 6,
                "coverage_reference_count": 6.0,
                "coverage_min_required_count": 5,
                "coverage_rejected": False,
                "comp_index": 0,
            },
        ],
    }

    fit_ranked_comparison_calibration_candidates(
        times,
        jd_times,
        airmass,
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={"midT": 0.5, "pPer": 1.0, "rprs": 0.1, "aRs": 10.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0},
        comparison_calibration=comparison_calibration,
        psf_data={},
        aper_data=aper_data,
        target_psf_flux=np.full(6, 100.0, dtype=float),
    )

    assert any("Starting comparison-star target-fit evaluation for Comp 1" in message for message in logged)
    assert any("Preparing comparison-candidate light curve for the full reduction." in message for message in logged)
    assert any("Full reduction starting. Optional out-of-transit baseline detrending is enabled." in message for message in logged)
    assert any("Completed comparison-star target-fit evaluation for Comp 1" in message and "transit_qc=PASS" in message for message in logged)
    assert any("Rp/R* posterior retry note: Applied 1 automatic Rp/R* posterior range refit(s)." in message for message in logged)
    assert any("OOT baseline detrending note: Skipped; need out-of-transit coverage on both sides of transit to fit a linear baseline." in message for message in logged)


def test_evaluate_lightcurve_candidate_requests_nested_fit(monkeypatch):
    def fake_diagnostics(*args, **kwargs):
        return {"usable_point_count": 6}

    def fake_fit_lightcurve(
        times,
        tflux,
        cflux,
        airmass,
        ld,
        p_dict,
        jd_times=None,
        **kwargs,
    ):
        assert kwargs.get("final_fit_mode") == "ns"
        fit = types.SimpleNamespace(
            residuals=np.full(6, 0.01, dtype=float),
            data=np.ones(6, dtype=float),
            parameters={"tmid": 0.5, "rprs": 0.1, "inc": 89.0, "a0": 1.0, "a2": 0.0},
            errors={"tmid": 0.001, "rprs": 0.001, "inc": 0.1, "a0": 0.01, "a2": 0.01},
            ns_type="ultranest",
            transit_qc={"status": "pass", "delta_bic": 12.0, "ktmf_metric": 3.8, "ktmf_contributions": []},
            transit_qc_status="pass",
            transit_qc_summary="ok",
            transit_qc_delta_bic=12.0,
            transit_qc_ktmf_metric=3.8,
        )
        return fit, np.asarray(tflux, dtype=float), np.asarray(cflux, dtype=float)

    monkeypatch.setattr("exotic.exotic.diagnose_lightcurve_fit_inputs", fake_diagnostics)
    monkeypatch.setattr("exotic.exotic.fit_lightcurve", fake_fit_lightcurve)

    result, tflux_fit, cflux_fit = evaluate_lightcurve_candidate(
        (
            np.linspace(0.0, 0.05, 6),
            np.full(6, 20.0),
            np.full(6, 10.0),
            np.linspace(1.0, 1.2, 6),
            [0.1, 0.1, 0.1, 0.1],
            {"midT": 0.5, "pPer": 1.0, "rprs": 0.1, "aRs": 10.0, "inc": 89.0, "ecc": 0.0, "omega": 0.0},
            2460000.0 + np.linspace(0.0, 0.05, 6),
            None,
            False,
            True,
            True,
            True,
        )
    )

    assert result["accepted"] is True
    assert result["ktmf_metric"] == pytest.approx(3.8)
    assert tflux_fit.shape == (6,)
    assert cflux_fit.shape == (6,)


def test_fit_lightcurve_refines_nested_tmid_bounds_from_two_sided_lm_fit(monkeypatch):
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
        duration_prior=None,
    ):
        captured_calls.append({
            "mode": mode,
            "bounds": {
                key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value
                for key, value in bounds.items()
            },
        })
        if mode == "lm":
            transit = np.ones_like(times, dtype=float)
            transit[(times >= 0.018) & (times <= 0.032)] = 0.98
            return types.SimpleNamespace(
                transit=transit,
                parameters={"tmid": 0.025, "rprs": 0.1, "inc": 89.0, "a2": 0.0, "per": 1.0},
                duration_expected=0.014,
            )
        return types.SimpleNamespace(parameters={"tmid": 0.025, "rprs": 0.1, "inc": 89.0, "a2": 0.0})

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

    times = np.linspace(0.0, 0.05, 21)
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
    assert captured_calls[0]["mode"] == "lm"
    assert captured_calls[1]["mode"] == "ns"
    assert captured_calls[0]["bounds"]["tmid"] == pytest.approx([0.011347361950458953, 0.03865263804954105])
    assert captured_calls[1]["bounds"]["tmid"] == pytest.approx([0.0175, 0.0325])
    assert myfit.nested_tmid_refinement_applied is True
    assert "recenter nested-sampling Tmid bounds" in myfit.nested_tmid_refinement_note


def test_fit_lightcurve_skips_nested_tmid_refinement_for_one_sided_lm_fit(monkeypatch):
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
        duration_prior=None,
    ):
        captured_calls.append({
            "mode": mode,
            "bounds": {
                key: list(value) if isinstance(value, (list, tuple, np.ndarray)) else value
                for key, value in bounds.items()
            },
        })
        if mode == "lm":
            transit = np.ones_like(times, dtype=float)
            transit[times >= 0.025] = 0.98
            return types.SimpleNamespace(
                transit=transit,
                parameters={"tmid": 0.025, "rprs": 0.1, "inc": 89.0, "a2": 0.0, "per": 1.0},
                duration_expected=0.014,
            )
        return types.SimpleNamespace(parameters={"tmid": 0.025, "rprs": 0.1, "inc": 89.0, "a2": 0.0})

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

    times = np.linspace(0.0, 0.05, 21)
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
    assert captured_calls[0]["mode"] == "lm"
    assert captured_calls[1]["mode"] == "ns"
    assert captured_calls[0]["bounds"]["tmid"] == pytest.approx([0.011347361950458953, 0.03865263804954105])
    assert captured_calls[1]["bounds"]["tmid"] == pytest.approx([0.011347361950458953, 0.03865263804954105])
    assert myfit.nested_tmid_refinement_applied is False
    assert "one side of the modeled transit" in myfit.nested_tmid_refinement_note


def test_run_target_driven_photometry_search_selects_best_method_across_psf_and_aperture(monkeypatch):
    evaluated = []

    class DummyFit:
        def __init__(self, residual_level, delta_bic, ktmf_metric):
            self.residuals = np.full(6, residual_level)
            self.data = np.ones(6)
            self.transit_qc_delta_bic = delta_bic
            self.transit_qc_ktmf_metric = ktmf_metric

    def fake_evaluate(task):
        _, tflux, cflux, *_ = task
        evaluated.append(np.asarray(cflux))
        cflux = np.asarray(cflux)
        tflux = np.asarray(tflux)
        if np.allclose(cflux, 20.0):
            return {
                "myfit": DummyFit(0.02, 9.0, 2.80),
                "res_std": 0.02,
                "transit_delta_bic": 9.0,
                "ktmf_metric": 2.80,
            }, tflux, cflux
        if np.allclose(cflux, 40.0):
            return {
                "myfit": DummyFit(0.01, 18.0, 4.85),
                "res_std": 0.01,
                "transit_delta_bic": 18.0,
                "ktmf_metric": 4.85,
            }, tflux, cflux
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
    assert result["selection_metric"] == "ktmf"
    assert result["best_candidate"]["method"] == "aperture"
    assert result["best_candidate"]["comp_index"] == 0
    assert result["selected_ktmf_metric"] == pytest.approx(4.85)
    assert result["selected_transit_delta_bic"] == pytest.approx(18.0)


def test_apply_raw_target_photometry_selection_sets_no_comparison_aperture_sentinel():
    fit = types.SimpleNamespace(time=np.linspace(0.0, 0.05, 6))
    target_flux = np.linspace(1000.0, 1010.0, 6)
    target_driven_search = {
        "best_candidate": {
            "method": "aperture",
            "a": 1,
            "an": 2,
            "aper": 5.0,
            "annulus": 12.0,
            "comp_index": None,
        },
        "best_fit_lc": fit,
        "selected_ktmf_metric": 3.5,
        "selected_transit_delta_bic": 12.0,
        "selection_metric": "ktmf",
        "selected_eebls_snr": 7.0,
        "flux_tar": target_flux,
        "flux_ref": np.ones(6),
        "selected_source_indices": np.arange(6),
        "candidate_summaries": [{"selected": True, "fit_point_count": 6}],
    }
    photometry_info = {"min_aperture": None, "comp_star_num": None}
    flux_values = {}
    centroid_positions = {}
    psf_data = {
        "target": np.column_stack([
            np.linspace(10.0, 15.0, 6),
            np.linspace(20.0, 25.0, 6),
        ])
    }

    applied = apply_raw_target_photometry_selection(
        target_driven_search,
        photometry_info,
        flux_values,
        centroid_positions,
        psf_data,
    )

    assert applied is True
    assert photometry_info["best_fit_lc"] is fit
    assert photometry_info["comp_star_num"] is None
    assert photometry_info["min_aperture"] == pytest.approx(-5.0)
    assert photometry_info["min_annulus"] == pytest.approx(12.0)
    assert photometry_info["selection_basis"] == "raw_target_flux_fallback"
    assert flux_values["flux_tar"] == pytest.approx(target_flux)
    assert flux_values["flux_ref"] == pytest.approx(np.ones(6))
    assert flux_values["flux_unc_ref"] == pytest.approx(np.zeros(6))
    assert np.isnan(centroid_positions["x_ref"]).all()
    assert np.isnan(centroid_positions["y_ref"]).all()


def test_run_target_driven_photometry_search_can_select_raw_target_without_comparison(monkeypatch):
    times = np.linspace(0.0, 0.05, 6)
    target_flux = np.linspace(1000.0, 1010.0, 6)

    def fake_evaluate(task):
        candidate_times, candidate_target_flux, candidate_reference_flux, *_ = task
        fit = types.SimpleNamespace(
            time=np.asarray(candidate_times, dtype=float),
            residuals=np.full(6, 0.01),
            data=np.ones(6),
        )
        return {
            "myfit": fit,
            "accepted": True,
            "ktmf_metric": 3.0,
            "fit_point_count": 6,
        }, np.asarray(candidate_target_flux), np.asarray(candidate_reference_flux)

    monkeypatch.setattr("exotic.exotic.evaluate_lightcurve_candidate", fake_evaluate)
    psf_data = {
        "target": np.column_stack([
            np.linspace(10.0, 15.0, 6),
            np.linspace(20.0, 25.0, 6),
        ])
    }
    aper_data = {"target": target_flux.reshape(6, 1, 1)}

    result = run_target_driven_photometry_search(
        times,
        2460000.0 + times,
        np.linspace(1.0, 1.5, 6),
        ld=[0.1, 0.1, 0.1, 0.1],
        p_dict={},
        comp_stars=[],
        psf_data=psf_data,
        aper_data=aper_data,
        apers=np.array([5.0]),
        annuli=np.array([12.0]),
        sigma=1.0,
        require_comp_star=False,
        use_psf_photometry=False,
        use_aperture_photometry=True,
    )

    assert result["best_candidate"]["comp_index"] is None
    assert result["best_candidate"]["method"] == "aperture"
    assert result["flux_tar"] == pytest.approx(target_flux)
    assert result["flux_ref"] == pytest.approx(np.ones(6))
    assert result["selected_source_indices"] == pytest.approx(np.arange(6))


def test_run_target_driven_photometry_search_can_prefer_highest_eebls_snr(monkeypatch):
    class DummyFit:
        def __init__(self, residual_level, eebls_snr, delta_bic):
            self.residuals = np.full(6, residual_level)
            self.data = np.ones(6)
            self.eebls_diagnostic_depth_snr = eebls_snr
            self.transit_qc_delta_bic = delta_bic

    def fake_evaluate(task):
        _, tflux, cflux, *_ = task
        cflux = np.asarray(cflux, dtype=float)
        tflux = np.asarray(tflux, dtype=float)
        if np.allclose(cflux, 20.0):
            return {"myfit": DummyFit(0.01, 4.0, 20.0), "res_std": 0.01, "eebls_snr": 4.0, "transit_delta_bic": 20.0}, tflux, cflux
        if np.allclose(cflux, 40.0):
            return {"myfit": DummyFit(0.02, 9.0, 12.0), "res_std": 0.02, "eebls_snr": 9.0, "transit_delta_bic": 12.0}, tflux, cflux
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
        pick_comparison_by_eebls_snr=True,
    )

    assert result["selection_metric"] == "eebls_snr"
    assert result["best_candidate"]["method"] == "aperture"
    assert result["selected_eebls_snr"] == pytest.approx(9.0)
    assert result["selected_transit_delta_bic"] == pytest.approx(12.0)


def test_fit_ranked_comparison_calibration_candidates_retries_next_best_candidate(monkeypatch):
    class DummyFit:
        def __init__(self):
            self.residuals = np.full(6, 0.01)
            self.data = np.ones(6)

    def fake_finalize(times, tflux, cflux, airmass, ld, p_dict, jd_times=None, **kwargs):
        cflux = np.asarray(cflux, dtype=float)
        if np.allclose(cflux, 0.0):
            return {
                "applied": False,
                "fit": None,
                "failure_reason": "the raw comparison-candidate photometry did not yield a usable light curve.",
                "note": "test full reduction",
            }
        return {
            "applied": True,
            "fit": DummyFit(),
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": cflux,
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

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
            {"comp_index": 0, "key": "comp1", "aggregate_score": 0.01, "coverage_rejected": False},
            {"comp_index": 1, "key": "comp2", "aggregate_score": 0.02, "coverage_rejected": False},
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
    assert result["attempts"][0]["fit"] is None
    assert result["attempts"][0]["fit_diagnostics"]["failure_reason"] is not None
    assert result["attempts"][1]["fit"] is not None


def test_fit_ranked_comparison_calibration_candidates_archives_qc_failed_run_and_tries_next(monkeypatch, tmp_path):
    class DummyFit:
        def __init__(self, qc_status):
            self.residuals = np.full(6, 0.01)
            self.data = np.ones(6)
            self.transit_qc_status = qc_status
            self.transit_qc_summary = (
                "Transit detection not supported strongly enough against a flat/null model (Delta BIC=2.50, Delta chi2=1.10)."
                if qc_status == "fail"
                else "Transit model strongly preferred over flat/null model (Delta BIC=18.40, Delta chi2=27.10)."
            )
            self.transit_qc = {"status": qc_status, "summary": self.transit_qc_summary}

    def fake_finalize(times, tflux, cflux, airmass, ld, p_dict, jd_times=None, **kwargs):
        cflux = np.asarray(cflux, dtype=float)
        fit = DummyFit("fail" if np.allclose(cflux, 8.0) else "pass")
        return {
            "applied": True,
            "fit": fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": cflux,
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    def fake_save(save_dir, provisional_fit, final_fit, p_dict, observation_date, comp_index, **kwargs):
        candidate_dir = Path(save_dir) / f"comp{comp_index + 1}"
        candidate_dir.mkdir(parents=True, exist_ok=True)
        return candidate_dir

    monkeypatch.setattr(
        "exotic.exotic.save_comparison_candidate_full_reduction_outputs",
        fake_save,
    )

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.5, 6)
    comparison_calibration = {
        "method": "aperture",
        "method_label": "Aperture photometry (aper=5.00px, annulus=12.00px)",
        "a": 0,
        "an": 0,
        "aper": 5.0,
        "annulus": 12.0,
        "best_comp_index": 0,
        "comp_summaries": [
            {"comp_index": 0, "key": "comp1", "label": "Comp 1", "position": [100.0, 200.0], "aggregate_score": 0.01, "coverage_rejected": False},
            {"comp_index": 1, "key": "comp2", "label": "Comp 2", "position": [300.0, 400.0], "aggregate_score": 0.02, "coverage_rejected": False},
        ],
    }
    aper_data = {
        "target": np.full((6, 1, 1), 10.0),
        "comp1": np.full((6, 1, 1), 8.0),
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
        save_dir=tmp_path,
        planet_name="HAT-P-32 b",
        observation_date="2026-04-28",
    )

    assert result["selected_result"]["comp_index"] == 1
    assert result["attempts"][0]["rejected_by_transit_qc"] is True
    assert result["attempts"][0]["fit_diagnostics"]["failed_stage"] == "transit_qc"
    failed_run_dir = result["attempts"][0]["failed_run_dir"]
    assert failed_run_dir is not None
    assert (
        tmp_path
        / "Diagnostics"
        / "comp_1_failed"
        / "working_artifacts"
        / "FailedFitSummary_HAT-P-32b_2026-04-28.json"
    ).exists()
    assert Path(failed_run_dir).exists()


def test_fit_ranked_comparison_calibration_candidates_falls_back_to_best_qc_rejected_fit(
    monkeypatch,
):
    class DummyFit:
        def __init__(self, ktmf, delta_bic):
            self.residuals = np.full(6, 0.01)
            self.data = np.ones(6)
            self.transit_qc_status = "fail"
            self.transit_qc_summary = (
                "Transit model is preferred over the flat/null model, but QC rejected the fit because "
                "the fit deviates too far from the expected published Rp/R* value "
                f"(Delta BIC={delta_bic:.2f}, Delta chi2=27.10)."
            )
            self.transit_qc = {
                "status": "fail",
                "summary": self.transit_qc_summary,
                "ktmf_metric": ktmf,
                "delta_bic": delta_bic,
            }
            self.transit_qc_ktmf_metric = ktmf
            self.transit_qc_delta_bic = delta_bic

    def fake_finalize(times, tflux, cflux, airmass, ld, p_dict, jd_times=None, **kwargs):
        cflux = np.asarray(cflux, dtype=float)
        comp_marker = int(np.nanmedian(cflux))
        fit = DummyFit(
            ktmf={8: 3.10, 5: 4.80}[comp_marker],
            delta_bic={8: 18.0, 5: 30.0}[comp_marker],
        )
        return {
            "applied": True,
            "fit": fit,
            "good_target_flux": np.asarray(tflux, dtype=float),
            "good_comp_flux": cflux,
            "source_indices": np.arange(len(times), dtype=int),
            "duration_samples": np.array([], dtype=float),
            "data_highres": None,
            "note": "test full reduction",
        }

    monkeypatch.setattr("exotic.exotic.finalize_comparison_candidate_full_reduction", fake_finalize)

    times = np.linspace(0.0, 0.05, 6)
    jd_times = 2460000.0 + times
    airmass = np.linspace(1.0, 1.5, 6)
    comparison_calibration = {
        "method": "aperture",
        "method_label": "Aperture photometry (aper=5.00px, annulus=12.00px)",
        "a": 0,
        "an": 0,
        "aper": 5.0,
        "annulus": 12.0,
        "best_comp_index": 0,
        "comp_summaries": [
            {"comp_index": 0, "key": "comp1", "label": "Comp 1", "aggregate_score": 0.01, "coverage_rejected": False},
            {"comp_index": 1, "key": "comp2", "label": "Comp 2", "aggregate_score": 0.02, "coverage_rejected": False},
        ],
    }
    aper_data = {
        "target": np.full((6, 1, 1), 10.0),
        "comp1": np.full((6, 1, 1), 8.0),
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

    assert [attempt["rejected_by_transit_qc"] for attempt in result["attempts"]] == [True, True]
    assert result["selection_metric"] == "ktmf_combined_quality"
    assert result["selected_result"]["comp_index"] == 1
    assert result["selected_result"]["selected_despite_transit_qc"] is True
    assert result["selected_result"]["ktmf_metric"] == pytest.approx(4.80)
    assert "best available fallback" in result["selected_result"]["selection_reason"]


def test_run_target_driven_photometry_search_skips_qc_failed_candidate(monkeypatch):
    class DummyFit:
        def __init__(self, residual_level, delta_bic=np.nan):
            self.residuals = np.full(6, residual_level)
            self.data = np.ones(6)
            self.transit_qc_delta_bic = delta_bic

    def fake_evaluate(task):
        _, tflux, cflux, *_ = task
        cflux = np.asarray(cflux, dtype=float)
        tflux = np.asarray(tflux, dtype=float)
        if np.allclose(cflux, 20.0):
            return {
                "myfit": DummyFit(0.005, 3.0),
                "accepted": False,
                "res_std": 0.005,
                "eebls_snr": 7.0,
                "transit_delta_bic": 3.0,
                "transit_qc_status": "fail",
                "transit_qc_summary": "Transit detection not supported strongly enough against a flat/null model.",
                "rejected_by_transit_qc": True,
                "fit_diagnostics": {"failed_stage": "transit_qc", "usable_point_count": 6},
                "failure_reason": "Transit detection not supported strongly enough against a flat/null model.",
                "fit_point_count": 6,
            }, tflux, cflux
        if np.allclose(cflux, 40.0):
            return {
                "myfit": DummyFit(0.02, 14.0),
                "accepted": True,
                "res_std": 0.02,
                "eebls_snr": 4.0,
                "transit_delta_bic": 14.0,
                "transit_qc_status": "pass",
                "transit_qc_summary": "Transit model strongly preferred over flat/null model.",
                "rejected_by_transit_qc": False,
                "fit_diagnostics": {"usable_point_count": 6},
                "failure_reason": None,
                "fit_point_count": 6,
            }, tflux, cflux
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

    assert len(result["candidate_summaries"]) == 2
    assert result["candidate_summaries"][0]["rejected_by_transit_qc"] is True
    assert result["best_candidate"]["method"] == "aperture"
    assert result["selected_transit_delta_bic"] == pytest.approx(14.0)


def test_diagnose_lightcurve_fit_inputs_allows_large_ratios():
    diagnostics = diagnose_lightcurve_fit_inputs(
        np.linspace(0.0, 0.05, 6),
        np.full(6, 30.0),
        np.full(6, 10.0),
        np.linspace(1.0, 1.5, 6),
    )

    assert diagnostics["failure_reason"] is None
    assert diagnostics["relative_flux_point_count"] == 6
    assert diagnostics["usable_point_count"] >= 5


def test_prepare_lightcurve_fit_input_series_normalizes_ratio_around_unity():
    times = np.linspace(0.0, 0.05, 6)
    prepared = prepare_lightcurve_fit_input_series(
        times,
        np.full(6, 30.0),
        np.full(6, 10.0),
        np.linspace(1.0, 1.5, 6),
    )

    assert prepared["applied"] is True
    assert np.nanmedian(prepared["debug_raw_ratio"]) == pytest.approx(3.0)
    assert prepared["approximate_baseline_level"] == pytest.approx(3.0)
    assert np.nanmedian(prepared["flux"]) == pytest.approx(1.0)


def test_prepare_lightcurve_fit_input_series_rejects_fewer_than_five_usable_points(monkeypatch):
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.array(
            [False, False, False, False, True, True],
            dtype=bool,
        ),
    )

    prepared = prepare_lightcurve_fit_input_series(
        np.linspace(0.0, 0.05, 6),
        np.full(6, 30.0),
        np.full(6, 10.0),
        np.linspace(1.0, 1.5, 6),
    )

    assert prepared["applied"] is False
    assert prepared["failure_reason"] == (
        "only 4 usable point(s) remained after filtering; "
        "need at least 5 for a lightcurve fit."
    )
    assert np.isnan(prepared["approximate_baseline_level"])


def test_prepare_lightcurve_fit_input_series_uses_per_star_flux_errors(monkeypatch):
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

    times = np.linspace(0.0, 0.05, 6)
    target_flux = np.full(6, 400.0)
    comp_flux = np.full(6, 100.0)
    target_error = np.full(6, 20.0)
    comp_error = np.full(6, 5.0)

    prepared = prepare_lightcurve_fit_input_series(
        times,
        target_flux,
        comp_flux,
        np.linspace(1.0, 1.5, 6),
        target_flux_error=target_error,
        comp_flux_error=comp_error,
    )

    propagated_relative_error = np.sqrt((20.0 / 100.0) ** 2 + (5.0 * 400.0 / 100.0 ** 2) ** 2)
    assert prepared["applied"] is True
    assert np.nanmedian(prepared["debug_relative_flux_error"]) == pytest.approx(propagated_relative_error)
    assert np.nanmedian(prepared["unc"]) == pytest.approx(propagated_relative_error / 4.0)
    assert np.allclose(prepared["target_flux_error"], target_error)
    assert np.allclose(prepared["comp_flux_error"], comp_error)


def test_prepare_lightcurve_fit_input_series_scales_target_only_counts_to_max_exposure(monkeypatch):
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

    times = np.linspace(0.0, 0.05, 6)
    target_flux = np.array([100.0, 200.0, 200.0, 200.0, 300.0, 300.0])
    comp_flux = np.ones(6)
    exposure_times = np.array([30.0, 60.0, 60.0, 60.0, 60.0, 60.0])

    prepared = prepare_lightcurve_fit_input_series(
        times,
        target_flux,
        comp_flux,
        np.linspace(1.0, 1.5, 6),
        exposure_times_seconds=exposure_times,
        gain_e_per_adu=2.0,
    )

    expected_flux = np.array([200.0, 200.0, 200.0, 200.0, 300.0, 300.0])
    expected_error = np.sqrt(target_flux / 2.0) * np.array([2.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    assert prepared["applied"] is True
    assert prepared["debug_target_flux"] == pytest.approx(expected_flux)
    assert prepared["target_flux"] == pytest.approx(expected_flux)
    assert prepared["target_flux_error"] == pytest.approx(expected_error)
    assert prepared["debug_relative_flux_error"] == pytest.approx(expected_error)


def test_compute_photometry_noise_budget_includes_optional_terms():
    config = {
        "gain_e_per_adu": 2.0,
        "read_noise_electrons": 4.0,
        "dark_current_electrons_per_second_per_pixel": 0.1,
        "flat_field_fractional_error": 0.01,
        "telescope_aperture_m": 0.3,
        "scintillation_coefficient": 0.09,
        "elevation_m": 100.0,
        "enabled_terms": (
            "source",
            "sky_aperture",
            "sky_estimate",
            "read",
            "dark",
            "flat",
            "scintillation",
        ),
    }

    budget = compute_photometry_noise_budget(
        10000.0,
        3.0,
        50.0,
        200.0,
        exposure_s=60.0,
        airmass=1.2,
        noise_config=config,
    )

    assert budget["source"] == pytest.approx(np.sqrt(10000.0 / 2.0))
    assert budget["read"] == pytest.approx(np.sqrt(50.0 * (4.0 / 2.0) ** 2))
    assert budget["dark"] == pytest.approx(np.sqrt(50.0 * 0.1 * 60.0 / 2.0 ** 2))
    assert budget["flat"] == pytest.approx(100.0)
    assert budget["scintillation"] > 0
    assert budget["total"] > budget["flat"]


def test_noise_budget_config_reads_inits_and_header_values():
    header = {
        "GAIN": 99.0,
        "EGAIN": 1.5,
        "RDNOISE": 7.0,
        "DARKCURR": 0.02,
        "FLATERR": 0.003,
        "APR-DIA": 250.0,
    }
    config = noise_budget_config_from_info(
        {
            "read_noise_electrons": 5.0,
        },
        header=header,
    )

    assert config["gain_e_per_adu"] == pytest.approx(1.5)
    assert config["read_noise_electrons"] == pytest.approx(5.0)
    assert config["dark_current_electrons_per_second_per_pixel"] == pytest.approx(0.02)
    assert config["flat_field_fractional_error"] == pytest.approx(0.003)
    assert config["telescope_aperture_m"] == pytest.approx(0.25)
    assert "read" in config["enabled_terms"]
    assert "flat" in config["enabled_terms"]


def test_prepare_lightcurve_fit_input_series_clips_prefit_raw_ratio_outliers(monkeypatch):
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

    times = np.linspace(0.0, 0.08, 21)
    comp_flux = np.full(times.shape, 1000.0, dtype=float)
    raw_ratio = np.ones(times.shape, dtype=float)
    raw_ratio[8:13] = 0.98
    raw_ratio[15] = 1.55
    raw_ratio[16] = 0.72
    target_flux = raw_ratio * comp_flux

    prepared = prepare_lightcurve_fit_input_series(
        times,
        target_flux,
        comp_flux,
        np.linspace(1.0, 1.4, times.shape[0]),
        expected_transit_depth=0.02,
    )

    assert prepared["applied"] is True
    assert prepared["initial_sigma_keep_mask"].all()
    assert prepared["prefit_raw_ratio_keep_mask"].tolist()[15:17] == [False, False]
    assert np.any(np.isclose(prepared["time"], times[10]))
    assert not np.any(np.isclose(prepared["time"], times[15]))
    assert any(
        diagnostic["stage"] == "Pre-fit raw-ratio outlier clip"
        and diagnostic["dropped_point_count"] == 2
        for diagnostic in prepared["filter_diagnostics"]
    )


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
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

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


def test_fit_lightcurve_forwards_exposure_times_to_fitter(monkeypatch):
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
        exposure_times_seconds=None,
    ):
        captured["times"] = np.asarray(times, dtype=float)
        captured["exposure_times_seconds"] = None if exposure_times_seconds is None else np.asarray(
            exposure_times_seconds,
            dtype=float,
        )
        return types.SimpleNamespace()

    monkeypatch.setattr("exotic.exotic.lc_fitter", fake_lc_fitter)
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

    times = np.linspace(0.0, 0.05, 6)
    exposure_times = np.array([60.0, 60.0, 90.0, 90.0, 120.0, 120.0])
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
        exposure_times_seconds=exposure_times,
    )

    assert captured["times"] == pytest.approx(times)
    assert captured["exposure_times_seconds"] == pytest.approx(exposure_times)


def test_build_initial_ars_bounds_prefers_published_uncertainty_when_available(monkeypatch):
    import exotic.exotic as exotic_module

    monkeypatch.setattr(exotic_module, "ARS_RANGE_RESTRICTION_ENABLED", True)
    monkeypatch.setattr(exotic_module, "ARS_RANGE_RESTRICTION_PERCENTAGE", 10.0)

    assert build_initial_ars_bounds(15.0, 0.1) == pytest.approx([13.5, 16.5])
    assert build_initial_ars_bounds(15.0, None) == pytest.approx([11.25, 18.75])

    monkeypatch.setattr(exotic_module, "ARS_RANGE_RESTRICTION_ENABLED", False)
    assert build_initial_ars_bounds(15.0, None) == pytest.approx([11.25, 18.75])


def test_build_single_transit_duration_prior_uses_published_geometry_uncertainties():
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

    assert duration_prior["applied"] is True
    assert duration_prior["expected_duration"] > 0
    assert duration_prior["sigma_log_duration"] > 0
    assert duration_prior["relative_sigma"] >= 0.049
    assert duration_prior["source"] == "published geometry uncertainties"
    assert "published geometry uncertainties" in duration_prior["note"]


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
    monkeypatch.setattr(
        "exotic.exotic.sigma_clip",
        lambda data, sigma=3, dt=21, po=2, times=None: np.zeros(len(data), dtype=bool),
    )

    times = np.linspace(0.0, 0.05, 6)
    tflux = np.full(times.shape[0], 2.0)
    cflux = np.full(times.shape[0], 2.0)
    airmass = np.array([1.10, 1.11, 1.12, 1.13, 1.14, 1.15])
    jd_times = 2460000.0 + times
    ld = [0.1, 0.1, 0.1, 0.1]
    p_dict = {
        "rprs": 0.1,
        "aRs": 15.0,
        "aRsUnc": 0.1,
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
    assert list(captured["bounds"])[:4] == ["rprs", "tmid", "ars", "inc"]
    assert captured["bounds"]["ars"] == pytest.approx([13.5, 16.5])
    assert "a2" not in captured["bounds"]
    assert myfit.airmass_fit_skipped is True


def _run_main_until_vertical_flux_bound(
        monkeypatch,
        tmp_path,
         disable_vertical_flux_normalization=Ellipsis,
         random_seed=123,
         override=True,
         nasa_result=None,
         target_ra=10.0,
         target_dec=20.0,
         ephemeris_overrides=None,
         expected_ephemeris=None,
         expected_error=None):
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
        "ra": target_ra,
        "dec": target_dec,
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
    if ephemeris_overrides:
        user_pdict.update(ephemeris_overrides)
    exotic_info = {
        "save": tmp_path,
        "prered_file": prered_file,
        "file_time": "BJD_TDB",
        "file_units": "flux",
        "airmass_already_corrected": False,
        "random_seed": random_seed,
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
        override=override,
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
    if nasa_result is not None:
        class FakeNASAExoplanetArchive:
            def __init__(self, planet, non_interactive=False):
                self.planet = planet
                self.non_interactive = non_interactive

            def planet_info(self):
                return nasa_result

        monkeypatch.setattr(exotic_module, "NASAExoplanetArchive", FakeNASAExoplanetArchive)
    monkeypatch.setattr(
        exotic_module,
        "get_ld_values",
        lambda *_args, **_kwargs: ([0.1, 0.1, 0.1, 0.1], [0.1], [0.1], [0.1], [0.1]),
    )

    def fake_apply_vertical_flux_normalization_bound(prior, bounds, flux_values, disabled):
        captured["disabled"] = disabled
        if expected_ephemeris is not None:
            assert prior['per'] == pytest.approx(expected_ephemeris['pPer'])
            assert prior['tmid'] == pytest.approx(expected_ephemeris['midT'])
        raise BoundReached()

    monkeypatch.setattr(
        exotic_module,
        "apply_vertical_flux_normalization_bound",
        fake_apply_vertical_flux_normalization_bound,
    )

    if expected_error is not None:
        with pytest.raises(ValueError, match=expected_error):
            exotic_module.main()
        return None

    with pytest.raises(BoundReached):
        exotic_module.main()

    return captured["disabled"]


def test_main_prereduced_defaults_vertical_flux_normalization_to_enabled(monkeypatch, tmp_path):
    disabled = _run_main_until_vertical_flux_bound(monkeypatch, tmp_path)

    assert disabled is False


def test_main_prereduced_respects_disable_vertical_flux_normalization_option(monkeypatch, tmp_path):
    disabled = _run_main_until_vertical_flux_bound(monkeypatch, tmp_path, disable_vertical_flux_normalization=True)

    assert disabled is True


def test_main_prereduced_override_invalid_coordinates_use_nasa_fallback_without_prompt(monkeypatch, tmp_path):
    monkeypatch.setattr(
        'builtins.input',
        lambda prompt: pytest.fail("non-interactive coordinate resolution must not prompt"),
    )

    disabled = _run_main_until_vertical_flux_bound(
        monkeypatch,
        tmp_path,
        override=True,
        nasa_result=("Test Planet b", False, {"ra": 123.456, "dec": -45.678}),
        target_ra="not-an-ra",
        target_dec="not-a-dec",
    )

    assert disabled is False


def test_main_prereduced_override_missing_ephemeris_uses_nasa_fallback(monkeypatch, tmp_path):
    archive_parameters = {
        'pPer': 2.5,
        'pPerUnc': 0.001,
        'midT': 2450000.25,
        'midTUnc': 0.002,
    }

    disabled = _run_main_until_vertical_flux_bound(
        monkeypatch,
        tmp_path,
        override=True,
        nasa_result=("Test Planet b", False, archive_parameters),
        ephemeris_overrides={'pPer': None, 'midT': 0.0},
        expected_ephemeris=archive_parameters,
    )

    assert disabled is False


def test_main_prereduced_stops_before_fitting_when_required_ephemeris_cannot_be_resolved(
        monkeypatch, tmp_path):
    _run_main_until_vertical_flux_bound(
        monkeypatch,
        tmp_path,
        override=True,
        nasa_result=("Test Planet b", False, {'pPer': np.nan, 'midT': None}),
        ephemeris_overrides={'pPer': None, 'midT': 0.0},
        expected_error=r"Cannot start EXOTIC reduction.*pPer.*midT",
    )


def test_main_prereduced_generates_seed_after_candidate_falls_back_to_inits(monkeypatch, tmp_path):
    disabled = _run_main_until_vertical_flux_bound(
        monkeypatch,
        tmp_path,
        random_seed=None,
        override=False,
        nasa_result=("TOI-3514.01", True, None),
    )

    assert disabled is False


def test_cli_logs_unhandled_exception_once(monkeypatch):
    import exotic.exotic as exotic_module

    logged = []

    monkeypatch.setattr(exotic_module, "configure_runtime_logging", lambda *args, **kwargs: None)
    monkeypatch.setattr(exotic_module, "install_exception_hooks", lambda: None)
    monkeypatch.setattr(exotic_module, "main", lambda: (_ for _ in ()).throw(RuntimeError("boom")))

    def fake_log_exception(message, exc_type, exc_value, exc_traceback):
        logged.append((message, exc_type, str(exc_value), exc_traceback is not None))

    monkeypatch.setattr(exotic_module, "_log_exception_with_fallback", fake_log_exception)

    with pytest.raises(RuntimeError, match="boom"):
        exotic_module.cli()

    assert logged == [("Unhandled exception during EXOTIC run", RuntimeError, "boom", True)]


def test_package_init_exports_lazy_main_and_cli(monkeypatch):
    import exotic

    monkeypatch.setattr(exotic, "_load_runtime_callable", lambda name: lambda: name)

    assert exotic.main() == "main"
    assert exotic.cli() == "cli"


def test_package_init_loads_nested_runtime_for_archive_layout(monkeypatch):
    import exotic

    def fake_import_module(module_name):
        if module_name == "exotic.exotic.exotic":
            return types.SimpleNamespace(main=lambda: "nested-main")
        raise AssertionError(f"unexpected import: {module_name}")

    monkeypatch.setitem(exotic.__dict__, "__name__", "exotic.exotic")
    monkeypatch.setattr(exotic, "import_module", fake_import_module)

    assert exotic._load_runtime_callable("main")() == "nested-main"


def test_configure_runtime_logging_rebinds_console_handler_to_current_stdout(monkeypatch, tmp_path):
    import io
    import exotic.exotic as exotic_module

    original_handlers = list(exotic_module.log.handlers)
    original_configured = exotic_module._RUNTIME_LOGGING_CONFIGURED
    original_basename = exotic_module._RUNTIME_LOG_BASENAME
    original_path = exotic_module._RUNTIME_LOG_PATH

    try:
        exotic_module.log.handlers = []
        exotic_module._RUNTIME_LOGGING_CONFIGURED = False
        exotic_module._RUNTIME_LOG_BASENAME = None
        exotic_module._RUNTIME_LOG_PATH = None
        monkeypatch.setattr(exotic_module, "_reset_runtime_traceback_watchdog", lambda: None)

        first_stdout = io.StringIO()
        monkeypatch.setattr(exotic_module.sys, "stdout", first_stdout)
        exotic_module.configure_runtime_logging(output_dir=tmp_path, start_new_run=True)
        handler = exotic_module._find_runtime_handler(exotic_module._RUNTIME_CONSOLE_HANDLER_NAME)
        assert handler.stream is first_stdout

        second_stdout = io.StringIO()
        monkeypatch.setattr(exotic_module.sys, "stdout", second_stdout)
        exotic_module.configure_runtime_logging(output_dir=tmp_path)
        assert handler.stream is second_stdout
    finally:
        exotic_module._close_runtime_file_handler()
        exotic_module.log.handlers = original_handlers
        exotic_module._RUNTIME_LOGGING_CONFIGURED = original_configured
        exotic_module._RUNTIME_LOG_BASENAME = original_basename
        exotic_module._RUNTIME_LOG_PATH = original_path


def test_configure_runtime_logging_does_not_use_environment_root_handlers(monkeypatch, tmp_path, capsys):
    import io
    import logging
    import exotic.exotic as exotic_module

    class DisconnectedColabStream(io.StringIO):
        def write(self, _value):
            raise OSError(107, "Transport endpoint is not connected")

        def flush(self):
            raise OSError(107, "Transport endpoint is not connected")

    original_handlers = list(exotic_module.log.handlers)
    original_propagate = exotic_module.log.propagate
    original_configured = exotic_module._RUNTIME_LOGGING_CONFIGURED
    original_basename = exotic_module._RUNTIME_LOG_BASENAME
    original_path = exotic_module._RUNTIME_LOG_PATH
    root_logger = logging.getLogger()
    original_root_handlers = list(root_logger.handlers)
    original_root_level = root_logger.level

    try:
        exotic_module.log.handlers = []
        exotic_module.log.propagate = True
        exotic_module._RUNTIME_LOGGING_CONFIGURED = False
        exotic_module._RUNTIME_LOG_BASENAME = None
        exotic_module._RUNTIME_LOG_PATH = None
        root_logger.handlers = [logging.StreamHandler(DisconnectedColabStream())]
        root_logger.setLevel(logging.WARNING)
        monkeypatch.setattr(exotic_module, "_reset_runtime_traceback_watchdog", lambda: None)

        exotic_module.configure_runtime_logging(output_dir=tmp_path, start_new_run=True)
        exotic_module.log.debug("frame progress written only to EXOTIC's file handler")

        assert exotic_module.log.propagate is False
        assert root_logger.level == logging.WARNING
        assert "Logging error" not in capsys.readouterr().err
    finally:
        exotic_module._close_runtime_file_handler()
        exotic_module.log.handlers = original_handlers
        exotic_module.log.propagate = original_propagate
        exotic_module._RUNTIME_LOGGING_CONFIGURED = original_configured
        exotic_module._RUNTIME_LOG_BASENAME = original_basename
        exotic_module._RUNTIME_LOG_PATH = original_path
        root_logger.handlers = original_root_handlers
        root_logger.setLevel(original_root_level)


def test_runtime_file_handler_suppresses_disconnected_mount_and_reopens(monkeypatch, tmp_path, capsys):
    import io
    import logging
    import exotic.exotic as exotic_module

    class DisconnectedDriveStream(io.StringIO):
        def write(self, _value):
            raise OSError(107, "Transport endpoint is not connected")

        def flush(self):
            raise OSError(107, "Transport endpoint is not connected")

        def close(self):
            pass

    log_path = tmp_path / "EXOTIC_RunLog_test.log"
    handler = exotic_module.FailSoftRuntimeFileHandler(log_path, mode="a", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(message)s"))
    handler.stream = DisconnectedDriveStream()
    recovered_stream = io.StringIO()
    monkeypatch.setattr(handler, "_open", lambda: recovered_stream)

    try:
        handler.emit(logging.LogRecord("exotic", logging.DEBUG, __file__, 1, "frame 18", (), None))
        first_output = capsys.readouterr()
        assert "Logging error" not in first_output.err
        assert "run log stream disconnected" in first_output.out
        assert handler.stream is None

        handler.emit(logging.LogRecord("exotic", logging.DEBUG, __file__, 1, "frame 19", (), None))
        second_output = capsys.readouterr()
        assert "Logging error" not in second_output.err
        assert "run log stream disconnected" not in second_output.out
        assert recovered_stream.getvalue() == "frame 19\n"
    finally:
        handler.stream = None
        handler.close()


def test_runtime_output_directory_is_read_from_command_line_init_file(tmp_path):
    import exotic.exotic as exotic_module

    output_dir = tmp_path / "run output"
    init_path = tmp_path / "inits.json"
    init_path.write_text(json.dumps({
        "user_info": {"Directory to Save Plots": str(output_dir)},
    }), encoding="utf-8")

    assert exotic_module._runtime_output_directory_from_command_line(
        ["-red", str(init_path), "-ov"]
    ) == str(output_dir)
    assert exotic_module._runtime_output_directory_from_command_line(
        [f"--reduce={init_path}"]
    ) == str(output_dir)


def test_runtime_logging_relocates_startup_content_and_keeps_runs_unique(monkeypatch, tmp_path):
    import exotic.exotic as exotic_module

    original_handlers = list(exotic_module.log.handlers)
    original_configured = exotic_module._RUNTIME_LOGGING_CONFIGURED
    original_basename = exotic_module._RUNTIME_LOG_BASENAME
    original_path = exotic_module._RUNTIME_LOG_PATH

    try:
        exotic_module.log.handlers = []
        exotic_module._RUNTIME_LOGGING_CONFIGURED = False
        exotic_module._RUNTIME_LOG_BASENAME = None
        exotic_module._RUNTIME_LOG_PATH = None
        monkeypatch.setattr(exotic_module.tempfile, "gettempdir", lambda: str(tmp_path / "staging"))
        monkeypatch.setattr(exotic_module, "_reset_runtime_traceback_watchdog", lambda: None)

        exotic_module.configure_runtime_logging(start_new_run=True)
        staged_log = Path(exotic_module._RUNTIME_LOG_PATH)
        exotic_module.log_info("startup message before the save directory was known")

        output_dir = tmp_path / "output"
        exotic_module.configure_runtime_logging(output_dir=output_dir)
        first_log = Path(exotic_module._RUNTIME_LOG_PATH)
        exotic_module.log_info("message after the save directory was known")
        exotic_module.close_runtime_logging()

        assert not staged_log.exists()
        assert first_log.parent == output_dir.resolve() / "Diagnostics"
        first_content = first_log.read_text(encoding="utf-8")
        assert "startup message before the save directory was known" in first_content
        assert "message after the save directory was known" in first_content

        exotic_module.configure_runtime_logging(output_dir=output_dir, start_new_run=True)
        second_log = Path(exotic_module._RUNTIME_LOG_PATH)
        exotic_module.log_info("second run message")
        exotic_module.close_runtime_logging()

        assert second_log != first_log
        assert len(list((output_dir / "Diagnostics").glob("EXOTIC_RunLog_*.log"))) == 2
        assert "second run message" in second_log.read_text(encoding="utf-8")
    finally:
        exotic_module._close_runtime_file_handler()
        exotic_module.log.handlers = original_handlers
        exotic_module._RUNTIME_LOGGING_CONFIGURED = original_configured
        exotic_module._RUNTIME_LOG_BASENAME = original_basename
        exotic_module._RUNTIME_LOG_PATH = original_path


def test_log_exception_with_fallback_writes_traceback_to_current_stdout(monkeypatch, capsys):
    import exotic.exotic as exotic_module

    monkeypatch.setattr(exotic_module, "_logger_has_current_stdout_handler", lambda logger: False)

    try:
        raise RuntimeError("boom")
    except RuntimeError as exc:
        exotic_module._log_exception_with_fallback(
            "Unhandled exception during EXOTIC run",
            type(exc),
            exc,
            exc.__traceback__,
        )

    output = capsys.readouterr().out
    assert "Unhandled exception during EXOTIC run" in output
    assert "Traceback" in output
    assert "RuntimeError: boom" in output


def test_main_logs_direct_call_exceptions_to_current_stdout(monkeypatch, capsys):
    import exotic.exotic as exotic_module

    monkeypatch.setattr(exotic_module, "configure_runtime_logging", lambda *args, **kwargs: None)
    monkeypatch.setattr(exotic_module, "install_exception_hooks", lambda: None)
    monkeypatch.setattr(exotic_module, "_logger_has_current_stdout_handler", lambda logger: False)
    monkeypatch.setattr(exotic_module, "_main_impl", lambda: (_ for _ in ()).throw(RuntimeError("boom")))

    with pytest.raises(RuntimeError, match="boom"):
        exotic_module.main()

    output = capsys.readouterr().out
    assert "Unhandled exception during EXOTIC run" in output
    assert "RuntimeError: boom" in output


def test_main_suppresses_all_internal_logging_error_tracebacks(monkeypatch, capsys):
    import io
    import logging
    import exotic.exotic as exotic_module

    class DisconnectedColabStream(io.StringIO):
        def write(self, _value):
            raise OSError(107, "Transport endpoint is not connected")

        def flush(self):
            raise OSError(107, "Transport endpoint is not connected")

    environment_logger = logging.getLogger("test.disconnected_colab_handler")
    environment_logger.handlers = [logging.StreamHandler(DisconnectedColabStream())]
    environment_logger.propagate = False
    original_raise_exceptions = logging.raiseExceptions

    monkeypatch.setattr(exotic_module, "configure_runtime_logging", lambda *args, **kwargs: None)
    monkeypatch.setattr(exotic_module, "install_exception_hooks", lambda: None)
    monkeypatch.setattr(exotic_module, "cancel_runtime_traceback_watchdog", lambda: None)
    monkeypatch.setattr(exotic_module, "close_runtime_logging", lambda: None)

    def report_through_disconnected_handler():
        environment_logger.error("frame progress")
        return "completed"

    monkeypatch.setattr(exotic_module, "_main_impl", report_through_disconnected_handler)

    try:
        logging.raiseExceptions = True
        assert exotic_module.main() == "completed"
        output = capsys.readouterr()
        assert "--- Logging error ---" not in output.err
        assert "Transport endpoint is not connected" not in output.err
        assert logging.raiseExceptions is True
    finally:
        environment_logger.handlers = []
        logging.raiseExceptions = original_raise_exceptions
