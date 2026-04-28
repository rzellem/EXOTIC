import importlib
import importlib.util
import sys
import types
from pathlib import Path
import numpy as np
import pytest


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
    adaptive_aperture_outlier_mask,
    annotate_transit_qc_expected_values,
    auto_tune_aperture_sigma_grid,
    build_initial_ars_bounds,
    build_single_transit_duration_prior,
    build_target_fit_candidate_jobs,
    build_time_rejection_diagnostic,
    check_coordinates,
    cheap_lightcurve_prescore,
    centroid_offset_matches_reference,
    choose_centroid_seed_position,
    compute_transit_qc_ktmf,
    apply_comparison_star_suitability_outlier_rejection,
    comparison_calibration_selection_reason,
    comparison_candidate_fit_selection_reason,
    comparison_star_coverage_summary,
    comparison_star_stability_summary,
    deduplicate_comparison_star_coords,
    diagnose_lightcurve_fit_inputs,
    detrend_flux_on_out_of_transit_baseline,
    ensure_lightcurve_fit_failure_reason,
    evaluate_lightcurve_candidate,
    evaluate_transit_detection_qc,
    finalize_comparison_candidate_full_reduction,
    fit_lightcurve,
    fit_final_lightcurve_with_oot_baseline_detrending,
    fit_lightcurve_to_every_comparison_candidate,
    fit_ranked_comparison_calibration_candidates,
    get_final_fit_baseline_duration_multiplier,
    estimate_ephemeris_tmid_and_bounds,
    estimate_tmid_and_bounds_with_eebls,
    is_adaptive_aperture_mode_enabled,
    is_comp_star_required,
    is_out_of_transit_baseline_detrending_enabled,
    is_target_driven_comp_selection_enabled,
    log_comparison_calibration_fit_attempt_summaries,
    log_comparison_candidate_fit_summaries,
    log_target_fit_candidate_summaries,
    normalize_flux_series_to_approximate_unity,
    phase_bin_sigma_clip,
    parse_deviation_from_expected_transit_in_qc_sigma,
    prepare_final_fit_lightcurve_series,
    prepare_lightcurve_fit_input_series,
    representative_psf_sigma,
    ranked_comparison_calibration_summaries,
    resolve_sky_annulus_geometry,
    run_target_driven_photometry_search,
    resolve_frame_aperture_radii,
    robust_flux_floor_mask,
    robust_target_reference_flux_mask,
    save_selected_photometry_debug_series,
    should_keep_header_wcs_alignment,
    sigma_clip,
    summarize_adaptive_aperture_usage,
    summarize_prior_transit_coverage,
    should_skip_airmass_fit,
    should_use_eebls_to_initialize_tmid_and_bounds,
    should_fit_lightcurve_to_every_comparison_candidate,
    should_detect_bad_pixels_before_photometry,
    should_use_aperture_photometry,
    should_pick_comparison_by_eebls_snr,
    should_use_psf_photometry,
    should_skip_low_comparison_coverage_rejection,
    should_assess_all_comparisons_before_selecting_best,
    should_use_fast_target_centroid,
    should_use_deviation_from_expected_transit_in_qc,
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
            "phase_clip_keep_mask_on_sigma_filtered": np.array([True, False], dtype=bool),
        }
    )

    output_path = save_selected_photometry_debug_series(tmp_path, "Qatar-10 b", "20260420", fit)

    assert output_path is not None
    assert output_path.exists()

    rows = np.loadtxt(output_path, delimiter=",", skiprows=1)
    assert rows.shape == (3, 6)
    assert rows[:, 4].astype(int).tolist() == [1, 0, 1]
    assert rows[:, 5].astype(int).tolist() == [1, 0, 0]


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

    modeled = detrend_flux_on_out_of_transit_baseline(times, flux, unc, fit)
    fallback = detrend_flux_on_out_of_transit_baseline(times, flux, unc, fit, prior=prior)
    prior_coverage = summarize_prior_transit_coverage(times, prior, flux_values=flux, flux_errors=unc)

    assert modeled["applied"] is False
    assert prior_coverage["valid"] is True
    assert prior_coverage["has_two_sided_oot"] is True
    assert fallback["applied"] is True
    assert fallback["used_prior_ephemeris"] is True
    assert "ephemeris-centered transit window" in fallback["note"]


def test_deduplicate_comparison_star_coords_merges_nearby_duplicates():
    unique_coords, duplicate_messages = deduplicate_comparison_star_coords(
        [
            [1826.0, 1499.0],
            [1827.0, 1511.0],
            [1828.0, 1487.0],
            [842.0, 1810.0],
        ],
        min_separation_pixels=15.0,
    )

    assert unique_coords == [[1826.0, 1499.0], [842.0, 1810.0]]
    assert len(duplicate_messages) == 2
    assert "Merged comparison star #2" in duplicate_messages[0]


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


def test_build_target_fit_candidate_jobs_masks_psf_target_and_comp_dropouts():
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
    target_amplitudes[7] = 1.0
    comp_amplitudes[13] = 1.0

    psf_data = {
        "target": build_psf_rows(target_amplitudes),
        "comp1": build_psf_rows(comp_amplitudes),
    }

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
    assert not candidate_jobs[0]["mask"][7]
    assert not candidate_jobs[0]["mask"][13]
    assert candidate_jobs[0]["coverage_count"] == 29


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


def test_should_assess_all_comparisons_before_selecting_best_parses_values():
    assert should_assess_all_comparisons_before_selecting_best(None) is True
    assert should_assess_all_comparisons_before_selecting_best("y") is True
    assert should_assess_all_comparisons_before_selecting_best("n") is False
    assert should_assess_all_comparisons_before_selecting_best(True) is True


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


def test_resolve_sky_annulus_geometry_enforces_fwhm_floor_and_min_sky_pixels():
    geometry = resolve_sky_annulus_geometry(aperture_radius=1.5, annulus_width=2.0, psf_sigma=1.0)

    assert geometry["inner_radius"] == pytest.approx(2.0 * 2.355)
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
                    "label": "Delta BIC",
                    "available": True,
                    "points": 1.25,
                    "max_points": 1.40,
                    "score": 0.89,
                    "detail": "Delta BIC=18.40",
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
    assert any("KTMF contribution: Delta BIC +1.25/1.40" in message for message in logged)
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
        "rprs_sigma": 6.0,
        "duration_ratio": 1.0,
        "eebls_depth_snr": 8.0,
    }

    ktmf_metric, contributions = compute_transit_qc_ktmf(summary)
    contributions_by_label = {contribution["label"]: contribution for contribution in contributions}

    assert "Model Evidence" in contributions_by_label
    assert "Delta BIC" not in contributions_by_label
    assert "Delta chi2" not in contributions_by_label
    assert contributions_by_label["Model Evidence"]["max_points"] == pytest.approx(0.8)
    assert contributions_by_label["Deviation From Expected Value"]["max_points"] == pytest.approx(1.5)
    assert contributions_by_label["Residual Scatter Around Full Model Fit"]["max_points"] == pytest.approx(0.7)
    assert contributions_by_label["Duration Consistency"]["max_points"] == pytest.approx(0.75)
    assert contributions_by_label["EEBLS Depth SNR"]["max_points"] == pytest.approx(0.75)

    model_evidence_score = ((1.0 - np.exp(-1.0)) + (1.0 - np.exp(-2.0))) / 2.0
    expected_ktmf = (
        0.8 * model_evidence_score
        + 1.5 * 0.6
        + 0.7 * 0.5
        + 0.5 * (1.0 - np.exp(-2.0))
        + 0.75 * 1.0
        + 0.75 * (1.0 - np.exp(-2.0))
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
    assert summary["image_outlier_valid_pair_counts"][-1] == 2
    assert summary["image_outlier_outlier_pair_counts"][-1] == 2


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

    np.random.seed(0)
    old_mask = sigma_clip(values, sigma=3, dt=37, times=None)
    np.random.seed(0)
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
    assert captured["calls"][1]["bounds"]["rprs"] == pytest.approx([0.108, 0.208])
    assert fit.rprs_posterior_refit_applied is True
    assert fit.rprs_posterior_refit_count == 1
    assert fit.rprs_posterior_refit_edge == "upper"
    assert fit.rprs_posterior_refit_bounds == pytest.approx([0.108, 0.208])


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
    assert captured["prior"]["a0"] == pytest.approx(0.05)
    assert captured["prior"]["a1"] == pytest.approx(0.05)
    assert captured["bounds"]["a0"] == pytest.approx([0.0375, 0.0625])


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
    assert captured_duration_priors[1]["expected_duration"] > 0


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
        "Finite/positive photometry filter",
    ]
    assert diagnostics[0]["dropped_point_count"] == 1
    assert diagnostics[0]["first_dropped_time"] == pytest.approx(1.0)
    assert diagnostics[1]["dropped_point_count"] == 1
    assert diagnostics[1]["first_dropped_time"] == pytest.approx(11.0)
    assert diagnostics[2]["dropped_point_count"] == 0


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


def test_evaluate_transit_detection_qc_rejects_large_expected_value_deviation():
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
        parameters={"rprs": 0.18, "tmid": 0.5, "inc": 89.0, "a2": 0.0},
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
    assert summary["rprs_deviation_sigma"] == pytest.approx(8.0)
    assert summary["deviation_from_expected_value"] == pytest.approx(0.0)
    assert summary["ktmf_metric"] <= 5.0


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


def test_evaluate_transit_detection_qc_failure_summary_reflects_expected_value_rejection():
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
        parameters={"rprs": 0.18, "tmid": 0.5, "inc": 89.0, "a2": 0.0, "per": 2.0},
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

    assert summary["status"] == "fail"
    assert "QC rejected the fit because" in summary["summary"]
    assert "expected published Tmid and/or Rp/R*" in summary["summary"]
    assert "not supported strongly enough against a flat/null model" not in summary["summary"]


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
    assert result["selection_metric"] == "ktmf"
    assert result["selected_result"]["comp_index"] == 1
    assert result["selected_result"]["rank"] == 1
    assert result["selected_result"]["selected"] is True
    assert result["selected_result"]["ktmf_metric"] == pytest.approx(4.70)
    assert "highest KTMF" in result["selected_result"]["selection_reason"]
    assert result["attempts"][0]["selection_reason"].startswith("not selected: KTMF")


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


def test_fit_ranked_comparison_calibration_candidates_evaluates_all_candidates_even_when_flag_disabled(
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
        assess_all_comparisons_before_selecting_best=False,
        save_dir=tmp_path,
        planet_name="HAT-P-32 b",
        observation_date="2026-04-28",
    )

    assert call_markers == [50, 40]
    assert len(result["attempts"]) == 2
    assert result["selection_metric"] == "ktmf"
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
            residual_level = 0.02
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
    assert "highest EEBLS SNR" in result["selected_result"]["selection_reason"]
    assert result["attempts"][0]["selection_reason"].startswith("not selected: EEBLS SNR")


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
    assert (tmp_path / "comp_1_failed" / "temp" / "FailedFitSummary_HAT-P-32 b_2026-04-28.json").exists()
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
                "the fit deviates too far from the expected published Tmid and/or Rp/R* values "
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
    assert result["selection_metric"] == "ktmf"
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


def test_build_initial_ars_bounds_prefers_published_uncertainty_when_available():
    assert build_initial_ars_bounds(15.0, 0.1) == pytest.approx([14.5, 15.5])
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
    assert captured["bounds"]["ars"] == pytest.approx([14.5, 15.5])
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


def test_cli_logs_unhandled_exception_once(monkeypatch):
    import exotic.exotic as exotic_module

    logged = []

    monkeypatch.setattr(exotic_module, "configure_runtime_logging", lambda: None)
    monkeypatch.setattr(exotic_module, "install_exception_hooks", lambda: None)
    monkeypatch.setattr(exotic_module, "main", lambda: (_ for _ in ()).throw(RuntimeError("boom")))

    def fake_log_exception(message, exc_type, exc_value, exc_traceback):
        logged.append((message, exc_type, str(exc_value), exc_traceback is not None))

    monkeypatch.setattr(exotic_module, "_log_exception_with_fallback", fake_log_exception)

    with pytest.raises(RuntimeError, match="boom"):
        exotic_module.cli()

    assert logged == [("Unhandled exception during EXOTIC run", RuntimeError, "boom", True)]
