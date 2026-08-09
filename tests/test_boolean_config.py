import pytest

import exotic.exotic as exotic_module


BOOLEAN_RUNTIME_PARSERS = (
    exotic_module.is_fast_aperture_mask_enabled,
    exotic_module.is_comp_star_required,
    exotic_module.is_target_driven_comp_selection_enabled,
    exotic_module.should_skip_low_comparison_coverage_rejection,
    exotic_module.should_fit_lightcurve_to_every_comparison_candidate,
    exotic_module.should_use_automatic_optimal_calibration_selector,
    exotic_module.should_use_ensemble_photometry_rather_than_single_comp,
    exotic_module.should_use_ensemble_photometry_for_stellar_variability,
    exotic_module.should_photometer_fortuitous_variables,
    exotic_module.should_use_single_comparison_for_fortuitous_variables,
    exotic_module.should_use_nextastro_vsx_cache_first,
    exotic_module.should_use_sparse_posterior_live_point_retry,
    exotic_module.should_run_fast_ultranest_before_final_run,
    exotic_module.should_run_final_residual_rejection,
    exotic_module.should_use_legacy_psf_flux_mode,
    exotic_module.should_run_final_fit_phase_residual_clip,
    exotic_module.should_pick_comparison_by_eebls_snr,
    exotic_module.should_use_deviation_from_expected_transit_in_qc,
    exotic_module.should_exit_at_first_qc_pass_solution,
    exotic_module.should_restrict_rprs_range,
    exotic_module.should_use_prior_rprs_when_posterior_pinned,
    exotic_module.should_restrict_ars_range,
    exotic_module.should_use_psf_photometry,
    exotic_module.should_use_aperture_photometry,
    exotic_module.should_use_eebls_to_initialize_tmid_and_bounds,
    exotic_module.should_detect_bad_pixels_before_photometry,
    exotic_module.is_adaptive_aperture_mode_enabled,
    exotic_module.should_use_aperture_corrections_and_full_image_fwhm,
    exotic_module.should_reject_overexposed_stars,
    exotic_module.should_ignore_header_wcs,
    exotic_module.should_prefer_pixel_values_over_wcs_for_target,
    exotic_module.is_vertical_flux_normalization_disabled,
    exotic_module.should_run_stellar_variability_only,
    exotic_module.is_out_of_transit_baseline_detrending_enabled,
    exotic_module.should_use_impactparameter_rather_than_inclination_to_fit,
    exotic_module.should_require_apparent_magnitudes,
    exotic_module.should_use_exactly_the_comps_provided,
)


TRUE_VARIANTS = (True, 1, "1", "y", "Y", "yes", "TRUE", "on")
FALSE_VARIANTS = (False, 0, "0", "n", "N", "no", "FALSE", "off")


@pytest.mark.parametrize("parser", BOOLEAN_RUNTIME_PARSERS, ids=lambda parser: parser.__name__)
def test_runtime_boolean_parser_accepts_every_supported_form(parser):
    for value in TRUE_VARIANTS:
        assert parser(value) is True
    for value in FALSE_VARIANTS:
        assert parser(value) is False


def test_multiprocess_bad_pixel_boolean_forms_are_consistent():
    for value in TRUE_VARIANTS:
        assert exotic_module.get_multiprocess_bad_pixel_precheck_processes(value) >= 1
    for value in FALSE_VARIANTS:
        assert exotic_module.get_multiprocess_bad_pixel_precheck_processes(value) is None
