import json
import requests
import pytest

import exotic.inputs as inputs_module
from exotic.inputs import Inputs, camera, parse_aavso_prereduced_overrides


def test_camera_accepts_cmos_as_ccd_without_prompt():
    assert camera("CMOS") == "CCD"


def test_camera_defaults_to_ccd_when_missing_or_unrecognized():
    assert camera(None) == "CCD"
    assert camera("") == "CCD"
    assert camera("mirrorless") == "CCD"


def test_camera_keeps_dslr_as_dslr():
    assert camera("DSLR") == "DSLR"
    assert camera("canon dslr") == "DSLR"


def test_comp_params_accepts_verbose_camera_key(tmp_path):
    init_data = {
        "user_info": {
            "Camera Type (e.g., CCD or DSLR; Note: if you are using a CMOS, please enter CCD here and then note your actual camera type in \"Observing Notes\")": "CCD"
        },
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["camera"] == "CCD"


def test_comp_params_defaults_require_comp_star_to_yes(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["require_comp_star"] == "y"


def test_comp_params_defaults_aavso_comp_to_no(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["aavso_comp"] == "n"


def test_comp_params_defaults_ignore_header_wcs_to_no(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["ignore_header_wcs"] == "n"


def test_comp_params_defaults_bad_wcs_threshold_percent_to_three(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["bad_wcs_threshold_percent"] == 3.0


def test_comp_params_defaults_pointing_rejection_sigma_to_four(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["pointing_rejection_sigma"] == pytest.approx(4.0)


def test_comp_params_defaults_skip_low_comparison_coverage_rejection_to_no(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["skip_low_comparison_coverage_rejection"] == "n"


def test_comp_params_defaults_fit_lightcurve_to_every_comparison_candidate_to_no(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["fit_lightcurve_to_every_comparison_candidate"] == "n"


def test_comp_params_defaults_ultranest_live_points_to_200(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["ultranest_min_num_live_points"] == 200


def test_comp_params_defaults_rprs_search_bound_max_to_half(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["rprs_search_bound_max"] == 0.5


def test_comp_params_defaults_sparse_posterior_live_point_retry_to_yes(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["use_sparse_posterior_live_point_retry"] == "y"


def test_comp_params_defaults_exit_at_first_qc_pass_solution_to_yes(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["exit_at_first_qc_pass_solution"] == "y"


def test_comp_params_defaults_disable_vertical_flux_normalization_to_false(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["disable_vertical_flux_normalization"] is False


def test_comp_params_defaults_detect_bad_pixels_before_photometry_to_yes(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["detect_bad_pixels_before_photometry"] == "y"


def test_comp_params_defaults_multiprocess_bad_pixel_precheck_to_no(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["multiprocess_bad_pixel_precheck"] == "n"


def test_comp_params_defaults_detrend_on_outoftransit_baseline_to_true(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["detrend_on_outoftransit_baseline"] is True


def test_comp_params_defaults_final_fit_baseline_duration_multiplier_to_one(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["final_fit_baseline_duration_multiplier"] == pytest.approx(1.0)


def test_comp_params_defaults_use_eebls_tmid_initializer_to_yes(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["use_eebls_to_initialize_tmid_and_bounds"] == "y"


def test_comp_params_defaults_use_impactparameter_fit_to_yes(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["use_impactparameter_rather_than_inclination_to_fit"] == "y"


def test_comp_params_defaults_use_psf_photometry_to_yes(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["use_psf_photometry"] == "y"


def test_comp_params_defaults_use_aperture_photometry_to_yes(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["use_aperture_photometry"] == "y"


def test_comp_params_defaults_fast_aperture_mask_to_false(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["fast_aperture_mask"] is False


def test_comp_params_defaults_use_adaptive_apertures_to_false(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["use_adaptive_apertures"] is False


def test_comp_params_reads_observatory_full_title_from_user_info(tmp_path):
    init_data = {
        "user_info": {"Observatory Full Title": "Whipple Observatory"},
        "optional_info": {},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["obs_name"] == "Whipple Observatory"


def test_comp_params_reads_require_comp_star_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"require_comp_star": "n"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["require_comp_star"] == "n"


def test_comp_params_reads_ignore_header_wcs_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"Ignore WCS in Header and Do Manual Alignment? (y/n)": "y"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["ignore_header_wcs"] == "y"


def test_comp_params_reads_bad_wcs_threshold_percent_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"bad_wcs_threshold_percent": 5.5},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["bad_wcs_threshold_percent"] == 5.5


def test_comp_params_reads_pointing_rejection_sigma_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"pointing_rejection_sigma": 3.5},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["pointing_rejection_sigma"] == 3.5


def test_comp_params_reads_skip_low_comparison_coverage_rejection_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"skip_low_comparison_coverage_rejection": "y"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["skip_low_comparison_coverage_rejection"] == "y"


def test_comp_params_reads_fit_lightcurve_to_every_comparison_candidate_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"fit_lightcurve_to_every_comparison_candidate": "y"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["fit_lightcurve_to_every_comparison_candidate"] == "y"


def test_comp_params_reads_ultranest_live_points_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"minimum number of live points for ultranest": 275},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["ultranest_min_num_live_points"] == 275


def test_comp_params_reads_rprs_search_bound_max_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"maximum Rp/Rs search bound": 0.35},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["rprs_search_bound_max"] == 0.35


def test_comp_params_reads_fast_ultranest_before_final_run_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"run fast ultranest before final run": "n"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["run_fast_ultranest_before_final_run"] == "n"


def test_comp_params_reads_sparse_posterior_live_point_retry_off_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"use_sparse_posterior_live_point_retry": "n"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["use_sparse_posterior_live_point_retry"] == "n"


def test_comp_params_reads_exit_at_first_qc_pass_solution_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"exit at first QC PASS solution": "n"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["exit_at_first_qc_pass_solution"] == "n"


def test_comp_params_reads_disable_vertical_flux_normalization_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"disable vertical flux normalization": True},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["disable_vertical_flux_normalization"] is True


def test_comp_params_reads_detect_bad_pixels_before_photometry_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"detect_bad_pixels_before_photometry": "n"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["detect_bad_pixels_before_photometry"] == "n"


def test_comp_params_reads_multiprocess_bad_pixel_precheck_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"multiprocess_bad_pixel_precheck": "y"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["multiprocess_bad_pixel_precheck"] == "y"


def test_comp_params_reads_detrend_on_outoftransit_baseline_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"detrend_on_outoftransit_baseline": True},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["detrend_on_outoftransit_baseline"] is True


def test_comp_params_reads_detrend_on_outoftransit_baseline_false_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"detrend_on_outoftransit_baseline": False},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["detrend_on_outoftransit_baseline"] is False


def test_comp_params_reads_final_fit_baseline_duration_multiplier_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"final_fit_baseline_duration_multiplier": 1.75},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["final_fit_baseline_duration_multiplier"] == pytest.approx(1.75)


def test_comp_params_reads_use_eebls_tmid_initializer_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"Use EEBLS to Initialize Tmid and Bounds? (y/n)": "n"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["use_eebls_to_initialize_tmid_and_bounds"] == "n"


def test_comp_params_reads_use_impactparameter_fit_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"use_impactparameter_rather_than_inclination_to_fit": "n"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["use_impactparameter_rather_than_inclination_to_fit"] == "n"


def test_comp_params_reads_use_psf_photometry_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"use_psf_photometry": "n"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["use_psf_photometry"] == "n"


def test_comp_params_reads_use_aperture_photometry_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"use_aperture_photometry": "n"},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["use_aperture_photometry"] == "n"


def test_comp_params_reads_use_adaptive_apertures_from_optional_info(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {"use_adaptive_apertures": True},
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["use_adaptive_apertures"] is True


class DummyResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


def test_comp_params_fetches_missing_gaia_astrometry_from_nextastro(tmp_path, monkeypatch):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {
            "Target Star RA": "01:02:03",
            "Target Star Dec": "+04:05:06",
            "Star Distance (pc)": None,
            "Star Proper Motion RA (mas/yr)": "",
            "Star Proper Motion DEC (mas/yr)": "null",
        },
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    called = {}

    def fake_get(url, params, timeout):
        called["url"] = url
        called["params"] = params
        called["timeout"] = timeout
        return DummyResponse({
            "gaia": {
                "distance_pc": 200.0,
                "pmra_mas_per_year": 10.0,
                "pmdec_mas_per_year": -20.0,
            }
        })

    monkeypatch.setattr(inputs_module.requests, "get", fake_get)

    inputs = Inputs(init_opt="y")
    planet_dict = inputs.comp_params(init_file, {})

    assert called["url"] == "https://archive.nextastro.org/single_star_gaia_distpm"
    assert called["timeout"] == 30
    assert called["params"]["ra"] == pytest.approx(15.5125)
    assert called["params"]["dec"] == pytest.approx(4.085)
    assert planet_dict["dist"] == 200.0
    assert planet_dict["pm_ra"] == 10.0
    assert planet_dict["pm_dec"] == -20.0


def test_comp_params_only_backfills_missing_gaia_fields(tmp_path, monkeypatch):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {
            "Target Star RA": 123.4501,
            "Target Star Dec": -12.3402,
            "Star Distance (pc)": 111.0,
            "Star Proper Motion RA (mas/yr)": None,
            "Star Proper Motion DEC (mas/yr)": "",
        },
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    monkeypatch.setattr(
        inputs_module.requests,
        "get",
        lambda url, params, timeout: DummyResponse({
            "gaia": {
                "distance_pc": 222.0,
                "pmra_mas_per_year": 8.5,
                "pmdec_mas_per_year": -4.25,
            }
        }),
    )

    inputs = Inputs(init_opt="y")
    planet_dict = inputs.comp_params(init_file, {})

    assert planet_dict["dist"] == 111.0
    assert planet_dict["pm_ra"] == 8.5
    assert planet_dict["pm_dec"] == -4.25


def test_comp_params_continues_when_nextastro_gaia_lookup_fails(tmp_path, monkeypatch):
    init_data = {
        "user_info": {},
        "optional_info": {},
        "planetary_parameters": {
            "Target Star RA": 123.4501,
            "Target Star Dec": -12.3402,
            "Star Distance (pc)": None,
            "Star Proper Motion RA (mas/yr)": "",
            "Star Proper Motion DEC (mas/yr)": None,
        },
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    def raise_request_exception(*args, **kwargs):
        raise requests.exceptions.RequestException("service unavailable")

    monkeypatch.setattr(inputs_module.requests, "get", raise_request_exception)

    inputs = Inputs(init_opt="y")
    planet_dict = inputs.comp_params(init_file, {})

    assert planet_dict["dist"] is None
    assert planet_dict["pm_ra"] == ""
    assert planet_dict["pm_dec"] is None


def test_prereduced_mode_forces_aavso_comp_to_no(tmp_path):
    pre_reduced_file = tmp_path / "prereduced.txt"
    pre_reduced_file.write_text("time flux uncertainty\n")

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "2020-01-01",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": {"ra": "", "dec": "", "x": "", "y": ""},
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["aavso_comp"] == "n"


def test_prereduced_allows_blank_observatory_location_for_bjd_tdb(tmp_path):
    pre_reduced_file = tmp_path / "prereduced.txt"
    pre_reduced_file.write_text("time flux uncertainty\n")

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "2020-01-01",
        "lat": "",
        "long": "",
        "elev": "",
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["lat"] is None
    assert info_dict["long"] is None
    assert info_dict["elev"] is None


def test_comp_params_accepts_blank_if_none_phot_comp_star_key(tmp_path):
    init_data = {
        "user_info": {},
        "optional_info": {
            "Comparison Star used in Photometry (blank if none)": {
                "ra": "",
                "dec": "",
                "x": "493",
                "y": "202",
            }
        },
        "planetary_parameters": {},
    }
    init_file = tmp_path / "inits.json"
    init_file.write_text(json.dumps(init_data))

    inputs = Inputs(init_opt="y")
    inputs.comp_params(init_file, {})

    assert inputs.info_dict["phot_comp_star"] == {"ra": "", "dec": "", "x": "493", "y": "202"}


def test_prereduced_uses_aavso_comp_star_metadata_without_prompt(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#COMP_STAR-XC={\"ra\": null, \"dec\": null, \"x\": \"493\", \"y\": \"202\"}\n"
        "#DATE,DIFF,ERR,DETREND_1\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "2020-01-01",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["phot_comp_star"] == {"ra": "", "dec": "", "x": "493", "y": "202"}


def test_prereduced_uses_aavso_observatory_metadata_without_prompt(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#OBSLAT=+32.41638889\n"
        "#OBSLON=-110.73444444\n"
        "#OBSELEV=2616\n"
        "#DATE,DIFF,ERR,DETREND_1\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "2020-01-01",
        "lat": "",
        "long": "",
        "elev": "",
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["lat"] == 32.41638889
    assert info_dict["long"] == -110.73444444
    assert info_dict["elev"] == 2616.0


def test_prereduced_uses_aavso_obsdate_metadata_without_prompt(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#OBSDATE=2026-03-08\n"
        "#DATE,DIFF,ERR,DETREND_1\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["date"] == "2026-03-08"


def test_prereduced_uses_aavso_filter_and_observing_metadata_without_prompt(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#OBSCODE=\n"
        "#SECONDARY_OBSCODES=\n"
        "#OBSNAME=Backyard Dome\n"
        "#OBSDATE=20260303\n"
        "#OBSTYPE=CCD\n"
        "#BINNING=1x1\n"
        "#EXPOSURE_TIME=30.0\n"
        "#OBSLAT=35.554298\n"
        "#OBSLON=-105.870197\n"
        "#OBSELEV=2194.0\n"
        "#GAIADIST=512.4\n"
        "#GAIAPMRA=13.25\n"
        "#GAIAPMDEC=-7.5\n"
        "#NOTES=na\n"
        "#DATE_TYPE=BJD_TDB\n"
        "#MEASUREMENT_TYPE=Rnflux\n"
        "#EXOPLANET_NAME=TOI-1259 A b\n"
        "#FILTER=CBB\n"
        "#FILTER-XC={\"name\": \"CBB\", \"desc\": \"Astrodon ExoPlanet-BB\", \"fwhm\": [{\"value\": \"500.0\", \"units\": \"nm\"}, {\"value\": \"1000.0\", \"units\": \"nm\"}]}\n"
        "#DATE,DIFF,ERR,DETREND_1\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": None,
        "second_obs": None,
        "date": "",
        "lat": "",
        "long": "",
        "elev": "",
        "camera": None,
        "pixel_bin": None,
        "filter": None,
        "notes": None,
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": None,
        "file_units": None,
        "file_time": None,
        "phot_comp_star": None,
        "wl_min": None,
        "wl_max": None,
    })

    info_dict, planet = inputs.prereduced(None)

    assert planet == "TOI-1259 A b"
    assert info_dict["aavso_num"] == ""
    assert info_dict["second_obs"] == ""
    assert info_dict["obs_name"] == "Backyard Dome"
    assert info_dict["date"] == "2026-03-03"
    assert info_dict["lat"] == 35.554298
    assert info_dict["long"] == -105.870197
    assert info_dict["elev"] == 2194.0
    assert info_dict["camera"] == "CCD"
    assert info_dict["pixel_bin"] == "1x1"
    assert info_dict["notes"] == "na"
    assert info_dict["file_time"] == "BJD_TDB"
    assert info_dict["file_units"] == "flux"
    assert info_dict["exposure"] == 30.0
    assert info_dict["dist"] == "512.4"
    assert info_dict["pm_ra"] == "13.25"
    assert info_dict["pm_dec"] == "-7.5"
    assert info_dict["filter"] == "CBB"
    assert info_dict["wl_min"] == "500.0"
    assert info_dict["wl_max"] == "1000.0"


def test_parse_aavso_prereduced_overrides_uses_known_filter_lookup_when_filter_xc_missing(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#FILTER=CBB\n"
        "#DATE,DIFF,ERR\n"
        "2461102.76092732,0.979108,0.0386426\n"
    )

    overrides = parse_aavso_prereduced_overrides(pre_reduced_file)

    assert overrides["filter"] == "CBB"
    assert overrides["filter_desc"] == "Astrodon ExoPlanet-BB"
    assert overrides["wl_min"] == "500.0"
    assert overrides["wl_max"] == "1000.0"


def test_parse_aavso_prereduced_overrides_uses_astrodon_exo_alias_lookup(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#FILTER=Astrodon-Exo\n"
        "#DATE,DIFF,ERR\n"
        "2461102.76092732,0.979108,0.0386426\n"
    )

    overrides = parse_aavso_prereduced_overrides(pre_reduced_file)

    assert overrides["filter"] == "Astrodon-Exo"
    assert overrides["filter_desc"] == "Astrodon ExoPlanet-BB"
    assert overrides["wl_min"] == "500.0"
    assert overrides["wl_max"] == "1000.0"


def test_lookup_aavso_filter_metadata_uses_c_alias_for_cv_filter() -> None:
    filter_metadata = inputs_module.lookup_aavso_filter_metadata("C")

    assert filter_metadata["name"] == "CV"
    assert filter_metadata["desc"] == "MObs CV"
    assert filter_metadata["fwhm"] == ("350.0", "850.0")


def test_lookup_aavso_filter_metadata_uses_luminosity_aliases_for_clearv_filter() -> None:
    for alias in ("lum", "Lum", "Luminosity", "luminosity"):
        filter_metadata = inputs_module.lookup_aavso_filter_metadata(alias)

        assert filter_metadata["name"] == "CV"
        assert filter_metadata["desc"] == "ClearV"
        assert filter_metadata["fwhm"] == ("350.0", "1000.0")


def test_parse_aavso_prereduced_overrides_uses_osc_split_filter_alias_lookup(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#FILTER=G2\n"
        "#DATE,DIFF,ERR\n"
        "2461102.76092732,0.979108,0.0386426\n"
    )

    overrides = parse_aavso_prereduced_overrides(pre_reduced_file)

    assert overrides["filter"] == "G2"
    assert overrides["filter_desc"] == "Photographic G"
    assert overrides["wl_min"] == "502.8"
    assert overrides["wl_max"] == "586.8"


def test_parse_aavso_prereduced_overrides_uses_c_alias_for_cv_filter_lookup(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#FILTER=C\n"
        "#DATE,DIFF,ERR\n"
        "2461102.76092732,0.979108,0.0386426\n"
    )

    overrides = parse_aavso_prereduced_overrides(pre_reduced_file)

    assert overrides["filter"] == "C"
    assert overrides["filter_desc"] == "MObs CV"
    assert overrides["wl_min"] == "350.0"
    assert overrides["wl_max"] == "850.0"


def test_parse_aavso_prereduced_overrides_uses_luminosity_alias_for_clearv_lookup(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#FILTER=Luminosity\n"
        "#DATE,DIFF,ERR\n"
        "2461102.76092732,0.979108,0.0386426\n"
    )

    overrides = parse_aavso_prereduced_overrides(pre_reduced_file)

    assert overrides["filter"] == "Luminosity"
    assert overrides["filter_desc"] == "ClearV"
    assert overrides["wl_min"] == "350.0"
    assert overrides["wl_max"] == "1000.0"


def test_parse_aavso_prereduced_overrides_marks_airmass_as_already_corrected(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#DETREND_PARAMETERS=AIRMASS, AIRMASS CORRECTION FUNCTION\n"
        "#DATE,DIFF,ERR,DETREND_1,DETREND_2\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172,0.998\n"
    )

    overrides = parse_aavso_prereduced_overrides(pre_reduced_file)

    assert overrides["airmass_already_corrected"] is True


def test_prereduced_prefers_aavso_obsdate_metadata_over_init_date(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#OBSDATE=2026-03-08\n"
        "#DATE,DIFF,ERR,DETREND_1\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "1999-01-01",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["date"] == "2026-03-08"


def test_prereduced_derives_obsdate_from_first_data_row_without_prompt(tmp_path):
    pre_reduced_file = tmp_path / "prereduced.txt"
    pre_reduced_file.write_text(
        "time,flux,uncertainty\n"
        "2458849.5,0.979108,0.0386426\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "JD_UTC",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["date"] == "2020-01-01"


def test_prereduced_leaves_phot_comp_star_blank_when_missing_from_aavso_metadata(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#DATE,DIFF,ERR,DETREND_1\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "2020-01-01",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["phot_comp_star"] == {"ra": "", "dec": "", "x": "", "y": ""}


def test_prereduced_carries_aavso_airmass_corrected_flag(tmp_path):
    pre_reduced_file = tmp_path / "aavso_prereduced.txt"
    pre_reduced_file.write_text(
        "#TYPE=EXOPLANET\n"
        "#DETREND_PARAMETERS=AIRMASS, AIRMASS CORRECTION FUNCTION\n"
        "#DATE,DIFF,ERR,DETREND_1,DETREND_2\n"
        "2461102.76092732,0.979108,0.0386426,1.3811172,0.998\n"
    )

    inputs = Inputs(init_opt="y")
    inputs.info_dict.update({
        "save": str(tmp_path),
        "aavso_num": "RTZ",
        "second_obs": "",
        "date": "2020-01-01",
        "lat": "+0.0",
        "long": "+0.0",
        "elev": 1.0,
        "camera": "CCD",
        "pixel_bin": "1x1",
        "notes": "na",
        "aavso_comp": "y",
        "prered_file": str(pre_reduced_file),
        "exposure": 60.0,
        "file_units": "flux",
        "file_time": "BJD_TDB",
        "phot_comp_star": None,
    })

    info_dict, _ = inputs.prereduced("HAT-P-32 b")

    assert info_dict["airmass_already_corrected"] is True
