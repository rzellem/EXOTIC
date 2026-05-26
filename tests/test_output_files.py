import json

import numpy as np
import pytest

from exotic.output_files import AIDOutputFiles, OutputFiles, fit_impact_parameter_value_error, save_comp_star_calibration_summary


class DummyFit:
    def __init__(self):
        self.parameters = {
            "tmid": 2450000.123456,
            "rprs": 0.1234,
            "ars": 12.0,
            "inc": 88.5,
            "ecc": 0.0,
            "omega": 90.0,
            "a1": 1.0,
            "a2": 0.0,
        }
        self.errors = {
            "tmid": 0.0001,
            "rprs": 0.001,
            "ars": 0.4,
            "inc": 0.2,
            "a1": 0.1,
            "a2": 0.1,
        }
        self.time = [2450000.123456]
        self.data = [1.0]
        self.dataerr = [0.01]
        self.residuals = 0.01
        self.airmass_model = [1.0]


def aavso_json_header(output_text, header_name):
    prefix = f"#{header_name}="
    for line in output_text.splitlines():
        if line.startswith(prefix):
            return json.loads(line[len(prefix):])
    raise AssertionError(f"Missing {header_name} header")


def test_aavso_output_includes_observatory_location_headers(tmp_path):
    fit = DummyFit()
    p_dict = {
        "pName": "HAT-P-32 b",
        "sName": "HAT-P-32",
        "pPer": 2.1500082,
        "pPerUnc": 1.3e-07,
        "rprs": 0.1488623525,
        "rprsUnc": 0.0005539487,
        "aRs": 5.344,
        "aRsUnc": 0.03949,
        "inc": 88.98,
        "incUnc": 0.7602,
        "ecc": 0.159,
        "dist": 245.7,
        "pm_ra": 14.25,
        "pm_dec": -9.5,
    }
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "aavso_num": "RTZ",
        "second_obs": "",
        "obs_name": "Whipple Observatory",
        "camera": "CCD",
        "pixel_bin": "1x1",
        "exposure": 60.0,
        "lat": "+32.41638889",
        "long": "-110.73444444",
        "elev": 2616,
        "notes": "na",
        "filter": "CV",
        "filter_desc": "Clear with V zero-point",
        "wl_min": None,
        "wl_max": None,
    }

    OutputFiles(fit, p_dict, i_dict, [0.1]).aavso(
        {"ra": "", "dec": "", "x": "493", "y": "202"},
        [1.0],
        (0.1, 0.01),
        (0.2, 0.01),
        (0.3, 0.01),
        (0.4, 0.01),
        None,
    )

    output_file = tmp_path / "AAVSO_HAT-P-32 b_2020-01-01.txt"
    output_text = output_file.read_text(encoding="utf-8")

    assert "#OBSDATE=2020-01-01" in output_text
    assert "#OBSNAME=Whipple Observatory" in output_text
    assert "#OBSLAT=+32.41638889" in output_text
    assert "#OBSLON=-110.73444444" in output_text
    assert "#OBSELEV=2616" in output_text
    assert "#GAIADIST=245.7" in output_text
    assert "#GAIAPMRA=14.25" in output_text
    assert "#GAIAPMDEC=-9.5" in output_text


def test_aavso_output_omits_obsname_header_when_blank(tmp_path):
    fit = DummyFit()
    p_dict = {
        "pName": "HAT-P-32 b",
        "sName": "HAT-P-32",
        "pPer": 2.1500082,
        "pPerUnc": 1.3e-07,
        "rprs": 0.1488623525,
        "rprsUnc": 0.0005539487,
        "aRs": 5.344,
        "aRsUnc": 0.03949,
        "inc": 88.98,
        "incUnc": 0.7602,
        "ecc": 0.159,
        "dist": None,
        "pm_ra": None,
        "pm_dec": None,
    }
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "aavso_num": "RTZ",
        "second_obs": "",
        "obs_name": "",
        "camera": "CCD",
        "pixel_bin": "1x1",
        "exposure": 60.0,
        "lat": "+32.41638889",
        "long": "-110.73444444",
        "elev": 2616,
        "notes": "na",
        "filter": "CV",
        "filter_desc": "Clear with V zero-point",
        "wl_min": None,
        "wl_max": None,
    }

    OutputFiles(fit, p_dict, i_dict, [0.1]).aavso(
        {"ra": "", "dec": "", "x": "493", "y": "202"},
        [1.0],
        (0.1, 0.01),
        (0.2, 0.01),
        (0.3, 0.01),
        (0.4, 0.01),
        None,
    )

    output_file = tmp_path / "AAVSO_HAT-P-32 b_2020-01-01.txt"
    output_text = output_file.read_text(encoding="utf-8")

    assert "#OBSNAME=" not in output_text
    assert "#GAIADIST=" not in output_text
    assert "#GAIAPMRA=" not in output_text
    assert "#GAIAPMDEC=" not in output_text


def test_aid_output_includes_nextastro_comparison_metadata(tmp_path):
    fit = DummyFit()
    p_dict = {
        "pName": "HAT-P-32 b",
        "sName": "HAT-P-32",
    }
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "aavso_num": "RTZ",
        "camera": "CCD",
        "filter": "V",
        "lat": "+32.41638889",
        "long": "-110.73444444",
        "elev": 2616,
    }
    vsp_params = [{
        "time": 2450000.12345,
        "mag": 12.34,
        "mag_err": 0.05,
        "airmass": 1.234,
        "cname": "RA=10.1000000 Dec=-20.2000000",
        "cmag": 12.1,
        "cmag_err": 0.03,
        "pos": [493, 202],
        "comp_ra": 10.1,
        "comp_dec": -20.2,
        "catalog_ra": 10.10001,
        "catalog_dec": -20.20001,
        "catalog_source": "NextAstro photometry catalog",
        "is_aavso_vsp": False,
        "mag_band": "V",
        "source_id": 12345,
        "separation_arcsec": 0.2,
    }]

    AIDOutputFiles(fit, p_dict, i_dict, auid=None, chart_id=None, vsp_params=vsp_params).aavso()

    output_text = (tmp_path / "AID_AAVSO_HAT-P-32_2020-01-01.txt").read_text(encoding="utf-8")
    metadata = aavso_json_header(output_text, "COMPARISON-CATALOG-XC")

    assert metadata["source"] == "NextAstro photometry catalog"
    assert metadata["is_aavso_vsp"] is False
    assert metadata["comparison_ra_deg"] == pytest.approx(10.1)
    assert metadata["comparison_dec_deg"] == pytest.approx(-20.2)
    assert metadata["apparent_magnitude"] == pytest.approx(12.1)
    assert metadata["apparent_magnitude_error"] == pytest.approx(0.03)
    assert "HAT-P-32,2450000.12345,12.340,0.050,V,NO,STD" in output_text


def test_aid_output_floors_reported_magnitude_errors(tmp_path):
    fit = DummyFit()
    p_dict = {
        "pName": "HAT-P-32 b",
        "sName": "HAT-P-32",
    }
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "aavso_num": "RTZ",
        "camera": "CCD",
        "filter": "V",
        "lat": "+32.41638889",
        "long": "-110.73444444",
        "elev": 2616,
    }
    vsp_params = [{
        "time": 2450000.12345,
        "mag": 12.34,
        "mag_err": 0.0,
        "airmass": 1.234,
        "cname": "RA=10.1000000 Dec=-20.2000000",
        "cmag": 12.1,
        "cmag_err": 0.0,
        "pos": [493, 202],
        "catalog_source": "NextAstro photometry catalog",
        "is_aavso_vsp": False,
        "mag_band": "V",
    }]

    AIDOutputFiles(fit, p_dict, i_dict, auid=None, chart_id=None, vsp_params=vsp_params).aavso()

    output_text = (tmp_path / "AID_AAVSO_HAT-P-32_2020-01-01.txt").read_text(encoding="utf-8")
    metadata = aavso_json_header(output_text, "COMPARISON-CATALOG-XC")

    assert metadata["apparent_magnitude_error"] == pytest.approx(0.001)
    assert "HAT-P-32,2450000.12345,12.340,0.001,V,NO,STD" in output_text


def test_aid_output_skips_over_30_magnitude_rows(tmp_path):
    fit = DummyFit()
    p_dict = {
        "pName": "HAT-P-32 b",
        "sName": "HAT-P-32",
    }
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "aavso_num": "RTZ",
        "camera": "CCD",
        "filter": "V",
        "lat": "+32.41638889",
        "long": "-110.73444444",
        "elev": 2616,
    }
    vsp_params = [{
        "time": 2450000.12345,
        "mag": 99.99,
        "mag_err": 0.05,
        "airmass": 1.234,
        "cname": "Comp",
        "cmag": 12.1,
        "cmag_err": 0.03,
        "pos": [493, 202],
    }]

    AIDOutputFiles(fit, p_dict, i_dict, auid=None, chart_id=None, vsp_params=vsp_params).aavso()

    output_text = (tmp_path / "AID_AAVSO_HAT-P-32_2020-01-01.txt").read_text(encoding="utf-8")

    assert "HAT-P-32,2450000.12345" not in output_text


def test_save_comp_star_calibration_summary_writes_selected_star(tmp_path):
    summary_path = save_comp_star_calibration_summary(
        tmp_path,
        "HAT-P-32 b",
        "2026-03-09",
        "PSF photometry",
        0.0012,
        [
            {
                "label": "Comp 1",
                "position": [101, 202],
                "selected": True,
                "aggregate_score": 0.0012,
                "ensemble_score": 0.0010,
                "pairwise_median_score": 0.0011,
                "pairwise_max_score": 0.0014,
                "self_score": 0.0009,
                "valid_pair_count": 2,
                "coverage_rejected": False,
                "suitability_outlier_rejected": False,
            },
            {
                "label": "Comp 2",
                "position": [303, 404],
                "selected": False,
                "aggregate_score": 0.0031,
                "ensemble_score": 0.0028,
                "pairwise_median_score": 0.0030,
                "pairwise_max_score": 0.0035,
                "self_score": 0.0012,
                "valid_pair_count": 2,
                "coverage_rejected": False,
                "suitability_outlier_rejected": True,
            },
        ],
        0,
    )

    text = summary_path.read_text()
    assert "# Selected comparison star,1" in text
    assert "suitability_outlier_rejected" in text
    assert "Comp 1,101,202,true" in text


def test_final_planetary_params_reports_skipped_airmass_correction(tmp_path):
    fit = DummyFit()
    fit.airmass_fit_skipped = True
    fit.airmass_correction_note = "Skipped (airmass span 0.0400 <= 0.05); no airmass correction applied."
    (tmp_path / "temp").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "temp" / "FinalParams_HAT-P-32 b_2020-01-01.json"
    output_text = output_file.read_text(encoding="utf-8")

    assert "Airmass correction" in output_text
    assert "no airmass correction applied" in output_text
    assert "Airmass coefficient 1 (a1)" not in output_text


def test_final_planetary_params_reports_nextastro_variability_reference(tmp_path):
    fit = DummyFit()
    (tmp_path / "temp").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}
    vsp_params = [{
        "cname": "RA=10.1000000 Dec=-20.2000000",
        "cmag": 12.345,
        "cmag_err": 0.067,
        "pos": [493, 202],
        "comp_ra": 10.1,
        "comp_dec": -20.2,
        "catalog_source": "NextAstro photometry catalog",
        "is_aavso_vsp": False,
        "mag_band": "V",
    }]

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=vsp_params,
    )

    output_file = tmp_path / "temp" / "FinalParams_HAT-P-32 b_2020-01-01.json"
    final_params = json.loads(output_file.read_text(encoding="utf-8"))["FINAL PLANETARY PARAMETERS"]

    reference = final_params["Variable Reference Star"]
    assert "NextAstro photometry catalog" in reference
    assert "RA=10.1000000" in reference
    assert "Dec=-20.2000000" in reference
    assert "V=12.345 +/- 0.067" in reference


def test_final_planetary_params_reports_ars_and_impact_parameter_under_inclination(tmp_path):
    fit = DummyFit()
    (tmp_path / "temp").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "temp" / "FinalParams_HAT-P-32 b_2020-01-01.json"
    output_data = json.loads(output_file.read_text(encoding="utf-8"))
    final_params = output_data["FINAL PLANETARY PARAMETERS"]
    keys = list(final_params)
    inclination_index = keys.index("Orbital Inclination (inc)")

    assert keys[inclination_index + 1] == "Ratio of Distance to Stellar Radius (a/Rs)"
    assert keys[inclination_index + 2] == "Impact Parameter (b)"
    assert final_params["Ratio of Distance to Stellar Radius (a/Rs)"] == "12.0 +/- 0.4"

    expected_b, expected_b_error = fit_impact_parameter_value_error(fit)
    assert expected_b == pytest.approx(12.0 * np.cos(np.deg2rad(88.5)))
    assert final_params["Impact Parameter (b)"] == "0.314 +/- 0.043"


def test_final_planetary_params_reports_fit_uncertainties_not_prior_uncertainties(tmp_path):
    fit = DummyFit()
    (tmp_path / "temp").mkdir()

    p_dict = {
        "pName": "HAT-P-32 b",
        "midTUnc": 9.9,
        "rprsUnc": 8.8,
        "aRsUnc": 7.7,
        "incUnc": 6.6,
    }
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "temp" / "FinalParams_HAT-P-32 b_2020-01-01.json"
    final_params = json.loads(output_file.read_text(encoding="utf-8"))["FINAL PLANETARY PARAMETERS"]

    assert final_params["Mid-Transit Time (Tmid)"].endswith("+/- 0.0001 BJD_TDB")
    assert final_params["Ratio of Planet to Stellar Radius (Rp/R*)"] == "0.1234 +/- 0.001"
    assert final_params["Orbital Inclination (inc)"] == "88.5 +/- 0.2 "
    assert final_params["Ratio of Distance to Stellar Radius (a/Rs)"] == "12.0 +/- 0.4"
    assert final_params["Impact Parameter (b)"] == "0.314 +/- 0.043"


def test_final_planetary_params_can_publish_accepted_copy_to_root(tmp_path):
    fit = DummyFit()
    (tmp_path / "temp").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
        publish_to_root=True,
    )

    temp_file = tmp_path / "temp" / "FinalParams_HAT-P-32 b_2020-01-01.json"
    root_file = tmp_path / "FinalParams_HAT-P-32 b_2020-01-01.json"

    assert temp_file.exists()
    assert root_file.exists()
    assert root_file.read_text(encoding="utf-8") == temp_file.read_text(encoding="utf-8")


def test_final_planetary_params_reports_adaptive_aperture_summary(tmp_path):
    fit = DummyFit()
    (tmp_path / "temp").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}
    adaptive_summary = {
        "aperture_sigma": 2.62,
        "annulus_sigma": 9.00,
        "aperture_median": 7.98,
        "aperture_std": 0.41,
        "aperture_min": 7.12,
        "aperture_max": 8.76,
        "annulus_median": 27.43,
        "annulus_std": 1.39,
        "annulus_min": 25.11,
        "annulus_max": 30.08,
    }

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=True,
        vsp_params=[],
        comp_star=9,
        comp_coords=[1446.0, 2399.0],
        min_aper=7.98,
        min_annul=27.43,
        adaptive_summary=adaptive_summary,
    )

    output_file = tmp_path / "temp" / "FinalParams_HAT-P-32 b_2020-01-01.json"
    output_text = output_file.read_text(encoding="utf-8")

    assert "Adaptive Aperture Scale" in output_text
    assert "2.62 sigma" in output_text
    assert "Optimal Aperture" in output_text
    assert "7.98 +/- 0.41 px" in output_text
    assert "Aperture Range" in output_text
    assert "7.12 to 8.76 px" in output_text


def test_final_planetary_params_reports_transit_qc_summary(tmp_path):
    fit = DummyFit()
    fit.transit_qc = {
        "status": "pass",
        "summary": "Transit model strongly preferred over flat/null model (Delta BIC=18.40, Delta chi2=27.10).",
        "delta_bic": 18.4,
        "delta_chi2": 27.1,
        "rprs_sigma": 6.2,
        "duration_ratio": 1.05,
        "eebls_depth_snr": 5.8,
        "residual_scatter": 0.0032,
        "deviation_from_expected_value": 0.91,
        "tmid_deviation_sigma": 1.1,
        "tmid_deviation_minutes": 3.2,
        "tmid_deviation_threshold_minutes": 14.4,
        "expected_tmid_unc_minutes": 2.88,
        "rprs_deviation_fit_unc": 0.0046,
        "rprs_deviation_sigma": 0.8,
        "deviation_sigma_threshold": 5.0,
        "ktmf_metric": 4.63,
        "ktmf_contributions": [
            {
                "label": "Delta BIC",
                "available": True,
                "points": 1.25,
                "max_points": 1.40,
                "score": 0.89,
                "detail": "Delta BIC=18.40",
            },
            {
                "label": "Deviation From Expected Value",
                "available": True,
                "points": 0.91,
                "max_points": 1.00,
                "score": 0.91,
                "detail": "score=0.91, Rp/R* sigma=0.80, fit uncertainty=0.004600",
            },
        ],
        "notes": ["The transit model is strongly preferred over the flat/null model."],
    }
    (tmp_path / "temp").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "temp" / "FinalParams_HAT-P-32 b_2020-01-01.json"
    output_text = output_file.read_text(encoding="utf-8")

    assert "Transit detection QC" in output_text
    assert "Transit vs flat model" in output_text
    assert "PASS" in output_text
    assert "Delta BIC=18.40" in output_text
    assert "Residual scatter around full model fit" in output_text
    assert "Deviation From Expected Value" in output_text
    assert "3.20 minutes" not in output_text
    assert "Expected-value Tmid offset" not in output_text
    assert "Expected-value Tmid QC window" not in output_text
    assert "KTMF" in output_text
    assert "KTMF contribution 1" in output_text


def test_final_planetary_params_reports_ktmf_decision_details(tmp_path):
    fit = DummyFit()
    fit.transit_qc = {
        "status": "pass",
        "summary": "Transit model strongly preferred over flat/null model.",
        "ktmf_metric": 4.63,
        "delta_bic": 18.4,
        "delta_chi2": 27.1,
        "ktmf_contributions": [
            {
                "label": "Model Evidence",
                "available": True,
                "points": 0.74,
                "max_points": 0.80,
                "score": 0.93,
                "detail": "Delta BIC=18.40",
            }
        ],
    }
    (tmp_path / "temp").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}
    photometry_info = {
        "selection_basis": "comparison_field_retry",
        "selection_metric": "ktmf",
        "comp_star_num": 2,
        "comparison_ktmf_metric": 4.60,
        "comparison_eebls_snr": 5.2,
        "comparison_transit_delta_bic": 18.4,
        "selected_comparison_selection_reason": "selected: highest KTMF among candidates",
        "selected_comparison_ktmf_contributions": [
            {
                "label": "Residual Scatter Around Full Model Fit",
                "available": True,
                "points": 0.63,
                "max_points": 0.70,
                "score": 0.90,
                "detail": "0.3500%",
            }
        ],
        "comparison_fit_attempt_summaries": [
            {
                "rank": 1,
                "comp_index": 0,
                "label": "Comp 1",
                "selected": False,
                "selection_reason": "not selected: KTMF 3.20/5.00 was lower than the selected 4.60/5.00",
                "ktmf_metric": 3.2,
                "transit_delta_bic": 8.1,
                "eebls_snr": 4.2,
                "transit_qc_status": "marginal",
            },
            {
                "rank": 2,
                "comp_index": 1,
                "label": "Comp 2",
                "selected": True,
                "selection_reason": "selected: highest KTMF among candidates",
                "ktmf_metric": 4.6,
                "transit_delta_bic": 18.4,
                "eebls_snr": 5.2,
                "transit_qc_status": "pass",
            },
        ],
    }

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=True,
        vsp_params=[],
        comp_star=2,
        comp_coords=[300.5, 400.5],
        min_aper=7.5,
        min_annul=22.5,
        photometry_info=photometry_info,
    )

    output_file = tmp_path / "temp" / "FinalParams_HAT-P-32 b_2020-01-01.json"
    output_data = json.loads(output_file.read_text(encoding="utf-8"))
    final_params = output_data["FINAL PLANETARY PARAMETERS"]

    assert final_params["KTMF target-fit decision"] == "PASS: KTMF=4.63 / 5.00"
    assert final_params["KTMF comparison selection mode"] == "basis=comparison_field_retry, metric=ktmf"
    assert "selected: highest KTMF" in final_params["KTMF selected comparison decision"]
    assert "Comp 1" in final_params["KTMF comparison candidate 1"]
    assert "not selected: KTMF" in final_params["KTMF comparison candidate 1"]
    assert "Residual Scatter Around Full Model Fit" in final_params["KTMF selected comparison contribution 1"]


def test_final_planetary_params_reports_absolute_fit_quality(tmp_path):
    fit = DummyFit()
    fit.data = np.array([1.0, 1.02, 0.98, 1.01, 0.99, 1.0])
    fit.model = np.ones(6, dtype=float)
    fit.residuals = fit.data - fit.model
    fit.dataerr = np.full(6, 0.01, dtype=float)
    fit.time = np.arange(6, dtype=float)
    fit.airmass_model = np.ones(6, dtype=float)
    fit.bounds = {"tmid": [0, 1], "rprs": [0, 1], "a1": [0, 2]}
    (tmp_path / "temp").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "temp" / "FinalParams_HAT-P-32 b_2020-01-01.json"
    output_data = json.loads(output_file.read_text(encoding="utf-8"))
    final_params = output_data["FINAL PLANETARY PARAMETERS"]

    assert final_params["Fit quality reduced chi-square"] == "3.333"
    assert final_params["Fit quality chi-square"] == "10.00"
    assert final_params["Fit quality degrees of freedom"] == "3"
    assert final_params["Fit quality RMS residual"] == "1.2910 %"
    assert final_params["Fit quality median absolute normalized residual"] == "1.00 sigma"
    assert final_params["Fit quality RMS residual / median uncertainty"] == "1.29"
    assert final_params["Fit quality point count"] == "6"


def test_aavso_output_writes_zero_airmass_terms_when_correction_is_skipped(tmp_path):
    fit = DummyFit()
    fit.airmass_fit_skipped = True
    fit.airmass_correction_note = "Skipped (input AAVSO file already reports AIRMASS, AIRMASS CORRECTION FUNCTION); no airmass correction applied."

    p_dict = {
        "pName": "HAT-P-32 b",
        "sName": "HAT-P-32",
        "pPer": 2.1500082,
        "pPerUnc": 1.3e-07,
        "rprs": 0.1488623525,
        "rprsUnc": 0.0005539487,
        "aRs": 5.344,
        "aRsUnc": 0.03949,
        "inc": 88.98,
        "incUnc": 0.7602,
        "ecc": 0.159,
        "dist": None,
        "pm_ra": None,
        "pm_dec": None,
    }
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "aavso_num": "RTZ",
        "second_obs": "",
        "obs_name": "",
        "camera": "CCD",
        "pixel_bin": "1x1",
        "exposure": 60.0,
        "lat": "+32.41638889",
        "long": "-110.73444444",
        "elev": 2616,
        "notes": "na",
        "filter": "CV",
        "filter_desc": "Clear with V zero-point",
        "wl_min": None,
        "wl_max": None,
    }

    OutputFiles(fit, p_dict, i_dict, [0.1]).aavso(
        {"ra": "", "dec": "", "x": "493", "y": "202"},
        [1.0],
        (0.1, 0.01),
        (0.2, 0.01),
        (0.3, 0.01),
        (0.4, 0.01),
        None,
    )

    output_file = tmp_path / "AAVSO_HAT-P-32 b_2020-01-01.txt"
    output_text = output_file.read_text(encoding="utf-8")

    assert "Am1=0 +/- 0" in output_text
    assert "Am2=0 +/- 0" in output_text
    assert output_text.strip().endswith("1.0")


def test_aavso_output_includes_extended_diagnostic_comment_headers(tmp_path):
    fit = DummyFit()
    fit.data = np.array([1.0, 1.02, 0.98, 1.01, 0.99, 1.0])
    fit.model = np.ones(6, dtype=float)
    fit.residuals = fit.data - fit.model
    fit.dataerr = np.full(6, 0.01, dtype=float)
    fit.time = np.arange(6, dtype=float) + 2450000.0
    fit.airmass_model = np.ones(6, dtype=float)
    fit.bounds = {"tmid": [0, 1], "rprs": [0, 1], "a1": [0, 2]}
    fit.transit_qc = {
        "computed": True,
        "status": "pass",
        "summary": "Transit model strongly preferred over flat/null model.",
        "delta_bic": 18.4,
        "delta_chi2": 27.1,
        "residual_scatter": 0.0032,
        "rprs_sigma": 6.2,
        "duration_ratio": 1.05,
        "eebls_depth_snr": 5.8,
        "deviation_from_expected_value": 0.91,
        "tmid_deviation_minutes": 3.2,
        "tmid_deviation_sigma": 1.1,
        "rprs_deviation_sigma": 0.8,
        "ktmf_metric": 4.63,
        "ktmf_contributions": [
            {
                "label": "Model Evidence",
                "available": True,
                "points": 0.74,
                "max_points": 0.80,
                "score": 0.93,
                "detail": "Delta BIC=18.40",
            }
        ],
    }
    fit.frame_filter_diagnostics = [
        {
            "stage": "Final-fit phase residual clip",
            "input_point_count": 4,
            "kept_point_count": 3,
            "dropped_point_count": 1,
            "dropped_ranges": [{"start": 2450000.2, "end": 2450000.2, "count": 1}],
        }
    ]
    p_dict = {
        "pName": "HAT-P-32 b",
        "sName": "HAT-P-32",
        "pPer": 2.1500082,
        "pPerUnc": 1.3e-07,
        "rprs": 0.1488623525,
        "rprsUnc": 0.0005539487,
        "aRs": 5.344,
        "aRsUnc": 0.03949,
        "inc": 88.98,
        "incUnc": 0.7602,
        "ecc": 0.159,
        "dist": 245.7,
        "pm_ra": 14.25,
        "pm_dec": -9.5,
    }
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "aavso_num": "RTZ",
        "second_obs": "",
        "obs_name": "",
        "camera": "CCD",
        "pixel_bin": "1x1",
        "exposure": 60.0,
        "lat": "+32.41638889",
        "long": "-110.73444444",
        "elev": 2616,
        "notes": "na",
        "filter": "CV",
        "filter_desc": "Clear with V zero-point",
        "wl_min": None,
        "wl_max": None,
    }
    photometry_info = {
        "comp_star_num": 2,
        "comp_star_coords": [300.5, 400.5],
        "min_aperture": 7.5,
        "min_annulus": 22.5,
        "aperture_index": 1,
        "annulus_index": 2,
        "calibration_field_score": 0.0042,
        "selection_basis": "comparison_field",
        "selection_metric": "ktmf",
        "comparison_ktmf_metric": 4.6,
        "comparison_eebls_snr": 5.2,
        "comparison_transit_delta_bic": 18.4,
        "selected_comparison_selection_reason": "selected: highest KTMF among candidates",
        "selected_comparison_ktmf_contributions": [
            {
                "label": "Residual Scatter Around Full Model Fit",
                "available": True,
                "points": 0.63,
                "max_points": 0.70,
                "score": 0.90,
                "detail": "0.3500%",
            }
        ],
        "comparison_fit_attempt_summaries": [
            {
                "rank": 1,
                "comp_index": 0,
                "label": "Comp 1",
                "selected": False,
                "selection_reason": "not selected: KTMF 3.20/5.00 was lower than the selected 4.60/5.00",
                "ktmf_metric": 3.2,
                "transit_delta_bic": 8.1,
                "eebls_snr": 4.2,
                "transit_qc_status": "marginal",
                "ktmf_contributions": [],
            },
            {
                "rank": 2,
                "comp_index": 1,
                "label": "Comp 2",
                "selected": True,
                "selection_reason": "selected: highest KTMF among candidates",
                "ktmf_metric": 4.6,
                "transit_delta_bic": 18.4,
                "eebls_snr": 5.2,
                "transit_qc_status": "pass",
                "ktmf_contributions": [
                    {
                        "label": "Residual Scatter Around Full Model Fit",
                        "available": True,
                        "points": 0.63,
                        "max_points": 0.70,
                        "score": 0.90,
                        "detail": "0.3500%",
                    }
                ],
            },
        ],
        "reuse_selected_full_reduction_fit": True,
        "selected_source_indices": np.array([0, 2, 3]),
        "selected_fit_good_times": np.array([2450000.0, 2450000.1, 2450000.2]),
        "adaptive_summary": {
            "aperture_sigma": 2.62,
            "annulus_sigma": 9.00,
            "frame_sigma": np.array([2.0, 2.1, 2.2]),
            "fwhm_series": np.array([4.7, 4.8, 4.9]),
            "sky_inner_series": np.array([12.0, 12.1, 12.2]),
            "sky_outer_series": np.array([18.0, 18.1, 18.2]),
            "sky_pixel_series": np.array([200.0, 201.0, 202.0]),
            "aperture_median": 7.98,
            "aperture_std": 0.41,
            "aperture_min": 7.12,
            "aperture_max": 8.76,
            "annulus_median": 27.43,
            "annulus_std": 1.39,
            "annulus_min": 25.11,
            "annulus_max": 30.08,
        },
    }
    comp_star_header = {"ra": "10.1", "dec": "-20.2", "x": "493", "y": "202"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).aavso(
        comp_star_header,
        np.ones(6, dtype=float),
        (0.1, 0.01),
        (0.2, 0.01),
        (0.3, 0.01),
        (0.4, 0.01),
        "abc123",
        photometry_info=photometry_info,
        frame_filtering_info={
            "initial_frame_count": 5,
            "after_missing_wcs_filter_frame_count": 4,
            "final_prephotometry_frame_count": 3,
            "ignore_header_wcs": False,
            "bad_wcs_threshold_percent": 3.0,
            "pointing_rejection_sigma": 3.0,
            "dropped_missing_wcs_files": [tmp_path / "missing_wcs.fits"],
            "dropped_pointing_files": [tmp_path / "bad_pointing.fits"],
        },
        astrometry_info={
            "wcs_file": tmp_path / "wcs.fits",
            "coordinate_source": "wcs",
            "target_ra_dec_deg": [10.0, -20.0],
            "comparison_ra_dec_deg": [[10.1, -20.2]],
        },
        bad_pixel_info={
            "enabled": True,
            "detected": True,
            "bad_pixel_count": 3,
            "frame_count": 10,
            "required_count": 4,
            "minimum_fraction": 0.3,
            "counts_path": tmp_path / "temp" / "BadPixelDetectionCounts.fits",
            "mask_path": tmp_path / "temp" / "BadPixelMask.fits",
        },
    )

    output_text = (tmp_path / "AAVSO_HAT-P-32 b_2020-01-01.txt").read_text(encoding="utf-8")

    results = aavso_json_header(output_text, "RESULTS-XC")
    assert "a/R*" in results
    assert "Impact Parameter (b)" in results
    assert results["Transit depth (Rp/R*)^2"]["units"] == "percent"
    assert results["Residual scatter around full model fit"]["value"] == "0.32"

    qc = aavso_json_header(output_text, "QC-XC")
    assert qc["status"] == "pass"
    assert qc["ktmf_metric"] == pytest.approx(4.63)
    assert qc["ktmf_contributions"][0]["label"] == "Model Evidence"

    fit_quality = aavso_json_header(output_text, "FIT_QUALITY-XC")
    assert fit_quality["reduced_chi_square"] == pytest.approx(10.0 / 3.0)
    assert fit_quality["chi_square"] == pytest.approx(10.0)
    assert fit_quality["degrees_of_freedom"] == 3
    assert fit_quality["median_absolute_normalized_residual"] == pytest.approx(1.0)

    ktmf_decision = aavso_json_header(output_text, "KTMF_DECISION-XC")
    assert ktmf_decision["target_fit"]["ktmf_metric"] == pytest.approx(4.63)
    assert ktmf_decision["comparison_selection"]["basis"] == "comparison_field"
    assert ktmf_decision["comparison_selection"]["metric"] == "ktmf"
    assert ktmf_decision["comparison_selection"]["selected"]["selection_reason"] == "selected: highest KTMF among candidates"
    assert ktmf_decision["comparison_selection"]["candidate_count"] == 2
    assert ktmf_decision["comparison_selection"]["candidates"][0]["selection_reason"].startswith("not selected: KTMF")

    photometry = aavso_json_header(output_text, "PHOTOMETRY-XC")
    assert photometry["selected_comparison_star"] == 2
    assert photometry["comparison_field_score_percent"] == pytest.approx(0.42)
    assert photometry["reused_selected_full_reduction_fit"] is True

    aperture = aavso_json_header(output_text, "APERTURE-XC")
    assert aperture["adaptive"] is True
    assert aperture["aperture_sigma"] == pytest.approx(2.62)
    assert aperture["fwhm_px"]["median"] == pytest.approx(4.8)

    frame_filtering = aavso_json_header(output_text, "FRAME_FILTERING-XC")
    assert frame_filtering["missing_wcs_rejections"]["files"] == ["missing_wcs.fits"]
    assert frame_filtering["pointing_rejections"]["files"] == ["bad_pointing.fits"]
    assert frame_filtering["lightcurve_dropped_point_count"] == 1

    astrometry = aavso_json_header(output_text, "ASTROMETRY-XC")
    assert astrometry["wcs_file"] == "wcs.fits"
    assert astrometry["comparison_star_aavso_header"] == comp_star_header

    bad_pixel = aavso_json_header(output_text, "BAD_PIXEL-XC")
    assert bad_pixel["enabled"] is True
    assert bad_pixel["bad_pixel_count"] == 3
    assert bad_pixel["counts_path"] == "BadPixelDetectionCounts.fits"
