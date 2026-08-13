import json
from types import SimpleNamespace

import numpy as np
import pytest

from exotic.output_files import (
    AREA_DEPTH_LABEL,
    OBSERVABLE_DEPTH_DELTA_LABEL,
    OBSERVABLE_DEPTH_LABEL,
    PRIOR_OBSERVABLE_DEPTH_LABEL,
    AIDOutputFiles,
    OutputFiles,
    aid_comparison_coordinate_headers,
    aavso_dicts,
    build_aavso_qc_metadata,
    fit_empirical_transit_uncertainty,
    fit_impact_parameter_value_error,
    differential_magnitude_series_from_fit,
    magnitude_series_from_fit,
    save_comp_star_calibration_summary,
    write_differential_magnitude_csv,
)
from exotic.transit_depth import (
    fit_transit_depth_summary,
    observable_depth_percent,
    radius_ratio_area_depth_percent,
)


def test_stellar_variability_differential_magnitudes_never_apply_airmass_correction():
    fit = SimpleNamespace(
        stellar_variability_only=True,
        time=np.array([2460000.1, 2460000.2, 2460000.3]),
        data=np.ones(3),
        dataerr=np.full(3, 0.01),
        airmass=np.array([1.1, 1.3, 1.5]),
        airmass_model=np.array([0.8, 1.0, 1.2]),
        transit=np.ones(3),
        stellar_variability_target_flux=np.array([80.0, 100.0, 120.0]),
        stellar_variability_comp_flux=np.full(3, 100.0),
        stellar_variability_target_flux_error=np.ones(3),
        stellar_variability_comp_flux_error=np.ones(3),
    )

    series = differential_magnitude_series_from_fit(fit)

    np.testing.assert_allclose(
        series['magnitude'],
        -2.5 * np.log10(fit.stellar_variability_target_flux / fit.stellar_variability_comp_flux),
    )
    assert series['airmass_corrected'] is False


def test_differential_csv_does_not_require_apparent_magnitude_calibration(tmp_path):
    fit = SimpleNamespace(
        stellar_variability_only=True,
        time=np.array([2460000.1, 2460000.2]),
        data=np.ones(2),
        dataerr=np.full(2, 0.01),
        airmass=np.array([1.1, 1.2]),
        airmass_model=np.array([0.9, 1.1]),
        transit=np.ones(2),
        stellar_variability_target_flux=np.array([500.0, 550.0]),
        stellar_variability_comp_flux=np.array([1000.0, 1000.0]),
        stellar_variability_target_flux_error=np.full(2, 2.0),
        stellar_variability_comp_flux_error=np.full(2, 3.0),
    )

    output_path = write_differential_magnitude_csv(
        fit,
        tmp_path,
        'Variable Star',
        observation_date='2026-08-02',
        observed_filter='V',
    )

    output_text = output_path.read_text(encoding='utf-8')
    assert '# AIRMASS_CORRECTION=NO' in output_text
    assert 'Differential Magnitude' in output_text
    assert 'Apparent' not in output_text
    assert ', 0.7526, 0.0054, V, ' in output_text
    assert ', 0.6491, 0.0051, V, ' in output_text


class DummyFit:
    def __init__(self):
        self.parameters = {
            "tmid": 2450000.123456,
            "rprs": 0.1234,
            "ars": 12.0,
            "per": 2.15,
            "inc": 88.5,
            "ecc": 0.0,
            "omega": 90.0,
            "u0": 0.0,
            "u1": 0.0,
            "u2": 0.0,
            "u3": 0.0,
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
        self.time_upsample = np.linspace(2450000.0, 2450000.2, 128)
        self.data = [1.0]
        self.dataerr = [0.01]
        self.residuals = 0.01
        self.airmass_model = [1.0]
        self.transit = [1.0 - self.parameters["rprs"] ** 2]
        self.transit_upsample = np.ones_like(self.time_upsample)
        self.transit_upsample[64] = 1.0 - self.parameters["rprs"] ** 2
        self.prior = {
            "tmid": 2450000.123456,
            "rprs": 0.1,
            "ars": 12.0,
            "per": 2.15,
            "inc": 88.5,
            "ecc": 0.0,
            "omega": 90.0,
            "u0": 0.0,
            "u1": 0.0,
            "u2": 0.0,
            "u3": 0.0,
        }


def test_prior_depth_uses_available_gj436_geometry():
    fit = DummyFit()
    prior = {
        "Published Mid-Transit Time": 2454510.80162,
        "Rp/Rs": 0.0822,
        "a/Rs": 13.73,
        "Orbital Period (days)": 2.64388312,
        "Orbital Inclination (deg)": 86.44,
        "Orbital Eccentricity": 0.13827,
        "Argument of Periastron (deg)": 351.0,
        "u0": 0.0,
        "u1": 0.0,
        "u2": 0.0,
        "u3": 0.0,
    }

    summary = fit_transit_depth_summary(
        fit,
        prior_parameters=prior,
        prior_errors={
            "Rp/Rs Uncertainty": 0.001,
            "a/Rs Uncertainty": 0.46,
            "Orbital Inclination Uncertainty": 0.17,
        },
    )

    assert summary["prior_observable_depth"] == pytest.approx(0.675684, abs=1.0e-5)
    assert summary["prior_observable_depth_error"] == pytest.approx(0.01644, abs=1.0e-6)


def test_prior_depth_respects_inclination_for_non_transiting_geometry():
    fit = DummyFit()
    prior = {
        "tmid": 2454510.80162,
        "rprs": 0.0822,
        "ars": 13.73,
        "per": 2.64388312,
        "inc": 0.0,
        "ecc": 0.13827,
        "omega": 351.0,
        "u0": 0.0,
        "u1": 0.0,
        "u2": 0.0,
        "u3": 0.0,
    }

    summary = fit_transit_depth_summary(fit, prior_parameters=prior)

    assert summary["prior_observable_depth"] == pytest.approx(0.0)


def test_final_params_writes_stellar_variability_only_payload(tmp_path):
    (tmp_path / "working_artifacts").mkdir()
    fit = SimpleNamespace(
        stellar_variability_only=True,
        time=np.arange(6, dtype=float),
        stellar_variability_scatter=0.00123,
        stellar_variability_transit_exclusion={
            'rejected_point_count': 2,
            'duration_days': 0.083,
            'note': 'Excluded synthetic transit-window points.',
        },
        airmass_fit_skipped=True,
        airmass_correction_note=(
            "Skipped in stellar-variability-only mode; no transit/systematics model was fit."
        ),
    )
    p_dict = {'pName': 'Syntheticb'}
    i_dict = {'save': str(tmp_path), 'date': '2020-01-01'}

    OutputFiles(fit, p_dict, i_dict, [0.083]).final_planetary_params(
        phot_opt=True,
        vsp_params=None,
        comp_star=2,
        comp_coords=[10, 20],
        min_aper=0,
        min_annul=15,
        photometry_info={'noise_budget_summary': 'gain only'},
        publish_to_root=True,
    )

    temp_file = next((tmp_path / "working_artifacts").glob("FinalParams_Syntheticb_2020-01-01.json"))
    root_file = tmp_path / temp_file.name
    payload = json.loads(temp_file.read_text())

    params = payload["FINAL STELLAR VARIABILITY PARAMETERS"]
    assert params["Analysis Mode"] == "Stellar variability only"
    assert params["Transit model fitting"] == "Skipped"
    assert params["Predicted in-transit points excluded"] == "2"
    assert params["Residual scatter around flat stellar-variability model"] == "0.1230 %"
    assert params["Stellar Variability Reference Star"] == "#2 - [10, 20]"
    assert params["Optimal Method"] == "PSF photometry"
    assert root_file.exists()


def test_final_params_describes_calibrated_stellar_variability_ensemble(tmp_path):
    (tmp_path / "working_artifacts").mkdir()
    fit = SimpleNamespace(
        stellar_variability_only=True,
        time=np.arange(3, dtype=float),
        stellar_variability_scatter=0.001,
        stellar_variability_transit_exclusion={'rejected_point_count': 0},
    )
    vsp_params = [{
        'time': 2460000.1,
        'mag': 12.3,
        'mag_err': 0.01,
        'cname': 'ENSEMBLE (2 stars)',
        'ensemble_reference': True,
        'ensemble_member_count': 2,
        'ensemble_member_labels': ['C1', 'C2'],
        'mag_band': 'V',
    }]
    p_dict = {'pName': 'Syntheticb'}
    i_dict = {'save': str(tmp_path), 'date': '2020-01-01'}

    OutputFiles(fit, p_dict, i_dict, []).final_planetary_params(
        phot_opt=True,
        vsp_params=vsp_params,
        comp_star='ensemble',
        comp_coords=None,
        min_aper=0,
        min_annul=15,
    )

    output_path = next((tmp_path / "working_artifacts").glob("FinalParams_Syntheticb_2020-01-01.json"))
    params = json.loads(output_path.read_text())["FINAL STELLAR VARIABILITY PARAMETERS"]
    assert params["Stellar Variability Reference Star"] == "ensemble"
    assert params["Variable Reference Star"] == (
        "Calibrated comparison-star ensemble (2 stars): C1, C2"
    )
    assert "calibrated comparison-star ensemble" in params["Variable Reference Measurement"]


def test_final_lightcurve_writes_stellar_variability_magnitudes(tmp_path):
    (tmp_path / "working_artifacts").mkdir()
    fit = SimpleNamespace(
        stellar_variability_only=True,
        stellar_variability_params=[
            {
                "time": 2461229.89899,
                "mag": 13.7378,
                "mag_err": 0.0042,
                "differential_mag": 1.2378,
                "differential_mag_err": 0.0021,
                "mag_band": "r",
                "airmass": 1.193135,
            },
            {
                "time": 2461229.90109,
                "mag": 13.7401,
                "mag_err": 0.0044,
                "differential_mag": 1.2401,
                "differential_mag_err": 0.0022,
                "mag_band": "r",
                "airmass": 1.1984942,
            },
        ],
    )
    p_dict = {'pName': 'WASP-194 b', 'sName': 'WASP-194'}
    i_dict = {'save': str(tmp_path), 'date': '2026-07-08', 'filter': 'SR'}

    OutputFiles(fit, p_dict, i_dict, []).final_lightcurve(np.array([]))

    output_text = next((tmp_path / "working_artifacts").glob("FinalLightCurve_WASP-194b_2026-07-08.csv")).read_text()

    assert "# FINAL STELLAR VARIABILITY TIMESERIES OF WASP-194" in output_text
    assert "Apparent Magnitude,Apparent Magnitude Uncertainty" in output_text
    assert "Differential Magnitude,Differential Magnitude Uncertainty" in output_text
    assert "2461229.89899, 13.7378, 0.0042, 1.2378, 0.0021, r, 1.193135" in output_text
    assert "Flux" not in output_text


def test_final_lightcurve_adds_transit_apparent_magnitude_columns_when_calibrated(tmp_path):
    (tmp_path / "working_artifacts").mkdir()
    fit = SimpleNamespace(
        time=np.array([2461229.9, 2461229.91]),
        detrended=np.array([1.0, 0.99]),
        dataerr=np.array([0.001, 0.001]),
        airmass_model=np.ones(2),
        transit=np.array([1.0, 0.99]),
        stellar_variability_params=[
            {"time": 2461229.9, "mag": 13.739, "mag_err": 0.001, "mag_band": "r"},
            {"time": 2461229.91, "mag": 13.741, "mag_err": 0.002, "mag_band": "r"},
        ],
    )
    p_dict = {'pName': 'WASP-194 b', 'sName': 'WASP-194'}
    i_dict = {'save': str(tmp_path), 'date': '2026-07-08', 'filter': 'SR'}

    OutputFiles(fit, p_dict, i_dict, []).final_lightcurve(np.array([0.1, 0.2]))

    output_text = next((tmp_path / "working_artifacts").glob("FinalLightCurve_WASP-194b_2026-07-08.csv")).read_text()

    assert "Differential Magnitude,Differential Magnitude Uncertainty" in output_text
    assert "Apparent Magnitude,Apparent Magnitude Uncertainty,Band" in output_text
    assert "2461229.9, 0.1, 1.0, 0.001, 1.0, 1.0, -0.0000, 0.0011, 13.7400" in output_text
    assert output_text.rstrip().endswith(", r")


def test_final_lightcurve_keeps_differential_magnitude_when_apparent_calibration_is_unavailable(tmp_path):
    (tmp_path / "working_artifacts").mkdir()
    fit = SimpleNamespace(
        time=np.array([2461229.9]),
        data=np.array([0.8]),
        detrended=np.array([0.8]),
        dataerr=np.array([0.008]),
        airmass_model=np.ones(1),
        transit=np.ones(1),
    )
    p_dict = {'pName': 'Uncalibrated b', 'sName': 'Uncalibrated'}
    i_dict = {'save': str(tmp_path), 'date': '2026-07-08', 'filter': 'V'}

    OutputFiles(fit, p_dict, i_dict, []).final_lightcurve(np.array([0.1]))

    output_text = next(
        (tmp_path / "working_artifacts").glob("FinalLightCurve_Uncalibratedb_2026-07-08.csv")
    ).read_text()
    expected_differential = -2.5 * np.log10(0.8)
    assert f"{expected_differential:.4f}" in output_text
    assert ", na, na, V" in output_text


def test_magnitude_series_preserves_raw_ratio_for_later_apparent_recalibration():
    target_flux = np.array([500.0, 550.0])
    reference_flux = np.full(2, 1000.0)
    target_error = np.full(2, 2.0)
    reference_error = np.full(2, 3.0)
    differential_mag = -2.5 * np.log10(target_flux / reference_flux)
    magnitude_factor = 2.5 / np.log(10.0)
    differential_error = magnitude_factor * np.sqrt(
        (target_error / target_flux) ** 2
        + (reference_error / reference_flux) ** 2
    )
    fit = SimpleNamespace(
        time=np.array([2461229.9, 2461229.91]),
        data=np.array([1.0, 1.1]),
        dataerr=np.full(2, 0.001),
        detrended=np.array([1.0, 1.1]),
        airmass=np.array([1.1, 1.2]),
        airmass_model=np.ones(2),
        transit=np.ones(2),
        stellar_variability_target_flux=target_flux,
        stellar_variability_comp_flux=reference_flux,
        stellar_variability_target_flux_error=target_error,
        stellar_variability_comp_flux_error=reference_error,
        stellar_variability_params=[{
            "time": 2461229.9,
            "mag": 12.0 + differential_mag[0],
            "mag_err": np.hypot(0.02, differential_error[0]),
            "differential_mag": differential_mag[0],
            "differential_mag_err": differential_error[0],
            "cmag": 12.0,
            "cmag_err": 0.02,
            "mag_band": "V",
        }],
    )

    series = magnitude_series_from_fit(fit, apply_airmass_correction=False)

    np.testing.assert_allclose(series['differential_magnitude'], differential_mag)
    np.testing.assert_allclose(series['differential_magnitude_error'], differential_error)
    np.testing.assert_allclose(series['apparent_magnitude'], 12.0 + differential_mag)
    np.testing.assert_allclose(
        series['apparent_magnitude_error'],
        np.hypot(0.02, differential_error),
    )


def test_calibrated_ensemble_keeps_apparent_magnitudes_independent_of_raw_differential():
    target_flux = np.array([1000.0, 1010.0])
    raw_reference_flux = np.full(2, 375.0)
    calibrated_magnitude = np.array([11.25, 11.27])
    calibrated_error = np.array([0.02, 0.021])
    expected_differential = -2.5 * np.log10(target_flux / raw_reference_flux)
    fit = SimpleNamespace(
        stellar_variability_only=True,
        time=np.array([2461229.9, 2461229.91]),
        data=np.ones(2),
        detrended=np.ones(2),
        dataerr=np.full(2, 0.001),
        airmass=np.array([1.1, 1.2]),
        airmass_model=np.ones(2),
        transit=np.ones(2),
        # The calibrated ensemble's normalized fitting reference remains
        # separate from the raw instrumental reference used for DIFFMAG.
        stellar_variability_target_flux=target_flux,
        stellar_variability_comp_flux=target_flux.copy(),
        stellar_variability_target_flux_error=np.ones(2),
        stellar_variability_comp_flux_error=np.ones(2),
        differential_magnitude_target_flux=target_flux,
        differential_magnitude_reference_flux=raw_reference_flux,
        differential_magnitude_target_flux_error=np.ones(2),
        differential_magnitude_reference_flux_error=np.ones(2),
        stellar_variability_ensemble_magnitudes=calibrated_magnitude,
        stellar_variability_ensemble_magnitude_errors=calibrated_error,
        stellar_variability_params=[
            {
                'time': 2461229.9,
                'mag': calibrated_magnitude[0],
                'mag_err': calibrated_error[0],
                'differential_mag': expected_differential[0],
                'differential_mag_err': 0.002,
                'mag_band': 'V',
            },
            {
                'time': 2461229.91,
                'mag': calibrated_magnitude[1],
                'mag_err': calibrated_error[1],
                'differential_mag': expected_differential[1],
                'differential_mag_err': 0.002,
                'mag_band': 'V',
            },
        ],
    )

    series = magnitude_series_from_fit(fit, apply_airmass_correction=False)

    np.testing.assert_allclose(series['differential_magnitude'], expected_differential)
    np.testing.assert_allclose(series['apparent_magnitude'], calibrated_magnitude)
    np.testing.assert_allclose(series['apparent_magnitude_error'], calibrated_error)
    assert series['apparent_calibrated'] is True


def aavso_json_header(output_text, header_name):
    prefix = f"#{header_name}="
    for line in output_text.splitlines():
        if line.startswith(prefix):
            return json.loads(line[len(prefix):])
    raise AssertionError(f"Missing {header_name} header")


def test_observable_depth_is_separate_from_area_depth_for_grazing_geometry():
    parameters = {
        "tmid": 0.0,
        "rprs": 0.2,
        "per": 3.0,
        "ars": 10.0,
        "inc": np.degrees(np.arccos(1.1 / 10.0)),
        "ecc": 0.0,
        "omega": 90.0,
        "u0": 0.0,
        "u1": 0.0,
        "u2": 0.0,
        "u3": 0.0,
    }
    errors = {"rprs": 0.01, "ars": 0.1, "inc": 0.1}

    area_depth, area_error = radius_ratio_area_depth_percent(parameters["rprs"], errors["rprs"])
    observable_depth, observable_error = observable_depth_percent(parameters, errors)

    assert area_depth == pytest.approx(4.0)
    assert area_error == pytest.approx(0.4)
    assert 0.0 < observable_depth < area_depth
    assert observable_error > 0.0


def test_aavso_output_includes_observatory_location_headers(tmp_path):
    fit = DummyFit()
    fit.stellar_variability_target_flux = np.array([500.0])
    fit.stellar_variability_comp_flux = np.array([1000.0])
    fit.stellar_variability_target_flux_error = np.array([2.0])
    fit.stellar_variability_comp_flux_error = np.array([3.0])
    differential_mag = float(-2.5 * np.log10(0.5))
    differential_error = float(
        (2.5 / np.log(10.0)) * np.hypot(2.0 / 500.0, 3.0 / 1000.0)
    )
    fit.stellar_variability_params = [{
        "time": fit.time[0],
        "mag": 12.0 + differential_mag,
        "mag_err": np.hypot(0.02, differential_error),
        "differential_mag": differential_mag,
        "differential_mag_err": differential_error,
        "cmag": 12.0,
        "cmag_err": 0.02,
        "mag_band": "V",
    }]
    p_dict = {
        "pName": "HAT-P-32b",
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
    final_plot_source = tmp_path / "FinalLightCurve_HAT-P-32b_2020-01-01.png"
    final_plot_source.write_bytes(b"final lightcurve")
    diagnostics_dir = tmp_path / "Diagnostics"
    diagnostics_dir.mkdir()
    diagnostic_sources = [
        diagnostics_dir / filename
        for filename in (
            "FinalTriangle_HAT-P-32b_2020-01-01.png",
            "Triangle_HAT-P-32b_2020-01-01.png",
            "ZoomedTrianglePlot_HAT-P-32b_2020-01-01.png",
            "KTMF_QC_HAT-P-32b_2020-01-01.png",
            "KTMF_QC_HAT-P-32b_2020-01-01.pdf",
            "PriorPosteriorComparison_HAT-P-32b_2020-01-01.png",
            "PriorPosteriorComparison_HAT-P-32b_2020-01-01.pdf",
        )
    ]
    for diagnostic_source in diagnostic_sources:
        diagnostic_source.write_bytes(diagnostic_source.name.encode("utf-8"))

    OutputFiles(fit, p_dict, i_dict, [0.1]).aavso(
        {"ra": "", "dec": "", "x": "493", "y": "202"},
        [1.0],
        (0.1, 0.01),
        (0.2, 0.01),
        (0.3, 0.01),
        (0.4, 0.01),
        None,
    )

    output_file = tmp_path / "AAVSO_Files" / "AAVSO_HAT-P-32b_2020-01-01.txt"
    output_text = output_file.read_text(encoding="utf-8")
    assert (
        tmp_path / "AAVSO_Files" / final_plot_source.name
    ).read_bytes() == b"final lightcurve"
    for diagnostic_source in diagnostic_sources:
        assert (
            tmp_path / "AAVSO_Files" / diagnostic_source.name
        ).read_bytes() == diagnostic_source.name.encode("utf-8")

    assert "#OBSDATE=2020-01-01" in output_text
    assert "#EXOPLANET_NAME=HAT-P-32 b" in output_text
    assert "#OBSNAME=Whipple Observatory" in output_text
    assert "#OBSLAT=+32.41638889" in output_text
    assert "#OBSLON=-110.73444444" in output_text
    assert "#OBSELEV=2616" in output_text
    assert "#GAIADIST=245.7" in output_text
    assert "#GAIAPMRA=14.25" in output_text
    assert "#GAIAPMDEC=-9.5" in output_text
    magnitude_fields = aavso_json_header(output_text, "MAGNITUDE_FIELDS-XC")
    magnitude_row = aavso_json_header(output_text, "MAGNITUDE-XC")
    assert magnitude_fields["apparent_calibrated"] is True
    assert magnitude_fields["differential_magnitude"].startswith("target minus")
    assert magnitude_row["differential_magnitude"] == round(differential_mag, 4)
    assert magnitude_row["differential_magnitude_error"] == round(differential_error, 4)
    assert magnitude_row["apparent_magnitude"] == round(12.0 + differential_mag, 4)
    assert magnitude_row["apparent_magnitude_error"] == round(
        np.hypot(0.02, differential_error),
        4,
    )


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

    output_file = tmp_path / "AAVSO_Files" / "AAVSO_HAT-P-32b_2020-01-01.txt"
    output_text = output_file.read_text(encoding="utf-8")

    assert "#OBSNAME=" not in output_text
    assert "#GAIADIST=" not in output_text
    assert "#GAIAPMRA=" not in output_text
    assert "#GAIAPMDEC=" not in output_text
    magnitude_fields = aavso_json_header(output_text, "MAGNITUDE_FIELDS-XC")
    magnitude_row = aavso_json_header(output_text, "MAGNITUDE-XC")
    assert magnitude_fields["apparent_calibrated"] is False
    assert magnitude_row["differential_magnitude"] == pytest.approx(0.0)
    assert magnitude_row["apparent_magnitude"] is None
    assert magnitude_row["apparent_magnitude_error"] is None


def test_aid_comparison_coordinate_headers_index_unique_comparisons_on_separate_lines():
    headers = aid_comparison_coordinate_headers(
        [
            {"cname": "Comp A", "comp_ra": 10.1, "comp_dec": -20.2},
            {"cname": "Comp A", "comp_ra": 10.1, "comp_dec": -20.2},
            {"cname": "Comp B", "comp_ra": 11.3, "comp_dec": -21.4},
        ],
        indexed=True,
    )

    assert headers.splitlines() == [
        "#COMPARISON_1_NAME=Comp A",
        "#COMPARISON_1_RA=10.1000000",
        "#COMPARISON_1_DEC=-20.2000000",
        "#COMPARISON_2_NAME=Comp B",
        "#COMPARISON_2_RA=11.3000000",
        "#COMPARISON_2_DEC=-21.4000000",
    ]


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
        "mag": 12.34567,
        "mag_err": 0.012345,
        "differential_mag": 0.24567,
        "differential_mag_err": 0.006789,
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
        "mag_band": "ClearV",
        "catalog_mag_band": "V",
        "source_id": 12345,
        "separation_arcsec": 0.2,
    }]
    working_artifacts_dir = tmp_path / "working_artifacts"
    working_artifacts_dir.mkdir()
    finder_source = (
        working_artifacts_dir / "FOV_HAT-P-32b_LinearStretch_2020-01-01.png"
    )
    finder_source.write_bytes(b"finder chart")

    AIDOutputFiles(fit, p_dict, i_dict, auid=None, chart_id=None, vsp_params=vsp_params).aavso()

    output_text = (
        tmp_path / "AAVSO_Files" / "AID_AAVSO_HAT-P-32_2020-01-01.txt"
    ).read_text(encoding="utf-8")
    assert (
        tmp_path / "AAVSO_Files" / finder_source.name
    ).read_bytes() == b"finder chart"
    metadata = aavso_json_header(output_text, "COMPARISON-CATALOG-XC")

    assert metadata["source"] == "NextAstro photometry catalog"
    assert metadata["is_aavso_vsp"] is False
    assert metadata["comparison_ra_deg"] == pytest.approx(10.1)
    assert metadata["comparison_dec_deg"] == pytest.approx(-20.2)
    assert metadata["apparent_magnitude"] == pytest.approx(12.1)
    assert metadata["apparent_magnitude_error"] == pytest.approx(0.03)
    assert metadata["magnitude_band"] == "V"
    assert metadata["reported_measurement_band"] == "ClearV"
    assert "#COMPARISON_RA=10.1000000\n#COMPARISON_DEC=-20.2000000\n" in output_text
    assert "#DATE=BJD_TDB" in output_text
    assert "HAT-P-32,2450000.12345,12.3457,0.0123,V,NO,STD" in output_text
    assert (
        "#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,"
        "GROUP,CHART,NOTES,DIFFMAG,DIFFERR\n"
    ) in output_text
    aid_header = next(line for line in output_text.splitlines() if line.startswith("#NAME,"))
    aid_data_row = next(line for line in output_text.splitlines() if not line.startswith("#"))
    assert len(aid_header.split(",")) == len(aid_data_row.split(",")) == 17
    assert aid_data_row.split(",")[-3:] == ["na", "0.2457", "0.0068"]
    assert "|DIFFMAG=" not in output_text
    assert "|DIFFERR=" not in output_text
    magnitude_fields = aavso_json_header(output_text, "MAGNITUDE_FIELDS-XC")
    assert magnitude_fields["apparent_magnitude"] == "MAG"
    assert magnitude_fields["differential_magnitude"] == "DIFFMAG"
    assert magnitude_fields["differential_magnitude_error"] == "DIFFERR"


def test_aid_output_records_calibrated_ensemble_members(tmp_path):
    fit = DummyFit()
    p_dict = {"pName": "Target b", "sName": "Target"}
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "aavso_num": "RTZ",
        "camera": "CCD",
        "filter": "V",
        "lat": "+32.4",
        "long": "-110.7",
        "elev": 2600,
    }
    vsp_params = [{
        "time": 2450000.12345,
        "mag": 12.34,
        "mag_err": 0.02,
        "differential_mag": 1.234567,
        "differential_mag_err": 0.00789,
        "airmass": 1.234,
        "cname": "ENSEMBLE (2 stars)",
        "cmag": None,
        "cmag_err": None,
        "catalog_source": "Calibrated comparison-star ensemble",
        "is_aavso_vsp": False,
        "mag_band": "V",
        "ensemble_reference": True,
        "ensemble_member_count": 2,
        "ensemble_member_labels": ["C1", "C2"],
        "ensemble_member_positions": [[10, 20], [30, 40]],
        "ensemble_member_catalog_magnitudes": [12.0, 12.5],
        "ensemble_member_catalog_errors": [0.01, 0.011],
        "ensemble_member_catalog_sources": ["Catalog A", "Catalog B"],
        "ensemble_member_ra_degs": [10.1, 10.2],
        "ensemble_member_dec_degs": [-20.1, -20.2],
        "ensemble_member_catalog_colors": [0.5, 0.6],
        "ensemble_member_catalog_color_labels": ["B-V", "B-V"],
        "ensemble_member_color_deltas": [0.02, 0.08],
        "ensemble_member_magnitude_deltas": [0.1, 0.4],
        "ensemble_member_similarity_scores": [0.102, 0.408],
    }]

    AIDOutputFiles(fit, p_dict, i_dict, auid=None, chart_id=None, vsp_params=vsp_params).aavso()

    output_text = (
        tmp_path / "AAVSO_Files" / "AID_AAVSO_Target_2020-01-01.txt"
    ).read_text(encoding="utf-8")
    metadata = aavso_json_header(output_text, "COMPARISON-CATALOG-XC")
    ensemble_metadata = aavso_json_header(output_text, "ENSEMBLE-COMPARISONS-XC")
    assert metadata["ensemble_reference"] is True
    assert metadata["ensemble_member_count"] == 2
    assert metadata["ensemble_member_labels"] == ["C1", "C2"]
    assert metadata["ensemble_member_ra_degs"] == [10.1, 10.2]
    assert metadata["ensemble_member_dec_degs"] == [-20.1, -20.2]
    assert metadata["ensemble_member_catalog_colors"] == [0.5, 0.6]
    assert metadata["ensemble_member_color_deltas"] == [0.02, 0.08]
    assert metadata["ensemble_member_magnitude_deltas"] == [0.1, 0.4]
    assert ensemble_metadata["member_count"] == 2
    assert ensemble_metadata["members"][0]["label"] == "C1"
    assert ensemble_metadata["members"][0]["ra_deg"] == pytest.approx(10.1)
    assert ensemble_metadata["members"][0]["dec_deg"] == pytest.approx(-20.1)
    assert ensemble_metadata["members"][1]["label"] == "C2"
    assert ensemble_metadata["members"][1]["ra_deg"] == pytest.approx(10.2)
    assert ensemble_metadata["members"][1]["dec_deg"] == pytest.approx(-20.2)
    assert "Target,2450000.12345,12.3400,0.0200,V,NO,STD,ENSEMBLE (2 stars),na" in output_text
    assert output_text.rstrip().endswith(",na,1.2346,0.0079")


def test_aid_output_samples_large_derived_anchor_label_lists(tmp_path):
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
    anchor_labels = [f"NextAstro-{index}" for index in range(20)]
    vsp_params = [{
        "time": 2450000.12345,
        "mag": 12.34,
        "mag_err": 0.05,
        "airmass": 1.234,
        "cname": "RA=10.1000000 Dec=-20.2000000",
        "cmag": 12.1,
        "cmag_err": 0.03,
        "pos": [493, 202],
        "catalog_source": "Derived from full-field catalog-calibrated stars",
        "is_aavso_vsp": False,
        "derived_catalog_reference": True,
        "derived_reference_anchor_count": len(anchor_labels),
        "derived_reference_anchor_labels": anchor_labels,
        "mag_band": "V",
    }]

    AIDOutputFiles(fit, p_dict, i_dict, auid=None, chart_id=None, vsp_params=vsp_params).aavso()

    output_text = (
        tmp_path / "AAVSO_Files" / "AID_AAVSO_HAT-P-32_2020-01-01.txt"
    ).read_text(encoding="utf-8")
    metadata = aavso_json_header(output_text, "COMPARISON-CATALOG-XC")

    assert metadata["derived_reference_anchor_count"] == 20
    assert "derived_reference_anchor_labels" not in metadata
    assert metadata["derived_reference_anchor_label_sample"] == anchor_labels[:10]


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

    output_text = (
        tmp_path / "AAVSO_Files" / "AID_AAVSO_HAT-P-32_2020-01-01.txt"
    ).read_text(encoding="utf-8")
    metadata = aavso_json_header(output_text, "COMPARISON-CATALOG-XC")

    assert metadata["apparent_magnitude_error"] == pytest.approx(0.001)
    assert "HAT-P-32,2450000.12345,12.3400,0.0010,V,NO,STD" in output_text


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

    output_text = (
        tmp_path / "AAVSO_Files" / "AID_AAVSO_HAT-P-32_2020-01-01.txt"
    ).read_text(encoding="utf-8")

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
    assert "intercomparison_score" in text
    assert "intercomparison_frame_rejected_count" in text
    assert "ensemble_score" not in text
    assert "suitability_outlier_rejected" in text
    assert "overexposure_rejected_count" in text
    assert "Comp 1,101,202,true" in text


def test_final_planetary_params_reports_skipped_airmass_correction(tmp_path):
    fit = DummyFit()
    fit.airmass_fit_skipped = True
    fit.airmass_correction_note = "Skipped (airmass span 0.0400 <= 0.05); no airmass correction applied."
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    output_text = output_file.read_text(encoding="utf-8")

    assert "Airmass correction" in output_text
    assert "no airmass correction applied" in output_text
    assert "Airmass coefficient 1 (a1)" not in output_text


def test_final_planetary_params_reports_nextastro_variability_reference(tmp_path):
    fit = DummyFit()
    (tmp_path / "working_artifacts").mkdir()

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

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    final_params = json.loads(output_file.read_text(encoding="utf-8"))["FINAL PLANETARY PARAMETERS"]

    reference = final_params["Variable Reference Star"]
    assert "NextAstro photometry catalog" in reference
    assert "RA=10.1000000" in reference
    assert "Dec=-20.2000000" in reference
    assert "V=12.3450 +/- 0.0670" in reference


def test_transit_outputs_use_rprs_fallback_uncertainty_when_model_error_missing(tmp_path):
    fit = DummyFit()
    fit.errors.pop("rprs")
    fit.rprs_prior_fallback_applied = True
    fit.rprs_prior_fallback_data_uncertainty = 0.005
    fit.rprs_prior_fallback_note = "Rp/R* fixed to prior."
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {
        "pName": "HAT-P-32 b",
        "pPer": 2.15,
        "pPerUnc": 0.001,
        "rprs": 0.1,
        "rprsUnc": 0.001,
        "aRs": 12.0,
        "aRsUnc": 0.4,
        "inc": 88.5,
        "incUnc": 0.2,
        "ecc": 0.0,
    }
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "filter": "V",
        "filter_desc": "Johnson V",
        "wl_min": None,
        "wl_max": None,
    }

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )
    final_params = json.loads(
        (tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json").read_text(encoding="utf-8")
    )["FINAL PLANETARY PARAMETERS"]

    assert "0.005" in final_params["Ratio of Planet to Stellar Radius (Rp/R*)"]

    _, _, results = aavso_dicts(
        p_dict,
        fit,
        i_dict,
        [0.1],
        (0.1, 0.01),
        (0.2, 0.02),
        (0.3, 0.03),
        (0.4, 0.04),
    )

    assert results["Rp/R*"]["uncertainty"] == "0.0050"


def test_final_planetary_params_reports_transit_comparison_catalog_reference(tmp_path):
    fit = DummyFit()
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}
    vsp_params = [
        {
            "cname": "000-BJX-718",
            "cmag": 9.751,
            "cmag_err": 0.018,
            "pos": [616, 113],
            "catalog_source": "AAVSO VSP",
            "is_aavso_vsp": True,
            "mag_band": "V",
        },
        {
            "cname": "000-BJX-718",
            "cmag": 9.751,
            "cmag_err": 0.018,
            "pos": [616, 113],
            "catalog_source": "AAVSO VSP",
            "is_aavso_vsp": True,
            "mag_band": "V",
        },
    ]

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=True,
        vsp_params=vsp_params,
        comp_star=1,
        comp_coords=[616, 113],
        min_aper=2.7,
        min_annul=10.15,
    )

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    final_params = json.loads(output_file.read_text(encoding="utf-8"))["FINAL PLANETARY PARAMETERS"]

    assert final_params["Transit Fit Comparison Star"] == "#1 - [616, 113]"
    assert "Best Comparison Star" not in final_params
    assert final_params["Variable Reference Star"] == "AAVSO Label: 000-BJX-718, Position: [616, 113]"
    assert "Remeasured 2 out-of-transit target/reference point(s)" in final_params["Variable Reference Measurement"]
    assert "AID rows list the BJD_TDB timestamps used" in final_params["Variable Reference Measurement"]
    assert "transit-fit catalog reference" in final_params["Variable Reference Measurement"]


def test_final_planetary_params_suppresses_variable_reference_without_transit_comparison(tmp_path):
    fit = DummyFit()
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}
    vsp_params = [{
        "cname": "000-BJX-718",
        "cmag": 9.751,
        "cmag_err": 0.018,
        "pos": [616, 113],
        "catalog_source": "AAVSO VSP",
        "is_aavso_vsp": True,
        "mag_band": "V",
    }]

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=True,
        vsp_params=vsp_params,
        comp_star=None,
        comp_coords=None,
        min_aper=-2.7,
        min_annul=10.15,
    )

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    final_params = json.loads(output_file.read_text(encoding="utf-8"))["FINAL PLANETARY PARAMETERS"]

    assert final_params["Transit Fit Comparison Star"] == "None"
    assert "Variable Reference Star" not in final_params
    assert "Variable Reference Measurement" not in final_params


def test_final_planetary_params_reports_ars_and_impact_parameter_under_inclination(tmp_path):
    fit = DummyFit()
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    output_data = json.loads(output_file.read_text(encoding="utf-8"))
    final_params = output_data["FINAL PLANETARY PARAMETERS"]
    keys = list(final_params)
    inclination_index = keys.index("Orbital Inclination (inc)")

    assert keys[inclination_index + 1] == "Ratio of Distance to Stellar Radius (a/Rs)"
    assert keys[inclination_index + 2] == "Impact Parameter (b)"
    assert final_params["Ratio of Distance to Stellar Radius (a/Rs)"] == "12.00 +/- 0.40"

    expected_b, expected_b_error = fit_impact_parameter_value_error(fit)
    assert expected_b == pytest.approx(12.0 * np.cos(np.deg2rad(88.5)))
    assert final_params["Impact Parameter (b)"] == "0.314 +/- 0.043"


def test_final_planetary_params_matches_values_to_two_sigfig_uncertainties(tmp_path):
    fit = DummyFit()
    fit.errors["a1"] = 0.00023
    fit.errors["a2"] = 0.0031
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "filter": "V",
        "filter_desc": "Johnson V",
        "wl_min": None,
        "wl_max": None,
    }

    OutputFiles(fit, p_dict, i_dict, [0.063, 0.083]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    final_params = json.loads(output_file.read_text(encoding="utf-8"))["FINAL PLANETARY PARAMETERS"]

    assert final_params["Flux normalization (a1)"] == "1.00000 +/- 0.00023"
    assert final_params["Airmass coefficient 2 (a2)"] == "0.0000 +/- 0.0031"
    assert final_params["Transit Duration (day)"] == "0.073 +/- 0.010"


def test_detrended_fixed_baseline_does_not_report_inherited_errors_as_fitted(tmp_path):
    fit = DummyFit()
    fit.errors["a1"] = 0.00023
    fit.errors["a2"] = 0.0031
    fit.oot_baseline_detrending_applied = True
    fit.pre_detrending_baseline_source = "test out-of-transit baseline fit"
    fit.pre_detrending_baseline_scale_parameter = "a1"
    fit.pre_detrending_baseline_scale_value = 1.004321
    fit.pre_detrending_baseline_scale_error = 0.00023
    fit.pre_detrending_baseline_a2_value = -0.01234
    fit.pre_detrending_baseline_a2_error = 0.0031
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {
        "save": str(tmp_path),
        "date": "2020-01-01",
        "filter": "V",
        "filter_desc": "Johnson V",
        "wl_min": None,
        "wl_max": None,
    }

    OutputFiles(fit, p_dict, i_dict, [0.063, 0.083]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    final_params = json.loads(output_file.read_text(encoding="utf-8"))["FINAL PLANETARY PARAMETERS"]
    _, _, results = aavso_dicts(
        {
            **p_dict,
            "pPer": 2.15,
            "pPerUnc": 0.001,
            "rprs": 0.1,
            "rprsUnc": 0.001,
            "aRs": 12.0,
            "aRsUnc": 0.4,
            "inc": 88.5,
            "incUnc": 0.2,
            "ecc": 0.0,
        },
        fit,
        i_dict,
        [0.063, 0.083],
        (0.1, 0.01),
        (0.2, 0.02),
        (0.3, 0.03),
        (0.4, 0.04),
    )

    assert final_params["Flux normalization (a1)"] == (
        "1.0 (fixed after out-of-transit baseline detrending)"
    )
    assert final_params["Airmass coefficient 2 (a2)"] == (
        "0.0 (fixed after out-of-transit baseline detrending)"
    )
    assert final_params["Pre-detrending baseline source"] == "test out-of-transit baseline fit"
    assert final_params["Pre-detrending airmass coefficient 1 (a1)"] == (
        "1.00432 +/- 0.00023"
    )
    assert final_params["Pre-detrending airmass coefficient 2 (a2)"] == (
        "-0.0123 +/- 0.0031"
    )
    assert results["Am1"] == {"value": "1.0", "uncertainty": "0"}
    assert results["Am2"] == {"value": "0.0", "uncertainty": "0"}


def test_final_planetary_params_reports_fit_uncertainties_not_prior_uncertainties(tmp_path):
    fit = DummyFit()
    (tmp_path / "working_artifacts").mkdir()

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

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    final_params = json.loads(output_file.read_text(encoding="utf-8"))["FINAL PLANETARY PARAMETERS"]

    assert final_params["Mid-Transit Time (Tmid)"].endswith("+/- 0.00010 BJD_TDB")
    assert final_params["Ratio of Planet to Stellar Radius (Rp/R*)"] == "0.1234 +/- 0.0010"
    assert "Transit depth (Rp/Rs)^2" not in final_params
    assert AREA_DEPTH_LABEL in final_params
    assert OBSERVABLE_DEPTH_LABEL in final_params
    assert PRIOR_OBSERVABLE_DEPTH_LABEL in final_params
    assert OBSERVABLE_DEPTH_DELTA_LABEL in final_params
    assert final_params["Orbital Inclination (inc)"] == "88.50 +/- 0.20 "
    assert final_params["Ratio of Distance to Stellar Radius (a/Rs)"] == "12.00 +/- 0.40"
    assert final_params["Impact Parameter (b)"] == "0.314 +/- 0.043"


def test_fit_empirical_transit_uncertainty_uses_residual_scatter_and_point_counts():
    fit = DummyFit()
    fit.parameters["rprs"] = 0.1
    fit.errors["rprs"] = 0.002
    fit.transit = np.array([1.0, 1.0, 0.99, 0.99, 1.0, 1.0])
    fit.model = np.array(fit.transit)
    fit.data = fit.model + np.array([0.0, 0.01, -0.01, 0.01, -0.01, 0.0])
    fit.residuals = fit.data - fit.model
    fit.dataerr = np.full_like(fit.model, 0.01)
    fit.airmass_model = np.ones_like(fit.model)

    empirical = fit_empirical_transit_uncertainty(fit)

    assert empirical["available"] is True
    assert empirical["in_transit_point_count"] == 2
    assert empirical["out_of_transit_point_count"] == 4
    assert empirical["data_rprs_uncertainty"] == pytest.approx(
        empirical["depth_uncertainty_fraction"] / 0.2
    )
    assert empirical["depth_flux_scatter_fraction"] == pytest.approx(
        empirical["residual_scatter"]
    )
    assert empirical["data_rprs_standard_error"] == pytest.approx(
        empirical["depth_standard_error_fraction"] / 0.2
    )
    assert empirical["data_rprs_flux_scatter_uncertainty"] == pytest.approx(
        empirical["depth_flux_scatter_fraction"] / 0.2
    )
    assert empirical["red_noise_beta_factor"] >= 1.0
    assert empirical["data_rprs_uncertainty"] >= empirical["data_rprs_standard_error"]
    assert empirical["data_rprs_flux_scatter_uncertainty"] > empirical["data_rprs_standard_error"]
    assert empirical["combined_rprs_uncertainty"] > empirical["model_rprs_uncertainty"]
    assert empirical["baseline_red_noise_uncertainty_fraction"] >= (
        empirical["baseline_standard_error_fraction"]
    )
    assert empirical["depth_uncertainty_fraction"] >= empirical["baseline_red_noise_uncertainty_fraction"]


def test_fit_empirical_transit_uncertainty_uses_data_only_for_prior_fallback():
    fit = DummyFit()
    fit.parameters["rprs"] = 0.1
    fit.errors["rprs"] = 0.5
    fit.rprs_prior_fallback_applied = True
    fit.rprs_prior_fallback_note = "Applied Rp/R* prior fallback."
    fit.transit = np.array([1.0, 1.0, 0.99, 0.99, 1.0, 1.0])
    fit.model = np.array(fit.transit)
    fit.data = fit.model + np.array([0.0, 0.01, -0.01, 0.01, -0.01, 0.0])
    fit.residuals = fit.data - fit.model
    fit.dataerr = np.full_like(fit.model, 0.01)
    fit.airmass_model = np.ones_like(fit.model)

    empirical = fit_empirical_transit_uncertainty(fit)

    assert empirical["rprs_uncertainty_basis"] == "prior_assumed_data_only"
    assert np.isnan(empirical["model_rprs_uncertainty"])
    assert empirical["combined_rprs_uncertainty"] == pytest.approx(
        empirical["data_rprs_uncertainty"]
    )
    assert empirical["conservative_rprs_uncertainty"] == pytest.approx(
        empirical["data_rprs_uncertainty"]
    )


def test_final_planetary_params_reports_model_and_red_noise_uncertainties(tmp_path):
    fit = DummyFit()
    fit.parameters["rprs"] = 0.1
    fit.errors["rprs"] = 0.002
    fit.transit = np.array([1.0, 1.0, 0.99, 0.99, 1.0, 1.0])
    fit.model = np.array(fit.transit)
    fit.data = fit.model + np.array([0.0, 0.01, -0.01, 0.01, -0.01, 0.0])
    fit.residuals = fit.data - fit.model
    fit.dataerr = np.full_like(fit.model, 0.01)
    fit.airmass_model = np.ones_like(fit.model)
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    final_params = json.loads(output_file.read_text(encoding="utf-8"))["FINAL PLANETARY PARAMETERS"]

    assert final_params["Ratio of Planet to Stellar Radius (Rp/R*)"] == (
        final_params["Ratio of Planet to Stellar Radius (Rp/R*) model+red-noise uncertainty"]
    )
    assert final_params["Ratio of Planet to Stellar Radius (Rp/R*) model-fit uncertainty"] == (
        "0.1000 +/- 0.0020"
    )
    assert "Ratio of Planet to Stellar Radius (Rp/R*) data-fit red-noise uncertainty" in final_params
    assert "Ratio of Planet to Stellar Radius (Rp/R*) model+red-noise uncertainty" in final_params
    assert "Ratio of Planet to Stellar Radius (Rp/R*) data-fit standard-error estimate" in final_params
    assert "Ratio of Planet to Stellar Radius (Rp/R*) flux-scatter equivalent" in final_params
    assert "Transit depth red-noise uncertainty" in final_params
    assert "Transit depth data-fit standard-error estimate" in final_params
    assert "Transit depth flux-scatter equivalent" in final_params
    assert final_params[AREA_DEPTH_LABEL] == (
        final_params[f"{AREA_DEPTH_LABEL} model+red-noise uncertainty"]
    )
    assert final_params["Mid-Transit Time (Tmid)"] == (
        final_params["Mid-Transit Time (Tmid) model+red-noise uncertainty"]
    )
    assert "Mid-Transit Time (Tmid) model-fit uncertainty" in final_params
    assert final_params["Orbital Inclination (inc)"] == (
        final_params["Orbital Inclination (inc) model+red-noise uncertainty"]
    )
    assert "Orbital Inclination (inc) model-fit uncertainty" in final_params
    assert final_params["Ratio of Distance to Stellar Radius (a/Rs)"] == (
        final_params["Ratio of Distance to Stellar Radius (a/Rs) model+red-noise uncertainty"]
    )
    assert "Ratio of Distance to Stellar Radius (a/Rs) model-fit uncertainty" in final_params
    assert final_params["Impact Parameter (b)"] == (
        final_params["Impact Parameter (b) model+red-noise uncertainty"]
    )
    assert "Impact Parameter (b) model-fit uncertainty" in final_params
    assert f"{AREA_DEPTH_LABEL} model-fit uncertainty" in final_params
    assert f"{AREA_DEPTH_LABEL} data-fit red-noise uncertainty" in final_params
    assert "Flux baseline red-noise uncertainty" in final_params
    assert "Flux baseline standard-error estimate" in final_params
    assert "Red-noise beta factor" in final_params
    assert final_params["Data-fit uncertainty point counts"] == "2 in transit, 4 out of transit"
    assert "primary Rp/R*" in final_params["Uncertainty interpretation note"]
    assert "baseline component" in final_params["Uncertainty interpretation note"]
    assert "time-binning" in final_params["Uncertainty interpretation note"]


def test_final_planetary_params_reports_prior_fallback_data_only_uncertainty(tmp_path):
    fit = DummyFit()
    fit.parameters["rprs"] = 0.1
    fit.errors["rprs"] = 0.5
    fit.rprs_prior_fallback_applied = True
    fit.rprs_prior_fallback_prior_value = 0.1
    fit.rprs_prior_fallback_original_fit_value = 0.11
    fit.rprs_prior_fallback_data_uncertainty = 0.02
    fit.rprs_prior_fallback_note = "Applied Rp/R* prior fallback."
    fit.transit = np.array([1.0, 1.0, 0.99, 0.99, 1.0, 1.0])
    fit.model = np.array(fit.transit)
    fit.data = fit.model + np.array([0.0, 0.01, -0.01, 0.01, -0.01, 0.0])
    fit.residuals = fit.data - fit.model
    fit.dataerr = np.full_like(fit.model, 0.01)
    fit.airmass_model = np.ones_like(fit.model)
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    final_params = json.loads(output_file.read_text(encoding="utf-8"))["FINAL PLANETARY PARAMETERS"]

    assert final_params["Rp/R* uncertainty basis"] == "prior_assumed_data_only"
    assert not any("Rp/R*) model-fit uncertainty" in key for key in final_params)
    assert not any("Rp/R*) model+standard-error" in key for key in final_params)
    assert "input prior Rp/R* value with a data-only" in final_params["Uncertainty interpretation note"]
    assert "Rp/R* prior fallback note" in final_params
    assert any("prior-assumed data-only uncertainty" in key for key in final_params)


def test_final_planetary_params_can_publish_accepted_copy_to_root(tmp_path):
    fit = DummyFit()
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
        publish_to_root=True,
    )

    temp_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    root_file = tmp_path / "FinalParams_HAT-P-32b_2020-01-01.json"

    assert temp_file.exists()
    assert root_file.exists()
    assert root_file.read_text(encoding="utf-8") == temp_file.read_text(encoding="utf-8")


def test_final_planetary_params_reports_adaptive_aperture_summary(tmp_path):
    fit = DummyFit()
    (tmp_path / "working_artifacts").mkdir()

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

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
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
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
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
        "tmid_gaussianity_score": 0.94,
        "tmid_gaussianity_score_uncertainty": 0.03,
        "tmid_gaussianity_effective_sample_count": 1840.0,
        "tmid_gaussianity_detail": "strongly Gaussian-like",
        "ktmf_contributions": [
            {
                "label": "EEBLS Depth SNR",
                "available": True,
                "points": 0.74,
                "max_points": 0.80,
                "score": 0.93,
                "detail": "5.80",
            },
            {
                "label": "Tmid Posterior Gaussianity",
                "available": True,
                "points": 0.94,
                "max_points": 1.00,
                "score": 0.94,
                "score_uncertainty": 0.03,
                "detail": "strongly Gaussian-like",
            },
        ],
    }
    (tmp_path / "working_artifacts").mkdir()

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

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    output_data = json.loads(output_file.read_text(encoding="utf-8"))
    final_params = output_data["FINAL PLANETARY PARAMETERS"]
    final_param_keys = list(final_params)

    assert final_param_keys[:2] == ["Transit detection QC", "KTMF"]
    assert final_params["Transit detection QC"] == "PASS"
    assert final_params["KTMF"] == "4.63 / 5.00"
    assert final_params["KTMF target-fit decision"] == "PASS: KTMF=4.63 / 5.00"
    assert final_params["KTMF comparison selection mode"] == "basis=comparison_field_retry, metric=ktmf"
    assert "selected: highest KTMF" in final_params["KTMF selected comparison decision"]
    assert "Comp 1" in final_params["KTMF comparison candidate 1"]
    assert "not selected: KTMF" in final_params["KTMF comparison candidate 1"]
    assert "Residual Scatter Around Full Model Fit" in final_params["KTMF selected comparison contribution 1"]
    assert "Tmid Posterior Gaussianity" in final_params["KTMF target contribution 2"]

    qc_metadata = build_aavso_qc_metadata(fit)
    assert qc_metadata["tmid_gaussianity_score"] == pytest.approx(0.94)
    assert qc_metadata["tmid_gaussianity_score_uncertainty"] == pytest.approx(0.03)
    assert qc_metadata["tmid_gaussianity_effective_sample_count"] == pytest.approx(1840.0)


def test_final_planetary_params_reports_absolute_fit_quality(tmp_path):
    fit = DummyFit()
    fit.data = np.array([1.0, 1.02, 0.98, 1.01, 0.99, 1.0])
    fit.model = np.ones(6, dtype=float)
    fit.residuals = fit.data - fit.model
    fit.dataerr = np.full(6, 0.01, dtype=float)
    fit.time = np.arange(6, dtype=float)
    fit.airmass_model = np.ones(6, dtype=float)
    fit.bounds = {"tmid": [0, 1], "rprs": [0, 1], "a1": [0, 2]}
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    output_data = json.loads(output_file.read_text(encoding="utf-8"))
    final_params = output_data["FINAL PLANETARY PARAMETERS"]

    assert final_params["Fit quality reduced chi-square"] == "3.333"
    assert final_params["Fit quality chi-square"] == "10.00"
    assert final_params["Fit quality degrees of freedom"] == "3"
    assert final_params["Fit quality RMS residual"] == "1.2910 %"
    assert final_params["Fit quality median absolute normalized residual"] == "1.00 sigma"
    assert final_params["Fit quality RMS residual / median uncertainty"] == "1.29"
    assert final_params["Fit quality point count"] == "6"


def test_final_planetary_params_reports_prior_assumed_geometry_note(tmp_path):
    fit = DummyFit()
    fit.partial_transit_geometry_prior_assumption_note = (
        "Applied prior-assumed transit geometry for a one-sided partial light curve."
    )
    (tmp_path / "working_artifacts").mkdir()

    p_dict = {"pName": "HAT-P-32 b"}
    i_dict = {"save": str(tmp_path), "date": "2020-01-01"}

    OutputFiles(fit, p_dict, i_dict, [0.1]).final_planetary_params(
        phot_opt=False,
        vsp_params=[],
    )

    output_file = tmp_path / "working_artifacts" / "FinalParams_HAT-P-32b_2020-01-01.json"
    final_params = json.loads(output_file.read_text(encoding="utf-8"))["FINAL PLANETARY PARAMETERS"]

    assert (
        final_params["Prior-assumed partial-transit geometry note"]
        == "Applied prior-assumed transit geometry for a one-sided partial light curve."
    )


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

    output_file = tmp_path / "AAVSO_Files" / "AAVSO_HAT-P-32b_2020-01-01.txt"
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
                "label": "EEBLS Depth SNR",
                "available": True,
                "points": 0.74,
                "max_points": 0.80,
                "score": 0.93,
                "detail": "5.80",
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
            "after_target_wcs_filter_frame_count": 3,
            "final_prephotometry_frame_count": 3,
            "ignore_header_wcs": False,
            "bad_wcs_threshold_percent": 3.0,
            "pointing_rejection_sigma": 3.0,
            "dropped_missing_wcs_files": [tmp_path / "missing_wcs.fits"],
            "dropped_target_wcs_files": [tmp_path / "target_off_frame.fits"],
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
            "counts_path": tmp_path / "working_artifacts" / "BadPixelDetectionCounts.fits",
            "mask_path": tmp_path / "working_artifacts" / "BadPixelMask.fits",
        },
    )

    output_text = (
        tmp_path / "AAVSO_Files" / "AAVSO_HAT-P-32b_2020-01-01.txt"
    ).read_text(encoding="utf-8")

    results = aavso_json_header(output_text, "RESULTS-XC")
    assert "a/R*" in results
    assert "Impact Parameter (b)" in results
    assert "Transit depth (Rp/R*)^2" not in results
    assert results[AREA_DEPTH_LABEL]["units"] == "percent"
    assert results[OBSERVABLE_DEPTH_LABEL]["units"] == "percent"
    assert results[PRIOR_OBSERVABLE_DEPTH_LABEL]["units"] == "percent"
    assert results[OBSERVABLE_DEPTH_DELTA_LABEL]["units"] == "percent"
    assert results["Residual scatter around full model fit"]["value"] == "0.32"

    qc = aavso_json_header(output_text, "QC-XC")
    assert qc["status"] == "pass"
    assert qc["ktmf_metric"] == pytest.approx(4.63)
    assert qc["ktmf_contributions"][0]["label"] == "EEBLS Depth SNR"

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
    assert frame_filtering["target_wcs_rejections"]["files"] == ["target_off_frame.fits"]
    assert frame_filtering["pointing_rejections"]["files"] == ["bad_pointing.fits"]
    assert frame_filtering["lightcurve_dropped_point_count"] == 1

    astrometry = aavso_json_header(output_text, "ASTROMETRY-XC")
    assert astrometry["wcs_file"] == "wcs.fits"
    assert astrometry["comparison_star_aavso_header"] == comp_star_header

    bad_pixel = aavso_json_header(output_text, "BAD_PIXEL-XC")
    assert bad_pixel["enabled"] is True
    assert bad_pixel["bad_pixel_count"] == 3
    assert bad_pixel["counts_path"] == "BadPixelDetectionCounts.fits"
