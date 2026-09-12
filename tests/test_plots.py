import matplotlib
matplotlib.use("Agg")

from types import SimpleNamespace

import numpy as np
import matplotlib.pyplot as plt
import pytest
from matplotlib.axes import Axes
import exotic.plots as plots_module

from exotic.plots import (
    _set_canonical_lightcurve_limits,
    _format_parameter_value,
    _short_ktmf_label,
    plot_fov,
    plot_adaptive_aperture_diagnostics,
    plot_comp_star_candidate_lightcurve_fits,
    plot_final_lightcurve,
    plot_individual_comp_star_calibration_series,
    plot_ktmf_qc_metrics,
    plot_obs_stats,
    plot_prior_posterior_comparison,
    plot_differential_magnitude,
    plot_stellar_variability,
    save_figure_formats,
)


class DummyFit:
    def __init__(self):
        self.time = np.array([1.0, 2.0, 3.0])
        self.airmass = np.array([1.1, 1.2, 1.3])


def test_format_parameter_value_uses_uncertainty_precision_without_scientific_notation():
    assert (
        _format_parameter_value(2461197.8645824, 0.0005037355680314821, split_error=True)
        == "2461197.86458\n+/- 0.00050"
    )
    assert _format_parameter_value(89.3511, 2.16, unit="deg") == "89.4 +/- 2.2 deg"


def test_save_figure_formats_keeps_native_square_png_and_adds_square_high_res_png(tmp_path):
    figure, axis = plt.subplots(figsize=(4, 4))
    axis.plot([0.0, 1.0], [0.0, 1.0])
    png_path = tmp_path / "FinalTriangle_Target_2026-03-09.png"

    save_figure_formats(figure, png_path)
    plt.close(figure)

    assert plt.imread(png_path).shape[:2] != (1920, 1920)
    high_res_path = tmp_path / "FinalTriangle_Target_2026-03-09_HighRes.png"
    assert plt.imread(high_res_path).shape[:2] == (1920, 1920)
    assert png_path.with_suffix(".pdf").exists()
    assert png_path.with_suffix(".eps").exists()


def test_plot_obs_stats_applies_relative_flux_mask(tmp_path, monkeypatch):
    fit = DummyFit()
    psf_rows = np.arange(35, dtype=float).reshape(5, 7)
    psf = {"target": psf_rows}
    si = np.array([2, 0, 4, 1, 3])
    gi = np.array([True, False, True, True, True])
    relative_flux_mask = np.array([True, False, True, True])
    captured = []

    original_plot = Axes.plot

    def spy_plot(self, x, y, *args, **kwargs):
        captured.append((np.asarray(x), np.asarray(y)))
        return original_plot(self, x, y, *args, **kwargs)

    monkeypatch.setattr(Axes, "plot", spy_plot)

    plot_obs_stats(
        fit,
        [],
        psf,
        si,
        gi,
        "Target",
        str(tmp_path),
        "2026-03-09",
        relative_flux_mask=relative_flux_mask,
    )

    assert captured
    np.testing.assert_array_equal(captured[0][0], fit.time)
    np.testing.assert_array_equal(captured[0][1], np.array([14.0, 7.0, 21.0]))
    assert (tmp_path / "working_artifacts" / "Observing_Statistics_target_2026-03-09.png").exists()


def test_stellar_variability_final_lightcurve_plots_calibrated_magnitude_by_time(tmp_path, monkeypatch):
    fit = SimpleNamespace(
        stellar_variability_only=True,
        time=np.array([2461229.5, 2461229.6, 2461229.8]),
        jd_times=np.array([2461229.4, 2461229.5, 2461229.7]),
        detrended=np.array([1.0, 1.01, 0.99]),
        detrendederr=np.array([0.001, 0.001, 0.001]),
        time_upsample=np.array([2461229.5, 2461229.8]),
        transit_upsample=np.ones(2),
        stellar_variability_params=[
            {
                "time": 2461229.5,
                "jd_time": 2461229.4,
                "mag": 13.738,
                "mag_err": 0.004,
                "cmag": 13.739,
                "cmag_err": 0.001,
                "mag_band": "r",
                "observed_filter": "SR",
                "comp_ra": 295.3085,
                "comp_dec": 56.1606,
                "cname": "RA=295.3085000 Dec=56.1606000",
            },
            {
                "time": 2461229.6,
                "jd_time": 2461229.5,
                "mag": 13.740,
                "mag_err": 0.004,
                "cmag": 13.739,
                "cmag_err": 0.001,
                "mag_band": "r",
                "observed_filter": "SR",
                "comp_ra": 295.3085,
                "comp_dec": 56.1606,
                "cname": "RA=295.3085000 Dec=56.1606000",
            },
            {
                "time": 2461229.8,
                "jd_time": 2461229.7,
                "mag": 13.735,
                "mag_err": 0.004,
                "cmag": 13.739,
                "cmag_err": 0.001,
                "mag_band": "r",
                "observed_filter": "SR",
                "comp_ra": 295.3085,
                "comp_dec": 56.1606,
                "cname": "RA=295.3085000 Dec=56.1606000",
            },
        ],
        stellar_variability_target_name="WASP-194",
    )
    captured_errorbar_x = []
    captured_errorbar_y = []
    captured_xlabels = []
    captured_ylabels = []
    inverted_axes = []

    original_errorbar = Axes.errorbar
    original_set_xlabel = Axes.set_xlabel
    original_set_ylabel = Axes.set_ylabel
    original_invert_yaxis = Axes.invert_yaxis

    def spy_errorbar(self, x, y, *args, **kwargs):
        captured_errorbar_x.append(np.asarray(x, dtype=float))
        captured_errorbar_y.append(np.asarray(y, dtype=float))
        return original_errorbar(self, x, y, *args, **kwargs)

    def spy_set_xlabel(self, xlabel, *args, **kwargs):
        captured_xlabels.append(xlabel)
        return original_set_xlabel(self, xlabel, *args, **kwargs)

    def spy_set_ylabel(self, ylabel, *args, **kwargs):
        captured_ylabels.append(ylabel)
        return original_set_ylabel(self, ylabel, *args, **kwargs)

    def spy_invert_yaxis(self, *args, **kwargs):
        inverted_axes.append(self)
        return original_invert_yaxis(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "errorbar", spy_errorbar)
    monkeypatch.setattr(Axes, "set_xlabel", spy_set_xlabel)
    monkeypatch.setattr(Axes, "set_ylabel", spy_set_ylabel)
    monkeypatch.setattr(Axes, "invert_yaxis", spy_invert_yaxis)

    plot_final_lightcurve(fit, np.ones(2), "Target", str(tmp_path), "2026-07-08")

    np.testing.assert_allclose(captured_errorbar_x[-1], np.array([2461229.4, 2461229.5, 2461229.7]))
    np.testing.assert_allclose(captured_errorbar_y[-1], np.array([13.738, 13.740, 13.735]))
    assert captured_xlabels[-1] == "Time [JD]"
    assert captured_xlabels[-1] != "Orbital Phase"
    assert captured_ylabels[-1] == "Magnitude (r)"
    assert len(inverted_axes) == 2
    assert "O-C [%]" not in captured_ylabels
    png_path = tmp_path / "FinalLightCurve_Target_2026-07-08.png"
    high_res_png_path = tmp_path / "FinalLightCurve_Target_2026-07-08_HighRes.png"
    assert png_path.exists()
    assert (tmp_path / "FinalLightCurve_Target_2026-07-08.pdf").exists()
    assert (tmp_path / "FinalLightCurve_Target_2026-07-08.eps").exists()
    assert high_res_png_path.exists()
    assert plt.imread(high_res_png_path).shape[:2] == (1080, 1920)


def test_stellar_variability_differential_plot_survives_without_apparent_magnitudes(
        tmp_path, monkeypatch):
    from matplotlib.axes import Axes

    captured_titles = []
    captured_errorbar_x = []
    captured_xlabels = []
    original_set_title = Axes.set_title
    original_errorbar = Axes.errorbar
    original_set_xlabel = Axes.set_xlabel

    def spy_set_title(self, title, *args, **kwargs):
        captured_titles.append(title)
        return original_set_title(self, title, *args, **kwargs)

    def spy_errorbar(self, x, y, *args, **kwargs):
        captured_errorbar_x.append(np.asarray(x, dtype=float))
        return original_errorbar(self, x, y, *args, **kwargs)

    def spy_set_xlabel(self, xlabel, *args, **kwargs):
        captured_xlabels.append(xlabel)
        return original_set_xlabel(self, xlabel, *args, **kwargs)

    monkeypatch.setattr(Axes, "set_title", spy_set_title)
    monkeypatch.setattr(Axes, "errorbar", spy_errorbar)
    monkeypatch.setattr(Axes, "set_xlabel", spy_set_xlabel)
    fit = SimpleNamespace(
        stellar_variability_only=True,
        time=np.array([2461229.5, 2461229.6, 2461229.8]),
        jd_times=np.array([2461229.4, 2461229.5, 2461229.7]),
        data=np.ones(3),
        dataerr=np.full(3, 0.001),
        airmass=np.array([1.1, 1.2, 1.3]),
        airmass_model=np.array([0.8, 1.0, 1.2]),
        transit=np.ones(3),
        stellar_variability_target_flux=np.array([900.0, 1000.0, 1100.0]),
        stellar_variability_comp_flux=np.full(3, 1000.0),
        stellar_variability_target_flux_error=np.ones(3),
        stellar_variability_comp_flux_error=np.ones(3),
        stellar_variability_params=[],
    )

    output_path = plot_differential_magnitude(
        fit,
        'Variable Star',
        tmp_path,
        '2026-08-02',
        observed_filter='V',
    )

    assert output_path.exists()
    assert (tmp_path / 'Stellar_Variability_DifferentialMagnitude.png').exists()
    assert (
        tmp_path / 'working_artifacts' / 'Stellar_Variability_DifferentialMagnitude.png'
    ).exists()
    assert captured_titles[-1] == 'Variable Star'
    np.testing.assert_allclose(captured_errorbar_x[-1], fit.jd_times)
    assert captured_xlabels[-1] == "Time [JD]"


def test_plot_obs_stats_uses_supplied_background_series(tmp_path, monkeypatch):
    fit = DummyFit()
    psf_rows = np.arange(35, dtype=float).reshape(5, 7)
    psf = {"target": psf_rows}
    si = np.array([2, 0, 4, 1, 3])
    gi = np.array([True, False, True, True, True])
    relative_flux_mask = np.array([True, False, True, True])
    background_series = {"target": np.array([100.0, 200.0, 300.0, 400.0, 500.0])}
    captured = []

    original_plot = Axes.plot

    def spy_plot(self, x, y, *args, **kwargs):
        captured.append((np.asarray(x), np.asarray(y)))
        return original_plot(self, x, y, *args, **kwargs)

    monkeypatch.setattr(Axes, "plot", spy_plot)

    plot_obs_stats(
        fit,
        [],
        psf,
        si,
        gi,
        "Target",
        str(tmp_path),
        "2026-03-09",
        relative_flux_mask=relative_flux_mask,
        background_series=background_series,
    )

    assert len(captured) >= 6
    np.testing.assert_array_equal(captured[5][0], fit.time)
    np.testing.assert_array_equal(captured[5][1], np.array([300.0, 200.0, 400.0]))


def test_plot_adaptive_aperture_diagnostics_writes_outputs(tmp_path):
    plot_adaptive_aperture_diagnostics(
        times=np.array([1.0, 2.0, 3.0]),
        aperture_series=np.array([7.5, 8.0, 8.5]),
        annulus_series=np.array([25.0, 26.0, 27.0]),
        fwhm_series=np.array([3.0, 3.2, 3.4]),
        airmass=np.array([1.1, 1.2, 1.3]),
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
        aperture_sigma=2.5,
        annulus_sigma=9.0,
    )

    assert (tmp_path / "working_artifacts" / "AdaptiveApertureDiagnostics_Target_2026-03-09.png").exists()
    assert (tmp_path / "working_artifacts" / "AdaptiveApertureDiagnostics_Target_2026-03-09.pdf").exists()


def test_plot_fov_psf_legend_omits_aperture_annulus_text(tmp_path, monkeypatch):
    labels = []

    original_legend = plt.legend

    def spy_legend(*args, **kwargs):
        legend = original_legend(*args, **kwargs)
        labels.extend(text.get_text() for text in legend.get_texts())
        return legend

    monkeypatch.setattr(plt, "legend", spy_legend)

    plot_fov(
        aper=20.0,
        annulus=60.0,
        sigma=4.0,
        x_targ=50.0,
        y_targ=60.0,
        x_ref=90.0,
        y_ref=100.0,
        image=np.ones((200, 200)),
        image_scale="Image scale in arcsec/pixel: 0.53",
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
        opt_method="PSF",
        min_aper_fov=20.44,
        min_annulus_fov=61.31,
    )

    assert labels
    assert set(labels) == {"PSF Photometry"}


def test_plot_fov_marks_every_ensemble_comparison(tmp_path, monkeypatch):
    plotted_labels = []
    original_text = Axes.text

    def spy_text(self, x, y, text, *args, **kwargs):
        plotted_labels.append(text)
        return original_text(self, x, y, text, *args, **kwargs)

    monkeypatch.setattr(Axes, "text", spy_text)

    plot_fov(
        aper=8.0,
        annulus=20.0,
        sigma=2.0,
        x_targ=50.0,
        y_targ=60.0,
        x_ref=90.0,
        y_ref=100.0,
        image=np.ones((220, 220)),
        image_scale="Image scale in arcsec/pixel: 0.53",
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
        opt_method="Aperture",
        min_aper_fov=8.0,
        min_annulus_fov=20.0,
        comparison_positions=[[90.0, 100.0], [130.0, 140.0], [170.0, 180.0]],
        comparison_labels=["Comp 1", "Comp 3", "Comp 4"],
    )

    assert {"Target", "Comp 1", "Comp 3", "Comp 4"}.issubset(plotted_labels)
    assert (
        tmp_path / "working_artifacts" / "FOV_Target_LinearStretch_2026-03-09.png"
    ).is_file()


def test_plot_individual_comp_star_calibration_series_writes_outputs(tmp_path):
    plot_individual_comp_star_calibration_series(
        times=np.array([1.0, 2.0, 3.0]),
        comp_summaries=[
            {
                "label": "Comp 1",
                "selected": True,
                "aggregate_score": 0.0012,
                "pairwise_ratio_series": {"vs 2": np.array([1.0, 1.01, 0.99])},
                "ensemble_ratio_series": np.array([1.0, 1.005, 0.995]),
            },
            {
                "label": "Comp 2",
                "selected": False,
                "aggregate_score": 0.0025,
                "pairwise_ratio_series": {"vs 1": np.array([0.99, 1.0, 1.01])},
                "ensemble_ratio_series": np.array([0.995, 1.0, 1.005]),
            },
        ],
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
        method_label="PSF photometry",
    )

    assert (tmp_path / "working_artifacts" / "CompStarCalibrationCurve_Comp1_Target_2026-03-09.png").exists()
    assert (tmp_path / "working_artifacts" / "CompStarCalibrationCurve_Comp1_Target_2026-03-09.pdf").exists()
    assert (tmp_path / "working_artifacts" / "CompStarCalibrationCurve_Comp2_Target_2026-03-09.png").exists()
    assert (tmp_path / "working_artifacts" / "CompStarCalibrationCurve_Comp2_Target_2026-03-09.pdf").exists()


def test_plot_individual_comp_star_calibration_series_masks_rejected_frame_lines(tmp_path, monkeypatch):
    captured_lines = {}
    original_plot = Axes.plot

    def spy_plot(self, x, y, *args, **kwargs):
        label = kwargs.get("label")
        if label in {"vs 2", "Intercomparison"}:
            captured_lines[label] = (np.asarray(x), np.asarray(y))
        return original_plot(self, x, y, *args, **kwargs)

    monkeypatch.setattr(Axes, "plot", spy_plot)

    plot_individual_comp_star_calibration_series(
        times=np.array([1.0, 2.0, 3.0, 4.0]),
        comp_summaries=[
            {
                "label": "Comp 1",
                "selected": False,
                "aggregate_score": 0.01,
                "pairwise_ratio_series": {"vs 2": np.array([1.0, 0.05, 1.01, 0.99])},
                "ensemble_ratio_series": np.array([1.0, 0.02, 1.005, 0.995]),
                "ensemble_frame_keep_mask": np.array([True, False, True, True]),
            },
        ],
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
        method_label="Aperture photometry",
    )

    assert np.isnan(captured_lines["vs 2"][1][1])
    assert np.isnan(captured_lines["Intercomparison"][1][1])


def test_plot_stellar_variability_labels_reference_coordinates(tmp_path, monkeypatch):
    titles = []
    errorbar_x = []
    xlabels = []
    ylabels = []
    inverted_axes = []
    original_set_title = Axes.set_title
    original_errorbar = Axes.errorbar
    original_set_xlabel = Axes.set_xlabel
    original_set_ylabel = Axes.set_ylabel
    original_invert_yaxis = Axes.invert_yaxis

    def spy_set_title(self, label, *args, **kwargs):
        titles.append(label)
        return original_set_title(self, label, *args, **kwargs)

    def spy_errorbar(self, x, y, *args, **kwargs):
        errorbar_x.append(float(x))
        return original_errorbar(self, x, y, *args, **kwargs)

    def spy_set_xlabel(self, label, *args, **kwargs):
        xlabels.append(label)
        return original_set_xlabel(self, label, *args, **kwargs)

    def spy_set_ylabel(self, label, *args, **kwargs):
        ylabels.append(label)
        return original_set_ylabel(self, label, *args, **kwargs)

    def spy_invert_yaxis(self, *args, **kwargs):
        inverted_axes.append(self)
        return original_invert_yaxis(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "set_title", spy_set_title)
    monkeypatch.setattr(Axes, "errorbar", spy_errorbar)
    monkeypatch.setattr(Axes, "set_xlabel", spy_set_xlabel)
    monkeypatch.setattr(Axes, "set_ylabel", spy_set_ylabel)
    monkeypatch.setattr(Axes, "invert_yaxis", spy_invert_yaxis)

    plot_stellar_variability(
        [
            {
                "time": 2450000.1,
                "jd_time": 2449999.9,
                "mag": 12.34,
                "mag_err": 0.05,
                "cmag": 12.345,
                "cmag_err": 0.067,
                "comp_ra": 10.1,
                "comp_dec": -20.2,
                "mag_band": "r",
                "observed_filter": "CV",
                "is_aavso_vsp": False,
            }
        ],
        str(tmp_path),
        "Host Star",
        "NextAstro-123",
    )

    assert titles[-1] == (
        "Host Star\n"
        "Label: NextAstro-123\n"
        "Comparison RA=10.100000 | Dec=-20.200000\n"
        "Original filter: CV | Comparison mag: r=12.3450 +/- 0.0670"
    )
    assert "No airmass correction applied to stellar variability" not in titles[-1]
    assert errorbar_x == pytest.approx([2449999.9])
    assert xlabels[-1] == "Time [JD]"
    assert ylabels[-1] == "Magnitude (r)"
    assert len(inverted_axes) == 1
    assert (tmp_path / "working_artifacts" / "Stellar_Variability.png").exists()
    assert (tmp_path / "Stellar_Variability.png").exists()


def test_plot_stellar_variability_labels_aavso_filter_and_assumed_comparison(tmp_path, monkeypatch):
    titles = []
    ylabels = []
    original_set_title = Axes.set_title
    original_set_ylabel = Axes.set_ylabel

    def spy_set_title(self, label, *args, **kwargs):
        titles.append(label)
        return original_set_title(self, label, *args, **kwargs)

    def spy_set_ylabel(self, label, *args, **kwargs):
        ylabels.append(label)
        return original_set_ylabel(self, label, *args, **kwargs)

    monkeypatch.setattr(Axes, "set_title", spy_set_title)
    monkeypatch.setattr(Axes, "set_ylabel", spy_set_ylabel)

    plot_stellar_variability(
        [
            {
                "time": 2450000.1,
                "mag": 12.34,
                "mag_err": 0.05,
                "cmag": 12.345,
                "cmag_err": 0.067,
                "comp_ra": 10.1,
                "comp_dec": -20.2,
                "mag_band": "ClearV",
                "catalog_mag_band": "V",
                "observed_filter": "CV",
                "is_aavso_vsp": True,
            }
        ],
        str(tmp_path),
        "Host Star",
        "000-BJX-718",
    )

    assert (
        "Label: 000-BJX-718\n"
        "Comparison RA=10.100000 | Dec=-20.200000"
    ) in titles[-1]
    assert "Original filter: CV" in titles[-1]
    assert "Comparison mag: V=12.3450 +/- 0.0670" in titles[-1]
    assert ylabels[-1] == "Magnitude (ClearV)"


def test_plot_stellar_variability_omits_invalid_reference_magnitudes(tmp_path, monkeypatch):
    titles = []
    original_set_title = Axes.set_title

    def spy_set_title(self, label, *args, **kwargs):
        titles.append(label)
        return original_set_title(self, label, *args, **kwargs)

    monkeypatch.setattr(Axes, "set_title", spy_set_title)

    plot_stellar_variability(
        [
            {
                "time": 2450000.1,
                "mag": 12.34,
                "mag_err": 0.05,
                "cmag": 99.99,
                "cmag_err": 99.99,
                "comp_ra": 10.1,
                "comp_dec": -20.2,
                "mag_band": "V",
                "observed_filter": "MObs CV",
                "is_aavso_vsp": False,
            }
        ],
        str(tmp_path),
        "Host Star",
        "NextAstro-123",
    )

    assert titles[-1] == (
        "Host Star\n"
        "Label: NextAstro-123\nComparison RA=10.100000 | Dec=-20.200000\n"
        "Original filter: MObs CV"
    )
    assert "99.99" not in titles[-1]
    assert "V=" not in titles[-1]


def test_plot_stellar_variability_skips_over_30_measurements(tmp_path):
    plot_stellar_variability(
        [
            {
                "time": 2450000.1,
                "mag": 99.99,
                "mag_err": 0.05,
                "cmag": 12.0,
                "cmag_err": 0.05,
                "mag_band": "V",
            }
        ],
        str(tmp_path),
        "Host Star",
        "Comp",
    )

    assert not (tmp_path / "working_artifacts" / "Stellar_Variability.png").exists()
    assert not (tmp_path / "Stellar_Variability.png").exists()


def test_plot_comp_star_candidate_lightcurve_fits_writes_outputs(tmp_path):
    class DummyCandidateFit:
        def __init__(self):
            self.kwargs = None

        def plot_bestfit(self, phase=False, show_flux_baseline_label=True):
            self.kwargs = {
                "phase": phase,
                "show_flux_baseline_label": show_flux_baseline_label,
            }
            fig, axes = plt.subplots(2, 1)
            return fig, axes
    selected_fit = DummyCandidateFit()
    other_fit = DummyCandidateFit()

    plot_comp_star_candidate_lightcurve_fits(
        candidate_fit_summaries=[
            {"label": "Comp 1", "selected": True, "fit": selected_fit, "res_std": 0.0012},
            {"label": "Comp 2", "selected": False, "fit": other_fit, "res_std": 0.0025},
            {"label": "Comp 3", "selected": False, "fit": None, "res_std": np.inf},
        ],
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
        method_label="Aperture photometry (aper=5.00px, annulus=12.00px)",
    )

    assert selected_fit.kwargs == {"phase": False, "show_flux_baseline_label": False}
    assert other_fit.kwargs == {"phase": False, "show_flux_baseline_label": False}
    assert (tmp_path / "working_artifacts" / "CompStarLightCurveFit_Comp1_Target_2026-03-09.png").exists()
    assert (tmp_path / "working_artifacts" / "CompStarLightCurveFit_Comp1_Target_2026-03-09.pdf").exists()
    assert (tmp_path / "working_artifacts" / "CompStarLightCurveFit_Comp2_Target_2026-03-09.png").exists()
    assert (tmp_path / "working_artifacts" / "CompStarLightCurveFit_Comp2_Target_2026-03-09.pdf").exists()
    assert not (tmp_path / "working_artifacts" / "CompStarLightCurveFit_Comp3_Target_2026-03-09.png").exists()


def test_plot_final_lightcurve_requests_uncertainty_bands_without_baseline_label(tmp_path):
    class DummyFinalFit:
        def __init__(self):
            self.kwargs = None
            self.phase_upsample = np.linspace(-0.05, 0.05, 5)
            self.transit_upsample = np.ones(5)

        def plot_bestfit(self, show_flux_baseline_label=True, show_model_uncertainty=False,
                         show_baseline_uncertainty=False):
            self.kwargs = {
                "show_flux_baseline_label": show_flux_baseline_label,
                "show_model_uncertainty": show_model_uncertainty,
                "show_baseline_uncertainty": show_baseline_uncertainty,
            }
            fig, axes = plt.subplots(2, 1)
            return fig, axes

    fit = DummyFinalFit()

    plot_final_lightcurve(
        fit,
        high_res=np.ones(5),
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
    )

    assert fit.kwargs == {
        "show_flux_baseline_label": False,
        "show_model_uncertainty": True,
        "show_baseline_uncertainty": True,
    }
    png_path = tmp_path / "FinalLightCurve_Target_2026-03-09.png"
    high_res_png_path = tmp_path / "FinalLightCurve_Target_2026-03-09_HighRes.png"
    assert png_path.exists()
    assert (tmp_path / "FinalLightCurve_Target_2026-03-09.pdf").exists()
    assert (tmp_path / "FinalLightCurve_Target_2026-03-09.eps").exists()
    assert high_res_png_path.exists()
    assert plt.imread(high_res_png_path).shape[:2] == (1080, 1920)


def test_plot_final_lightcurve_splits_full_data_diagnostic_from_canonical_plot(
    tmp_path,
    monkeypatch,
):
    class DummyFinalFit:
        def __init__(self):
            self.calls = []
            self.phase_upsample = np.linspace(-0.05, 0.05, 5)
            self.transit_upsample = np.ones(5)
            self.restricted_baseline_points = {
                "point_count": 2,
                "times": np.array([-0.08, 0.08]),
                "flux": np.array([1.0, 1.0]),
                "unc": np.array([0.01, 0.01]),
            }

        def plot_bestfit(self, **kwargs):
            self.calls.append(dict(kwargs))
            fig, axes = plt.subplots(2, 1)
            return fig, axes

    rejected_point_draws = []
    monkeypatch.setattr(
        plots_module,
        "_plot_final_residual_rejected_points",
        lambda *args: rejected_point_draws.append(args),
    )

    fit = DummyFinalFit()
    plot_final_lightcurve(
        fit,
        high_res=np.ones(5),
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
    )

    assert len(fit.calls) == 2
    assert fit.calls[0]["show_restricted_baseline_points"] is True
    assert fit.calls[0]["show_binned_points"] is True
    assert fit.calls[1]["show_restricted_baseline_points"] is False
    assert fit.calls[1]["show_binned_points"] is False
    assert len(rejected_point_draws) == 2
    assert (
        tmp_path
        / "Diagnostics"
        / "FullDataFullLightCurve_Target_2026-03-09.png"
    ).exists()
    assert (tmp_path / "FinalLightCurve_Target_2026-03-09.png").exists()


def test_canonical_lightcurve_limits_follow_visible_points_not_errorbar_extents():
    figure, (ax_lc, ax_res) = plt.subplots(2, 1)
    fit = SimpleNamespace(
        phase=np.array([-0.05, -0.04, -0.03]),
        detrended=np.array([0.98, 1.00, 1.04]),
        final_residual_rejection={
            "applied": True,
            "rejected_phase": [0.02],
            "rejected_flux": [0.50],
        },
    )

    _set_canonical_lightcurve_limits(fit, ax_lc, ax_res)

    # Median cadence is 0.01 phase, so the visible x range gets a half-cadence
    # margin; the rejected red marker is included in the x extent.
    np.testing.assert_allclose(ax_lc.get_xlim(), [-0.055, 0.025])
    np.testing.assert_allclose(ax_res.get_xlim(), [-0.055, 0.025])
    # Y limits come from black measured points +/- 0.01, not the red cross or
    # any large uncertainty bars attached to the points.
    np.testing.assert_allclose(ax_lc.get_ylim(), [0.97, 1.05])
    plt.close(figure)


def test_plot_final_lightcurve_adds_apparent_magnitude_axis_when_calibrated(tmp_path, monkeypatch):
    class DummyFinalFit:
        def __init__(self):
            self.kwargs = None
            self.phase_upsample = np.linspace(-0.05, 0.05, 5)
            self.transit_upsample = np.ones(5)
            self.stellar_variability_params = [
                {"time": 1.0, "mag": 13.739, "mag_err": 0.001, "mag_band": "r"},
                {"time": 2.0, "mag": 13.741, "mag_err": 0.002, "mag_band": "r"},
            ]

        def plot_bestfit(self, show_flux_baseline_label=True, show_model_uncertainty=False,
                         show_baseline_uncertainty=False):
            self.kwargs = {
                "show_flux_baseline_label": show_flux_baseline_label,
                "show_model_uncertainty": show_model_uncertainty,
                "show_baseline_uncertainty": show_baseline_uncertainty,
            }
            fig, axes = plt.subplots(2, 1)
            return fig, axes

    secondary_calls = []
    secondary_labels = []

    class FakeSecondaryAxis:
        def set_ylabel(self, label):
            secondary_labels.append(label)

    def spy_secondary_yaxis(self, location, functions=None, *args, **kwargs):
        secondary_calls.append((location, functions))
        return FakeSecondaryAxis()

    monkeypatch.setattr(Axes, "secondary_yaxis", spy_secondary_yaxis)

    plot_final_lightcurve(
        DummyFinalFit(),
        high_res=np.ones(5),
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
    )

    assert secondary_labels == ["Apparent Magnitude (r)"]
    assert secondary_calls[0][0] == "right"
    flux_to_mag, mag_to_flux = secondary_calls[0][1]
    assert flux_to_mag(np.array([1.0])) == pytest.approx(np.array([13.740]))
    assert mag_to_flux(np.array([13.740])) == pytest.approx(np.array([1.0]))


def test_plot_final_lightcurve_draws_data_scatter_uncertainty_band(tmp_path, monkeypatch):
    captured = []
    original_fill_between = Axes.fill_between

    def spy_fill_between(self, x, y1, y2=0, *args, **kwargs):
        captured.append({
            "x": np.asarray(x, dtype=float),
            "y1": np.asarray(y1, dtype=float),
            "y2": np.asarray(y2, dtype=float),
            "color": kwargs.get("color"),
            "alpha": kwargs.get("alpha"),
            "label": kwargs.get("label"),
        })
        return original_fill_between(self, x, y1, y2, *args, **kwargs)

    monkeypatch.setattr(Axes, "fill_between", spy_fill_between)

    class DummyFinalFit:
        def __init__(self):
            self.phase_upsample = np.linspace(-0.05, 0.05, 41)
            depth_shape = np.exp(-0.5 * (self.phase_upsample / 0.015) ** 2)
            self.transit_upsample = 1.0 - 0.01 * depth_shape
            self.time_upsample = self.phase_upsample.copy()
            self.phase = self.phase_upsample.copy()
            self.transit = self.transit_upsample.copy()
            self.model = self.transit.copy()
            residual_pattern = 0.02 * np.sin(np.linspace(0, 6 * np.pi, self.model.size))
            self.data = self.model + residual_pattern
            self.residuals = self.data - self.model
            self.dataerr = np.full_like(self.model, 0.02)
            self.parameters = {"rprs": 0.1}
            self.errors = {"rprs": 0.001}

        def transit_model_uncertainty(self, times):
            return self.transit_upsample - 0.001, self.transit_upsample + 0.001

        def plot_bestfit(self, show_flux_baseline_label=True, show_model_uncertainty=False,
                         show_baseline_uncertainty=False):
            fig, axes = plt.subplots(2, 1)
            axes[0].plot(self.phase_upsample, self.transit_upsample, 'r-', label='model')
            axes[0].legend(loc='best')
            return fig, axes

    plot_final_lightcurve(
        DummyFinalFit(),
        high_res=np.ones(41),
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
    )

    purple_bands = [item for item in captured if item["color"] == "#6a1b9a"]
    assert purple_bands
    assert all(item["label"] == "_nolegend_" for item in purple_bands)
    assert all(item["alpha"] <= 0.16 for item in purple_bands)
    assert any(np.nanmax(np.abs(item["y2"] - item["y1"])) > 0.001 for item in purple_bands)


def test_plot_final_lightcurve_marks_final_residual_rejections(tmp_path, monkeypatch):
    captured = []
    original_scatter = Axes.scatter

    def spy_scatter(self, x, y, *args, **kwargs):
        captured.append({
            "x": np.asarray(x, dtype=float),
            "y": np.asarray(y, dtype=float),
            "label": kwargs.get("label"),
            "color": kwargs.get("color"),
        })
        return original_scatter(self, x, y, *args, **kwargs)

    monkeypatch.setattr(Axes, "scatter", spy_scatter)

    class DummyFinalFit:
        def __init__(self):
            self.phase_upsample = np.linspace(-0.05, 0.05, 5)
            self.transit_upsample = np.ones(5)
            self.final_residual_rejection = {
                "applied": True,
                "rejected_phase": [0.01],
                "rejected_flux": [0.92],
                "rejected_residual_percent": [-8.0],
            }

        def plot_bestfit(self, show_flux_baseline_label=True, show_model_uncertainty=False,
                         show_baseline_uncertainty=False):
            fig, axes = plt.subplots(2, 1)
            return fig, axes

    plot_final_lightcurve(
        DummyFinalFit(),
        high_res=np.ones(5),
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
    )

    residual_rejection_points = [
        item for item in captured
        if item["label"] == "_nolegend_" and item["color"] == "red"
    ]
    assert len(residual_rejection_points) == 2
    np.testing.assert_allclose(residual_rejection_points[0]["x"], [0.01])
    np.testing.assert_allclose(residual_rejection_points[1]["y"], [-8.0])


def test_plot_prior_posterior_comparison_omits_prior_fallback_rprs(tmp_path, monkeypatch):
    captured_text = []
    original_text = Axes.text

    def spy_text(self, x, y, s, *args, **kwargs):
        captured_text.append(str(s))
        return original_text(self, x, y, s, *args, **kwargs)

    monkeypatch.setattr(Axes, "text", spy_text)

    class DummyFit:
        def __init__(self):
            self.parameters = {
                "tmid": 1.012,
                "rprs": 0.1,
                "ars": 10.6,
                "inc": 88.2,
            }
            self.errors = {
                "tmid": 0.002,
                "rprs": 0.005,
                "ars": 0.4,
                "inc": 0.3,
            }
            self.rprs_prior_fallback_applied = True
            self.empirical_transit_uncertainty = {
                "available": True,
                "combined_rprs_uncertainty": 0.02,
                "rprs_uncertainty_basis": "prior_assumed_data_only",
            }

    planet_dict = {
        "midT": 1.0,
        "midTUnc": 0.001,
        "rprs": 0.1,
        "rprsUnc": 0.003,
        "aRs": 10.0,
        "aRsUnc": 0.2,
        "inc": 89.0,
        "incUnc": 0.4,
    }

    output = plot_prior_posterior_comparison(
        DummyFit(),
        planet_dict,
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
    )

    assert output == tmp_path / "PriorPosteriorComparison_Target_2026-03-09.png"
    assert output.exists()
    assert (tmp_path / "PriorPosteriorComparison_Target_2026-03-09.pdf").exists()
    assert any("Rp/R* omitted: prior value assumed, not measured" in text for text in captured_text)
    assert any("Prior\n" in text and "Posterior\n" in text for text in captured_text)
    assert not any("BJD_TDB" in text for text in captured_text)
    assert not any("Prior epoch" in text for text in captured_text)
    assert not any("Posterior (Prior)" in text for text in captured_text)


def test_plot_ktmf_qc_metrics_writes_outputs_and_annotations(tmp_path, monkeypatch):
    captured_text = []
    original_text = Axes.text

    def spy_text(self, x, y, s, *args, **kwargs):
        captured_text.append(str(s))
        return original_text(self, x, y, s, *args, **kwargs)

    monkeypatch.setattr(Axes, "text", spy_text)

    class DummyKTMFFit:
        def __init__(self):
            self.transit_qc = {
                "status": "fail",
                "ktmf_metric": 3.07,
                "ktmf_contributions": [
                    {
                        "label": "EEBLS Depth SNR",
                        "available": True,
                        "points": 0.11,
                        "max_points": 0.89,
                        "score": 0.12,
                        "detail": "2.00",
                    },
                    {
                        "label": "Residual Scatter Around Full Model Fit",
                        "available": True,
                        "points": 0.14,
                        "max_points": 0.78,
                        "score": 0.18,
                        "detail": "2.3437%",
                    },
                ],
            }

    output = plot_ktmf_qc_metrics(
        DummyKTMFFit(),
        targ_name="Target",
        save=str(tmp_path),
        date="2026-03-09",
    )

    assert output == tmp_path / "KTMF_QC_Target_2026-03-09.png"
    assert output.exists()
    assert (tmp_path / "KTMF_QC_Target_2026-03-09.pdf").exists()
    assert any("KTMF\n3.07 / 5.00\nMARGINAL" in text for text in captured_text)
    assert any("0.11 / 0.89" in text for text in captured_text)
    assert not any("2.00" in text for text in captured_text)
    assert not any("2.3437" in text for text in captured_text)


def test_ktmf_plot_shortens_tmid_posterior_gaussianity_label():
    assert _short_ktmf_label("Tmid Posterior Gaussianity") == "Tmid Gaussianity"


def test_final_lightcurve_title_carries_low_coverage_verdict():
    class DummyFinalFit:
        def __init__(self):
            self.phase_upsample = np.linspace(-0.05, 0.05, 5)
            self.transit_upsample = np.ones(5)
            self.pre_ultranest_transit_coverage_valid = True
            self.pre_ultranest_transit_coverage_expected_successful = False
            self.pre_ultranest_transit_coverage_status = "very low"
            self.pre_ultranest_transit_coverage = {
                "valid": True,
                "observed_segment": "pre-transit baseline only",
                "transit_fraction_observed": 0.0,
                "success_label": "very low",
            }

        def plot_bestfit(self, **kwargs):
            fig, axes = plt.subplots(2, 1)
            return fig, axes

    fit = DummyFinalFit()
    figure = plots_module._build_final_lightcurve_figure(
        fit,
        np.ones(5),
        "Target",
        show_restricted_baseline_points=False,
        show_binned_points=False,
    )
    title = figure.axes[0].get_title()
    plt.close(figure)

    assert title.splitlines()[0] == "Target"
    assert title.splitlines()[1:] == [
        "COVERAGE 0%: pre-transit baseline only",
        "fit success VERY LOW",
    ]


def test_final_lightcurve_title_is_unchanged_when_coverage_was_expected_to_succeed():
    class DummyFinalFit:
        def __init__(self):
            self.phase_upsample = np.linspace(-0.05, 0.05, 5)
            self.transit_upsample = np.ones(5)
            self.pre_ultranest_transit_coverage_valid = True
            self.pre_ultranest_transit_coverage_expected_successful = True
            self.pre_ultranest_transit_coverage = {"observed_segment": "full transit"}

        def plot_bestfit(self, **kwargs):
            fig, axes = plt.subplots(2, 1)
            return fig, axes

    figure = plots_module._build_final_lightcurve_figure(
        DummyFinalFit(),
        np.ones(5),
        "Target",
        show_restricted_baseline_points=False,
        show_binned_points=False,
    )
    title = figure.axes[0].get_title()
    plt.close(figure)

    assert title == "Target"
    assert plots_module.transit_coverage_caption(SimpleNamespace()) is None
