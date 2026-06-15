import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

from exotic.plots import (
    plot_fov,
    plot_adaptive_aperture_diagnostics,
    plot_comp_star_candidate_lightcurve_fits,
    plot_final_lightcurve,
    plot_individual_comp_star_calibration_series,
    plot_obs_stats,
    plot_stellar_variability,
)


class DummyFit:
    def __init__(self):
        self.time = np.array([1.0, 2.0, 3.0])
        self.airmass = np.array([1.1, 1.2, 1.3])


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
    assert (tmp_path / "temp" / "Observing_Statistics_target_2026-03-09.png").exists()


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

    assert (tmp_path / "temp" / "AdaptiveApertureDiagnostics_Target_2026-03-09.png").exists()
    assert (tmp_path / "temp" / "AdaptiveApertureDiagnostics_Target_2026-03-09.pdf").exists()


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

    assert (tmp_path / "temp" / "CompStarCalibrationCurve_Comp1_Target_2026-03-09.png").exists()
    assert (tmp_path / "temp" / "CompStarCalibrationCurve_Comp1_Target_2026-03-09.pdf").exists()
    assert (tmp_path / "temp" / "CompStarCalibrationCurve_Comp2_Target_2026-03-09.png").exists()
    assert (tmp_path / "temp" / "CompStarCalibrationCurve_Comp2_Target_2026-03-09.pdf").exists()


def test_plot_individual_comp_star_calibration_series_masks_rejected_frame_lines(tmp_path, monkeypatch):
    captured_lines = {}
    original_plot = Axes.plot

    def spy_plot(self, x, y, *args, **kwargs):
        label = kwargs.get("label")
        if label in {"vs 2", "Ensemble"}:
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
    assert np.isnan(captured_lines["Ensemble"][1][1])


def test_plot_stellar_variability_labels_reference_coordinates(tmp_path, monkeypatch):
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
                "mag_band": "r",
                "observed_filter": "CV",
                "is_aavso_vsp": False,
            }
        ],
        str(tmp_path),
        "Host Star",
        "NextAstro-123",
    )

    assert titles[-1] == "Host Star\nRA=10.100000, Dec=-20.200000"
    assert "Comparison:" not in titles[-1]
    assert "Observed filter" not in titles[-1]
    assert "r=12.345 +/- 0.067" not in titles[-1]
    assert ylabels[-1] == "Magnitude (r)"
    assert (tmp_path / "temp" / "Stellar_Variability.png").exists()


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

    assert titles[-1] == "Host Star\nRA=10.100000, Dec=-20.200000"
    assert "Observed filter" not in titles[-1]
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

    assert not (tmp_path / "temp" / "Stellar_Variability.png").exists()


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
    assert (tmp_path / "temp" / "CompStarLightCurveFit_Comp1_Target_2026-03-09.png").exists()
    assert (tmp_path / "temp" / "CompStarLightCurveFit_Comp1_Target_2026-03-09.pdf").exists()
    assert (tmp_path / "temp" / "CompStarLightCurveFit_Comp2_Target_2026-03-09.png").exists()
    assert (tmp_path / "temp" / "CompStarLightCurveFit_Comp2_Target_2026-03-09.pdf").exists()
    assert not (tmp_path / "temp" / "CompStarLightCurveFit_Comp3_Target_2026-03-09.png").exists()


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
    assert (tmp_path / "FinalLightCurve_Target_2026-03-09.png").exists()
    assert (tmp_path / "FinalLightCurve_Target_2026-03-09.pdf").exists()
