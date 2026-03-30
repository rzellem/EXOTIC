import matplotlib
matplotlib.use("Agg")

import numpy as np
from matplotlib.axes import Axes

from exotic.plots import plot_adaptive_aperture_diagnostics, plot_obs_stats


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
