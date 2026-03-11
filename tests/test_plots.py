import matplotlib
matplotlib.use("Agg")

import numpy as np
from matplotlib.axes import Axes

from exotic.plots import plot_obs_stats


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
