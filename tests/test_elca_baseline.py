import importlib
import os
import sys
import types

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np
import pytest


def load_elca_with_stubs(monkeypatch, tmp_path):
    import exotic.api as exotic_api

    root = tmp_path / "stubdeps"
    package_dir = root / "pylightcurve"
    model_dir = package_dir / "models"
    model_dir.mkdir(parents=True)

    (package_dir / "__init__.py").write_text("", encoding="utf-8")
    (model_dir / "__init__.py").write_text("", encoding="utf-8")
    (model_dir / "exoplanet_lc.py").write_text(
        "import numpy as np\n"
        "def transit(ld, rprs, per, ars, ecc, inc, omega, tmid, times, method=None, precision=None):\n"
        "    times = np.asarray(times, dtype=float)\n"
        "    return 1.0 - (rprs ** 2) * np.exp(-0.5 * ((times - tmid) / 0.01) ** 2)\n",
        encoding="utf-8",
    )

    monkeypatch.syspath_prepend(str(root))

    # Track both import caches with monkeypatch so the real ELCA/PyLightCurve
    # modules are restored after each test.  A raw sys.modules.pop() leaves the
    # stubbed ELCA module behind and contaminates tests collected after this file.
    monkeypatch.delattr(exotic_api, "elca", raising=False)
    for name in list(sys.modules):
        if name == "exotic.api.elca" or name.startswith("pylightcurve"):
            monkeypatch.delitem(sys.modules, name, raising=False)

    fake_ultranest = types.ModuleType("ultranest")
    fake_ultranest.ReactiveNestedSampler = type("ReactiveNestedSampler", (), {})
    fake_plotting = types.ModuleType("plotting")
    fake_plotting.corner = lambda *args, **kwargs: None
    fake_ultranest_utils = types.ModuleType("ultranest_utils")
    fake_ultranest_utils.run_reactive_sampler = lambda *args, **kwargs: None

    monkeypatch.setitem(sys.modules, "ultranest", fake_ultranest)
    monkeypatch.setitem(sys.modules, "plotting", fake_plotting)
    monkeypatch.setitem(sys.modules, "exotic.api.plotting", fake_plotting)
    monkeypatch.setitem(sys.modules, "ultranest_utils", fake_ultranest_utils)
    monkeypatch.setitem(sys.modules, "exotic.api.ultranest_utils", fake_ultranest_utils)

    import exotic.api.elca as elca

    return importlib.reload(elca)


def make_prior():
    return {
        "rprs": 0.1,
        "ars": 12.0,
        "per": 3.0,
        "inc": 89.0,
        "u0": 0.0,
        "u1": 0.0,
        "u2": 0.0,
        "u3": 0.0,
        "ecc": 0.0,
        "omega": 90.0,
        "tmid": 0.0,
        "a0": 1.0,
        "a2": 0.0,
    }


def make_expanded_prior_warmstart_source(
    prior,
    time,
    data,
    dataerr,
    airmass,
    sample_points,
):
    sample_points = np.asarray(sample_points, dtype=float)
    sample_count = sample_points.shape[0]
    return types.SimpleNamespace(
        ns_type="ultranest",
        bounds={
            "rprs": [0.08, 0.12],
            "tmid": [-0.005, 0.005],
        },
        sampled_keys=["rprs", "tmid"],
        prior=prior.copy(),
        time=np.asarray(time, dtype=float),
        data=np.asarray(data, dtype=float),
        dataerr=np.asarray(dataerr, dtype=float),
        airmass=np.asarray(airmass, dtype=float),
        exposure_times_days=None,
        baseline_fit_mask=None,
        duration_prior=None,
        fixed_flux_baseline=False,
        use_impactparameter_rather_than_inclination_to_fit=False,
        results={
            "weighted_samples": {
                "points": sample_points,
                "weights": np.full(sample_count, 1.0 / sample_count),
                "logl": np.linspace(-5.0, -1.0, sample_count),
            },
        },
    )


def make_dummy_nested_result(sample_points, auxiliary=False):
    sample_points = np.asarray(sample_points, dtype=float)
    if auxiliary:
        result_points = np.column_stack([
            sample_points,
            np.zeros(sample_points.shape[0], dtype=float),
        ])
    else:
        result_points = sample_points
    parameter_count = result_points.shape[1]
    maximum_likelihood = np.zeros(parameter_count, dtype=float)
    maximum_likelihood[0] = 0.1
    return {
        "maximum_likelihood": {"point": maximum_likelihood},
        "posterior": {
            "stdev": np.full(parameter_count, 0.001),
            "errlo": np.full(parameter_count, -0.001),
            "errup": np.full(parameter_count, 0.001),
        },
        "weighted_samples": {
            "points": result_points,
            "weights": np.full(result_points.shape[0], 1.0 / result_points.shape[0]),
            "logl": np.linspace(-5.0, -1.0, result_points.shape[0]),
        },
        "samples": result_points.copy(),
    }


def test_lc_fitter_recovers_explicit_a0_baseline(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 301)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = 0.98 * elca.transit(time, prior)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]},
        mode="lm",
        verbose=False,
    )

    oot_mask = np.abs(fit.time - fit.parameters["tmid"]) > 0.02
    assert fit.parameters["a0"] == pytest.approx(0.98, abs=1e-4)
    assert fit.parameters["a1"] == pytest.approx(0.98, abs=1e-4)
    assert np.median(fit.detrended[oot_mask]) == pytest.approx(1.0, abs=5e-4)


def test_lc_fitter_lm_minimizes_normalized_residuals_and_reports_covariance(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    true_prior = prior.copy()
    true_prior["tmid"] = 0.0025
    time = np.linspace(-0.04, 0.04, 161)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 8e-4)
    rng = np.random.default_rng(2048)
    data = elca.transit(time, true_prior) + rng.normal(0.0, dataerr)

    captured = {}
    scipy_least_squares = elca.least_squares

    def capture_initial_residuals(fun, *args, **kwargs):
        x0 = np.asarray(kwargs.get("x0", args[0] if args else []), dtype=float)
        captured["residuals"] = np.asarray(fun(x0), dtype=float)
        return scipy_least_squares(fun, *args, **kwargs)

    monkeypatch.setattr(elca, "least_squares", capture_initial_residuals)
    initial_model = elca.transit(time, prior)
    initial_model *= elca.solve_flux_baseline(initial_model, data, dataerr)
    expected_residuals = (data - initial_model) / dataerr

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.01, 0.01]},
        mode="lm",
        verbose=False,
    )

    np.testing.assert_allclose(captured["residuals"], expected_residuals)
    assert fit.parameters["tmid"] == pytest.approx(true_prior["tmid"], abs=3e-4)
    assert np.isfinite(fit.errors["tmid"])
    assert fit.errors["tmid"] > 0


@pytest.mark.parametrize("time", [
    np.linspace(-0.035, -0.002, 90),
    np.linspace(0.007, 0.040, 90),
])
def test_lc_fitter_lm_recovers_tmid_from_a_single_transit_edge(monkeypatch, tmp_path, time):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    true_prior = prior.copy()
    true_prior["tmid"] = 0.0025
    dataerr = np.full_like(time, 2e-4)
    data = elca.transit(time, true_prior)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        np.zeros_like(time),
        prior.copy(),
        {"tmid": [-0.01, 0.01]},
        mode="lm",
        verbose=False,
    )

    assert fit.parameters["tmid"] == pytest.approx(true_prior["tmid"], abs=2e-5)
    assert np.isfinite(fit.errors["tmid"])


def test_lc_fitter_explicit_a0_tracks_mean_airmass_normalization(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    prior["a2"] = -0.35
    time = np.linspace(-0.03, 0.03, 301)
    airmass = np.linspace(1.15, 1.85, len(time))
    dataerr = np.full_like(time, 1e-3)
    true_a0 = 0.985
    data = true_a0 * elca.airmass_trend(prior["a2"], airmass) * elca.transit(time, prior)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]},
        mode="lm",
        verbose=False,
    )

    assert fit.airmass_reference == pytest.approx(np.mean(airmass))
    assert fit.parameters["a0"] == pytest.approx(true_a0, abs=1e-4)
    assert fit.parameters["a1"] == pytest.approx(true_a0, abs=1e-4)


def test_lc_fitter_auto_solves_baseline_when_a0_is_not_free(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 301)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = 1.03 * elca.transit(time, prior)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005]},
        mode="lm",
        verbose=False,
    )

    oot_mask = np.abs(fit.time - fit.parameters["tmid"]) > 0.02
    assert fit.parameters["a0"] == pytest.approx(1.03, abs=1e-4)
    assert fit.parameters["a1"] == pytest.approx(1.03, abs=1e-4)
    assert np.median(fit.detrended[oot_mask]) == pytest.approx(1.0, abs=5e-4)


def test_lc_fitter_falls_back_to_lm_after_nested_linalg_error(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    calls = []

    def fake_fit_nested(self):
        calls.append("ns")
        raise np.linalg.LinAlgError("Singular matrix")

    def fake_fit_LM(self):
        calls.append("lm")
        self.parameters = self.prior.copy()
        self.errors = {}
        self.quantiles = {}
        self.sampled_keys = []
        self.sample_bounds = {}
        self.sample_parameters = {}
        self.sample_errors = {}
        self.sample_quantiles = {}

    monkeypatch.setattr(elca.lc_fitter, "fit_nested", fake_fit_nested)
    monkeypatch.setattr(elca.lc_fitter, "fit_LM", fake_fit_LM)

    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 31)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = elca.transit(time, prior)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005]},
        mode="ns",
        verbose=False,
    )

    assert calls == ["ns", "lm"]
    assert fit.mode == "lm"
    assert fit.ns_type == "lm"
    assert fit.nested_fit_fallback is True
    assert fit.nested_fit_failure_reason == "LinAlgError: Singular matrix"


def test_lc_fitter_auto_solves_mean_airmass_normalization(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    prior["a2"] = 0.22
    time = np.linspace(-0.03, 0.03, 301)
    airmass = np.linspace(1.05, 1.75, len(time))
    dataerr = np.full_like(time, 1e-3)
    true_a0 = 1.018
    data = true_a0 * elca.airmass_trend(prior["a2"], airmass) * elca.transit(time, prior)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005]},
        mode="lm",
        verbose=False,
    )

    assert fit.parameters["a0"] == pytest.approx(true_a0, abs=1e-4)
    assert fit.parameters["a1"] == pytest.approx(true_a0, abs=1e-4)


def test_create_fit_variables_solves_baseline_from_out_of_transit_mask(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 301)
    dataerr = np.full_like(time, 1e-3)
    airmass = np.zeros_like(time)
    transit_model = elca.transit(time, prior)
    data = 1.02 * transit_model
    in_transit = np.abs(time - prior["tmid"]) < 0.012
    data[in_transit] *= 0.90

    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.time = time
    fit.data = data
    fit.dataerr = dataerr
    fit.airmass = airmass
    fit.airmass_reference = elca.get_airmass_reference(airmass)
    fit.prior = prior.copy()
    fit.bounds = {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005]}
    fit.mode = "ns"
    fit.parameters = prior.copy()
    fit.errors = {"rprs": 0.0, "tmid": 0.0, "a2": 0.0}
    fit.quantiles = {}
    fit.baseline_fit_mask = ~in_transit
    fit.fixed_parameter_errors = {}

    fit.create_fit_variables()

    assert fit.parameters["a0"] == pytest.approx(1.02, abs=1e-5)
    assert fit.parameters["a1"] == pytest.approx(1.02, abs=1e-5)


def test_lc_fitter_rejects_redundant_a0_and_a1_bounds(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    prior["a1"] = 1.0
    time = np.linspace(-0.03, 0.03, 31)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = elca.transit(time, prior)

    with pytest.raises(ValueError, match="Use only one of 'a0' or 'a1'"):
        elca.lc_fitter(
            time,
            data,
            dataerr,
            airmass,
            prior,
            {"a0": [0.95, 1.05], "a1": [0.95, 1.05]},
            mode="lm",
            verbose=False,
        )


def test_create_fit_variables_preserves_explicit_baseline_in_nested_mode(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    sampled = make_prior()
    sampled["rprs"] = 0.09
    sampled["a0"] = 0.98
    sampled["a1"] = 0.98

    truth = make_prior()
    truth["rprs"] = 0.12

    fit.time = np.linspace(-0.03, 0.03, 301)
    fit.data = elca.transit(fit.time, truth)
    fit.dataerr = np.full_like(fit.time, 1e-3)
    fit.airmass = np.zeros_like(fit.time)
    fit.prior = sampled.copy()
    fit.bounds = {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]}
    fit.mode = "ns"
    fit.parameters = sampled.copy()
    fit.errors = {"rprs": 1e-3, "tmid": 1e-4, "a0": 2e-3}

    fit.create_fit_variables()

    assert fit.parameters["a0"] == pytest.approx(0.98, abs=1e-9)
    assert fit.parameters["a1"] == pytest.approx(0.98, abs=1e-9)


def test_create_fit_variables_respects_plot_time_range(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    prior = make_prior()
    fit.time = np.array([-0.015, 0.010], dtype=float)
    fit.data = 0.99 * elca.transit(fit.time, prior)
    fit.dataerr = np.full_like(fit.time, 1e-3)
    fit.airmass = np.zeros_like(fit.time)
    fit.prior = prior.copy()
    fit.bounds = {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]}
    fit.mode = "ns"
    fit.parameters = prior.copy()
    fit.errors = {"rprs": 1e-3, "tmid": 1e-4, "a0": 2e-3}
    fit.plot_time_range = (-0.12, 0.18)

    fit.create_fit_variables()

    assert fit.time_upsample[0] == pytest.approx(-0.12, abs=1e-12)
    assert fit.time_upsample[-1] == pytest.approx(0.18, abs=1e-12)
    assert fit.phase_upsample[0] == pytest.approx(-0.04, abs=1e-12)
    assert fit.phase_upsample[-1] == pytest.approx(0.06, abs=1e-12)


def test_plot_bestfit_uses_full_plot_time_range_for_phase_xlim(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.015, 0.010, 51)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = 0.99 * elca.transit(time, prior)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]},
        mode="lm",
        verbose=False,
    )
    fit.plot_time_range = (-0.12, 0.18)
    fit._update_plot_geometry()

    fig, axes = fit.plot_bestfit()

    assert axes[0].get_xlim() == pytest.approx((-0.04, 0.06), abs=1e-6)
    assert axes[1].get_xlim() == pytest.approx((-0.04, 0.06), abs=1e-6)
    assert len(axes[1].child_axes) == 1
    assert axes[1].child_axes[0].get_xlabel() == "Time [hours]"
    plt.close(fig)


def test_plot_bestfit_adds_half_hour_offsets_between_flux_and_residuals(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.015, 0.010, 51)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = 0.99 * elca.transit(time, prior)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]},
        mode="lm",
        verbose=False,
    )
    fit.plot_time_range = (-0.08, 0.08)
    fit._update_plot_geometry()

    fig, axes = fit.plot_bestfit()
    fig.canvas.draw()

    assert axes[1].get_xlabel() == "Phase"
    assert len(axes[1].child_axes) == 1
    time_axis = axes[1].child_axes[0]
    assert time_axis.get_xlabel() == "Time [hours]"
    assert time_axis.get_xlim() == pytest.approx((-1.92, 1.92), abs=1e-5)
    visible_labels = {
        tick.get_text()
        for tick in time_axis.get_xticklabels()
        if tick.get_visible()
    }
    assert {"-1.5", "-1.0", "-0.5", "0", "+0.5", "+1.0", "+1.5"} <= visible_labels
    plt.close(fig)


def test_differential_hours_axis_skips_invalid_period(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    _, axis = plt.subplots()

    assert elca._add_differential_hours_axis(axis, np.nan) is None
    assert axis.child_axes == []
    plt.close(axis.figure)


def test_plot_bestfit_can_hide_flux_baseline_label(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.015, 0.010, 51)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = 0.99 * elca.transit(time, prior)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]},
        mode="lm",
        verbose=False,
    )

    fig, axes = fit.plot_bestfit(show_flux_baseline_label=False)
    legend_text = "\n".join(text.get_text() for text in axes[0].get_legend().get_texts())

    assert "$a_0$" not in legend_text
    plt.close(fig)


def test_format_value_error_for_plot_preserves_two_sigfig_uncertainty_places(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)

    assert elca.format_value_error_for_plot(0.027, 0.05) == ("0.027", "0.050")
    assert elca.format_value_error_for_plot(2461209.81, 0.087) == ("2461209.810", "0.087")
    assert elca.format_value_error_for_plot(89.3511, 2.16) == ("89.4", "2.2")


def test_plot_bestfit_marks_prior_rprs_fallback_uncertainty(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.015, 0.010, 51)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = 0.99 * elca.transit(time, prior)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"tmid": [-0.005, 0.005], "a0": [0.95, 1.05]},
        mode="lm",
        verbose=False,
        fixed_parameter_errors={"rprs": 0.02},
    )
    fit.rprs_prior_fallback_applied = True
    fit.empirical_transit_uncertainty = {
        "available": True,
        "combined_rprs_uncertainty": 0.02,
    }

    fig, axes = fit.plot_bestfit(show_flux_baseline_label=False)
    legend_text = "\n".join(text.get_text() for text in axes[0].get_legend().get_texts())

    assert "(Prior)" in legend_text
    assert "0.0100" in legend_text
    assert "0.0040" in legend_text
    plt.close(fig)


def test_plot_bestfit_can_draw_transit_model_uncertainty_band(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.015, 0.010, 51)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = 0.99 * elca.transit(time, prior)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]},
        mode="lm",
        verbose=False,
    )
    fit.errors["rprs"] = 0.01
    fit.errors["tmid"] = 0.001
    envelope = fit.transit_model_uncertainty(fit.time_upsample)

    fig, axes = fit.plot_bestfit(show_model_uncertainty=True)
    labels = [artist.get_label() for artist in axes[0].collections]
    legend_text = "\n".join(text.get_text() for text in axes[0].get_legend().get_texts())
    uncertainty_line_count = sum(1 for line in axes[0].lines if line.get_linestyle() == "--")

    assert envelope is not None
    assert np.nanmax(envelope[1] - envelope[0]) > 0
    assert "_nolegend_" in labels
    assert r'1-$\sigma$ model uncertainty' not in legend_text
    assert uncertainty_line_count >= 2
    plt.close(fig)


def test_plot_bestfit_draws_unbinned_points_black_with_grey_errorbars(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.015, 0.010, 51)
    airmass = np.zeros_like(time)
    dataerr = np.linspace(5e-4, 4e-3, time.size)
    data = 0.99 * elca.transit(time, prior)
    captured_errorbars = []

    original_errorbar = Axes.errorbar

    def spy_errorbar(self, *args, **kwargs):
        captured_errorbars.append(kwargs.copy())
        return original_errorbar(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "errorbar", spy_errorbar)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]},
        mode="lm",
        verbose=False,
    )

    fig, _ = fit.plot_bestfit()

    assert captured_errorbars[0]["color"] == "black"
    assert captured_errorbars[0]["ecolor"] == "0.72"
    assert captured_errorbars[0]["alpha"] == 1.0
    plot_yerr = np.asarray(captured_errorbars[0]["yerr"], dtype=float)
    assert np.ptp(plot_yerr) > 0
    np.testing.assert_allclose(plot_yerr, fit.detrendederr)
    plt.close(fig)


def test_plot_bestfit_draws_restricted_baseline_points_blue(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.015, 0.010, 51)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = 0.99 * elca.transit(time, prior)
    captured_errorbars = []

    original_errorbar = Axes.errorbar

    def spy_errorbar(self, *args, **kwargs):
        captured_errorbars.append(kwargs.copy())
        return original_errorbar(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "errorbar", spy_errorbar)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]},
        mode="lm",
        verbose=False,
    )
    fit.restricted_baseline_points = {
        "times": np.array([-0.03, 0.03]),
        "flux": np.array([1.0, 1.0]),
        "unc": np.array([1e-3, 1e-3]),
    }

    fig, _ = fit.plot_bestfit()

    assert any(
        errorbar.get("color") == "#1565c0"
        and errorbar.get("label") == "Excluded baseline points"
        for errorbar in captured_errorbars
    )
    plt.close(fig)


def test_plot_bestfit_can_hide_restricted_and_binned_blue_points(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.015, 0.010, 51)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = 0.99 * elca.transit(time, prior)
    captured_errorbars = []
    captured_scatter = []

    original_errorbar = Axes.errorbar

    def spy_errorbar(self, *args, **kwargs):
        captured_errorbars.append(kwargs.copy())
        return original_errorbar(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "errorbar", spy_errorbar)
    original_scatter = Axes.scatter

    def spy_scatter(self, *args, **kwargs):
        captured_scatter.append(kwargs.copy())
        return original_scatter(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "scatter", spy_scatter)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]},
        mode="lm",
        verbose=False,
    )
    fit.restricted_baseline_points = {
        "times": np.array([-0.03, 0.03]),
        "flux": np.array([1.0, 1.0]),
        "unc": np.array([1e-3, 1e-3]),
        "rejected_times": np.array([-0.04]),
        "rejected_flux": np.array([0.75]),
        "rejected_unc": np.array([1e-3]),
    }

    fig, _ = fit.plot_bestfit(
        show_restricted_baseline_points=False,
        show_binned_points=False,
    )

    assert not any(
        errorbar.get("color") in {"blue", "#1565c0"}
        for errorbar in captured_errorbars
    )
    assert any(
        scatter.get("color") == "#d62728" and scatter.get("marker") == "x"
        for scatter in captured_scatter
    )
    plt.close(fig)


def test_plot_bestfit_draws_prefit_rejected_baseline_points_as_red_crosses(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.015, 0.010, 51)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = 0.99 * elca.transit(time, prior)
    captured_scatter = []

    original_scatter = Axes.scatter

    def spy_scatter(self, *args, **kwargs):
        captured_scatter.append(kwargs.copy())
        return original_scatter(self, *args, **kwargs)

    monkeypatch.setattr(Axes, "scatter", spy_scatter)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]},
        mode="lm",
        verbose=False,
    )
    fit.restricted_baseline_points = {
        "times": np.array([-0.03]),
        "flux": np.array([1.0]),
        "unc": np.array([1e-3]),
        "rejected_times": np.array([-0.04]),
        "rejected_flux": np.array([0.75]),
        "rejected_unc": np.array([1e-3]),
    }

    fig, _ = fit.plot_bestfit()

    assert any(
        scatter.get("color") == "#d62728" and scatter.get("marker") == "x"
        for scatter in captured_scatter
    )
    plt.close(fig)


def test_transit_model_uncertainty_includes_baseline_terms(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 51)

    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.time = time
    fit.data = elca.transit(time, prior)
    fit.dataerr = np.full_like(time, 1e-3)
    fit.airmass = np.linspace(1.0, 1.5, time.size)
    fit.airmass_reference = elca.get_airmass_reference(fit.airmass)
    fit.prior = prior.copy()
    fit.bounds = {}
    fit.mode = "ns"
    fit.parameters = prior.copy()
    fit.parameters["a0"] = 1.0
    fit.parameters["a1"] = 1.0
    fit.parameters["a2"] = 0.1
    fit.errors = {"a0": 0.01, "a1": 0.01, "a2": 0.05}
    fit.quantiles = {}
    fit.results = None

    envelope = fit.transit_model_uncertainty(time)

    assert envelope is not None
    assert np.nanmax(envelope[1] - envelope[0]) > 0


def test_baseline_model_uncertainty_is_centered_on_unity_and_includes_a2(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 51)

    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.time = time
    fit.airmass = np.linspace(1.0, 2.0, time.size)
    fit.airmass_reference = elca.get_airmass_reference(fit.airmass)
    fit.parameters = prior.copy()
    fit.parameters["a0"] = 1.0
    fit.parameters["a1"] = 1.0
    fit.parameters["a2"] = 0.1
    fit.errors = {"a0": 0.01, "a1": 0.01, "a2": 0.05}

    lower, upper = fit.baseline_model_uncertainty(time)
    width = upper - lower

    np.testing.assert_allclose(0.5 * (lower + upper), np.ones_like(time), atol=1e-12)
    assert np.nanmin(lower) < 1.0
    assert np.nanmax(upper) > 1.0
    assert width[0] > width[len(width) // 2]


def test_baseline_model_uncertainty_does_not_double_count_analytic_a0(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 51)
    airmass = np.linspace(1.0, 2.0, time.size)

    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.time = time
    fit.data = elca.transit(time, prior)
    fit.dataerr = np.full_like(time, 1e-3)
    fit.airmass = airmass
    fit.airmass_reference = elca.get_airmass_reference(fit.airmass)
    fit.prior = prior.copy()
    fit.bounds = {"a2": [-1.0, 1.0]}
    fit.fixed_flux_baseline = False
    fit.mode = "ns"
    fit.parameters = prior.copy()
    fit.parameters["a0"] = 1.0
    fit.parameters["a1"] = 1.0
    fit.parameters["a2"] = 0.0
    fit.errors = {"a0": 0.5, "a1": 0.5, "a2": 0.02}
    fit.results = None

    lower, upper = fit.baseline_model_uncertainty(time)
    half_width = np.nanmax(np.maximum(1.0 - lower, upper - 1.0))

    assert half_width < 0.03


def test_baseline_model_uncertainty_includes_empirical_flux_floor(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 51)

    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.time = time
    fit.airmass = np.linspace(1.0, 2.0, time.size)
    fit.airmass_reference = elca.get_airmass_reference(fit.airmass)
    fit.parameters = prior.copy()
    fit.parameters["a0"] = 1.0
    fit.parameters["a1"] = 1.0
    fit.parameters["a2"] = 0.0
    fit.errors = {"a0": 0.001, "a1": 0.001, "a2": 0.001}
    fit.results = None
    fit.empirical_transit_uncertainty = {
        "available": True,
        "baseline_red_noise_uncertainty_fraction": 0.02,
    }

    lower, upper = fit.baseline_model_uncertainty(time)
    half_width = np.nanmax(np.maximum(1.0 - lower, upper - 1.0))

    assert half_width >= 0.02


def test_baseline_model_uncertainty_prefers_posterior_samples(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 51)
    sample_count = 41
    a0_samples = 1.0 + np.linspace(-0.004, 0.004, sample_count)
    a2_samples = np.linspace(-0.02, 0.02, sample_count)

    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.time = time
    fit.data = elca.transit(time, prior)
    fit.dataerr = np.full_like(time, 1e-3)
    fit.airmass = np.linspace(1.0, 2.0, time.size)
    fit.airmass_reference = elca.get_airmass_reference(fit.airmass)
    fit.prior = prior.copy()
    fit.bounds = {"a0": [0.5, 1.5], "a2": [-1.0, 1.0]}
    fit.sampled_keys = ["a0", "a2"]
    fit.mode = "ns"
    fit.ns_type = "ultranest"
    fit.parameters = prior.copy()
    fit.parameters["a0"] = 1.0
    fit.parameters["a1"] = 1.0
    fit.parameters["a2"] = 0.0
    fit.errors = {"a0": 0.5, "a1": 0.5, "a2": 0.5}
    fit.results = {
        "weighted_samples": {
            "points": np.column_stack([a0_samples, a2_samples]),
            "logl": np.zeros(sample_count, dtype=float),
            "weights": np.ones(sample_count, dtype=float),
        }
    }

    lower, upper = fit.baseline_model_uncertainty(time)
    half_width = np.nanmax(np.maximum(1.0 - lower, upper - 1.0))

    assert half_width < 0.02


def test_plot_bestfit_can_draw_baseline_uncertainty_band(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 51)
    airmass = np.linspace(1.0, 2.0, time.size)
    dataerr = np.full_like(time, 1e-3)
    data = 0.99 * elca.transit(time, prior)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "tmid": [-0.005, 0.005], "a0": [0.95, 1.05]},
        mode="lm",
        verbose=False,
    )
    fit.parameters["a2"] = 0.1
    fit.errors["a0"] = 0.01
    fit.errors["a2"] = 0.05

    fig, axes = fit.plot_bestfit(show_baseline_uncertainty=True)
    labels = [artist.get_label() for artist in axes[0].collections]
    legend_text = "\n".join(text.get_text() for text in axes[0].get_legend().get_texts())
    baseline_line_count = sum(
        1
        for line in axes[0].lines
        if line.get_linestyle() == "--" and line.get_color() == "gold"
    )

    assert "_nolegend_" in labels
    assert r'$a_0/a_2$ 1-$\sigma$ baseline uncertainty' not in legend_text
    assert baseline_line_count == 2
    plt.close(fig)


def test_posterior_model_uncertainty_recenters_on_best_fit_model(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 51)

    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.time = time
    fit.data = elca.transit(time, prior)
    fit.dataerr = np.full_like(time, 1e-3)
    fit.airmass = np.zeros_like(time)
    fit.airmass_reference = elca.get_airmass_reference(fit.airmass)
    fit.prior = prior.copy()
    fit.bounds = {"a0": [0.99, 1.03]}
    fit.sampled_keys = ["a0"]
    fit.sample_bounds = {"a0": [0.99, 1.03]}
    fit.mode = "ns"
    fit.ns_type = "ultranest"
    fit.parameters = prior.copy()
    fit.parameters["a0"] = 1.0
    fit.parameters["a1"] = 1.0
    fit.errors = {"a0": 0.002}
    fit.quantiles = {}
    fit.results = {
        "weighted_samples": {
            "points": np.linspace(1.008, 1.012, 41)[:, None],
            "logl": np.zeros(41, dtype=float),
            "weights": np.ones(41, dtype=float),
        }
    }

    lower, upper = fit.transit_model_uncertainty(time)
    center = 0.5 * (lower + upper)

    np.testing.assert_allclose(center, elca.transit(time, fit.parameters), atol=5e-5)
    assert np.nanmedian(center[:3]) < 1.001


def test_glc_plot_bestfit_median_limits_use_full_phase_span(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    phase = np.array([0.01, 0.02], dtype=float)
    phase_upsample = np.linspace(-0.04, 0.06, 100)
    residuals = np.array([1e-4, -1e-4], dtype=float)
    times = prior["tmid"] + phase * prior["per"]

    fit = elca.glc_fitter.__new__(elca.glc_fitter)
    fit.parameters = prior.copy()
    fit.errors = {"rprs": 1e-3, "tmid": 1e-4}
    fit.lc_data = [{
        "time": times,
        "flux": np.ones_like(times),
        "detrend": np.ones_like(times),
        "ferr": np.full_like(times, 1e-3),
        "residuals": residuals,
        "phase": phase,
        "phase_upsample": phase_upsample,
        "time_upsample": prior["tmid"] + phase_upsample * prior["per"],
        "transit_upsample": np.ones_like(phase_upsample),
        "priors": prior.copy(),
        "errors": {"rprs": 1e-3, "tmid": 1e-4},
        "name": "dataset",
    }]

    fig, axes = fit.plot_bestfit(phase_limits="median")

    assert axes[0].get_xlim() == pytest.approx((-0.04, 0.06), abs=1e-6)
    assert axes[1].get_xlim() == pytest.approx((-0.04, 0.06), abs=1e-6)
    assert len(axes[1].child_axes) == 1
    assert axes[1].child_axes[0].get_xlabel() == "Time [hours]"
    plt.close(fig)


def test_plot_triangle_clips_ranges_to_parameter_bounds(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    captured = {}

    def fake_corner(*args, **kwargs):
        captured["points"] = args[0]
        captured["labels"] = kwargs["labels"]
        captured["range"] = kwargs["range"]
        captured["titles"] = kwargs["titles"]
        captured["title_kwargs"] = kwargs["title_kwargs"]
        captured["label_kwargs"] = kwargs["label_kwargs"]
        return "figure"

    monkeypatch.setattr(elca, "corner", fake_corner)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.125],
        "inc": [84.0, 90.0],
        "a0": [0.95, 1.05],
    }
    fit.prior = make_prior()
    fit.quantiles = {"rprs": [], "inc": [], "a0": []}
    fit.parameters = {"rprs": 0.10, "inc": 88.42, "a0": 0.94962}
    fit.errors = {"rprs": 0.01, "inc": 0.75, "a0": 0.00394}

    points = np.array(
        [
            [0.099, 88.30, 0.9501],
            [0.101, 88.55, 0.9502],
            [0.102, 88.10, 0.9515],
            [0.098, 88.70, 0.9520],
            [0.100, 88.40, 0.9508],
        ]
    )
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.array([-5.0, -4.0, -4.5, -5.5, -4.2]),
        },
        "samples": points.copy(),
    }

    fig = fit.plot_triangle()

    assert fig == "figure"
    assert captured["labels"][1] == r"$\Delta i$"
    assert captured["range"][0][0] == pytest.approx(0.0)
    assert captured["range"][0][1] == pytest.approx(0.125)
    expected_inc_distance_limit = np.max(np.abs(np.array([84.0, 90.0]) - fit.parameters["inc"]))
    assert captured["range"][1][0] == pytest.approx(-expected_inc_distance_limit)
    assert captured["range"][1][1] == pytest.approx(expected_inc_distance_limit)
    assert captured["range"][2][0] == pytest.approx(0.95)
    assert captured["range"][2][1] == pytest.approx(1.05)
    assert captured["points"].shape == (10, 3)
    expected_inc_distance = np.abs(points[:, 1] - fit.parameters["inc"])
    np.testing.assert_allclose(captured["points"][:5, 1], expected_inc_distance)
    np.testing.assert_allclose(captured["points"][5:, 1], -expected_inc_distance)
    assert captured["titles"][1].startswith("b=")
    assert "\ni=" in captured["titles"][1]
    assert captured["title_kwargs"]["loc"] == "left"
    assert captured["label_kwargs"]["labelpad"] == 10


def test_internal_impact_parameter_transform_round_trips_inclination(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = True
    fit.prior = make_prior()
    fit.bounds = {
        "rprs": [0.08, 0.12],
        "inc": [87.0, 90.0],
        "tmid": [-0.005, 0.005],
    }

    sample_point = fit._sample_point_from_unit_cube(np.array([0.25, 0.4, 0.75]))
    physical = fit._physical_values_from_sample_point(sample_point)
    expected_rprs = 0.08 + 0.25 * (0.12 - 0.08)
    expected_b = 0.4 * (1.0 + expected_rprs)
    expected_inc = float(elca.inclination_from_impact_parameter(
        {**fit.prior, "rprs": expected_rprs},
        expected_b,
    ))

    assert fit._get_sampled_keys() == ["rprs", "b", "tmid"]
    assert sample_point[1] == pytest.approx(expected_b)
    assert physical["inc"] == pytest.approx(expected_inc)
    assert physical["b"] == pytest.approx(sample_point[1])


def test_internal_impact_parameter_samples_grazing_range_beyond_one(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = True
    fit.prior = make_prior()
    fit.bounds = {
        "rprs": [0.08, 0.12],
        "inc": [89.8, 90.0],
        "tmid": [-0.005, 0.005],
    }

    sample_point = fit._sample_point_from_unit_cube(np.array([1.0, 0.99, 0.5]))
    physical = fit._physical_values_from_sample_point(sample_point)

    assert sample_point[0] == pytest.approx(0.12)
    assert sample_point[1] == pytest.approx(0.99 * 1.12)
    assert sample_point[1] > 1.0
    assert physical["inc"] < 90.0
    assert fit._get_sample_bounds()["b"] == pytest.approx([0.0, 1.12])


def test_internal_impact_parameter_transform_does_not_deepcopy_prior(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = True
    fit.prior = make_prior()
    fit.bounds = {
        "rprs": [0.08, 0.12],
        "ars": [10.0, 14.0],
        "inc": [87.0, 90.0],
        "tmid": [-0.005, 0.005],
    }

    def fail_deepcopy(value, memo=None):
        raise AssertionError("sampling transforms should not deepcopy parameter dictionaries")

    monkeypatch.setattr(elca.copy, "deepcopy", fail_deepcopy)

    sample_point = fit._sample_point_from_unit_cube(np.array([0.25, 0.5, 0.4, 0.75]))
    unit_points = np.array([
        [0.25, 0.5, 0.4, 0.75],
        [1.0, 0.25, 0.9, 0.5],
    ])
    sample_points = fit._sample_point_from_unit_cube(unit_points)
    physical = fit._physical_values_from_sample_point(sample_point)
    sample_bounds = fit._get_sample_bounds()

    assert sample_point[0] == pytest.approx(0.09)
    np.testing.assert_allclose(
        sample_points,
        np.vstack([fit._sample_point_from_unit_cube(row) for row in unit_points]),
    )
    assert physical["b"] == pytest.approx(sample_point[2])
    assert "b" in sample_bounds


def test_unit_cube_transform_vectorizes_simple_bounds(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = False
    fit.prior = make_prior()
    fit.bounds = {
        "rprs": [0.08, 0.12],
        "inc": [87.0, 90.0],
        "tmid": [-0.005, 0.005],
    }

    unit_points = np.array([
        [0.0, 0.5, 1.0],
        [1.0, 0.25, 0.0],
    ])

    np.testing.assert_allclose(
        fit._sample_point_from_unit_cube(unit_points),
        np.array([
            [0.08, 88.5, 0.005],
            [0.12, 87.75, -0.005],
        ]),
    )


def test_unit_cube_inverse_maps_expanded_prior_samples_back_to_unit_cube(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = False
    fit.prior = make_prior()
    fit.bounds = {
        "rprs": [0.05, 0.15],
        "tmid": [-0.01, 0.01],
    }

    unit_points = np.array([
        [0.30, 0.25],
        [0.70, 0.75],
    ])
    sample_points = fit._sample_point_from_unit_cube(unit_points)

    np.testing.assert_allclose(
        fit._unit_cube_from_sample_points(sample_points),
        unit_points,
    )


def test_expanded_prior_warmstart_uses_corrected_guarded_auxiliary_problem(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 81)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = elca.transit(time, prior)
    source_points = np.column_stack([
        np.linspace(0.085, 0.115, 64),
        np.linspace(-0.004, 0.004, 64),
    ])
    source = make_expanded_prior_warmstart_source(
        prior,
        time,
        data,
        dataerr,
        airmass,
        source_points,
    )
    captured = {}

    class DummySampler:
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs
            captured["sampler"] = self

    def fake_run_reactive_sampler(sampler, *args, **kwargs):
        unit_values = np.array([
            [0.4, 0.6, 0.25],
            [0.4, 0.6, 0.75],
        ])
        transformed = np.asarray(sampler.args[2](unit_values), dtype=float)
        captured["transformed"] = transformed
        captured["physical_loglike"] = np.asarray(
            [
                sampler.args[1](np.append(row[:2], 0.0))
                for row in transformed
            ],
            dtype=float,
        )
        captured["corrected_loglike"] = np.asarray(
            sampler.args[1](transformed),
            dtype=float,
        )
        return make_dummy_nested_result(source_points, auxiliary=True)

    monkeypatch.setattr(elca, "ReactiveNestedSampler", DummySampler)
    monkeypatch.setattr(elca, "run_reactive_sampler", fake_run_reactive_sampler)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {
            "rprs": [0.05, 0.15],
            "tmid": [-0.005, 0.005],
        },
        mode="ns",
        verbose=False,
        use_impactparameter_rather_than_inclination_to_fit=False,
        ultranest_warmstart_source=source,
    )

    assert captured["sampler"].args[0] == ["rprs", "tmid", "aux_logweight"]
    np.testing.assert_allclose(captured["transformed"][0, :2], [0.09, 0.001])
    assert np.all(np.isfinite(captured["transformed"]))
    assert captured["transformed"][1, 0] > captured["transformed"][0, 0]
    np.testing.assert_allclose(
        captured["corrected_loglike"] - captured["physical_loglike"],
        captured["transformed"][:, 2],
    )
    assert fit.ultranest_expanded_prior_warmstart_attempted is True
    assert fit.ultranest_expanded_prior_warmstart_applied is True
    assert fit.ultranest_expanded_prior_warmstart_expanded_keys == ["rprs"]
    assert fit.ultranest_expanded_prior_warmstart_source_sample_count == 64
    assert fit._get_triangle_plot_samples()[0].shape[1] == 2


def test_expanded_prior_warmstart_failure_falls_back_to_clean_sampler(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 81)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = elca.transit(time, prior)
    source_points = np.column_stack([
        np.linspace(0.085, 0.115, 64),
        np.linspace(-0.004, 0.004, 64),
    ])
    source = make_expanded_prior_warmstart_source(
        prior,
        time,
        data,
        dataerr,
        airmass,
        source_points,
    )
    samplers = []

    class DummySampler:
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs
            samplers.append(self)

    def fake_run_reactive_sampler(sampler, *args, **kwargs):
        if len(sampler.args[0]) == 3:
            raise RuntimeError("synthetic corrected-warmstart failure")
        return make_dummy_nested_result(source_points, auxiliary=False)

    monkeypatch.setattr(elca, "ReactiveNestedSampler", DummySampler)
    monkeypatch.setattr(elca, "run_reactive_sampler", fake_run_reactive_sampler)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {
            "rprs": [0.05, 0.15],
            "tmid": [-0.005, 0.005],
        },
        mode="ns",
        verbose=False,
        use_impactparameter_rather_than_inclination_to_fit=False,
        ultranest_warmstart_source=source,
    )

    assert [sampler.args[0] for sampler in samplers] == [
        ["rprs", "tmid", "aux_logweight"],
        ["rprs", "tmid"],
    ]
    assert fit.ultranest_expanded_prior_warmstart_attempted is True
    assert fit.ultranest_expanded_prior_warmstart_applied is False
    assert "reran from the full expanded prior" in fit.ultranest_expanded_prior_warmstart_note
    assert fit.parameters["rprs"] == pytest.approx(0.1)


def test_expanded_prior_warmstart_restores_physical_likelihood_for_best_fit(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    results = {
        "maximum_likelihood": {
            "point": np.array([0.09, 0.0, 5.0]),
            "logl": 5.0,
        },
        "weighted_samples": {
            "points": np.array([
                [0.09, 0.0, 5.0],
                [0.11, 0.0, -5.0],
            ]),
            "logl": np.array([5.0, -5.0]),
        },
    }

    def physical_loglike(points):
        points = np.asarray(points, dtype=float)
        return -((points[:, 0] - 0.11) / 0.01) ** 2

    elca.lc_fitter._restore_expanded_prior_physical_likelihoods(
        results,
        physical_loglike,
        2,
    )

    np.testing.assert_allclose(
        results["weighted_samples"]["auxiliary_logl"],
        np.array([5.0, -5.0]),
    )
    np.testing.assert_allclose(
        results["weighted_samples"]["logl"],
        np.array([-4.0, 0.0]),
    )
    assert results["weighted_samples"]["points"].shape == (2, 2)
    assert results["weighted_samples"]["auxiliary_points"].shape == (2, 1)
    assert results["maximum_likelihood"]["point"].shape == (2,)
    assert results["maximum_likelihood"]["point"][0] == pytest.approx(0.11)
    assert results["maximum_likelihood"]["logl"] == pytest.approx(0.0)


def test_nested_fit_can_keep_inclination_parameterization_when_requested(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = False
    fit.prior = make_prior()
    fit.bounds = {
        "rprs": [0.08, 0.12],
        "inc": [87.0, 90.0],
        "tmid": [-0.005, 0.005],
    }

    assert fit._get_sampled_keys() == ["rprs", "inc", "tmid"]


def test_nested_fit_reports_inclination_from_internal_impact_parameter(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 101)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = elca.transit(time, prior)
    data[0] += 5e-5

    class DummySampler:
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs

    b_ml = float(elca.impact_parameter_from_inclination(prior, 88.8))
    sample_points = np.array(
        [
            [0.100, b_ml - 0.02, 0.0000],
            [0.101, b_ml - 0.01, 0.0002],
            [0.099, b_ml + 0.01, -0.0001],
            [0.100, b_ml + 0.02, 0.0001],
        ]
    )

    monkeypatch.setattr(elca, "ReactiveNestedSampler", DummySampler)
    monkeypatch.setattr(
        elca,
        "run_reactive_sampler",
        lambda *args, **kwargs: {
            "maximum_likelihood": {"point": np.array([0.100, b_ml, 0.0])},
            "posterior": {
                "stdev": np.array([0.005, 0.02, 0.0005]),
                "errlo": np.array([-0.005, -0.02, -0.0005]),
                "errup": np.array([0.005, 0.02, 0.0005]),
            },
            "weighted_samples": {
                "points": sample_points,
                "logl": np.array([-4.0, -3.0, -3.2, -3.8]),
            },
            "samples": sample_points.copy(),
        },
    )

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "inc": [87.0, 90.0], "tmid": [-0.005, 0.005]},
        mode="ns",
        verbose=False,
    )

    assert fit.sampled_keys == ["rprs", "b", "tmid"]
    assert fit.sample_parameters["b"] == pytest.approx(b_ml)
    assert fit.parameters["inc"] == pytest.approx(88.8, abs=1e-6)
    assert fit.errors["inc"] > 0


def test_nested_fit_tracks_free_ars_with_internal_impact_parameter(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 101)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = elca.transit(time, prior)

    class DummySampler:
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs

    ml_values = prior.copy()
    ml_values["ars"] = 12.3
    b_ml = float(elca.impact_parameter_from_inclination(ml_values, 88.8))
    sample_points = np.array(
        [
            [0.100, 12.10, b_ml - 0.02, 0.0000],
            [0.101, 12.20, b_ml - 0.01, 0.0002],
            [0.099, 12.40, b_ml + 0.01, -0.0001],
            [0.100, 12.50, b_ml + 0.02, 0.0001],
        ]
    )

    monkeypatch.setattr(elca, "ReactiveNestedSampler", DummySampler)
    monkeypatch.setattr(
        elca,
        "run_reactive_sampler",
        lambda *args, **kwargs: {
            "maximum_likelihood": {"point": np.array([0.100, 12.30, b_ml, 0.0])},
            "posterior": {
                "stdev": np.array([0.005, 0.1, 0.02, 0.0005]),
                "errlo": np.array([-0.005, -0.1, -0.02, -0.0005]),
                "errup": np.array([0.005, 0.1, 0.02, 0.0005]),
            },
            "weighted_samples": {
                "points": sample_points,
                "logl": np.array([-4.0, -3.0, -3.2, -3.8]),
            },
            "samples": sample_points.copy(),
        },
    )

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "ars": [11.5, 12.5], "inc": [87.0, 89.5], "tmid": [-0.005, 0.005]},
        mode="ns",
        verbose=False,
    )

    assert fit.sampled_keys == ["rprs", "ars", "b", "tmid"]
    assert fit.parameters["ars"] == pytest.approx(12.3, abs=1e-12)
    assert fit.parameters["inc"] == pytest.approx(88.8, abs=1e-6)
    assert fit.sample_bounds["b"] == pytest.approx([0.0, 1.12])


def test_nested_fit_replaces_degenerate_ultranest_errors_from_loglike_neighborhood(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 101)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = elca.transit(time, prior)

    class DummySampler:
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs

    sample_points = np.array(
        [
            [0.095, -0.0010],
            [0.097, -0.0008],
            [0.099, -0.0003],
            [0.100, 0.0000],
            [0.101, 0.0002],
            [0.103, 0.0005],
            [0.105, 0.0008],
            [0.106, 0.0010],
            [0.120, 0.0030],
            [0.080, -0.0030],
        ],
        dtype=float,
    )
    logl = np.array([-0.4, -0.3, -0.1, 0.0, -0.1, -0.2, -0.3, -0.4, -2.0, -3.0])

    monkeypatch.setattr(elca, "ReactiveNestedSampler", DummySampler)
    monkeypatch.setattr(
        elca,
        "run_reactive_sampler",
        lambda *args, **kwargs: {
            "maximum_likelihood": {"point": np.array([0.100, 0.0])},
            "posterior": {
                "stdev": np.array([1e-15, 1e-15]),
                "errlo": np.array([0.100, 0.0]),
                "errup": np.array([0.100, 0.0]),
            },
            "weighted_samples": {
                "points": sample_points,
                "logl": logl,
            },
            "samples": np.repeat(np.array([[0.100, 0.0]]), 10, axis=0),
        },
    )

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.0, 0.2], "tmid": [-0.01, 0.01]},
        mode="ns",
        verbose=False,
    )

    assert fit.errors["rprs"] > 1e-3
    assert fit.errors["tmid"] > 1e-4
    assert set(fit.ultranest_error_fallbacks) == {"rprs", "tmid"}
    assert fit.ultranest_error_fallbacks["rprs"]["sample_count"] == 8


def test_nested_fit_replaces_prior_width_like_error_with_local_likelihood_width(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.prior = make_prior()
    fit.bounds = {"rprs": [0.0, 0.3]}
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = True
    fit.fixed_parameter_errors = {}

    center = 0.152
    broad_points = np.linspace(0.0, 0.3, 40)
    local_points = center + np.linspace(-0.006, 0.006, 17)
    points = np.concatenate([broad_points, local_points])[:, None]
    broad_logl = np.full(broad_points.shape, -100.0)
    local_logl = -0.5 * ((local_points - center) / 0.0038) ** 2
    logl = np.concatenate([broad_logl, local_logl])
    fit.results = {
        "maximum_likelihood": {"point": np.array([center])},
        "posterior": {
            "stdev": np.array([0.082]),
            "errlo": np.array([-0.082]),
            "errup": np.array([0.082]),
        },
        "weighted_samples": {
            "points": points,
            "logl": logl,
        },
        "samples": points.copy(),
    }

    fit._finalize_ultranest_fit_results(
        ["rprs"],
        ["rprs"],
        lambda point: {"rprs": float(point[0])},
    )

    fallback = fit.ultranest_error_fallbacks["rprs"]
    assert fit.errors["rprs"] < 0.01
    assert fit.errors["rprs"] == pytest.approx(fallback["error"])
    assert fallback["reported_error"] == pytest.approx(0.082)
    assert fallback["reason"] == "posterior_summary_inflated_relative_to_local_fit"
    assert fallback["delta_chi2"] <= 1.0


def _make_gaussian_nested_fit(elca, sigmas, seed=1401, count=4000):
    """A healthy d-dimensional Gaussian posterior as UltraNest would report it:
    weighted_samples drawn from the posterior, logl = -chi2/2, stdev = the true sigmas."""
    keys = ["tmid", "rprs", "ars", "b"][: len(sigmas)]
    sigmas = np.asarray(sigmas, dtype=float)
    rng = np.random.default_rng(seed)
    z = rng.standard_normal((count, sigmas.size))
    points = z * sigmas
    logl = -0.5 * np.sum(z ** 2, axis=1)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.prior = make_prior()
    fit.bounds = {key: [-10.0 * sigma, 10.0 * sigma] for key, sigma in zip(keys, sigmas)}
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = True
    fit.fixed_parameter_errors = {}
    fit.results = {
        "maximum_likelihood": {"point": np.zeros(sigmas.size)},
        "posterior": {
            "stdev": sigmas.copy(),
            "errlo": -sigmas.copy(),
            "errup": sigmas.copy(),
        },
        "weighted_samples": {"points": points, "logl": logl},
        "samples": points.copy(),
    }
    return fit, keys, sigmas


def test_loglike_neighborhood_uncertainty_recovers_marginal_sigma_in_four_dimensions(monkeypatch, tmp_path):
    # Issue #1401: the std of the delta-chi2<=1 points is ~sigma/sqrt(d+2) in d
    # dimensions, so it under-reported a healthy 4-parameter Tmid bar by ~3x. The
    # profile half-range is the marginal sigma of a Gaussian posterior.
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit, keys, sigmas = _make_gaussian_nested_fit(elca, [0.005, 0.02, 1.0, 0.25])
    for index, sigma in enumerate(sigmas):
        local = fit._loglike_neighborhood_uncertainty(index, 0.0)
        assert local is not None
        assert local["delta_chi2"] <= 1.0
        assert local["interior_std"] < 0.6 * sigma  # the old reference, documented bias
        assert 0.6 * sigma < local["error"] < 1.3 * sigma  # the new reference


def test_nested_fit_keeps_healthy_posterior_summary_in_four_dimensions(monkeypatch, tmp_path):
    # With the biased reference and factor 3.0 this fired on real runs (Tmid 0.0056 -> 0.0018).
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit, keys, sigmas = _make_gaussian_nested_fit(elca, [0.005, 0.02, 1.0, 0.25])
    fit._finalize_ultranest_fit_results(
        keys,
        keys,
        lambda point: {key: float(value) for key, value in zip(keys, point)},
    )
    assert fit.ultranest_error_fallbacks == {}
    for key, sigma in zip(keys, sigmas):
        assert fit.errors[key] == pytest.approx(sigma)


def test_nested_fit_duration_prior_penalizes_wrong_transit_length(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(0.20, 0.30, 51)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = np.ones_like(time)
    data[0] += 1e-4

    class DummySampler:
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs

    captured = {}
    good_point = np.array([prior["rprs"], prior["ars"], prior["inc"], prior["tmid"]], dtype=float)
    bad_point = np.array([prior["rprs"], 30.0, prior["inc"], prior["tmid"]], dtype=float)

    monkeypatch.setattr(elca, "ReactiveNestedSampler", DummySampler)

    def fake_run_reactive_sampler(sampler, *args, **kwargs):
        loglike = sampler.args[1]
        captured["good"] = float(loglike(good_point))
        captured["bad"] = float(loglike(bad_point))
        return {
            "maximum_likelihood": {"point": good_point.copy()},
            "posterior": {
                "stdev": np.array([0.001, 0.1, 0.05, 0.0001]),
                "errlo": np.array([-0.001, -0.1, -0.05, -0.0001]),
                "errup": np.array([0.001, 0.1, 0.05, 0.0001]),
            },
            "weighted_samples": {
                "points": np.vstack([good_point, bad_point]),
                "logl": np.array([captured["good"], captured["bad"]]),
            },
            "samples": np.vstack([good_point, bad_point]),
        }

    monkeypatch.setattr(elca, "run_reactive_sampler", fake_run_reactive_sampler)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "ars": [10.0, 35.0], "inc": [88.5, 89.5], "tmid": [-0.005, 0.005]},
        mode="ns",
        verbose=False,
        use_impactparameter_rather_than_inclination_to_fit=False,
        duration_prior={
            "applied": True,
            "expected_duration": elca.transit_duration(prior),
            "sigma_log_duration": 0.05,
        },
    )

    expected_penalty = -0.5 * (
        np.log(elca.transit_duration({"per": prior["per"], "rprs": prior["rprs"], "ars": 30.0, "inc": prior["inc"], "ecc": prior["ecc"], "omega": prior["omega"]}) / elca.transit_duration(prior))
        / 0.05
    ) ** 2

    assert fit.parameters["ars"] == pytest.approx(prior["ars"], abs=1e-12)
    assert captured["good"] > captured["bad"]
    assert (captured["bad"] - captured["good"]) == pytest.approx(expected_penalty, rel=1e-6, abs=1e-6)


def test_nested_fit_duration_prior_returns_finite_floor_for_invalid_geometry(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    prior = make_prior()
    time = np.linspace(-0.03, 0.03, 51)
    airmass = np.zeros_like(time)
    dataerr = np.full_like(time, 1e-3)
    data = elca.transit(time, prior)

    class DummySampler:
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs

    captured = {}
    good_point = np.array([prior["rprs"], prior["ars"], prior["inc"], prior["tmid"]], dtype=float)
    invalid_point = np.array([prior["rprs"], 50.0, 80.0, prior["tmid"]], dtype=float)

    monkeypatch.setattr(elca, "ReactiveNestedSampler", DummySampler)

    def fake_run_reactive_sampler(sampler, *args, **kwargs):
        loglike = sampler.args[1]
        prior_transform = sampler.args[2]
        captured["invalid"] = float(loglike(invalid_point))
        captured["vector"] = np.asarray(loglike(np.vstack([good_point, invalid_point])), dtype=float)
        captured["transformed"] = prior_transform(np.full((2, 4), 0.5, dtype=float))
        return {
            "maximum_likelihood": {"point": good_point.copy()},
            "posterior": {
                "stdev": np.array([0.001, 0.1, 0.05, 0.0001]),
                "errlo": np.array([-0.001, -0.1, -0.05, -0.0001]),
                "errup": np.array([0.001, 0.1, 0.05, 0.0001]),
            },
            "weighted_samples": {
                "points": np.vstack([good_point, invalid_point]),
                "logl": captured["vector"],
            },
            "samples": np.vstack([good_point, invalid_point]),
        }

    monkeypatch.setattr(elca, "run_reactive_sampler", fake_run_reactive_sampler)

    fit = elca.lc_fitter(
        time,
        data,
        dataerr,
        airmass,
        prior.copy(),
        {"rprs": [0.08, 0.12], "ars": [10.0, 50.0], "inc": [80.0, 89.5], "tmid": [-0.005, 0.005]},
        mode="ns",
        verbose=False,
        use_impactparameter_rather_than_inclination_to_fit=False,
        duration_prior={
            "applied": True,
            "expected_duration": elca.transit_duration(prior),
            "sigma_log_duration": 0.05,
        },
    )

    assert fit.parameters["ars"] == pytest.approx(prior["ars"], abs=1e-12)
    assert np.isfinite(captured["invalid"])
    assert captured["invalid"] == pytest.approx(elca.BAD_LOG_LIKELIHOOD)
    assert captured["vector"].shape == (2,)
    assert np.all(np.isfinite(captured["vector"]))
    assert captured["vector"][1] == pytest.approx(elca.BAD_LOG_LIKELIHOOD)
    assert np.asarray(captured["transformed"]).shape == (2, 4)


def test_rprs_posterior_recenter_diagnostics_detect_upper_bound_clipping(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = True
    fit.prior = make_prior()
    fit.bounds = {"rprs": [0.0, 0.15], "tmid": [-0.005, 0.005]}
    fit.sampled_keys = ["rprs", "tmid"]
    fit.sample_bounds = {"rprs": [0.0, 0.15], "tmid": [-0.005, 0.005]}

    rprs_samples = np.concatenate([
        np.linspace(0.090, 0.120, 12),
        np.linspace(0.128, 0.149, 28),
    ])
    tmid_samples = np.linspace(-2e-4, 2e-4, rprs_samples.size)
    points = np.column_stack([rprs_samples, tmid_samples])
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-6.0, -3.0, rprs_samples.size),
        },
        "samples": points.copy(),
    }

    diagnostics = fit.get_parameter_posterior_recenter_diagnostics("rprs")

    assert diagnostics["clipped"] is True
    assert diagnostics["edge"] == "upper"
    assert diagnostics["mode"] > 0.13
    assert diagnostics["std"] > 0
    assert diagnostics["upper_edge_peak_fraction"] >= 0.20
    assert diagnostics["bounds"][0] >= 0.0
    assert diagnostics["bounds"][1] > 0.15


def test_ars_posterior_recenter_diagnostics_detect_lower_bound_clipping(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = True
    fit.prior = make_prior()
    fit.bounds = {"ars": [10.0, 15.0], "tmid": [-0.005, 0.005]}
    fit.sampled_keys = ["ars", "tmid"]
    fit.sample_bounds = {"ars": [10.0, 15.0], "tmid": [-0.005, 0.005]}

    ars_samples = np.concatenate([
        np.linspace(10.001, 10.040, 30),
        np.linspace(10.060, 10.800, 12),
    ])
    tmid_samples = np.linspace(-2e-4, 2e-4, ars_samples.size)
    points = np.column_stack([ars_samples, tmid_samples])
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-6.0, -3.0, ars_samples.size),
        },
        "samples": points.copy(),
    }

    diagnostics = fit.get_parameter_posterior_recenter_diagnostics("ars")

    assert diagnostics["clipped"] is True
    assert diagnostics["edge"] == "lower"
    assert diagnostics["mode"] < 10.5
    assert diagnostics["std"] > 0
    assert diagnostics["lower_edge_peak_fraction"] >= 0.20
    assert diagnostics["bounds"][0] < 10.0
    assert diagnostics["bounds"][0] >= 0.0


def test_rprs_posterior_recenter_diagnostics_ignores_upper_edge_below_twenty_percent(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = True
    fit.prior = make_prior()
    fit.bounds = {"rprs": [0.0, 0.15], "tmid": [-0.005, 0.005]}
    fit.sampled_keys = ["rprs", "tmid"]
    fit.sample_bounds = {"rprs": [0.0, 0.15], "tmid": [-0.005, 0.005]}

    rprs_samples = np.concatenate([
        np.linspace(0.106, 0.119, 40),
        np.linspace(0.120, 0.134, 15),
        np.linspace(0.145, 0.149, 5),
    ])
    tmid_samples = np.linspace(-2e-4, 2e-4, rprs_samples.size)
    points = np.column_stack([rprs_samples, tmid_samples])
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-6.0, -3.0, rprs_samples.size),
        },
        "samples": points.copy(),
    }

    diagnostics = fit.get_parameter_posterior_recenter_diagnostics("rprs")

    assert diagnostics["clipped"] is False
    assert diagnostics["edge"] is None
    assert diagnostics["upper_edge_peak_fraction"] < 0.20
    assert diagnostics["bounds"] == pytest.approx([0.0, 0.15])
    assert "not treated as truncated" in diagnostics["reason"]


def test_rprs_posterior_recenter_diagnostics_ignores_lower_edge_below_twenty_percent(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = True
    fit.prior = make_prior()
    fit.bounds = {"rprs": [0.0, 0.15], "tmid": [-0.005, 0.005]}
    fit.sampled_keys = ["rprs", "tmid"]
    fit.sample_bounds = {"rprs": [0.0, 0.15], "tmid": [-0.005, 0.005]}

    rprs_samples = np.concatenate([
        np.linspace(0.001, 0.005, 5),
        np.linspace(0.016, 0.029, 15),
        np.linspace(0.031, 0.044, 40),
    ])
    tmid_samples = np.linspace(-2e-4, 2e-4, rprs_samples.size)
    points = np.column_stack([rprs_samples, tmid_samples])
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-6.0, -3.0, rprs_samples.size),
        },
        "samples": points.copy(),
    }

    diagnostics = fit.get_parameter_posterior_recenter_diagnostics("rprs")

    assert diagnostics["clipped"] is False
    assert diagnostics["edge"] is None
    assert diagnostics["lower_edge_peak_fraction"] < 0.20
    assert diagnostics["bounds"] == pytest.approx([0.0, 0.15])
    assert "not treated as truncated" in diagnostics["reason"]


def test_plot_triangle_uses_direct_fitted_impact_parameter_axis(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    captured = {}

    def fake_corner(*args, **kwargs):
        captured["points"] = args[0]
        captured["labels"] = kwargs["labels"]
        captured["range"] = kwargs["range"]
        captured["titles"] = kwargs["titles"]
        captured["truths"] = kwargs["truths"]
        captured["label_kwargs"] = kwargs["label_kwargs"]
        return "figure"

    monkeypatch.setattr(elca, "corner", fake_corner)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.125],
        "inc": [84.0, 90.0],
        "a0": [0.95, 1.05],
    }
    fit.prior = make_prior()
    fit.sampled_keys = ["rprs", "b", "a0"]
    fit.sample_bounds = {
        "rprs": [0.0, 0.125],
        "b": [0.0, 1.25434156],
        "a0": [0.95, 1.05],
    }
    fit.sample_parameters = {"rprs": 0.10, "b": 0.314, "a0": 0.94962}
    fit.sample_errors = {"rprs": 0.01, "b": 0.05, "a0": 0.00394}
    fit.parameters = {"rprs": 0.10, "inc": 88.5, "a0": 0.94962}
    fit.errors = {"rprs": 0.01, "inc": 0.75, "a0": 0.00394}

    points = np.array(
        [
            [0.099, 0.300, 0.9501],
            [0.101, 0.330, 0.9502],
            [0.102, 0.290, 0.9515],
            [0.098, 0.360, 0.9520],
            [0.100, 0.314, 0.9508],
        ]
    )
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.array([-5.0, -4.0, -4.5, -5.5, -4.2]),
        },
        "samples": points.copy(),
    }

    fig = fit.plot_triangle()

    assert fig == "figure"
    assert captured["labels"][1] == r"Impact parameter $b$"
    assert captured["range"][1] == pytest.approx([0.0, 1.25434156])
    assert captured["points"].shape == (5, 3)
    np.testing.assert_allclose(captured["points"][:, 1], points[:, 1])
    assert captured["truths"][1] == pytest.approx(fit.sample_parameters["b"])
    assert captured["titles"][1].startswith("b=")
    assert "\ni=" in captured["titles"][1]
    assert captured["label_kwargs"]["labelpad"] == 10


def test_triangle_contour_levels_drop_duplicate_chi2_percentiles(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    chi2 = np.full(12, 42.0)
    mask = np.ones(chi2.size, dtype=bool)

    levels = fit._triangle_contour_levels(chi2, mask, mask, mask)

    assert levels == [pytest.approx(42.0)]


def test_triangle_plot_sigma_window_ranges_clip_to_solved_point_uncertainties(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    payload = {
        "ranges": [[0.0, 1.0], [0.0, 1.2], [0.95, 1.05]],
        "mask_centers": [0.20, 0.80, 1.0],
        "mask_errors": [0.02, 0.05, 0.001],
        "display_points": np.array(
            [
                [0.18, 0.75, 0.999],
                [0.20, 0.80, 1.000],
                [0.22, 0.85, 1.001],
            ]
        ),
    }

    ranges = fit._triangle_plot_sigma_window_ranges(payload, sigma=5.0)

    assert ranges[0] == pytest.approx([0.10, 0.30])
    assert ranges[1] == pytest.approx([0.55, 1.05])
    assert ranges[2] == pytest.approx([0.995, 1.005])


def test_plot_triangle_accepts_zoom_sigma(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    captured = {}

    def fake_corner(*args, **kwargs):
        captured["range"] = kwargs["range"]
        return "figure"

    monkeypatch.setattr(elca, "corner", fake_corner)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 1.0],
        "inc": [84.0, 90.0],
    }
    fit.sample_bounds = {
        "rprs": [0.0, 1.0],
        "b": [0.0, 1.2],
    }
    fit.sampled_keys = ["rprs", "b"]
    fit.prior = make_prior()
    fit.parameters = {"rprs": 0.20, "inc": 86.0}
    fit.errors = {"rprs": 0.02, "inc": 0.5}
    fit.sample_parameters = {"rprs": 0.20, "b": 0.80}
    fit.sample_errors = {"rprs": 0.02, "b": 0.05}
    points = np.column_stack([
        np.linspace(0.18, 0.22, 40),
        np.linspace(0.75, 0.85, 40),
    ])
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-4.0, -1.0, points.shape[0]),
        },
        "samples": points.copy(),
    }

    fig = fit.plot_triangle(zoom_sigma=5.0)

    assert fig == "figure"
    assert captured["range"][0][0] > 0.0
    assert captured["range"][0][1] < 1.0
    assert captured["range"][1][0] > 0.0
    assert captured["range"][1][1] < 1.2


def test_triangle_payload_recenter_uses_visible_zoom_peak(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    ars_values = np.concatenate([
        np.linspace(5.22, 5.30, 60),
        np.linspace(6.45, 6.55, 20),
        np.linspace(8.0, 9.0, 20),
    ])
    payload = {
        "sampled_keys": ["ars"],
        "display_points": ars_values[:, None],
        "display_weights": None,
        "ranges": [[5.0, 7.0]],
        "titles": ["6.0 +/- 1.0"],
        "truths": [6.0],
        "mask_centers": [6.0],
        "mask_errors": [1.0],
        "display_spec": None,
        "geometry_summary": {},
    }

    updated = fit._recenter_triangle_plot_payload_for_visible_ranges(payload)

    assert updated["truths"][0] == pytest.approx(5.3)
    assert updated["mask_centers"][0] == pytest.approx(5.3)
    assert updated["titles"][0].startswith("5.3 +/-")


def test_triangle_payload_expands_degenerate_error_ranges_to_sample_cloud(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.2],
        "tmid": [-0.01, 0.01],
    }
    fit.sample_bounds = dict(fit.bounds)
    fit.sampled_keys = ["rprs", "tmid"]
    fit.prior = make_prior()
    fit.parameters = {"rprs": 0.100, "tmid": 0.0}
    fit.errors = {"rprs": 1e-15, "tmid": 1e-15}
    fit.sample_parameters = dict(fit.parameters)
    fit.sample_errors = dict(fit.errors)
    points = np.column_stack(
        [
            np.linspace(0.050, 0.150, 50),
            np.linspace(-0.004, 0.004, 50),
        ]
    )
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-4.0, -1.0, points.shape[0]),
        },
        "samples": np.repeat(np.array([[0.100, 0.0]]), points.shape[0], axis=0),
    }

    payload = fit._get_triangle_plot_payload()

    assert payload["ranges"][0][0] <= 0.052
    assert payload["ranges"][0][1] >= 0.148
    assert payload["ranges"][1][0] <= -0.0038
    assert payload["ranges"][1][1] >= 0.0038


def test_triangle_payload_titles_follow_weighted_posterior_display_estimate(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.34],
        "a0": [0.95, 1.05],
    }
    fit.sample_bounds = dict(fit.bounds)
    fit.sampled_keys = ["rprs", "a0"]
    fit.prior = make_prior()
    fit.parameters = {"rprs": 0.33796, "a0": 1.0}
    fit.errors = {"rprs": 0.09150, "a0": 0.001}
    fit.sample_parameters = dict(fit.parameters)
    fit.sample_errors = dict(fit.errors)
    rprs_samples = np.concatenate([
        np.linspace(0.108, 0.122, 20),
        np.linspace(0.318, 0.338, 80),
    ])
    weights = np.concatenate([
        np.ones(20, dtype=float),
        np.full(80, 0.01, dtype=float),
    ])
    points = np.column_stack([rprs_samples, np.linspace(0.998, 1.002, rprs_samples.size)])
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-4.0, -1.0, points.shape[0]),
            "weights": weights,
        },
        "samples": points.copy(),
    }

    payload = fit._get_triangle_plot_payload()

    assert payload["titles"][0] == "0.1190 +/- 0.0089"
    assert payload["truths"][0] == pytest.approx(0.1186, abs=5e-4)
    np.testing.assert_allclose(payload["display_weights"], weights)


def test_plot_triangle_passes_ultranest_weights_to_visible_histograms(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    captured = {}

    def fake_corner(*args, **kwargs):
        captured["weights"] = kwargs["weights"]
        captured["truths"] = kwargs["truths"]
        captured["data_kwargs"] = kwargs["data_kwargs"]
        return "figure"

    monkeypatch.setattr(elca, "corner", fake_corner)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.2],
        "a0": [0.95, 1.05],
    }
    fit.sample_bounds = dict(fit.bounds)
    fit.sampled_keys = ["rprs", "a0"]
    fit.prior = make_prior()
    fit.parameters = {"rprs": 0.1, "a0": 1.0}
    fit.errors = {"rprs": 0.01, "a0": 0.001}
    fit.sample_parameters = dict(fit.parameters)
    fit.sample_errors = dict(fit.errors)
    points = np.column_stack([
        np.linspace(0.09, 0.11, 12),
        np.linspace(0.998, 1.002, 12),
    ])
    weights = np.linspace(1.0, 2.0, points.shape[0])
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-4.0, -1.0, points.shape[0]),
            "weights": weights,
        },
        "samples": points.copy(),
    }

    fig = fit.plot_triangle()

    assert fig == "figure"
    np.testing.assert_allclose(captured["weights"], weights)
    np.testing.assert_allclose(captured["truths"], [0.1018961, 1.00037922], rtol=1e-6)
    assert captured["data_kwargs"]["s"] == pytest.approx(1.6)
    assert captured["data_kwargs"]["alpha"] == pytest.approx(0.38)


def test_triangle_payload_expands_sparse_visible_ranges_to_sample_cloud(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.2],
        "a0": [0.95, 1.05],
    }
    fit.sample_bounds = dict(fit.bounds)
    fit.sampled_keys = ["rprs", "a0"]
    fit.prior = make_prior()
    fit.parameters = {"rprs": 0.100, "a0": 1.0}
    fit.errors = {"rprs": 0.01, "a0": 1e-4}
    fit.sample_parameters = dict(fit.parameters)
    fit.sample_errors = dict(fit.errors)
    points = np.column_stack(
        [
            np.linspace(0.090, 0.110, 100),
            np.linspace(0.980, 1.020, 100),
        ]
    )
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-4.0, -1.0, points.shape[0]),
        },
        "samples": points.copy(),
    }

    payload = fit._get_triangle_plot_payload()

    assert payload["ranges"][1][0] <= 0.981
    assert payload["ranges"][1][1] >= 1.019


def test_triangle_payload_uses_tested_rprs_range_when_posterior_is_narrow(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.2],
        "a0": [0.95, 1.05],
    }
    fit.sample_bounds = dict(fit.bounds)
    fit.sampled_keys = ["rprs", "a0"]
    fit.prior = make_prior()
    fit.parameters = {"rprs": 0.100, "a0": 1.0}
    fit.errors = {"rprs": 0.001, "a0": 0.001}
    fit.sample_parameters = dict(fit.parameters)
    fit.sample_errors = dict(fit.errors)
    points = np.column_stack(
        [
            np.linspace(0.090, 0.110, 120),
            np.linspace(0.998, 1.002, 120),
        ]
    )
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-4.0, -1.0, points.shape[0]),
        },
        "samples": points.copy(),
    }

    payload = fit._get_triangle_plot_payload()
    rprs_range = payload["ranges"][0]
    plot_bins = int(max(1, np.sqrt(points.shape[0])))
    lower_fraction, upper_fraction, _ = fit._histogram_edge_peak_fractions(
        points[:, 0],
        rprs_range,
        plot_bins,
    )

    assert rprs_range == pytest.approx([0.0, 0.2])
    assert lower_fraction < elca.TRIANGLE_PLOT_EDGE_PEAK_FRACTION_MAX
    assert upper_fraction < elca.TRIANGLE_PLOT_EDGE_PEAK_FRACTION_MAX


def test_triangle_payload_expands_rprs_lower_edge_past_narrow_recorded_sample_bounds(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.36],
        "a0": [0.95, 1.05],
    }
    fit.sample_bounds = {
        "rprs": [0.104, 0.184],
        "a0": [0.95, 1.05],
    }
    fit.sampled_keys = ["rprs", "a0"]
    fit.prior = make_prior()
    fit.parameters = {"rprs": 0.144, "a0": 1.0}
    fit.errors = {"rprs": 0.008, "a0": 0.001}
    fit.sample_parameters = dict(fit.parameters)
    fit.sample_errors = dict(fit.errors)
    rprs_samples = np.concatenate([
        np.linspace(0.104, 0.120, 80),
        np.linspace(0.120, 0.180, 20),
    ])
    points = np.column_stack([
        rprs_samples,
        np.linspace(0.998, 1.002, rprs_samples.size),
    ])
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-4.0, -1.0, points.shape[0]),
        },
        "samples": points.copy(),
    }

    payload = fit._get_triangle_plot_payload()

    assert payload["ranges"][0][0] < 0.08
    assert payload["truths"][0] < 0.13
    assert not payload["titles"][0].startswith("0.144")


def test_triangle_payload_keeps_direct_impact_parameter_full_sample_range(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.2],
        "inc": [84.0, 90.0],
        "a0": [0.95, 1.05],
    }
    fit.sampled_keys = ["rprs", "b", "a0"]
    fit.sample_bounds = {
        "rprs": [0.0, 0.2],
        "b": [0.0, 1.2],
        "a0": [0.95, 1.05],
    }
    fit.prior = make_prior()
    fit.sample_parameters = {"rprs": 0.10, "b": 0.30, "a0": 1.0}
    fit.sample_errors = {"rprs": 0.01, "b": 0.01, "a0": 0.001}
    fit.parameters = {"rprs": 0.10, "inc": 88.6, "a0": 1.0}
    fit.errors = {"rprs": 0.01, "inc": 0.75, "a0": 0.001}
    points = np.column_stack(
        [
            np.linspace(0.090, 0.110, 100),
            np.linspace(0.10, 0.80, 100),
            np.linspace(0.998, 1.002, 100),
        ]
    )
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-4.0, -1.0, points.shape[0]),
        },
        "samples": points.copy(),
    }

    payload = fit._get_triangle_plot_payload()

    assert payload["labels"][1] == r"Impact parameter $b$"
    assert payload["ranges"][1] == pytest.approx([0.0, 1.2])
    assert payload["display_points"].shape == points.shape
    np.testing.assert_allclose(payload["display_points"][:, 1], points[:, 1])
    assert payload["truths"][1] == pytest.approx(fit.sample_parameters["b"])
    assert payload["display_spec"]["mirror"] is False


def test_triangle_payload_uses_full_b_range_for_direct_impact_parameter(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.2],
        "inc": [84.0, 90.0],
        "a0": [0.95, 1.05],
    }
    fit.sampled_keys = ["rprs", "b", "a0"]
    fit.sample_bounds = {
        "rprs": [0.0, 0.2],
        "b": [0.0, 1.2],
        "a0": [0.95, 1.05],
    }
    fit.prior = make_prior()
    fit.sample_parameters = {"rprs": 0.10, "b": 0.856, "a0": 1.0}
    fit.sample_errors = {"rprs": 0.01, "b": 0.002, "a0": 0.001}
    fit.parameters = {"rprs": 0.10, "inc": 85.0, "a0": 1.0}
    fit.errors = {"rprs": 0.01, "inc": 2.7, "a0": 0.001}
    b_samples = np.linspace(0.846, 0.866, 100)
    points = np.column_stack([
        np.linspace(0.090, 0.110, b_samples.size),
        b_samples,
        np.linspace(0.998, 1.002, b_samples.size),
    ])
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.linspace(-4.0, -1.0, points.shape[0]),
        },
        "samples": points.copy(),
    }

    payload = fit._get_triangle_plot_payload()
    reference_values = [reference["value"] for reference in payload["display_spec"]["reference_lines"]]

    assert payload["labels"][1] == r"Impact parameter $b$"
    assert payload["ranges"][1] == pytest.approx([0.0, 1.2])
    assert payload["truths"][1] == pytest.approx(0.856)
    assert payload["display_spec"]["mirror"] is False
    assert reference_values == pytest.approx([1.0, 1.10])


def test_triangle_payload_tracks_left_and_right_geometry_branches_for_inclination(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.125],
        "inc": [84.0, 90.0],
        "a0": [0.95, 1.05],
    }
    fit.prior = make_prior()
    fit.parameters = {"rprs": 0.10, "inc": 88.42, "a0": 0.94962}
    fit.errors = {"rprs": 0.01, "inc": 0.75, "a0": 0.00394}
    points = np.array(
        [
            [0.099, 88.30, 0.9501],
            [0.101, 88.55, 0.9502],
            [0.102, 88.10, 0.9515],
            [0.098, 88.70, 0.9520],
            [0.100, 88.40, 0.9508],
        ]
    )
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.array([-5.0, -4.0, -4.5, -5.5, -4.2]),
        },
        "samples": points.copy(),
    }

    payload = fit._get_triangle_plot_payload()

    np.testing.assert_allclose(
        payload["geometry_overlay"]["left_mirrored"],
        np.array([-0.12, -0.32, -0.02, 0.12, 0.32, 0.02]),
        atol=1e-12,
    )
    np.testing.assert_allclose(
        payload["geometry_overlay"]["right_mirrored"],
        np.array([0.13, 0.28, -0.13, -0.28]),
        atol=1e-12,
    )


def test_triangle_payload_skips_mirrored_overlay_for_direct_impact_parameter(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.125],
        "inc": [84.0, 90.0],
        "a0": [0.95, 1.05],
    }
    fit.prior = make_prior()
    fit.sampled_keys = ["rprs", "b", "a0"]
    fit.sample_bounds = {
        "rprs": [0.0, 0.125],
        "b": [0.0, 1.25434156],
        "a0": [0.95, 1.05],
    }
    fit.sample_parameters = {"rprs": 0.10, "b": 0.314, "a0": 0.94962}
    fit.sample_errors = {"rprs": 0.01, "b": 0.05, "a0": 0.00394}
    fit.parameters = {"rprs": 0.10, "inc": 88.5, "a0": 0.94962}
    fit.errors = {"rprs": 0.01, "inc": 0.75, "a0": 0.00394}
    points = np.array(
        [
            [0.099, 0.300, 0.9501],
            [0.101, 0.330, 0.9502],
            [0.102, 0.290, 0.9515],
            [0.098, 0.360, 0.9520],
            [0.100, 0.314, 0.9508],
        ]
    )
    fit.results = {
        "weighted_samples": {
            "points": points,
            "logl": np.array([-5.0, -4.0, -4.5, -5.5, -4.2]),
        },
        "samples": points.copy(),
    }

    payload = fit._get_triangle_plot_payload()

    assert payload["geometry_overlay"] is None
    assert payload["display_spec"]["mirror"] is False
    assert payload["ranges"][1] == pytest.approx([0.0, 1.25434156])
    np.testing.assert_allclose(payload["display_points"][:, 1], points[:, 1])


def test_triangle_geometry_curves_fall_back_to_surviving_branch(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    curves = fit._build_triangle_plot_geometry_curves(
        {
            "left_count": 0,
            "right_count": 4,
            "left_mirrored": np.array([], dtype=float),
            "right_mirrored": np.array([0.05, 0.10, 0.15, -0.05, -0.10, -0.15], dtype=float),
        },
        [-0.3, 0.3],
        31,
    )

    np.testing.assert_allclose(curves["main_curve"], fit._smooth_triangle_plot_counts(curves["right_curve"]))
    assert np.allclose(curves["left_curve"], 0.0)


def test_triangle_geometry_overlay_reuses_shared_title_and_label_kwargs(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    fig, axes = plt.subplots(2, 2)
    payload = {
        "sampled_keys": ["rprs", "b"],
        "display_points": np.zeros((10, 2)),
        "display_spec": {
            "index": 1,
            "center": 0.92,
            "reference_lines": [
                {"value": 1.0, "linestyle": ":", "color": "#707070"},
                {"value": 1.10, "linestyle": "-.", "color": "#a35d00"},
            ],
        },
        "geometry_overlay": {
            "index": 1,
            "left_count": 2,
            "right_count": 2,
            "left_mirrored": np.array([-0.1, 0.1]),
            "right_mirrored": np.array([-0.2, 0.2]),
        },
        "ranges": [[0.0, 0.1], [-0.3, 0.3]],
        "titles": ["rprs", "b=0.32 +/- 0.18\ni=88.5 +/- 1.35 deg"],
        "labels": ["rprs", r"$\Delta b$"],
    }

    fit._overlay_triangle_plot_geometry_histograms(
        fig,
        payload,
        title_kwargs={"loc": "left", "pad": 4, "fontsize": 12},
        label_kwargs={"labelpad": 10},
    )

    ax = axes[1, 1]
    assert ax.title.get_fontsize() == pytest.approx(12.0)
    assert ax.xaxis.label.get_text() == r"$\Delta b$"
    assert ax.xaxis.labelpad == pytest.approx(10.0)
    reference_offsets = []
    for line in ax.lines:
        xdata = np.asarray(line.get_xdata(), dtype=float)
        if xdata.size == 2 and np.allclose(xdata, xdata[0]) and line.get_linestyle() in (":", "-."):
            reference_offsets.append(float(xdata[0]))
    assert sorted(reference_offsets) == pytest.approx([-0.18, -0.08, 0.08, 0.18])
    plt.close(fig)


def test_nested_fit_keeps_posterior_summary_on_real_ultranest_dead_points(monkeypatch, tmp_path):
    # Real UltraNest output from the WBoM final fit of an 83-point CoRoT-2 b light curve
    # (MicroObservatory, 2026-08-08; issue #1401), sampled keys ars, b, rprs, tmid, with
    # tmid stored relative to 2461261.8 BJD_TDB. The posterior stdev for tmid is 0.0056 d;
    # the previous guard (interior std, factor 3.0) replaced it with 0.0018 d while rprs,
    # at ratio 2.96, escaped. Neither may be replaced: the posterior is healthy.
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fixture = np.load(os.path.join(os.path.dirname(__file__), "data", "ultranest_dead_points_corot2_2026-08-08.npz"))
    keys = [str(key) for key in fixture["keys"]]
    points = np.asarray(fixture["points"], dtype=float)
    stdev = np.asarray(fixture["stdev"], dtype=float)
    ml = np.asarray(fixture["ml"], dtype=float)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)
    fit.prior = make_prior()
    fit.bounds = {key: [float(points[:, index].min()), float(points[:, index].max())] for index, key in enumerate(keys)}
    fit.mode = "ns"
    fit.use_impactparameter_rather_than_inclination_to_fit = True
    fit.fixed_parameter_errors = {}
    fit.results = {
        "maximum_likelihood": {"point": ml},
        "posterior": {"stdev": stdev, "errlo": ml - stdev, "errup": ml + stdev},
        "weighted_samples": {"points": points, "logl": np.asarray(fixture["logl"], dtype=float)},
        "samples": points.copy(),
    }

    fit._finalize_ultranest_fit_results(
        keys,
        keys,
        lambda point: {key: float(value) for key, value in zip(keys, point)},
    )
    tmid_index = keys.index("tmid")
    assert fit.ultranest_error_fallbacks == {}
    assert fit.errors["tmid"] == pytest.approx(stdev[tmid_index])
    assert fit.errors["rprs"] == pytest.approx(stdev[keys.index("rprs")])

    local = fit._loglike_neighborhood_uncertainty(tmid_index, ml[tmid_index])
    assert local["interior_std"] < 0.4 * stdev[tmid_index]  # what used to be quoted
    assert local["error"] > 0.5 * stdev[tmid_index]
