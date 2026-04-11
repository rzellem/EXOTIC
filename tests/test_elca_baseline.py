import importlib
import sys
import types

import numpy as np
import pytest


def load_elca_with_stubs(monkeypatch, tmp_path):
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

    for name in list(sys.modules):
        if name == "exotic.api.elca" or name.startswith("pylightcurve"):
            sys.modules.pop(name, None)

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


def test_plot_triangle_clips_ranges_to_parameter_bounds(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    captured = {}

    def fake_corner(*args, **kwargs):
        captured["points"] = args[0]
        captured["labels"] = kwargs["labels"]
        captured["range"] = kwargs["range"]
        return "figure"

    monkeypatch.setattr(elca, "corner", fake_corner)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.125],
        "inc": [84.0, 90.0],
        "a0": [0.95, 1.05],
    }
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
    assert captured["labels"][1] == "Distance from fitted Inc. [deg] (mirrored)"
    assert captured["range"][0][0] == pytest.approx(0.05)
    assert captured["range"][0][1] == pytest.approx(0.125)
    expected_inc_distance_limit = np.max(np.abs(np.array([84.67, 90.0]) - fit.parameters["inc"]))
    assert captured["range"][1][0] == pytest.approx(-expected_inc_distance_limit)
    assert captured["range"][1][1] == pytest.approx(expected_inc_distance_limit)
    assert captured["range"][2][0] == pytest.approx(0.95)
    assert captured["range"][2][1] == pytest.approx(0.96932)
    assert captured["points"].shape == (10, 3)
    expected_inc_distance = np.abs(points[:, 1] - fit.parameters["inc"])
    np.testing.assert_allclose(captured["points"][:5, 1], expected_inc_distance)
    np.testing.assert_allclose(captured["points"][5:, 1], -expected_inc_distance)


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
    expected_inc = 87.0 + 0.4 * (90.0 - 87.0)

    assert fit._get_sampled_keys() == ["rprs", "b", "tmid"]
    assert physical["inc"] == pytest.approx(expected_inc)
    assert physical["b"] == pytest.approx(sample_point[1])


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


def test_plot_triangle_uses_mirrored_distance_from_fitted_impact_parameter_axis(monkeypatch, tmp_path):
    elca = load_elca_with_stubs(monkeypatch, tmp_path)
    fit = elca.lc_fitter.__new__(elca.lc_fitter)

    captured = {}

    def fake_corner(*args, **kwargs):
        captured["points"] = args[0]
        captured["labels"] = kwargs["labels"]
        captured["range"] = kwargs["range"]
        return "figure"

    monkeypatch.setattr(elca, "corner", fake_corner)

    fit.ns_type = "ultranest"
    fit.bounds = {
        "rprs": [0.0, 0.125],
        "inc": [84.0, 90.0],
        "a0": [0.95, 1.05],
    }
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
    assert captured["labels"][1] == "Distance from fitted b (mirrored)"
    assert captured["range"][1][0] == pytest.approx(-0.25)
    assert captured["range"][1][1] == pytest.approx(0.25)
    assert captured["points"].shape == (10, 3)
    expected_b_distance = np.abs(points[:, 1] - fit.sample_parameters["b"])
    np.testing.assert_allclose(captured["points"][:5, 1], expected_b_distance)
    np.testing.assert_allclose(captured["points"][5:, 1], -expected_b_distance)
