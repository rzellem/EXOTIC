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
