import numpy as np
import pytest

from exotic import exotic as ex
from exotic.api.elca import transit


def make_prior():
    return {
        "rprs": 0.15,
        "ars": 7.0,
        "per": 2.15,
        "inc": 88.9,
        "u0": 0.4,
        "u1": 0.2,
        "u2": 0.0,
        "u3": 0.0,
        "ecc": 0.0,
        "omega": 90.0,
        "tmid": 0.0,
        "a1": 1.0,
        "a2": 0.0,
    }


def make_lightcurve(true_noise, quoted_noise, seed=7):
    prior = make_prior()
    time = np.linspace(-0.15, 0.15, 120)
    airmass = np.ones_like(time)
    rng = np.random.default_rng(seed)
    flux = transit(time, prior) + rng.normal(0.0, true_noise, time.size)
    errors = np.full_like(time, quoted_noise)
    bounds = {"rprs": [0.1, 0.2], "tmid": [-0.02, 0.02], "a1": [0.9, 1.1]}
    return time, flux, errors, airmass, prior, bounds


def test_underestimated_uncertainties_are_inflated_to_unit_reduced_chi2():
    time, flux, errors, airmass, prior, bounds = make_lightcurve(true_noise=0.004, quoted_noise=0.002)
    inflated, summary = ex.inflate_uncertainties_to_unit_reduced_chi2(
        time, flux, errors, airmass, prior, bounds,
    )
    assert summary["applied"] is True
    assert summary["factor"] == pytest.approx(2.0, rel=0.15)
    assert np.allclose(inflated, errors * summary["factor"])
    assert summary["reduced_chi2"] == pytest.approx(summary["factor"] ** 2)


def test_overestimated_uncertainties_are_never_deflated():
    time, flux, errors, airmass, prior, bounds = make_lightcurve(true_noise=0.001, quoted_noise=0.003)
    kept, summary = ex.inflate_uncertainties_to_unit_reduced_chi2(
        time, flux, errors, airmass, prior, bounds,
    )
    assert summary["applied"] is False
    assert summary["factor"] == 1.0
    assert summary["reduced_chi2"] < 1.0
    assert np.array_equal(kept, errors)


def test_misaligned_arrays_return_input_unchanged():
    time, flux, errors, airmass, prior, bounds = make_lightcurve(true_noise=0.004, quoted_noise=0.002)
    kept, summary = ex.inflate_uncertainties_to_unit_reduced_chi2(
        time, flux[:-1], errors, airmass, prior, bounds,
    )
    assert summary["applied"] is False
    assert np.array_equal(kept, errors)


@pytest.mark.parametrize("value, expected", [("y", True), ("n", False), (None, True), (False, False)])
def test_inflation_setting_parses(value, expected):
    assert ex.should_inflate_final_fit_uncertainties(value) is expected
