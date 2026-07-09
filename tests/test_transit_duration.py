"""Tests for the analytic transit-duration formula (issue #1383).

The eccentric duration must follow Winn (2010), "Transits and Occultations",
arXiv:1001.2010: the impact parameter uses the eccentricity factor
(1 - e^2)/(1 + e sin w) (eq. 7), the arcsin argument is normalized by the
plain a/R* (eq. 14), and the whole expression is multiplied by the velocity
factor sqrt(1 - e^2)/(1 + e sin w) (eq. 16). Putting the eq.-7 factor inside
the arcsin instead inverts the eccentricity dependence: periastron-at-transit
(fastest planet, shortest transit) comes out longest, and the error grows as
(1 + e sin w)^2 / (1 - e^2)^{3/2}, reaching ~60x for Kepler-1704 b.
"""

import math

import numpy as np
import pytest

from exotic.transit_depth import transit_duration_days


def winn_2010_duration_days(period, ars, inc_deg, rprs, ecc, omega_deg):
    """Reference implementation: Winn (2010) eqs. 7, 14, 16."""
    inc = math.radians(inc_deg)
    omega = math.radians(omega_deg)
    denom = 1.0 + ecc * math.sin(omega)
    b = ars * math.cos(inc) * (1.0 - ecc ** 2) / denom
    chord_sq = (1.0 + rprs) ** 2 - b ** 2
    if chord_sq <= 0:
        return float("nan")
    argument = min(math.sqrt(chord_sq) / (ars * math.sin(inc)), 1.0)
    velocity_factor = math.sqrt(1.0 - ecc ** 2) / denom
    return (period / math.pi) * math.asin(argument) * velocity_factor


def duration_params(period, ars, inc, rprs, ecc, omega):
    return {
        "per": period,
        "ars": ars,
        "inc": inc,
        "rprs": rprs,
        "ecc": ecc,
        "omega": omega,
    }


def test_circular_matches_winn_exactly():
    params = duration_params(3.0, 10.0, 90.0, 0.1, 0.0, 90.0)
    expected = winn_2010_duration_days(3.0, 10.0, 90.0, 0.1, 0.0, 90.0)
    assert transit_duration_days(params) == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("ecc", [0.1, 0.3, 0.5, 0.7, 0.9])
@pytest.mark.parametrize("omega", [0.0, 45.0, 90.0, 135.0, 180.0, 270.0])
@pytest.mark.parametrize("inc", [90.0, 88.0, 85.0])
def test_eccentric_matches_winn(ecc, omega, inc):
    params = duration_params(3.0, 10.0, inc, 0.1, ecc, omega)
    expected = winn_2010_duration_days(3.0, 10.0, inc, 0.1, ecc, omega)
    result = transit_duration_days(params)
    if math.isnan(expected):
        assert math.isnan(result)
    else:
        assert result == pytest.approx(expected, rel=1e-9)


def test_periastron_transit_is_shorter_and_apastron_longer():
    circular = transit_duration_days(duration_params(3.0, 10.0, 90.0, 0.1, 0.0, 90.0))
    periastron = transit_duration_days(duration_params(3.0, 10.0, 90.0, 0.1, 0.5, 90.0))
    apastron = transit_duration_days(duration_params(3.0, 10.0, 90.0, 0.1, 0.5, 270.0))
    assert periastron < circular < apastron


@pytest.mark.parametrize(
    "name, period, ars, inc, rprs, ecc, omega, published_hours",
    [
        # NASA Exoplanet Archive `ps` default rows, pl_trandur in hours.
        ("Kepler-1704 b", 988.88112, 256.4, 89.00, 0.0644, 0.920, 82.40, 6.007),
        ("HD 17156 b", 21.2164294, 23.11, 86.51, 0.07412, 0.6772, 122.06, 3.1505),
        ("HD 80606 b", 111.436765, 94.452, 89.24, 0.1009, 0.93183, -58.887, 11.98),
    ],
)
def test_reproduces_published_durations(name, period, ars, inc, rprs, ecc, omega, published_hours):
    params = duration_params(period, ars, inc, rprs, ecc, omega)
    hours = transit_duration_days(params) * 24.0
    assert hours == pytest.approx(published_hours, rel=0.05), name


def _exotic_main_module():
    return pytest.importorskip(
        "exotic.exotic", reason="exotic.exotic imports the full pipeline dependency stack"
    )


def test_qc_contact_duration_matches_winn():
    exotic_main = _exotic_main_module()
    params = {"per": 3.0, "ars": 10.0, "inc": 89.0, "ecc": 0.5, "omega": 90.0}
    result = exotic_main.transit_qc_geometry_contact_duration(params, 1.1)
    inc = math.radians(89.0)
    omega = math.radians(90.0)
    denom = 1.0 + 0.5 * math.sin(omega)
    b = 10.0 * math.cos(inc) * (1.0 - 0.25) / denom
    argument = min(math.sqrt(1.1 ** 2 - b ** 2) / (10.0 * math.sin(inc)), 1.0)
    expected = (3.0 / math.pi) * math.asin(argument) * (math.sqrt(0.75) / denom)
    assert result == pytest.approx(expected, rel=1e-9)


def test_prior_geometry_duration_matches_winn():
    exotic_main = _exotic_main_module()
    prior = {"per": 3.0, "ars": 10.0, "inc": 88.0, "rprs": 0.1, "ecc": 0.4, "omega": 120.0}
    result = exotic_main.estimate_transit_duration_from_prior_geometry(prior)
    expected = winn_2010_duration_days(3.0, 10.0, 88.0, 0.1, 0.4, 120.0)
    assert result == pytest.approx(expected, rel=1e-9)
