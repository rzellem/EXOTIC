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


def test_baseline_restriction_keeps_one_hour_before_ingress_and_after_egress():
    exotic_main = _exotic_main_module()
    planet = {
        "pPer": 3.0,
        "midT": 2450000.0,
        "rprs": 0.1,
        "aRs": 10.0,
        "inc": 88.0,
        "ecc": 0.0,
        "omega": 90.0,
    }
    duration = exotic_main.estimate_transit_duration_from_prior_geometry({
        "per": planet["pPer"],
        "rprs": planet["rprs"],
        "ars": planet["aRs"],
        "inc": planet["inc"],
        "ecc": planet["ecc"],
        "omega": planet["omega"],
    })
    edge = duration / 2.0 + 1.0 / 24.0
    times = np.array([
        planet["midT"] - edge - 1e-5,
        planet["midT"] - edge + 1e-5,
        planet["midT"],
        planet["midT"] + edge - 1e-5,
        planet["midT"] + edge + 1e-5,
    ])

    keep, summary = exotic_main.build_baseline_restriction_mask(times, planet)

    assert keep.tolist() == [False, True, True, True, False]
    assert summary["applied"] is True
    assert summary["excluded_point_count"] == 2
    assert summary["kept_point_count"] == 3


def test_baseline_restriction_disabled_keeps_all_points():
    exotic_main = _exotic_main_module()
    planet = {"pPer": 3.0, "midT": 2450000.0, "rprs": 0.1, "aRs": 10.0, "inc": 88.0}
    times = np.array([planet["midT"] - 10.0, planet["midT"], planet["midT"] + 10.0])

    keep, summary = exotic_main.build_baseline_restriction_mask(times, planet, enabled=False)

    assert np.all(keep)
    assert summary["applied"] is False
    assert summary["excluded_point_count"] == 0


def test_baseline_restriction_falls_back_when_window_leaves_too_few_points():
    exotic_main = _exotic_main_module()
    planet = {"pPer": 3.0, "midT": 2450000.0, "rprs": 0.1, "aRs": 10.0, "inc": 88.0}
    times = planet["midT"] + np.arange(6, dtype=float) * 0.01 + 0.5

    keep, summary = exotic_main.build_baseline_restriction_mask(times, planet)

    assert np.all(keep)
    assert summary["applied"] is False
    assert summary["excluded_point_count"] == 0
    assert "below EXOTIC's minimum" in summary["note"]


def test_comparison_preflight_preserves_baseline_level_for_excluded_plot_points():
    exotic_main = _exotic_main_module()
    times = np.linspace(2450000.0, 2450000.1, 20)
    prepared = exotic_main.prepare_comparison_candidate_full_reduction_series(
        times,
        np.full(times.shape, 200.0),
        np.full(times.shape, 100.0),
        np.ones(times.shape),
    )

    assert prepared["applied"] is True
    assert prepared["approximate_baseline_level"] == pytest.approx(2.0)
    assert np.nanmedian(prepared["flux"]) == pytest.approx(1.0)


def test_unrestricted_prefit_clip_masks_outside_window_points_before_plotting(monkeypatch):
    exotic_main = _exotic_main_module()
    times = np.arange(4.0)
    target_flux = np.full(times.shape, 2.0)
    comp_flux = np.ones(times.shape)
    airmass = np.ones(times.shape)

    def fake_prepare(*args, **kwargs):
        return {
            "applied": True,
            "source_indices": np.array([0, 2, 3], dtype=int),
            "filter_diagnostics": [{"name": "Initial sigma clip"}],
        }

    monkeypatch.setattr(exotic_main, "prepare_lightcurve_fit_input_series", fake_prepare)
    summary = exotic_main.build_unrestricted_candidate_prefit_clip_summary(
        times,
        target_flux,
        comp_flux,
        airmass,
        np.ones(times.shape, dtype=bool),
    )

    assert summary["applied"] is True
    assert summary["keep_mask"].tolist() == [True, False, True, True]
    assert summary["rejected_mask"].tolist() == [False, True, False, False]


def test_restricted_baseline_payload_separates_prefit_rejects_from_blue_points():
    exotic_main = _exotic_main_module()
    times = np.arange(4.0)
    payload = exotic_main.build_restricted_baseline_plot_payload(
        times,
        np.full(times.shape, 2.0),
        np.ones(times.shape),
        np.ones(times.shape),
        np.array([False, True, True, False]),
        prefit_keep_mask=np.array([True, True, False, True]),
        normalization_level=2.0,
    )

    assert payload["point_count"] == 1
    assert payload["rejected_point_count"] == 1
    assert payload["flux"].tolist() == [1.0]
    assert payload["rejected_flux"].tolist() == [1.0]


def test_restricted_baseline_payload_can_include_all_prefit_rejects():
    exotic_main = _exotic_main_module()
    times = np.arange(4.0)
    payload = exotic_main.build_restricted_baseline_plot_payload(
        times,
        np.full(times.shape, 2.0),
        np.ones(times.shape),
        np.ones(times.shape),
        np.array([False, True, True, False]),
        prefit_keep_mask=np.array([True, True, False, True]),
        prefit_rejected_mask=np.array([True, False, True, False]),
        normalization_level=2.0,
    )

    assert payload["point_count"] == 1
    assert payload["rejected_point_count"] == 2
    assert payload["rejected_times"].tolist() == [0.0, 2.0]


def test_restricted_baseline_payload_keeps_prefit_rejects_when_window_has_no_blue_points():
    exotic_main = _exotic_main_module()
    times = np.arange(3.0)
    payload = exotic_main.build_restricted_baseline_plot_payload(
        times,
        np.full(times.shape, 2.0),
        np.ones(times.shape),
        np.ones(times.shape),
        np.zeros(times.shape, dtype=bool),
        prefit_keep_mask=np.array([False, True, True]),
        prefit_rejected_mask=np.array([True, False, False]),
        normalization_level=2.0,
    )

    assert payload["point_count"] == 0
    assert payload["rejected_point_count"] == 1
    assert payload["rejected_times"].tolist() == [0.0]
