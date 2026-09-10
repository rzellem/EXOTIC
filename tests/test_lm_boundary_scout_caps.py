"""Issue #1406: the LM boundary scout must never expand past a physical envelope."""
import numpy as np
import pytest

import exotic.exotic as ex


def test_envelope_uses_prior_relative_caps():
    env = ex.lm_boundary_scout_physical_envelope({'rprs': 0.1667, 'ars': 6.70})
    assert env['ars'][0] == pytest.approx(max(1.5, 1.1667, 0.35 * 6.70))
    assert env['ars'][1] == pytest.approx(3.0 * 6.70)
    assert env['rprs'][1] == pytest.approx(min(ex.RPRS_SEARCH_BOUND_MAX, max(0.30, 2 * 0.1667)))


def test_envelope_without_prior_keeps_absolute_floors():
    env = ex.lm_boundary_scout_physical_envelope({})
    assert env['ars'][0] == pytest.approx(ex.ULTRANEST_LM_BOUNDARY_SCOUT_ARS_PHYSICAL_FLOOR)
    assert env['ars'][1] == pytest.approx(ex.ARS_SEARCH_BOUND_FALLBACK_MAX)
    assert env['rprs'] == (ex.RPRS_SEARCH_BOUND_MIN, ex.RPRS_SEARCH_BOUND_MAX)


def test_envelope_never_below_one_plus_rprs():
    env = ex.lm_boundary_scout_physical_envelope({'rprs': 0.6, 'ars': 2.0})
    assert env['ars'][0] >= 1.6


def test_generic_expansion_clamps_to_envelope_floor():
    env = ex.lm_boundary_scout_physical_envelope({'rprs': 0.1667, 'ars': 6.70})
    floor, ceiling = env['ars']
    # An LM solution sitting on the lower edge asks for room below; the floor must hold.
    bounds, adjusted = ex.lm_boundary_scout_expanded_bounds(
        4.70, 0.5, [4.69, 8.71], minimum_bound=floor, maximum_bound=ceiling,
    )
    assert bounds[0] >= floor
    assert bounds[1] <= ceiling
    # Repeated expansion cannot creep below the floor either.
    for _ in range(5):
        bounds, _ = ex.lm_boundary_scout_expanded_bounds(
            bounds[0] + 1e-3, 0.5, bounds, minimum_bound=floor, maximum_bound=ceiling,
        )
    assert bounds[0] >= floor > ex.ARS_SEARCH_BOUND_MIN


class _FakeFit:
    def __init__(self, parameters):
        self.parameters = parameters


def test_geometry_envelope_does_not_fall_back_to_unphysical_floor(monkeypatch):
    # Force the geometric lower limit to be uncomputable (rprs missing) while the
    # central-duration limit stays finite: the old code fell back to 1e-6 here.
    monkeypatch.setattr(ex, 'estimate_transit_duration_from_fit', lambda fit: 0.095)
    fit = _FakeFit({'per': 1.743, 'ecc': 0.0, 'omega': 0.0})  # no 'rprs'
    prior = {'per': 1.743, 'ars': 6.70}  # no 'rprs' in the prior either
    bounds = {'ars': [0.67, 12.73], 'rprs': [0.09, 0.305]}
    safe, adjustments = ex.lm_boundary_scout_geometry_safe_bounds(prior, fit, bounds, {'ars'})
    assert safe['ars'][0] >= ex.ULTRANEST_LM_BOUNDARY_SCOUT_ARS_PHYSICAL_FLOOR
    assert safe['ars'][0] > ex.ARS_SEARCH_BOUND_MIN * 1e3


def test_geometry_envelope_caps_ars_ceiling_to_prior_multiple(monkeypatch):
    monkeypatch.setattr(ex, 'estimate_transit_duration_from_fit', lambda fit: 0.095)
    fit = _FakeFit({'per': 1.743, 'rprs': 0.1667, 'ecc': 0.0, 'omega': 0.0})
    prior = {'per': 1.743, 'rprs': 0.1667, 'ars': 6.70}
    bounds = {'ars': [4.69, 8.71], 'rprs': [0.09, 0.305]}
    safe, _ = ex.lm_boundary_scout_geometry_safe_bounds(prior, fit, bounds, {'ars'})
    assert safe['ars'][1] <= 3.0 * 6.70 + 1e-9
    assert safe['ars'][0] >= 0.35 * 6.70 - 1e-9
