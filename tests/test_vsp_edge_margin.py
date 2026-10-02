"""Regression test for the VSP comparison-star edge margin.

Motivating failure (Exoplanet Watch, 2026-08-17): on a 650x500
MicroObservatory night that drifted ~73 px in y, the VSP fallback
auto-added a comparison star at [99, 1] - one pixel from the frame
edge. The old filter accepted anything more than 1 px inside the
frame; the comp left the frame almost immediately under drift and
the reduction failed after the star had already been accepted.
Comps must now sit at least 5% of each dimension from every edge.
"""
from exotic.exotic import (
    VSP_COMPARISON_EDGE_MARGIN_FRACTION,
    vsp_candidate_clear_of_edges,
)

AXIS = (650, 500)  # standard MicroObservatory frame


def test_margin_fraction_is_five_percent():
    assert VSP_COMPARISON_EDGE_MARGIN_FRACTION == 0.05


def test_motivating_case_rejected():
    # The star that broke the 2026-08-16 CoRoT-2 reduction.
    assert not vsp_candidate_clear_of_edges(99, 1, AXIS)


def test_frame_center_accepted():
    assert vsp_candidate_clear_of_edges(325, 250, AXIS)


def test_exactly_on_margin_accepted():
    assert vsp_candidate_clear_of_edges(32.5, 25.0, AXIS)
    assert vsp_candidate_clear_of_edges(650 - 32.5, 500 - 25.0, AXIS)


def test_just_inside_old_filter_now_rejected():
    # Would have passed the old (1 < x < axis) test.
    assert not vsp_candidate_clear_of_edges(2, 250, AXIS)
    assert not vsp_candidate_clear_of_edges(325, 499, AXIS)


def test_all_four_edges_enforced():
    assert not vsp_candidate_clear_of_edges(10, 250, AXIS)   # left
    assert not vsp_candidate_clear_of_edges(645, 250, AXIS)  # right
    assert not vsp_candidate_clear_of_edges(325, 5, AXIS)    # bottom
    assert not vsp_candidate_clear_of_edges(325, 495, AXIS)  # top


def test_garbage_inputs_rejected_not_raised():
    assert not vsp_candidate_clear_of_edges(None, 250, AXIS)
    assert not vsp_candidate_clear_of_edges(325, 250, None)
    assert not vsp_candidate_clear_of_edges("x", 250, AXIS)
    assert not vsp_candidate_clear_of_edges(325, 250, (650,))
