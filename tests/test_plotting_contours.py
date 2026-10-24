import numpy as np

from exotic.api import plotting


def test_contour_levels_inside_surface_are_preserved():
    levels = plotting._contour_levels_within_surface([0.2, 0.5, 0.8], 0.0, 1.0)

    np.testing.assert_allclose(levels, np.array([0.2, 0.5, 0.8]))


def test_contour_levels_outside_surface_are_clipped_into_drawable_range():
    levels = plotting._contour_levels_within_surface([10.0, 20.0, 30.0], 0.0, 1.0)

    assert levels.size == 1
    assert 0.0 < levels[0] < 1.0

