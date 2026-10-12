import numpy as np
import pytest
from scipy.ndimage import gaussian_filter

from exotic import exotic as ex

NY = NX = 300
YY, XX = np.mgrid[0:NY, 0:NX]


def donut(diameter, cx=150.3, cy=149.7, total=2e6, sky=500.0, seed=1):
    radius = diameter / 2
    r = np.hypot(XX - cx, YY - cy)
    mask = gaussian_filter(((r <= radius) & (r >= 0.35 * radius)).astype(float), 1.2)
    rng = np.random.default_rng(seed)
    return rng.poisson(sky + total * mask / mask.sum()).astype(float)


def gaussian_star(fwhm, cx=150.0, cy=150.0, amp=5e4, sky=500.0):
    sigma = fwhm / 2.355
    return sky + amp * np.exp(-((XX - cx) ** 2 + (YY - cy) ** 2) / (2 * sigma ** 2))


@pytest.mark.parametrize("fwhm", [3, 5, 8])
def test_ordinary_stars_keep_the_default_box(fwhm):
    assert ex.adaptive_centroid_box(gaussian_star(fwhm), [150, 150], 15) == 15


def test_empty_sky_returns_nan_and_keeps_box():
    sky = np.random.default_rng(3).poisson(np.full((NY, NX), 500.0)).astype(float)
    assert np.isnan(ex.estimate_star_outer_radius(sky, [150, 150]))
    assert ex.adaptive_centroid_box(sky, [150, 150], 15) == 15


def test_bright_neighbour_does_not_inflate_a_small_star():
    img = gaussian_star(5) + gaussian_star(5, cx=190.0, cy=150.0, amp=2e5) - 500.0
    assert ex.adaptive_centroid_box(img, [150, 150], 15) == 15


@pytest.mark.parametrize("diameter", [36, 60])
def test_defocused_donut_widens_box_and_fits(diameter):
    img = donut(diameter)
    box = ex.adaptive_centroid_box(img, [150, 150], 15)
    assert box > 15
    fit = ex.fit_centroid(img, [150, 150], 0)
    assert fit[0] == pytest.approx(150.3, abs=0.5)
    assert fit[1] == pytest.approx(149.7, abs=0.5)
    assert fit[2] > 0
    assert min(fit[3], fit[4]) > diameter / 6
