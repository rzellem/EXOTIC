"""Exercise compiled dependencies, including the bundled Windows ABI wheels."""

import numpy as np
import pytest
import sep
from ultranest.mlfriends import compute_mean_pair_distance, int_dtype


def test_ultranest_native_pair_distance():
    points = np.array([[0.0, 0.0], [3.0, 4.0]])
    # Match the native integer buffer supplied by UltraNest's sampler.
    clusters = np.ones(len(points), dtype=int_dtype)

    assert compute_mean_pair_distance(points, clusters) == pytest.approx(5.0)


def test_sep_native_source_extraction():
    image = np.zeros((20, 20), dtype=np.float64)
    image[7, 11] = 100.0

    sources = sep.extract(image, 10.0, minarea=1, filter_kernel=None)

    assert len(sources) == 1
    assert sources['x'][0] == pytest.approx(11.0)
    assert sources['y'][0] == pytest.approx(7.0)
