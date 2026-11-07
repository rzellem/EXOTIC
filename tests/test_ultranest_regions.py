import numpy as np
import pytest
from ultranest import ReactiveNestedSampler
from ultranest.mlfriends import ScalingLayer

from exotic.api.ultranest_regions import OneDimensionalRegion
from exotic.api.ultranest_utils import run_reactive_sampler


@pytest.fixture(autouse=True)
def _preserve_sampling_random_state():
    state = np.random.get_state()
    np.random.seed(74)
    try:
        yield
    finally:
        np.random.set_state(state)


class _BootstrapSelection:
    def __init__(self, indices):
        self.indices = np.array(indices)

    def randint(self, high, size):
        assert size == len(self.indices)
        return self.indices


@pytest.mark.parametrize(
    'values,indices',
    [
        ([0.2, 0.5, 0.8], [0, 1, 2]),  # no held-out points
        ([0.2, 0.5, 0.8], [0, 2, 0]),  # held-out point is at selected mean
        ([0.2, 0.5, 0.8], [1, 1, 1]),  # selected variance is zero
        ([0.5, 0.5, 0.5], [0, 0, 0]),  # coincident live points
        ([0.3], [0]),
    ],
)
@pytest.mark.parametrize('minvol', [0.0, 0.25])
def test_scalar_region_keeps_live_points_for_degenerate_bootstrap_rounds(values, indices, minvol):
    points = np.array(values)[:, None]
    region = OneDimensionalRegion(points, ScalingLayer())
    region.maxradiussq, region.enlarge = region.compute_enlargement(
        nbootstraps=1, minvol=minvol, rng=_BootstrapSelection(indices),
    )
    region.create_ellipsoid(minvol=minvol)

    assert np.isfinite(region.enlarge) and region.enlarge >= 1.0
    assert region.ellipsoid_cov.shape == (1, 1)
    assert np.isfinite(region.ellipsoid_invcov).all()
    assert region.inside(points).all()
    assert 2 * region.ellipsoid_axlens[0] >= minvol
    proposals = region.sample(nsamples=100)
    assert proposals.shape[1] == 1
    assert ((proposals > 0) & (proposals < 1)).all()
    assert np.unique(proposals).size > 1


@pytest.mark.parametrize('posterior', ['gaussian', 'bimodal'])
def test_scalar_sampler_recovers_known_posterior(monkeypatch, posterior):
    monkeypatch.setenv('EXOTIC_ULTRANEST_CPU_GATE', '0')
    monkeypatch.setenv('EXOTIC_ULTRANEST_WORKERS', '1')
    monkeypatch.setenv('EXOTIC_ULTRANEST_MAX_NCALLS', 'auto')
    random_state = np.random.get_state()
    np.random.seed(74)
    try:
        if posterior == 'gaussian':
            def loglike(points):
                return -0.5 * ((points[:, 0] - 0.6) / 0.04) ** 2
        else:
            def loglike(points):
                left = -0.5 * ((points[:, 0] - 0.3) / 0.025) ** 2
                right = -0.5 * ((points[:, 0] - 0.7) / 0.025) ** 2
                return np.logaddexp(left, right)

        sampler = ReactiveNestedSampler(['tmid'], loglike, lambda points: points, vectorized=True)
        result = run_reactive_sampler(sampler, verbose=False)
        samples = result['weighted_samples']['points'][:, 0]
        weights = result['weighted_samples']['weights']

        assert type(sampler.region) is OneDimensionalRegion
        if posterior == 'gaussian':
            assert result['posterior']['mean'][0] == pytest.approx(0.6, abs=0.008)
            assert result['posterior']['stdev'][0] == pytest.approx(0.04, abs=0.008)
        else:
            assert weights[samples < 0.5].sum() == pytest.approx(0.5, abs=0.1)
            assert result['posterior']['stdev'][0] == pytest.approx(np.hypot(0.2, 0.025), abs=0.025)
    finally:
        np.random.set_state(random_state)
