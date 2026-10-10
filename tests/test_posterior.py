import numpy as np
import pytest

from exotic.posterior import summarize_posterior, format_posterior_interval


def test_unequal_weights_do_not_turn_prior_tail_into_uncertainty():
    values = np.r_[np.linspace(0, 100, 1000), np.linspace(0.2, 0.4, 1000)]
    weights = np.r_[np.full(1000, 1e-12), np.ones(1000)]
    summary = summarize_posterior(values, weights)
    assert summary['median'] == pytest.approx(0.3, abs=1e-6)
    assert summary['lower'] == pytest.approx(0.232, abs=0.001)
    assert summary['upper_limit_95'] < 0.4
    assert summary['stdev'] < 0.06


@pytest.mark.parametrize('center,scale', [(0, 0.1), (1.02, 0.025), (0.5, 0.02)])
def test_central_grazing_and_resolved_geometries_remain_reportable(center, scale):
    draws = np.maximum(0, center + scale * np.linspace(-3, 3, 1001))
    summary = summarize_posterior(draws)
    assert 0 <= summary['lower'] <= summary['median'] <= summary['upper']
    assert summary['upper_limit_95'] >= summary['upper']
    text = format_posterior_interval(summary, include_limit='b')
    assert '68% CrI' in text and '(95%)' in text and '; <' in text


def test_posterior_interval_retains_asymmetry_and_inclination_lower_limit():
    summary = summarize_posterior(np.linspace(80, 90, 1001) ** 2 / 90)
    text = format_posterior_interval(summary, 'deg', include_limit='inc')
    assert ' +/- ' not in text
    assert '/+' in text and '; >' in text and 'deg (95%)' in text


def test_zero_weight_outliers_have_no_effect():
    base = summarize_posterior([1, 2, 3])
    assert summarize_posterior([1, 2, 3, 10000], [1, 1, 1, 0]) == base


def test_central_transit_limits_match_known_half_normal_posterior():
    from scipy.stats import norm
    u = (np.arange(10000) + 0.5) / 10000
    b = 0.1 * norm.ppf((1 + u) / 2)
    summary = summarize_posterior(b)
    assert summary['median'] == pytest.approx(0.067449, abs=1e-6)
    assert summary['upper_limit_95'] == pytest.approx(0.195996, abs=1e-6)
    inc = summarize_posterior(np.degrees(np.arccos(b / 10)))
    assert inc['lower_limit_95'] == pytest.approx(np.degrees(np.arccos(summary['upper_limit_95'] / 10)), abs=1e-6)


def test_gaussian_fit_uses_posterior_weights_and_preserves_bjd_precision():
    from exotic.posterior import posterior_distribution_estimates
    offset = 2459377.
    estimate = posterior_distribution_estimates(offset + np.array([0., .01, .02]), [1, 2, 1])
    gaussian = estimate['gaussian_fit']
    assert gaussian['center'] == pytest.approx(offset + .01, abs=1e-9)
    assert gaussian['sigma'] == pytest.approx(np.sqrt(.00005), abs=1e-9)
    assert gaussian['peak'] == gaussian['center']
    assert estimate['posterior_histogram_peak']['bin_width'] > 0


def test_skewed_posterior_peak_and_gaussian_centre_are_reported_separately():
    from scipy.stats import norm
    from exotic.posterior import posterior_distribution_estimates
    u = (np.arange(10000) + .5) / 10000
    samples = norm.ppf((1 + u) / 2)
    estimates = posterior_distribution_estimates(samples)
    peak = estimates['posterior_histogram_peak']
    assert peak['peak_bin'][0] == pytest.approx(samples.min())
    assert peak['value'] == pytest.approx(samples.min() + peak['bin_width'] / 2)
    assert estimates['gaussian_fit']['center'] == pytest.approx(np.sqrt(2 / np.pi), abs=.0001)
    assert estimates['gaussian_fit']['sigma'] == pytest.approx(np.sqrt(1 - 2 / np.pi), abs=.0002)
    # Zero-weight prior/dead-point tails must not define a second estimate.
    with_tail = posterior_distribution_estimates(np.r_[samples, 1e8], np.r_[np.ones(samples.size), 0])
    assert with_tail == estimates


def test_degenerate_posterior_reports_zero_variance_instead_of_fitting_a_curve():
    from exotic.posterior import posterior_distribution_estimates
    estimate = posterior_distribution_estimates([3, 3, 3])
    assert estimate['posterior_histogram_peak']['value'] == 3
    assert estimate['posterior_histogram_peak']['bin_width'] == 0
    assert estimate['gaussian_fit']['status'] == 'zero_variance'


def test_plot_titles_can_omit_repeated_confidence_labels():
    summary = summarize_posterior([1, 2, 3])
    assert 'CrI' not in format_posterior_interval(summary, include_probability=False)
    assert 'CrI' in format_posterior_interval(summary)


def test_mathtext_interval_stacks_upper_and_lower_errors_without_probability_label():
    summary = {'median': 2457495.79043, 'error_minus': .00077, 'error_plus': .00066}
    text = format_posterior_interval(
        summary, r'BJD$_{TDB}$', include_probability=False, mathtext=True,
    )
    assert text == r'$2457495.79043^{+0.00066}_{-0.00077}$ BJD$_{TDB}$'


def test_adaptive_bins_maximize_the_documented_weighted_objective():
    from scipy.special import gammaln
    from exotic.posterior import posterior_distribution_estimates
    values = np.r_[np.linspace(-2., -1., 51), np.linspace(1., 2., 51)]
    weights = np.r_[np.ones(51), np.full(51, .25)]
    peak = posterior_distribution_estimates(values, weights)['posterior_histogram_peak']
    weights /= weights.sum()
    neff = 1 / np.sum(weights ** 2)
    scores = []
    for bins in range(1, len(values) + 1):
        counts = np.histogram(values, bins=bins, weights=weights)[0] * neff
        scores.append(neff * np.log(bins) + gammaln(bins / 2) - bins * gammaln(.5)
                      - gammaln(neff + bins / 2) + np.sum(gammaln(counts + .5)))
    assert peak['bin_count'] == np.argmax(scores) + 1
    assert peak['binning']['selected_log_score'] == pytest.approx(max(scores))
    assert peak['binning']['candidate_bin_counts'] == [1, len(values)]
    other = posterior_distribution_estimates(np.r_[values, 1e7], np.r_[weights * 17, 0])
    assert other['posterior_histogram_peak']['bin_count'] == peak['bin_count']


@pytest.mark.parametrize('key', ['b', 'inc'])
def test_bounded_geometry_has_no_gaussian_even_if_resolved(key):
    from exotic.posterior import posterior_distribution_estimates
    estimate = posterior_distribution_estimates(np.linspace(.3, .5, 51), parameter_key=key)
    assert estimate['gaussian_fit']['status'] == 'not_fitted_bounded_geometry'
    assert 'center' not in estimate['gaussian_fit']


def test_legacy_cached_geometry_gaussians_are_suppressed_without_changing_the_fit():
    from types import SimpleNamespace
    from exotic.posterior import fit_posterior_distribution_estimates
    cached = {'inc': {'gaussian_fit': {'center': 89., 'sigma': 1.}},
              'tmid': {'gaussian_fit': {'center': 2., 'sigma': .1}}}
    fit = SimpleNamespace(posterior_distribution_estimates=cached)
    report = fit_posterior_distribution_estimates(fit)
    assert 'center' not in report['inc']['gaussian_fit']
    assert report['tmid'] == cached['tmid']
    assert cached['inc']['gaussian_fit']['center'] == 89.
