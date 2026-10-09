"""Posterior summaries and presentation, separate from best-fit model errors."""

import numpy as np
from scipy.special import gammaln


def optimal_weighted_histogram(values, weights):
    """Select uniform bins by an effective-count adaptation of Knuth's score.

    This is an explicitly approximate adaptation for Monte Carlo weights, not
    a claim that nested samples are independent observations. The search covers
    every integer from one to the number of distinct positive-weight draws.
    Full-posterior endpoints are fixed, so zooming cannot change the estimate.
    """
    effective_count = float(1 / np.sum(weights ** 2))
    lower, upper = float(np.min(values)), float(np.max(values))
    maximum_bins = int(np.unique(values).size)
    if lower == upper:
        return 1, {'method': 'degenerate_distribution', 'candidate_bin_counts': [1, 1]}
    # Work in offsets to avoid loss of BJD precision in histogram construction.
    scaled = (values - lower) / (upper - lower)
    best_score, best_bins = -np.inf, 1
    for bins in range(1, maximum_bins + 1):
        indexes = np.minimum((scaled * bins).astype(int), bins - 1)
        counts = np.bincount(indexes, weights=weights, minlength=bins) * effective_count
        score = (effective_count * np.log(bins) + gammaln(bins / 2)
                 - bins * gammaln(.5) - gammaln(effective_count + bins / 2)
                 + np.sum(gammaln(counts + .5)))
        if score > best_score:
            best_score, best_bins = float(score), bins
    return best_bins, {
        'method': 'effective_count_weighted_Knuth',
        'objective': 'maximum_effective_count_multinomial_Dirichlet_log_evidence',
        'candidate_bin_counts': [1, maximum_bins], 'selected_log_score': best_score,
        'effective_sample_count': effective_count,
        'description': 'Weighted adaptation; effective count is 1/sum(normalised weights squared). '
                       'Optimal within the stated uniform-bin candidates and fixed full-posterior range. '
                       'Bin origin and Monte Carlo resolution still affect a histogram peak.',
        'reference': 'https://arxiv.org/abs/physics/0605197',
    }


def summarize_posterior(values, weights=None):
    """Summarize a marginal with posterior weights (or equal-weight samples).

    Limits are always recorded for geometry, without a threshold that decides
    whether a transit is central or whether a result should be accepted.
    """
    values = np.asarray(values, dtype=float).reshape(-1)
    if weights is None:
        weights = np.ones(values.size, dtype=float)
    else:
        weights = np.asarray(weights, dtype=float).reshape(-1)
    if weights.shape != values.shape:
        raise ValueError("Posterior values and weights must have the same shape.")
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[valid], weights[valid]
    if not values.size:
        return None
    weights = weights / np.sum(weights)
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order]) - 0.5 * weights[order]
    q05, q16, q50, q84, q95 = np.interp(
        [0.05, 0.16, 0.5, 0.84, 0.95], cumulative, values[order],
    )
    mean = np.sum(weights * values)
    return {
        'median': float(q50),
        'lower': float(q16),
        'upper': float(q84),
        'error_minus': float(q50 - q16),
        'error_plus': float(q84 - q50),
        'stdev': float(np.sqrt(np.sum(weights * (values - mean) ** 2))),
        'credible_probability': 0.68,
        'interval_method': 'equal_tail',
        'lower_limit_95': float(q05),
        'upper_limit_95': float(q95),
        'limit_probability': 0.95,
    }


def format_posterior_interval(summary, unit='', include_limit=None, include_probability=True,
                              mathtext=False):
    """Format a median and interval, with optional stacked plot uncertainties."""
    error = max(summary['error_minus'], summary['error_plus'])
    digits = max(0, 1 - int(np.floor(np.log10(error)))) if error > 0 else 6
    median = summary['median']
    # Retain precision when an error is smaller than the displayed unit.
    if mathtext:
        text = (f"${median:.{digits}f}"
                f"^{{+{summary['error_plus']:.{digits}f}}}"
                f"_{{-{summary['error_minus']:.{digits}f}}}$")
    else:
        text = (f"{median:.{digits}f} -{summary['error_minus']:.{digits}f}"
                f"/+{summary['error_plus']:.{digits}f}")
    if unit:
        text += f" {unit}"
    if include_probability:
        text += " (68% CrI)"
    if include_limit in ('b', 'inc'):
        direction = '<' if include_limit == 'b' else '>'
        key = 'upper_limit_95' if include_limit == 'b' else 'lower_limit_95'
        text += f"; {direction} {summary[key]:.{digits}f}"
        if unit:
            text += f" {unit}"
        text += " (95%)"
    return text


def fit_posterior_summary(fit, key):
    return (getattr(fit, 'posterior_summaries', None) or {}).get(key)


def posterior_distribution_estimates(values, weights=None, parameter_key=None):
    """Report a histogram peak and the maximum-likelihood normal approximation.

    The Gaussian fit uses full posterior weights: its centre is the weighted
    mean and its sigma is the weighted population standard deviation. This is
    a probability-model fit, independent of binning and plot zoom, rather than
    fitting an arbitrary Gaussian to just the tallest histogram bins.
    """
    values = np.asarray(values, dtype=float).reshape(-1)
    weights = np.ones(values.size) if weights is None else np.asarray(weights, dtype=float).reshape(-1)
    if values.shape != weights.shape:
        raise ValueError('Posterior values and weights must have the same shape.')
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[valid], weights[valid]
    if not values.size:
        return None
    weights = weights / weights.sum()
    # Subtract an anchor first, preserving precision for full BJD timestamps.
    anchor = float(values[0])
    mean_offset = float(np.sum(weights * (values - anchor)))
    center = anchor + mean_offset
    sigma = float(np.sqrt(np.sum(weights * ((values - anchor) - mean_offset) ** 2)))
    bin_count, binning = optimal_weighted_histogram(values, weights)
    lo, hi = float(np.min(values)), float(np.max(values))
    if lo == hi:
        peak, peak_bin, width = lo, [lo, hi], 0.
    else:
        counts, edges = np.histogram(values, bins=bin_count, range=(lo, hi), weights=weights)
        index = int(np.argmax(counts))
        peak_bin = [float(edges[index]), float(edges[index + 1])]
        width = peak_bin[1] - peak_bin[0]
        peak = peak_bin[0] + width / 2
    gaussian = {'center': float(center), 'peak': float(center), 'sigma': sigma,
                'method': 'weighted_maximum_likelihood_normal_on_full_posterior',
                'status': 'available' if sigma > 0 else 'zero_variance'}
    if parameter_key in ('b', 'inc'):
        gaussian = {'status': 'not_fitted_bounded_geometry',
                    'reason': 'Use the weighted posterior, asymmetric intervals and one-sided limits; '
                              'no normal or truncated-normal shape is imposed on geometry.'}
    return {
        'effective_sample_count': float(1 / np.sum(weights ** 2)),
        'posterior_histogram_peak': {'value': float(peak), 'method': 'maximum_weighted_histogram_bin',
                                     'bin_count': bin_count, 'bin_width': float(width),
                                     'peak_bin': peak_bin, 'range': [lo, hi], 'binning': binning},
        'gaussian_fit': gaussian,
    }


def fit_posterior_sample_values(fit):
    """Return physical marginal arrays and their actual posterior weights."""
    results = getattr(fit, 'results', {}) or {}
    weighted = results.get('weighted_samples') or {}
    points = np.asarray(weighted.get('points', []), dtype=float)
    weights = np.asarray(weighted.get('weights', []), dtype=float)
    if (points.ndim == 2 and weights.shape == (points.shape[0],)
            and np.isfinite(weights).all() and np.all(weights >= 0) and weights.sum() > 0):
        source = 'weighted_samples'
    else:
        points, weights = np.asarray(results.get('samples', []), dtype=float), None
        source = 'equal_weight_posterior_samples'
    keys = list(getattr(fit, 'sampled_keys', []))
    if points.ndim != 2 or not points.shape[0] or not keys or points.shape[1] < len(keys):
        return {}, None, 'unavailable'
    values = {key: points[:, index] for index, key in enumerate(keys)}
    parameters = dict(getattr(fit, 'prior', {}) or {})
    parameters.update(getattr(fit, 'parameters', {}) or {})
    ars, ecc, omega = [values.get(key, parameters.get(key)) for key in ('ars', 'ecc', 'omega')]
    if ars is not None and ecc is not None and omega is not None:
        factor = (1 - np.asarray(ecc) ** 2) / (1 + np.asarray(ecc) * np.sin(np.radians(omega)))
        with np.errstate(invalid='ignore', divide='ignore'):
            if 'b' in values:
                values['inc'] = np.degrees(np.arccos(values['b'] / (ars * factor)))
            elif 'inc' in values:
                values['b'] = ars * factor * np.cos(np.radians(values['inc']))
    return values, weights, source


def fit_posterior_distribution_estimates(fit):
    """Summarise saved joint draws without modifying fit values or weights."""
    existing = getattr(fit, 'posterior_distribution_estimates', None)
    if existing is not None:
        # Older saved fit objects may contain normal approximations to geometry.
        # Suppress those in reporting without changing the saved object itself.
        return {
            key: dict(estimate, gaussian_fit={
                'status': 'not_fitted_bounded_geometry',
                'reason': 'Geometry uses the weighted posterior, asymmetric intervals and one-sided limits.',
            }) if key in ('b', 'inc') else estimate
            for key, estimate in existing.items()
        }
    values, weights, source = fit_posterior_sample_values(fit)
    estimates = {}
    for key, samples in values.items():
        estimate = posterior_distribution_estimates(samples, weights, parameter_key=key)
        if estimate is not None:
            estimate['sample_source'] = source
            estimates[key] = estimate
    return estimates
