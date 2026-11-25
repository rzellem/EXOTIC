"""Fit hybrid4 to the same PHOENIX intensity profiles used for Claret LD.

The flux implementation in limbdarkfast.py is vendored from David Kipping's
MIT-licensed limbdarkfast, revision 120caecc7d6d1d60c10cf0ed4ec27e4dede952db.
See limbdarkfast_LICENSE.txt. No external checkout is required.
"""

import numpy as np

try:
    from .limbdarkfast import intensity
except ImportError:
    from limbdarkfast import intensity


HYBRID4_KEYS = ('h0', 'h1', 'h2', 'h3')
LIMB_DARKENING_LAW_KEY = '_limb_darkening_law'


def fit_hybrid4_profiles(mu, profiles, uncertainties, mumin=0.1):
    """Weighted least squares in hybrid4's four shape weights.

    Use precisely the Claret fit's mu selection and uncertainty weighting.
    Physicality and residuals are diagnostics, never acceptance thresholds.
    Coefficient errors describe the atmosphere-profile fit; transit posteriors
    remain conditional on fixed LD, as in the existing Claret transit fit.
    """
    mu = np.asarray(mu, dtype=float)
    profiles = np.atleast_2d(np.asarray(profiles, dtype=float))
    uncertainties = np.atleast_2d(np.asarray(uncertainties, dtype=float))
    if mu.ndim != 1 or profiles.shape != uncertainties.shape or profiles.shape[1] != mu.size:
        raise ValueError('Atmosphere mu, intensity and uncertainty arrays have incompatible shapes.')
    selected = mu > mumin
    basis = np.column_stack([
        1.0 - intensity(mu[selected], row, 'hybrid4') for row in np.eye(4)
    ])
    results = []
    for profile, sigma in zip(profiles, uncertainties):
        weights_sigma = np.where(sigma[selected] == 0.0, 1e-10, sigma[selected])
        design = basis / weights_sigma[:, None]
        target = (1.0 - profile[selected]) / weights_sigma
        coefficients, _, rank, _ = np.linalg.lstsq(design, target, rcond=None)
        residual = profile[selected] - intensity(mu[selected], coefficients, 'hybrid4')
        chi2 = float(np.sum((residual / weights_sigma) ** 2))
        dof = int(np.count_nonzero(selected) - 4)
        inverse_design = np.linalg.pinv(design)
        covariance = inverse_design @ inverse_design.T
        if dof > 0:
            covariance *= chi2 / dof  # lmfit's existing scale_covar convention
        evaluated = intensity(mu, coefficients, 'hybrid4')
        order = np.argsort(mu)
        results.append({
            'law': 'hybrid4',
            'coefficient_names': list(HYBRID4_KEYS),
            'coefficients': coefficients.tolist(),
            'coefficient_errors': np.sqrt(np.maximum(np.diag(covariance), 0.0)).tolist(),
            'coefficient_covariance': covariance.tolist(),
            'mu': mu.tolist(),
            'intensity': profile.tolist(),
            'intensity_uncertainty': sigma.tolist(),
            'fit_mask': selected.tolist(),
            'mu_min_exclusive': float(mumin),
            'profile_fit_chi2': chi2,
            'profile_fit_dof': dof,
            'profile_design_rank': int(rank),
            'coefficient_covariance_note': (
                'Weighted least-squares covariance, scaled by reduced profile chi2 when dof > 0; '
                'a deficient design rank leaves coefficient combinations undetermined.'
            ),
            'profile_fit_rms': float(np.sqrt(np.mean(residual ** 2))),
            'profile_fit_max_abs_residual': float(np.max(np.abs(residual))),
            'minimum_model_intensity': float(np.min(evaluated)),
            'model_intensity_monotonic_in_mu': bool(np.all(np.diff(evaluated[order]) >= 0.0)),
            'source': 'Same LDTk PHOENIX profile means, uncertainties and mu mask as the Claret fit',
            'upstream_revision': '120caecc7d6d1d60c10cf0ed4ec27e4dede952db',
        })
    return results
