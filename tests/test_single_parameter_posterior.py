import matplotlib
matplotlib.use('Agg')
import numpy as np
import pytest
from ultranest import mlfriends

import exotic.exotic as exotic_module
from exotic.api.elca import lc_fitter, transit


@pytest.mark.parametrize('sparse_egress', [False, True], ids=['partial_transit', 'sparse_egress'])
def test_partial_transit_samples_tmid_and_saves_triangle_with_scalar_covariance_bug(
    monkeypatch, tmp_path, sparse_egress,
):
    monkeypatch.setenv('EXOTIC_ULTRANEST_CPU_GATE', '0')
    monkeypatch.setenv('EXOTIC_ULTRANEST_WORKERS', '1')
    monkeypatch.setenv('EXOTIC_ULTRANEST_MAX_NCALLS', 'auto')
    monkeypatch.setattr(exotic_module, 'lc_fitter', lc_fitter)

    original_ellipsoid = mlfriends.bounding_ellipsoid

    def legacy_ellipsoid(points, minvol=0.0):
        if points.shape[1] == 1:
            raise ValueError('Buffer has wrong number of dimensions (expected 2, got 0)')
        return original_ellipsoid(points, minvol=minvol)

    # Reproduce the affected native region path even on builds with the fix.
    monkeypatch.setattr(mlfriends, 'bounding_ellipsoid', legacy_ellipsoid)
    prior = {
        'rprs': 0.11091, 'ars': 10.05, 'inc': 87.24, 'per': 2.8997,
        'ecc': 0.0, 'omega': 90.0, 'tmid': 2460659.01263,
        'u0': 0.4, 'u1': 0.2, 'u2': 0.1, 'u3': 0.05,
        'a0': 1.0, 'a2': 0.0,
    }
    if sparse_egress:
        # CoRoT-1 b geometry and the 26-point time range from the crash log.
        prior.update({
            'rprs': 0.1388, 'ars': 4.92, 'inc': 85.1, 'per': 1.5089686,
            'tmid': 2461085.559863, 'u0': 2.1073246634873106,
            'u1': -3.8147790856466592, 'u2': 4.317206873960663,
            'u3': -1.609752451801314,
        })
        times = np.linspace(2461085.60958101, 2461085.65379324, 26)
        uncertainty = 0.009
        coverage = {
            'valid': True, 'transit_fraction_observed': 0.024,
            'in_transit_points': 2, 'pre_ingress_points': 0, 'post_egress_points': 24,
            'observed_segment': 'egress plus post-egress baseline',
        }
    else:
        times = prior['tmid'] + np.linspace(-0.07, 0.01, 31)
        uncertainty = 0.0008
        coverage = {
            'valid': True, 'transit_fraction_observed': 0.6,
            'in_transit_points': 21, 'pre_ingress_points': 10, 'post_egress_points': 0,
            'observed_segment': 'pre-ingress baseline plus ingress plus mid-transit',
        }
    truth = dict(prior, tmid=prior['tmid'] + 0.002)
    flux = transit(times, truth) + np.random.default_rng(74).normal(0, uncertainty, times.size)

    random_state = np.random.get_state()
    np.random.seed(74)
    try:
        fit = exotic_module.run_nested_lightcurve_fit_with_rprs_posterior_retry(
            times, flux, np.full(times.shape, uncertainty), np.ones(times.shape), prior,
            {
                'tmid': [prior['tmid'] - 0.14, prior['tmid'] + 0.14],
                'rprs': [0.7 * prior['rprs'], 1.3 * prior['rprs']],
                'ars': [0.9 * prior['ars'], 1.1 * prior['ars']], 'inc': [84.0, 90.0],
            },
            fixed_flux_baseline=True,
            keep_ultranest_sampler=True,
            pre_ultranest_coverage_assessment=coverage,
        )

        assert fit.mode == 'ns'
        assert fit.ns_type == 'ultranest'
        assert fit.nested_fit_fallback is False
        assert fit.sampled_keys == ['tmid']
        assert fit.partial_transit_geometry_prior_assumption_mode == 'tmid_only'
        if not sparse_egress:
            assert fit.parameters['tmid'] == pytest.approx(truth['tmid'], abs=0.001)
        assert abs(fit.parameters['tmid'] - prior['tmid']) <= 0.14
        assert np.isfinite(fit.errors['tmid']) and fit.errors['tmid'] > 0
        assert fit.results['weighted_samples']['points'].shape[1] == 1
        for key in ('rprs', 'ars', 'inc'):
            assert fit.parameters[key] == prior[key]

        # Extensions must select the same valid one-dimensional region.
        assert fit.extend_ultranest_fit(min_num_live_points=200)
        assert fit.ns_type == 'ultranest'
        output_path = exotic_module.save_final_triangle_plot(
            fit, tmp_path, 'CoRoT-1 b' if sparse_egress else 'HAT-P-3 b', '2024-12-14',
        )
        assert output_path.is_file()
        planet_token = 'CoRoT-1b' if sparse_egress else 'HAT-P-3b'
        assert (output_path.parent / f'ZoomedTrianglePlot_{planet_token}_2024-12-14.png').is_file()
    finally:
        np.random.set_state(random_state)
