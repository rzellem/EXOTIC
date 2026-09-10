import json
import sys
from types import SimpleNamespace

import numpy as np
import pytest

import exotic.exotic as exotic_module
from exotic.exotic_gui import gui_reduction_command


def test_quick_look_cli_parses_optional_init_file(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['exotic', '--quick-look', 'inits.json'])

    args = exotic_module.parse_args()

    assert args.quick_look == 'inits.json'
    assert args.realtime is None
    assert args.prereduced is None


def test_quick_look_runtime_log_isolated_under_configured_output(tmp_path):
    init_path = tmp_path / 'inits.json'
    init_path.write_text(json.dumps({
        'user_info': {'Directory to Save Plots': str(tmp_path / 'results')},
    }), encoding='utf-8')

    output = exotic_module._runtime_output_directory_from_command_line([
        '--quick-look', str(init_path),
    ])

    assert output == str(tmp_path / 'results' / 'QuickLook')


def test_quick_look_preset_is_in_memory_and_scientific_aperture_only(tmp_path):
    original = {
        'save': str(tmp_path),
        'use_psf_photometry': 'y',
        'use_adaptive_apertures': True,
        'photometer_fortuitous_variables': True,
    }

    quick = exotic_module.apply_quick_look_runtime_preset(original)

    assert original['save'] == str(tmp_path)
    assert original['use_psf_photometry'] == 'y'
    assert quick['save'] == str(tmp_path / 'QuickLook')
    assert quick['quick_look_mode'] is True
    assert quick['use_psf_photometry'] == 'n'
    assert quick['use_aperture_photometry'] == 'y'
    assert quick['use_adaptive_apertures'] is False
    assert quick['use_aperture_corrections_and_full_image_fwhm'] is False
    assert quick['fast_aperture_mask'] is False
    assert quick['photometer_fortuitous_variables'] is False
    assert quick['stellar_variability_only'] is False
    assert quick['use_ensemble_photometry_for_stellar_variability'] is False
    assert quick['fit_lightcurve_to_every_comparison_candidate'] is False


def test_quick_look_gui_command_supports_fits_and_prereduced_inputs():
    assert gui_reduction_command(3, 1) == '--quick-look'
    assert gui_reduction_command(3, 2) == '--quick-look'
    assert gui_reduction_command(2, 1) == '--reduce'
    assert gui_reduction_command(2, 2) == '--prereduced'


def test_quick_look_init_switch_is_opt_in_and_defaults_off():
    assert exotic_module.quick_look_mode_from_init_data({}) is False
    assert exotic_module.quick_look_mode_from_init_data({'optional_info': {}}) is False
    assert exotic_module.quick_look_mode_from_init_data({
        'optional_info': {'quick_look_mode': False},
    }) is False
    assert exotic_module.quick_look_mode_from_init_data({
        'optional_info': {'quick_look_mode': True},
    }) is True
    assert exotic_module.quick_look_mode_from_init_data({
        'optional_info': {'quick_look_mode': 'y'},
    }) is True


def test_quick_look_input_kind_detects_prereduced_init(tmp_path):
    init_path = tmp_path / 'inits.json'
    init_path.write_text(json.dumps({
        'user_info': {'Directory to Save Plots': str(tmp_path)},
        'optional_info': {
            'quick_look_mode': True,
            'Pre-reduced File:': str(tmp_path / 'lightcurve.txt'),
        },
    }), encoding='utf-8')

    assert exotic_module.quick_look_input_kind_from_init_file(init_path) == 2


def test_quick_look_runtime_log_honors_init_switch_for_normal_route(tmp_path):
    init_path = tmp_path / 'inits.json'
    init_path.write_text(json.dumps({
        'user_info': {'Directory to Save Plots': str(tmp_path / 'results')},
        'optional_info': {'quick_look_mode': True},
    }), encoding='utf-8')

    output = exotic_module._runtime_output_directory_from_command_line([
        '--prereduced', str(init_path),
    ])

    assert output == str(tmp_path / 'results' / 'QuickLook')


def test_missing_quick_look_init_switch_keeps_normal_output_route(tmp_path):
    init_path = tmp_path / 'inits.json'
    init_path.write_text(json.dumps({
        'user_info': {'Directory to Save Plots': str(tmp_path / 'results')},
        'optional_info': {},
    }), encoding='utf-8')

    output = exotic_module._runtime_output_directory_from_command_line([
        '--prereduced', str(init_path),
    ])

    assert output == str(tmp_path / 'results')


def test_prereduced_init_switch_enters_quick_look_before_mpi_validation(monkeypatch, tmp_path):
    class RouteReached(Exception):
        pass

    init_path = tmp_path / 'inits.json'
    init_path.write_text(json.dumps({
        'user_info': {'Directory to Save Plots': str(tmp_path / 'results')},
        'optional_info': {
            'quick_look_mode': True,
            'Pre-reduced File:': str(tmp_path / 'lightcurve.txt'),
        },
    }), encoding='utf-8')
    args = SimpleNamespace(
        multiprocess_transformations=None,
        multiprocess_lightcurve_fits=None,
        realtime=None,
        reduce=None,
        prereduced=str(init_path),
        photometry=None,
        quick_look=None,
    )

    monkeypatch.setattr(exotic_module, 'parse_args', lambda: args)
    monkeypatch.setattr(
        exotic_module,
        'validate_ultranest_mpi_runtime',
        lambda: (_ for _ in ()).throw(AssertionError('Quick Look must skip MPI validation')),
    )

    class FakeInputs:
        def __init__(self, init_opt):
            raise RouteReached

    monkeypatch.setattr(exotic_module, 'Inputs', FakeInputs)

    with pytest.raises(RouteReached):
        exotic_module._main_impl()


def test_quick_look_does_not_expand_tracking_pool_for_v_calibration():
    common = {
        'preferred_band': 'V',
        'usable_science_nextastro_v': False,
        'fortuitous_auto_scan_performed': False,
        'use_exactly_the_comps_provided': False,
    }

    assert not exotic_module.should_expand_full_field_v_calibration_pool(
        quick_look_mode=True,
        **common,
    )
    assert exotic_module.should_expand_full_field_v_calibration_pool(
        quick_look_mode=False,
        **common,
    )


def test_quick_look_reuses_selected_lm_fit_as_final_model():
    selected_attempt = {
        'fit': SimpleNamespace(inference_method='Least-squares (LM)'),
        'full_reduction_applied': False,
    }

    assert exotic_module.should_reuse_selected_comparison_fit(True, selected_attempt)
    assert not exotic_module.should_reuse_selected_comparison_fit(False, selected_attempt)
    selected_attempt['full_reduction_applied'] = True
    assert exotic_module.should_reuse_selected_comparison_fit(False, selected_attempt)


def test_quick_look_aperture_sample_spans_sequence_and_is_capped_at_twelve():
    indices = exotic_module.evenly_spaced_aperture_tuning_indices(
        524,
        max_frames=exotic_module.QUICK_LOOK_APERTURE_AUTOTUNE_MAX_FRAMES,
    )

    assert len(indices) == 12
    assert indices[0] == 0
    assert indices[-1] == 523
    assert np.all(np.diff(indices) > 0)


def test_quick_look_fit_uses_only_lm_and_labels_uncertainty(monkeypatch):
    calls = []

    def fake_lc_fitter(times, flux, unc, airmass, prior, bounds, mode=None, **kwargs):
        calls.append(mode)
        return SimpleNamespace(
            time=np.asarray(times),
            data=np.asarray(flux),
            dataerr=np.asarray(unc),
            airmass=np.asarray(airmass),
            parameters={'tmid': 1.05, 'rprs': 0.1, 'ars': 10.0, 'inc': 89.0, 'a1': 1.0},
            errors={'tmid': 0.001, 'rprs': 0.002, 'ars': 0.2, 'inc': 0.1, 'a1': 0.01},
        )

    monkeypatch.setattr(exotic_module, 'lc_fitter', fake_lc_fitter)
    monkeypatch.setattr(exotic_module, 'annotate_transit_detection_qc', lambda fit: None)
    monkeypatch.setattr(
        exotic_module,
        'run_nested_lightcurve_fit_with_rprs_posterior_retry',
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError('Quick Look must not enter posterior inference')
        ),
    )
    times = np.linspace(1.0, 1.1, 20)
    fit, fitted_flux, fitted_unc = exotic_module.fit_quick_look_lightcurve_least_squares(
        times,
        np.ones(20),
        np.full(20, 0.002),
        np.linspace(1.1, 1.5, 20),
        {'tmid': 1.05, 'rprs': 0.1, 'ars': 10.0, 'inc': 89.0, 'per': 1.0},
        {'tmid': [1.04, 1.06], 'rprs': [0.05, 0.15]},
        detrend_on_outoftransit_baseline=False,
    )

    assert calls == ['lm']
    assert fit.quick_look_mode is True
    assert fit.inference_method == 'Least-squares (LM)'
    assert fit.uncertainty_type == (
        'Local covariance with existing empirical red-noise scaling; not a posterior'
    )
    assert fit.submission_ready is False
    np.testing.assert_array_equal(fitted_flux, np.ones(20))
    np.testing.assert_array_equal(fitted_unc, np.full(20, 0.002))


def test_quick_look_partial_uses_observed_edge_as_timing_anchor(monkeypatch):
    captured = {}

    monkeypatch.setattr(
        exotic_module,
        'build_expected_transit_coverage_assessment',
        lambda *args, **kwargs: {
            'valid': True,
            'transit_fraction_observed': 0.45,
            'in_transit_points': 12,
            'pre_ingress_points': 18,
            'post_egress_points': 0,
            'covers_ingress': True,
            'covers_mid_transit': False,
            'covers_egress': False,
            'observed_segment': 'ingress-only partial transit',
        },
    )

    def fake_lc_fitter(times, flux, unc, airmass, prior, bounds, **kwargs):
        captured['bounds'] = dict(bounds)
        captured['fixed_parameter_errors'] = dict(kwargs.get('fixed_parameter_errors') or {})
        return SimpleNamespace(
            time=np.asarray(times),
            data=np.asarray(flux),
            dataerr=np.asarray(unc),
            airmass=np.asarray(airmass),
            parameters=dict(prior),
            errors={key: 0.001 for key in bounds},
        )

    monkeypatch.setattr(exotic_module, 'lc_fitter', fake_lc_fitter)
    monkeypatch.setattr(exotic_module, 'annotate_transit_detection_qc', lambda fit: None)
    times = np.linspace(1.0, 1.05, 30)
    fit, _, _ = exotic_module.fit_quick_look_lightcurve_least_squares(
        times,
        np.ones(times.size),
        np.full(times.size, 0.002),
        np.linspace(1.1, 1.4, times.size),
        {
            'tmid': 1.05, 'rprs': 0.1, 'rprs_unc': 0.002,
            'ars': 10.0, 'ars_unc': 0.2, 'inc': 88.0, 'inc_unc': 0.1,
            'per': 1.0, 'a0': 1.0, 'a2': 0.0,
        },
        {
            'tmid': [1.03, 1.07], 'rprs': [0.05, 0.15],
            'ars': [8.0, 12.0], 'inc': [85.0, 90.0],
            'a0': [0.95, 1.05], 'a2': [-1.0, 1.0],
        },
        detrend_on_outoftransit_baseline=False,
    )

    assert list(captured['bounds']) == ['tmid']
    assert set(captured['fixed_parameter_errors']) >= {'rprs', 'ars', 'inc'}
    assert fit.partial_transit_geometry_prior_assumption_applied is True
    assert fit.partial_transit_geometry_prior_assumption_mode == 'tmid_only'
    assert fit.partial_transit_geometry_prior_assumption_sampled_parameters == ['tmid']


def test_quick_look_ranked_comparisons_stop_at_first_qc_pass_without_posterior(monkeypatch):
    inference_methods = []

    monkeypatch.setattr(
        exotic_module,
        'diagnose_lightcurve_fit_inputs',
        lambda *args, **kwargs: {'usable_point_count': len(args[0])},
    )
    monkeypatch.setattr(
        exotic_module,
        'build_comparison_candidate_preflight',
        lambda *args, **kwargs: {
            'coverage_priority': 0,
            'scout': {'score': 0.0},
            'prepared_series': None,
        },
    )
    monkeypatch.setattr(
        exotic_module,
        'rank_comparison_candidate_preflight_plans',
        lambda plans: plans,
    )
    monkeypatch.setattr(
        exotic_module,
        'refit_selected_fast_comparison_on_full_lightcurve',
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError('Quick Look must not run the selected posterior refit')
        ),
    )
    monkeypatch.setattr(
        exotic_module,
        'extend_selected_comparison_live_points_if_needed',
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError('Quick Look must not extend posterior live points')
        ),
    )

    def fake_finalize(times, target_flux, comp_flux, airmass, *args, inference_method=None, **kwargs):
        inference_methods.append(inference_method)
        passed = np.nanmedian(comp_flux) < 45.0
        status = 'pass' if passed else 'fail'
        fit = SimpleNamespace(
            time=np.asarray(times),
            data=np.asarray(target_flux) / np.asarray(comp_flux),
            dataerr=np.full(len(times), 0.002),
            airmass=np.asarray(airmass),
            residuals=np.full(len(times), 0.001 if passed else 0.01),
            transit=np.ones(len(times)),
            parameters={'tmid': 1.05, 'rprs': 0.1, 'ars': 10.0, 'inc': 89.0},
            errors={'tmid': 0.001, 'rprs': 0.002, 'ars': 0.2, 'inc': 0.1},
            transit_qc={
                'status': status,
                'summary': f'synthetic {status}',
                'ktmf_metric': 4.5 if passed else 2.0,
                'ktmf_contributions': [],
            },
            transit_qc_status=status,
            transit_qc_summary=f'synthetic {status}',
        )
        return {
            'applied': True,
            'fit': fit,
            'good_times': np.asarray(times),
            'good_flux': fit.data,
            'good_unc': fit.dataerr,
            'good_airmass': np.asarray(airmass),
            'good_jd_times': 2460000.0 + np.asarray(times),
            'good_target_flux': np.asarray(target_flux),
            'good_comp_flux': np.asarray(comp_flux),
            'good_target_flux_error': np.full(len(times), 1.0),
            'good_comp_flux_error': np.full(len(times), 1.0),
            'source_indices': np.arange(len(times)),
            'duration_samples': np.array([]),
            'data_highres': None,
            'fast_ultranest_binning': {'applied': False},
            'note': 'synthetic Quick Look candidate',
        }

    monkeypatch.setattr(
        exotic_module,
        'finalize_comparison_candidate_full_reduction',
        fake_finalize,
    )
    times = np.linspace(1.0, 1.1, 20)
    calibration = {
        'method': 'aperture',
        'method_label': 'Aperture photometry',
        'a': 0,
        'an': 0,
        'comp_summaries': [
            {
                'label': 'Comp 1', 'position': (10.0, 10.0), 'aggregate_score': 0.01,
                'coverage_count': 20, 'coverage_total_frame_count': 20,
                'coverage_reference_count': 20.0, 'coverage_min_required_count': 5,
                'coverage_rejected': False, 'comp_index': 0,
            },
            {
                'label': 'Comp 2', 'position': (20.0, 20.0), 'aggregate_score': 0.02,
                'coverage_count': 20, 'coverage_total_frame_count': 20,
                'coverage_reference_count': 20.0, 'coverage_min_required_count': 5,
                'coverage_rejected': False, 'comp_index': 1,
            },
        ],
    }
    aperture_data = {
        'target': np.full((20, 1, 1), 100.0),
        'comp1': np.full((20, 1, 1), 50.0),
        'comp2': np.full((20, 1, 1), 40.0),
    }

    result = exotic_module.fit_ranked_comparison_calibration_candidates(
        times,
        2460000.0 + times,
        np.linspace(1.1, 1.5, 20),
        [0.1, 0.1, 0.1, 0.1],
        {'midT': 1.05, 'pPer': 1.0, 'rprs': 0.1, 'aRs': 10.0, 'inc': 89.0},
        calibration,
        {},
        aperture_data,
        np.full(20, 100.0),
        inference_method='lm',
        save_dir=None,
    )

    assert inference_methods == ['lm', 'lm']
    assert len(result['attempts']) == 2
    assert result['selected_result']['comp_index'] == 1
    assert result['selection_metric'] == 'first_qc_pass'
    assert result['stopped_after_first_qc_pass'] is True
