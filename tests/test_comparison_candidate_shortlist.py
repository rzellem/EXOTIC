from types import SimpleNamespace

import numpy as np
import pytest

import exotic.exotic as runtime


@pytest.fixture
def candidate_search(monkeypatch):
    """Exercise real candidate ordering/selection with recorded full-fit outcomes."""
    calls = []
    messages = []
    times = np.linspace(0.0, 0.05, 12)
    comp_indices = [7, 3, 1, 2, 0, 8, 5]  # Ivan's logged preflight order.
    monkeypatch.setattr(runtime, 'log_info', lambda message, **kwargs: messages.append(message))
    monkeypatch.setattr(runtime, 'diagnose_lightcurve_fit_inputs', lambda *args, **kwargs: {})

    def preflight(*args, **kwargs):
        marker = int(np.median(args[6])) - 100
        return {'coverage_priority': 0, 'scout': {'score': 10.0 - marker}}

    monkeypatch.setattr(runtime, 'build_comparison_candidate_preflight', preflight)

    def run(statuses, *, scores=None, failed_reductions=(), **options):
        calls.clear()
        messages.clear()
        active_indices = comp_indices[:len(statuses)]
        scores = scores or [3.5] * len(comp_indices)

        def finalize(fit_times, target, comp, airmass, ld, p_dict, **kwargs):
            marker = int(np.median(comp)) - 100
            calls.append(marker)
            if marker in failed_reductions:
                return {'applied': False, 'failure_reason': 'recorded processing failure'}
            status = statuses[marker]
            delta_bic = -1.81 if status == 'marginal' else -7.31
            residuals = np.tile([-0.01, 0.01], len(fit_times) // 2)
            fit = SimpleNamespace(
                residuals=residuals,
                data=np.ones(len(fit_times)),
                parameters={'tmid': 0.5, 'rprs': 0.1, 'inc': 89.0, 'a0': 1.0, 'a2': 0.0},
                errors={'tmid': 0.001, 'rprs': 0.001, 'inc': 0.1, 'a0': 0.01, 'a2': 0.01},
                transit_qc={'status': status, 'summary': 'recorded outcome',
                            'delta_bic': delta_bic, 'ktmf_metric': scores[marker]},
                transit_qc_status=status,
                transit_qc_ktmf_metric=scores[marker],
                transit_qc_delta_bic=delta_bic,
            )
            return {'applied': True, 'fit': fit, 'good_target_flux': target,
                    'good_comp_flux': comp, 'source_indices': np.arange(len(fit_times)),
                    'duration_samples': np.array([]), 'data_highres': None}

        monkeypatch.setattr(runtime, 'finalize_comparison_candidate_full_reduction', finalize)
        # Field order differs from preflight order; the real preflight ranker must be used.
        summaries = [
            {'label': f'Comp {index + 1}', 'comp_index': index,
             'aggregate_score': 0.01, 'coverage_rejected': False}
            for index in reversed(active_indices)
        ]
        apertures = {'target': np.full((len(times), 1, 1), 1000.0)}
        apertures.update({
            f'comp{index + 1}': np.full((len(times), 1, 1), 100.0 + marker)
            for marker, index in enumerate(active_indices)
        })
        return runtime.fit_ranked_comparison_calibration_candidates(
            times, 2460000.0 + times, np.linspace(1.0, 1.2, len(times)),
            [0.1] * 4,
            {'midT': 0.5, 'pPer': 1.0, 'rprs': 0.1, 'aRs': 10.0,
             'inc': 89.0, 'ecc': 0.0, 'omega': 0.0},
            {'method': 'aperture', 'a': 0, 'an': 0, 'comp_summaries': summaries},
            psf_data={}, aper_data=apertures, target_psf_flux=np.full(len(times), 1000.0),
            **options,
        )

    return run, calls, messages


def test_ivan_outcomes_select_same_comp_after_three_candidates(candidate_search):
    run, calls, messages = candidate_search
    result = run(['fail', 'fail', 'marginal', 'fail', 'fail', 'fail', 'fail'],
                 scores=[2.39, 2.45, 3.72, 1.95, 1.71, 1.48, 1.13])

    assert calls == [0, 1, 2]
    assert result['selected_result']['label'] == 'Comp 2'
    assert result['selected_result']['transit_qc_status'] == 'marginal'
    assert result['selected_result']['ktmf_metric'] == pytest.approx(3.72)
    assert result['selection_metric'] == 'ktmf_combined_quality'
    assert result['stopped_after_acceptable_shortlist']
    assert any('3/7 preflight-ranked candidates' in message for message in messages)


def test_shortlist_uses_final_fit_metrics_among_all_three(candidate_search):
    run, calls, _ = candidate_search
    result = run(['marginal'] * 7, scores=[3.1, 3.9, 3.6, 5.0, 5.0, 5.0, 5.0])

    assert calls == [0, 1, 2]
    assert result['selected_result']['label'] == 'Comp 4'
    assert result['selection_metric'] == 'ktmf_combined_quality'


def test_all_three_fail_extends_search_until_acceptable(candidate_search):
    run, calls, messages = candidate_search
    result = run(['fail', 'fail', 'fail', 'fail', 'marginal', 'marginal', 'marginal'])

    assert calls == [0, 1, 2, 3, 4]
    assert result['selected_result']['label'] == 'Comp 1'
    assert any('continuing with the remaining' in message for message in messages)


def test_unusable_reductions_do_not_count_as_acceptable(candidate_search):
    run, calls, _ = candidate_search
    result = run(['marginal'] * 7, failed_reductions=(0, 1, 2))

    assert calls == [0, 1, 2, 3]
    assert result['selected_result']['label'] == 'Comp 3'


def test_all_fail_evaluates_entire_pool_and_retains_failed_qc_reporting(candidate_search):
    run, calls, _ = candidate_search
    result = run(['fail'] * 7)

    assert calls == list(range(7))
    assert not result['stopped_after_acceptable_shortlist']
    assert result['selected_result']['selected_despite_transit_qc']
    assert result['selected_result']['transit_qc_status'] == 'fail'


@pytest.mark.parametrize('statuses', [['marginal'], ['fail', 'marginal'], ['fail'] * 3])
def test_short_candidate_pools_are_processed_without_a_three_star_prerequisite(candidate_search, statuses):
    run, calls, messages = candidate_search
    result = run(statuses)

    assert calls == list(range(len(statuses)))
    assert result['selected_result'] is not None
    assert not result['stopped_after_acceptable_shortlist']
    assert not any('continuing with the remaining' in message for message in messages)


@pytest.mark.parametrize('options', [
    {'exit_at_first_qc_pass_solution': False},
    {'inference_method': 'lm'},
])
def test_exhaustive_and_quick_look_search_still_evaluate_all_candidates(candidate_search, options):
    run, calls, _ = candidate_search
    result = run(['marginal'] * 7, **options)

    assert calls == list(range(7))
    assert not result['stopped_after_acceptable_shortlist']
