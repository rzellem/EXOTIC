import io
from types import SimpleNamespace

import numpy as np
import pytest

from exotic.api import ultranest_utils as utils


@pytest.fixture(autouse=True)
def sampler_environment(monkeypatch):
    monkeypatch.setenv(utils.ULTRANEST_MAX_NCALLS_ENV, '')
    monkeypatch.setenv('EXOTIC_ULTRANEST_CPU_GATE', '0')
    monkeypatch.setenv('EXOTIC_ULTRANEST_WORKERS', '1')
    monkeypatch.setenv('EXOTIC_ULTRANEST_PLAIN_PROGRESS', '1')


@pytest.mark.parametrize('requested,expected', [(None, 50), (200, 50), (20, 20), (1, 1)])
def test_ceiling_preserves_smaller_explicit_budgets(monkeypatch, requested, expected):
    monkeypatch.setenv(utils.ULTRANEST_MAX_NCALLS_ENV, '50')
    assert utils._apply_default_run_kwargs({'max_ncalls': requested})['max_ncalls'] == expected


def test_default_kwargs_preserve_api_ceiling_for_automatic_estimator():
    assert utils._apply_default_run_kwargs({'max_ncalls': 200000})['max_ncalls'] == 200000


def test_auto_can_replace_a_previous_explicit_environment_setting(monkeypatch):
    monkeypatch.setenv(utils.ULTRANEST_MAX_NCALLS_ENV, '50000')
    assert utils.configure_ultranest_max_ncalls('auto') is None
    assert utils.configured_ultranest_max_ncalls() is None


def test_scope_without_nested_sampling_has_no_budget_metadata():
    fit = SimpleNamespace()
    result = {'fit': fit}

    @utils.ultranest_budget_scope
    def quick_look():
        return result

    assert quick_look() == {'fit': fit}
    assert not hasattr(fit, 'ultranest_evaluation_budget')
    assert utils._EVALUATION_BUDGET.get() is None


@pytest.mark.parametrize('value', [0, -1, 1.5, True, 'bad', float('nan'), float('inf')])
def test_ceiling_rejects_values_ultranest_cannot_use(value):
    with pytest.raises(ValueError, match='positive integer'):
        utils.configure_ultranest_max_ncalls(value)


def test_single_evaluation_allowed_and_omitted_setting_preserves_environment(monkeypatch):
    monkeypatch.setenv(utils.ULTRANEST_MAX_NCALLS_ENV, '50')
    assert utils.configure_ultranest_max_ncalls() == 50
    assert utils.configure_ultranest_max_ncalls(1) == 1


@pytest.mark.parametrize('verbose', [True, False])
def test_budget_stop_retains_results_and_reports_incomplete_convergence(monkeypatch, verbose):
    stream = io.StringIO()
    monkeypatch.setenv(utils.ULTRANEST_MAX_NCALLS_ENV, '100')
    result = {'ncall': 105, 'samples': np.array([[0.2], [0.8]])}

    class Sampler:
        def run(self, **kwargs):
            assert kwargs['max_ncalls'] == 100
            if callable(kwargs['viz_callback']):
                kwargs['viz_callback'](info={'progress': 0.7, 'ncall': 95})
            return result

    actual = utils.run_reactive_sampler(Sampler(), {'max_ncalls': 200000}, verbose=verbose, stream=stream)
    assert actual is result
    np.testing.assert_array_equal(actual['samples'], [[0.2], [0.8]])
    assert actual['exotic_call_limit'] == {'max_ncalls': 100, 'ncall': 105, 'reached': True}
    if verbose:
        assert 'Evaluation budget reached' in stream.getvalue()
        assert '[ultranest] stopped' in stream.getvalue()
        assert '100.00%' not in stream.getvalue()
    else:
        assert stream.getvalue() == ''


def test_failure_does_not_print_success():
    stream = io.StringIO()

    class Sampler:
        def run(self, **kwargs):
            raise RuntimeError('fit failed')

    with pytest.raises(RuntimeError, match='fit failed'):
        utils.run_reactive_sampler(Sampler(), stream=stream)
    assert '[ultranest] failed' in stream.getvalue()
    assert '100.00%' not in stream.getvalue()


def test_candidate_shares_allowance_across_retries_and_retains_last_completed_fit(monkeypatch):
    monkeypatch.setenv(utils.ULTRANEST_MAX_NCALLS_ENV, '100')
    calls = []
    first_fit = SimpleNamespace(name='first')
    second_fit = SimpleNamespace(name='second')

    class Sampler:
        def __init__(self, expected, used):
            self.expected = expected
            self.used = used
            self.ncall = 0

        def run(self, **kwargs):
            assert kwargs['max_ncalls'] == self.expected
            calls.append(self.expected)
            self.ncall += self.used
            return {'ncall': self.ncall}

    @utils.ultranest_budget_scope(retain_fit=True)
    def logical_fit():
        utils.run_reactive_sampler(Sampler(100, 40), verbose=False)
        utils.record_ultranest_fit(first_fit)
        utils.run_reactive_sampler(Sampler(60, 60), verbose=False)
        utils.record_ultranest_fit(second_fit)
        utils.run_reactive_sampler(Sampler(0, 0), verbose=False)
        pytest.fail('Sampling must stop when the shared effort allowance is spent.')

    @utils.ultranest_budget_scope
    def candidate():
        fit = logical_fit()
        assert utils.ultranest_budget_exhausted()
        return {'fit': fit}

    for _ in range(2):
        result = candidate()
        assert result['fit'] is second_fit
        assert result['ultranest_evaluation_budget'] == {'max_ncalls': 100, 'ncall': 100, 'reached': True}
        assert not utils.ultranest_budget_exhausted()
    assert calls == [100, 60, 100, 60]


def test_shared_budget_preserves_actual_fit_errors(monkeypatch):
    monkeypatch.setenv(utils.ULTRANEST_MAX_NCALLS_ENV, '100')

    @utils.ultranest_budget_scope(retain_fit=True)
    def fit():
        utils.record_ultranest_fit(SimpleNamespace())
        raise ValueError('actual model failure')

    with pytest.raises(ValueError, match='actual model failure'):
        fit()
    assert not utils.ultranest_budget_exhausted()


def test_failed_sampler_consumes_its_actual_evaluations(monkeypatch):
    monkeypatch.setenv(utils.ULTRANEST_MAX_NCALLS_ENV, '100')

    class FailedSampler:
        ncall = 0

        def run(self, **kwargs):
            self.ncall = 25
            raise ValueError('actual model failure')

    class RetrySampler:
        ncall = 0

        def run(self, **kwargs):
            assert kwargs['max_ncalls'] == 75
            self.ncall = 30
            return {'ncall': 30}

    @utils.ultranest_budget_scope
    def candidate():
        with pytest.raises(ValueError, match='actual model failure'):
            utils.run_reactive_sampler(FailedSampler(), verbose=False)
        utils.run_reactive_sampler(RetrySampler(), verbose=False)
        return {'fit': SimpleNamespace()}

    assert candidate()['ultranest_evaluation_budget']['ncall'] == 55


def test_final_fit_tuple_keeps_the_matching_measurements():
    utils.configure_ultranest_max_ncalls(100)
    fit = SimpleNamespace()
    flux = np.array([1.0, 0.9])
    error = np.array([0.01, 0.01])

    @utils.ultranest_budget_scope
    def final_fit():
        return fit, flux, error

    actual = final_fit()
    assert actual[0] is fit
    assert actual[1] is flux
    assert actual[2] is error
    assert fit.ultranest_evaluation_budget == {'max_ncalls': 100, 'ncall': 0, 'reached': False}


class AutomaticSampler:
    """Deterministic sampler with a costly pilot and a resumable sampling tree."""
    x_dim = 6
    ncall = 0

    def __init__(self, tail=0.15, ess=800, error=0.4, finish_after_pilot=False):
        self.tail, self.ess, self.error = tail, ess, error
        self.finish_after_pilot = finish_after_pilot
        self.allowances = []
        self.points = object()

    def run(self, **kwargs):
        self.allowances.append(kwargs['max_ncalls'])
        self.ncall = kwargs['max_ncalls']
        if self.finish_after_pilot and len(self.allowances) > 1:
            self.ncall -= 1000
        return {
            'ncall': self.ncall, 'niter': 2200, 'ess': self.ess,
            'logzerr': self.error, 'logzerr_tail': np.log1p(self.tail),
            'samples': self.points,
        }


@pytest.mark.parametrize('dimensions,live_points,min_ess,expected', [
    (2, 200, 200, 19200), (6, 200, 200, 44800),
    (6, 400, 200, 89600), (6, 200, 400, 89600),
])
def test_automatic_pilot_scales_with_model_and_sampling_requirements(dimensions, live_points, min_ess, expected):
    sampler = AutomaticSampler(tail=0.01)
    sampler.x_dim = dimensions
    result = utils.run_reactive_sampler(sampler, {
        'max_ncalls': 200000, 'min_num_live_points': live_points, 'min_ess': min_ess,
    }, verbose=False)
    assert sampler.allowances == [expected]
    assert result['exotic_call_limit']['automatic_budget']['estimate']['pilot_max_ncalls'] == expected
    assert result['samples'] is sampler.points


def test_automatic_projection_uses_measured_cost_and_keeps_the_same_sampler():
    sampler = AutomaticSampler(finish_after_pilot=True)
    result = utils.run_reactive_sampler(sampler, {'max_ncalls': 200000}, verbose=False)
    pilot = 44800
    projected_work = 2 * (pilot / (2200 - 200)) * 200 * np.log(0.15 / 0.05)
    assert sampler.allowances == [pilot, pilot + int(np.ceil(projected_work))]
    assert result['samples'] is sampler.points
    assert result['exotic_call_limit']['reached'] is False
    estimate = result['exotic_call_limit']['automatic_budget']['estimate']
    assert estimate['pilot_tail_fraction'] == pytest.approx(0.15)
    assert estimate['projected'] is True


@pytest.mark.parametrize('tail,ess,error', [(0.8, 800, 0.4), (0.01, 20, 0.4), (0.01, 800, 2.0)])
def test_remaining_tail_ess_or_evidence_error_can_increase_automatic_effort(tail, ess, error):
    sampler = AutomaticSampler(tail=tail, ess=ess, error=error)
    result = utils.run_reactive_sampler(sampler, {'max_ncalls': 200000}, verbose=False)
    assert len(sampler.allowances) == 2
    assert 44800 < sampler.allowances[1] <= 200000
    assert result['exotic_call_limit']['reached'] is True


@pytest.mark.parametrize('settings', [{'frac_remain': 0}, {'dlogz': 0}])
def test_exact_integration_targets_remain_permitted(settings):
    sampler = AutomaticSampler()
    utils.run_reactive_sampler(sampler, {'max_ncalls': 60000, **settings}, verbose=False)
    assert sampler.allowances == [44800, 60000]


def test_automatic_budget_respects_a_smaller_api_ceiling():
    sampler = AutomaticSampler()
    result = utils.run_reactive_sampler(sampler, {'max_ncalls': 1000}, verbose=False)
    assert sampler.allowances == [1000]
    assert result['exotic_call_limit']['max_ncalls'] == 1000


def test_smaller_api_ceiling_on_a_later_refit_is_preserved():
    class NaturalSampler(AutomaticSampler):
        def run(self, **kwargs):
            result = super().run(**kwargs)
            self.ncall = result['ncall'] = 44000
            return result

    @utils.ultranest_budget_scope
    def candidate():
        utils.run_reactive_sampler(NaturalSampler(), {'max_ncalls': 200000}, verbose=False)
        retry = AutomaticSampler(ess=1)
        utils.run_reactive_sampler(retry, {'max_ncalls': 800}, verbose=False)
        assert retry.allowances == [800]
        return {'fit': SimpleNamespace()}

    candidate()


def test_natural_completion_preserves_remaining_allowance_for_other_refinements():
    class NaturalSampler(AutomaticSampler):
        def run(self, **kwargs):
            result = super().run(**kwargs)
            self.ncall = result['ncall'] = 1000
            return result

    @utils.ultranest_budget_scope
    def candidate():
        for _ in range(2):
            utils.run_reactive_sampler(NaturalSampler(), {'max_ncalls': 200000}, verbose=False)
            assert not utils.ultranest_budget_exhausted()
        return {'fit': SimpleNamespace()}

    summary = candidate()['ultranest_evaluation_budget']
    assert summary['ncall'] == 2000
    assert summary['reached'] is False
    assert summary['estimate']['projected'] is False


def test_automatic_projection_is_shared_and_cannot_be_renewed_by_retries():
    first_fit, retry_fit = SimpleNamespace(), SimpleNamespace()
    samplers = []

    @utils.ultranest_budget_scope(retain_fit=True)
    def logical_fit():
        sampler = AutomaticSampler(finish_after_pilot=True)
        samplers.append(sampler)
        utils.run_reactive_sampler(sampler, {'max_ncalls': 200000}, verbose=False)
        utils.record_ultranest_fit(first_fit)
        retry = AutomaticSampler(ess=1)
        samplers.append(retry)
        utils.run_reactive_sampler(retry, {'max_ncalls': 200000}, verbose=False)
        utils.record_ultranest_fit(retry_fit)
        utils.run_reactive_sampler(AutomaticSampler(), verbose=False)
        pytest.fail('A retry must not renew the automatically estimated allowance.')

    for _ in range(2):
        result = logical_fit()
        assert result is retry_fit
        summary = result.ultranest_evaluation_budget
        assert summary['mode'] == 'automatic'
        assert summary['ncall'] == summary['max_ncalls']
        assert samplers[-1].allowances == [1000]
        assert not utils.ultranest_budget_exhausted()


def test_missing_sampling_diagnostics_retains_results_and_reports_budget_stop():
    class MinimalSampler(AutomaticSampler):
        def run(self, **kwargs):
            return {'ncall': kwargs['max_ncalls'], 'samples': self.points}

    sampler = MinimalSampler()
    result = utils.run_reactive_sampler(sampler, verbose=False)
    assert result['samples'] is sampler.points
    assert result['exotic_call_limit']['reached'] is True
    assert 'unavailable' in result['exotic_call_limit']['automatic_budget']['estimate']['reason']


def test_automatic_budget_propagates_real_fit_failures_and_resets_context():
    @utils.ultranest_budget_scope(retain_fit=True)
    def fit():
        utils.record_ultranest_fit(SimpleNamespace())
        raise ValueError('actual model failure')

    with pytest.raises(ValueError, match='actual model failure'):
        fit()
    assert utils._EVALUATION_BUDGET.get() is None


def test_budget_spent_by_failed_initial_sampling_reports_an_ordinary_failure(monkeypatch):
    monkeypatch.setenv(utils.ULTRANEST_MAX_NCALLS_ENV, '100')

    class FailedSampler:
        ncall = 0

        def run(self, **kwargs):
            self.ncall = 100
            raise ValueError('failed initial sampling')

    @utils.ultranest_budget_scope(retain_fit=True)
    def logical_fit():
        with pytest.raises(ValueError, match='failed initial sampling'):
            utils.run_reactive_sampler(FailedSampler(), verbose=False)
        utils.run_reactive_sampler(AutomaticSampler(), verbose=False)

    with pytest.raises(RuntimeError, match='before a usable fit was available'):
        logical_fit()
    assert utils._EVALUATION_BUDGET.get() is None


def test_deferred_extension_resumes_the_original_allowance_instead_of_renewing_it():
    @utils.ultranest_budget_scope
    def initial_fit():
        sampler = AutomaticSampler(finish_after_pilot=True)
        result = utils.run_reactive_sampler(sampler, {'max_ncalls': 200000}, verbose=False)
        return SimpleNamespace(results=result)

    @utils.ultranest_budget_scope(resume_fit=True)
    def deferred_extension(fit):
        if not utils.ultranest_budget_exhausted():
            sampler = AutomaticSampler()
            utils.run_reactive_sampler(sampler, {'max_ncalls': 200000}, verbose=False)
            assert sampler.allowances == [1000]
        return fit

    fit = initial_fit()
    limit = fit.ultranest_evaluation_budget['max_ncalls']
    assert fit.ultranest_evaluation_budget['ncall'] == limit - 1000
    assert deferred_extension(fit) is fit
    assert fit.ultranest_evaluation_budget['ncall'] == limit
    assert deferred_extension(fit) is fit
    assert fit.ultranest_evaluation_budget['max_ncalls'] == limit
    assert utils._EVALUATION_BUDGET.get() is None


def test_deferred_direct_api_extension_can_read_the_recorded_automatic_budget():
    result = utils.run_reactive_sampler(AutomaticSampler(tail=0.01), verbose=False)
    fit = SimpleNamespace(results=result)

    @utils.ultranest_budget_scope(resume_fit=True)
    def deferred_extension(fit):
        assert utils.ultranest_budget_exhausted()
        return fit

    assert deferred_extension(fit) is fit
    assert fit.ultranest_evaluation_budget['reached'] is True
