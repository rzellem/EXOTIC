import io
from contextlib import contextmanager
from types import SimpleNamespace

import pytest

import exotic.api.ultranest_utils as utils


@pytest.fixture(autouse=True)
def clean_gate_environment(monkeypatch):
    monkeypatch.delenv(utils.ULTRANEST_CPU_GATE_ENV, raising=False)
    monkeypatch.delenv(utils.ULTRANEST_CPU_THRESHOLD_ENV, raising=False)
    monkeypatch.delenv(utils.RESOURCE_WAIT_EVENTS_ENV, raising=False)


def test_cpu_gate_waits_for_load_to_drop_before_workers_and_sampling(monkeypatch):
    events = []
    loads = iter([94.0, 51.0, 50.0])

    def sample(interval):
        events.append(('sample', interval))
        return next(loads)

    @contextmanager
    def workers(sampler):
        events.append('workers')
        yield 1, 'none', None

    def run(**kwargs):
        events.append('run')
        return {'fit': 'complete'}

    monkeypatch.setattr(utils.psutil, 'cpu_percent', sample)
    monkeypatch.setattr(utils.time, 'sleep', lambda seconds: events.append(('sleep', seconds)))
    monkeypatch.setattr(utils, '_parallel_vectorized_loglike', workers)
    monkeypatch.setattr(utils, 'is_mpi_worker_process', lambda: False)
    stream = io.StringIO()
    result = utils.run_reactive_sampler(SimpleNamespace(run=run), stream=stream)

    assert result == {'fit': 'complete'}
    assert events == [
        ('sample', 1.0), ('sleep', 5.0),
        ('sample', 1.0), ('sleep', 5.0),
        ('sample', 1.0), 'workers', 'run',
    ]
    assert 'system CPU 94.0% > 50%' in stream.getvalue()
    assert 'CPU gate ready' in stream.getvalue()


@pytest.mark.parametrize('disabled', ['0', 'off', 'false', 'no', 'n', False])
def test_disabled_cpu_gate_does_not_measure_or_sleep(monkeypatch, disabled):
    utils.configure_ultranest_cpu_gate(enabled=disabled, threshold='unused')

    def unexpected(*args, **kwargs):
        pytest.fail('disabled CPU gate must not sample or wait')

    monkeypatch.setattr(utils.psutil, 'cpu_percent', unexpected)
    monkeypatch.setattr(utils.time, 'sleep', unexpected)
    sampler = SimpleNamespace(run=lambda **kwargs: 'complete')
    assert utils.run_reactive_sampler(sampler, verbose=False) == 'complete'


@pytest.mark.parametrize('threshold,load', [(0, 0), (100, 100), (75, 74.5), (50, 50)])
def test_cpu_gate_permits_load_at_or_below_configured_threshold(monkeypatch, threshold, load):
    utils.configure_ultranest_cpu_gate(threshold=threshold)
    monkeypatch.setattr(utils.psutil, 'cpu_percent', lambda **kwargs: load)
    monkeypatch.setattr(utils.time, 'sleep', lambda seconds: pytest.fail('load should be admitted'))
    utils._wait_for_cpu_capacity()


def test_cpu_gate_defaults_and_environment_inheritance(monkeypatch):
    assert utils.configure_ultranest_cpu_gate() == (True, 50.0)
    monkeypatch.setenv(utils.ULTRANEST_CPU_THRESHOLD_ENV, '67.5')
    assert utils.configure_ultranest_cpu_gate(enabled='y') == (True, 67.5)


def test_silent_cpu_gate_still_waits(monkeypatch):
    loads = iter([80.0, 30.0])
    sleeps = []
    stream = io.StringIO()
    monkeypatch.setattr(utils.psutil, 'cpu_percent', lambda **kwargs: next(loads))
    monkeypatch.setattr(utils.time, 'sleep', sleeps.append)
    utils._wait_for_cpu_capacity(verbose=False, stream=stream)
    assert sleeps == [5.0]
    assert stream.getvalue() == ''


@pytest.mark.parametrize('loads', [[80.0, 30.0], [30.0]])
def test_gate_events_include_sampling_and_waits_even_in_silent_fits(monkeypatch, capsys, loads):
    monkeypatch.setenv(utils.RESOURCE_WAIT_EVENTS_ENV, '1')
    now = SimpleNamespace(value=100.0)
    readings = iter(loads)
    monkeypatch.setattr(utils.time, 'monotonic', lambda: now.value)

    def sample(interval):
        now.value += interval
        return next(readings)

    monkeypatch.setattr(utils.psutil, 'cpu_percent', sample)
    monkeypatch.setattr(utils.time, 'sleep', lambda seconds: setattr(now, 'value', now.value + seconds))
    stream = io.StringIO()
    worker = f'{utils.os.getpid()}:{utils.threading.get_ident()}'

    def run(**kwargs):
        # The gate has ended before the actual sampler is started.
        assert capsys.readouterr().out.splitlines() == [
            f'[exotic-runtime] gate_begin {worker} 100.000000000',
            f'[exotic-runtime] gate_end {worker} {now.value:.9f}',
            f'[exotic-runtime] processing_begin {worker} {now.value:.9f}',
        ]
        return 'complete'

    assert utils.run_reactive_sampler(SimpleNamespace(run=run), verbose=False, stream=stream) == 'complete'
    assert now.value == (107 if len(loads) == 2 else 101)
    assert stream.getvalue() == ''
    assert capsys.readouterr().out.splitlines() == [
        f'[exotic-runtime] processing_end {worker} {now.value:.9f}',
    ]


def test_gate_events_end_and_preserve_the_actual_measurement_failure(monkeypatch, capsys):
    monkeypatch.setenv(utils.RESOURCE_WAIT_EVENTS_ENV, '1')
    monkeypatch.setattr(utils.time, 'monotonic', lambda: 100.0)

    def sample(interval):
        raise OSError('CPU counters unavailable')

    monkeypatch.setattr(utils.psutil, 'cpu_percent', sample)
    with pytest.raises(OSError, match='CPU counters unavailable'):
        utils._wait_for_cpu_capacity(verbose=False)
    worker = f'{utils.os.getpid()}:{utils.threading.get_ident()}'
    assert capsys.readouterr().out.splitlines() == [
        f'[exotic-runtime] gate_begin {worker} 100.000000000',
        f'[exotic-runtime] gate_end {worker} 100.000000000',
    ]


def test_disabled_gate_has_no_events_or_excluded_time(monkeypatch, capsys):
    monkeypatch.setenv(utils.RESOURCE_WAIT_EVENTS_ENV, '1')
    utils.configure_ultranest_cpu_gate(enabled=False)
    monkeypatch.setattr(utils.psutil, 'cpu_percent', lambda **kw: pytest.fail('Gate is disabled'))
    utils._wait_for_cpu_capacity(verbose=False)
    assert capsys.readouterr().out == ''


def test_cpu_gate_has_no_wait_timeout(monkeypatch):
    loads = iter([90.0] * 100 + [49.0])
    sleeps = []
    monkeypatch.setattr(utils.psutil, 'cpu_percent', lambda **kwargs: next(loads))
    monkeypatch.setattr(utils.time, 'sleep', sleeps.append)
    utils._wait_for_cpu_capacity(verbose=False)
    assert len(sleeps) == 100


def test_cpu_measurement_failure_is_reported_without_starting_sampler(monkeypatch):
    def failed_sample(**kwargs):
        raise OSError('CPU counters unavailable')

    monkeypatch.setattr(utils.psutil, 'cpu_percent', failed_sample)
    sampler = SimpleNamespace(run=lambda **kwargs: pytest.fail('must not bypass measurement failure'))
    with pytest.raises(OSError, match='CPU counters unavailable'):
        utils.run_reactive_sampler(sampler, verbose=False)


@pytest.mark.parametrize('value', ['bad', 'nan', 'inf', -1, 101])
def test_cpu_gate_reports_invalid_percentage(value):
    with pytest.raises(ValueError, match='CPU percentage'):
        utils.configure_ultranest_cpu_gate(threshold=value)


def test_cpu_gate_reports_invalid_switch():
    with pytest.raises(ValueError, match='on/off boolean'):
        utils.configure_ultranest_cpu_gate(enabled='maybe')


@pytest.mark.parametrize('load', [float('nan'), float('inf'), -1, 101])
def test_cpu_gate_reports_invalid_measurement(monkeypatch, load):
    monkeypatch.setattr(utils.psutil, 'cpu_percent', lambda **kwargs: load)
    with pytest.raises(RuntimeError, match='Cannot measure system CPU usage'):
        utils._wait_for_cpu_capacity()


def test_mpi_cpu_gate_root_broadcasts_admission_and_workers_do_not_sample(monkeypatch):
    messages = []
    calls = []

    def broadcast(value, root):
        assert root == 0
        messages.append(value)
        return value

    comm = SimpleNamespace(bcast=broadcast)
    monkeypatch.setattr(utils, '_wait_for_cpu_capacity', lambda **kwargs: calls.append('sample'))
    utils._wait_for_sampler_cpu_capacity(SimpleNamespace(use_mpi=True, mpi_rank=0, comm=comm))
    assert messages == [None]
    utils._wait_for_sampler_cpu_capacity(SimpleNamespace(use_mpi=True, mpi_rank=1, comm=comm))
    assert calls == ['sample']


def test_mpi_cpu_gate_propagates_root_failure_to_workers(monkeypatch):
    messages = []

    def failed_wait(**kwargs):
        raise OSError('CPU counters unavailable')

    comm = SimpleNamespace(bcast=lambda value, root: messages.append(value))
    monkeypatch.setattr(utils, '_wait_for_cpu_capacity', failed_wait)
    with pytest.raises(OSError, match='CPU counters unavailable'):
        utils._wait_for_sampler_cpu_capacity(SimpleNamespace(use_mpi=True, mpi_rank=0, comm=comm))
    comm.bcast = lambda value, root: messages[0]
    with pytest.raises(RuntimeError, match='CPU counters unavailable'):
        utils._wait_for_sampler_cpu_capacity(SimpleNamespace(use_mpi=True, mpi_rank=1, comm=comm))
