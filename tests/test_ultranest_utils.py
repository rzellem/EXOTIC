import io
import logging
import sys
import types

import numpy as np

import exotic.api.ultranest_utils as ultranest_utils
from exotic.api.ultranest_utils import run_reactive_sampler
from exotic.api.ultranest_utils import supports_ultranest_live_status


_MPI_ENV_KEYS = (
    "OMPI_COMM_WORLD_SIZE",
    "PMI_SIZE",
    "PMIX_SIZE",
    "MV2_COMM_WORLD_SIZE",
    "OMPI_COMM_WORLD_RANK",
    "PMI_RANK",
    "PMIX_RANK",
    "MV2_COMM_WORLD_RANK",
)


class _FakeStream(io.StringIO):
    def __init__(self, tty):
        super().__init__()
        self._tty = tty

    def isatty(self):
        return self._tty


def _reset_ultranest_env(monkeypatch):
    monkeypatch.delenv("EXOTIC_ULTRANEST_PLAIN_PROGRESS", raising=False)
    monkeypatch.delenv("EXOTIC_ULTRANEST_RICH_PROGRESS", raising=False)
    monkeypatch.delenv("EXOTIC_ULTRANEST_MIN_NUM_LIVE_POINTS", raising=False)
    monkeypatch.delenv("EXOTIC_ULTRANEST_MIN_LIVE_POINTS", raising=False)
    monkeypatch.delenv("EXOTIC_ULTRANEST_WORKERS", raising=False)
    monkeypatch.delenv("EXOTIC_ULTRANEST_WORKER_BACKEND", raising=False)
    monkeypatch.delenv("NEXTASTRO_EXOTIC_ULTRANEST_WORKERS", raising=False)
    for env_key in _MPI_ENV_KEYS:
        monkeypatch.delenv(env_key, raising=False)
    monkeypatch.delenv("CI", raising=False)


def test_supports_ultranest_live_status_overrides(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("TERM", "xterm-256color")
    stream = _FakeStream(tty=True)

    monkeypatch.setenv("EXOTIC_ULTRANEST_RICH_PROGRESS", "1")
    assert supports_ultranest_live_status(stream=stream) is True

    monkeypatch.setenv("EXOTIC_ULTRANEST_PLAIN_PROGRESS", "1")
    assert supports_ultranest_live_status(stream=stream) is False


def test_run_reactive_sampler_compat_mode_emits_progress(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    stream = _FakeStream(tty=False)

    class FakeSampler:
        def __init__(self):
            self.kwargs = None

        def run(self, **kwargs):
            self.kwargs = kwargs
            callback = kwargs["viz_callback"]
            callback(None, {"it": 120, "ncall": 540, "logz": -151.9, "logz_remain": -150.0})
            callback(None, {"it": 360, "ncall": 907, "logz": -139.2, "logz_remain": -141.0})
            return {"status": "ok"}

    sampler = FakeSampler()
    result = run_reactive_sampler(
        sampler,
        run_kwargs={"max_ncalls": 1000},
        verbose=True,
        stream=stream,
        interval_seconds=0.0,
    )

    assert result == {"status": "ok"}
    assert sampler.kwargs["show_status"] is False
    assert callable(sampler.kwargs["viz_callback"])
    output = stream.getvalue()
    assert "Using simple progress updates" in output
    assert "[ultranest] running" in output
    assert "[ultranest] done 100.00%" in output


def test_run_reactive_sampler_silent_mode(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    stream = _FakeStream(tty=False)

    class FakeSampler:
        def __init__(self):
            self.kwargs = None

        def run(self, **kwargs):
            self.kwargs = kwargs
            return {"status": "ok"}

    sampler = FakeSampler()
    result = run_reactive_sampler(
        sampler,
        run_kwargs={"max_ncalls": 1000},
        verbose=False,
        stream=stream,
    )

    assert result == {"status": "ok"}
    assert sampler.kwargs["show_status"] is False
    assert sampler.kwargs["viz_callback"] is False
    assert stream.getvalue() == ""


def test_run_reactive_sampler_applies_fast_defaults(monkeypatch):
    _reset_ultranest_env(monkeypatch)

    class FakeSampler:
        def __init__(self):
            self.kwargs = None

        def run(self, **kwargs):
            self.kwargs = kwargs
            return {"status": "ok"}

    sampler = FakeSampler()
    run_reactive_sampler(
        sampler,
        run_kwargs={"max_ncalls": 1000},
        verbose=False,
    )

    assert sampler.kwargs["min_num_live_points"] == 200
    assert sampler.kwargs["min_ess"] == 200
    assert sampler.kwargs["dlogz"] == 1.0
    assert sampler.kwargs["dKL"] == 1.0
    assert sampler.kwargs["frac_remain"] == 0.05
    assert sampler.kwargs["max_num_improvement_loops"] == 1
    assert sampler.kwargs["max_ncalls"] == 1000


def test_run_reactive_sampler_uses_env_live_point_override(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("EXOTIC_ULTRANEST_MIN_NUM_LIVE_POINTS", "320")

    class FakeSampler:
        def __init__(self):
            self.kwargs = None

        def run(self, **kwargs):
            self.kwargs = kwargs
            return {"status": "ok"}

    sampler = FakeSampler()
    run_reactive_sampler(sampler, verbose=False)

    assert sampler.kwargs["min_num_live_points"] == 320


def test_run_reactive_sampler_auto_scales_draw_size_by_workers_and_ram(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("EXOTIC_ULTRANEST_WORKERS", "72")
    monkeypatch.setenv("EXOTIC_ULTRANEST_WORKER_BACKEND", "thread")
    monkeypatch.setattr(ultranest_utils, "_available_cpu_count", lambda: 72)
    monkeypatch.setattr(ultranest_utils, "_system_total_memory_bytes", lambda: 128 * 1024 ** 3)

    class FakeSampler:
        def __init__(self):
            self.ndraw_min = 128
            self.ndraw_max = 65536
            self.draw_multiple = True
            self.x_dim = 6
            self.num_params = 6
            self.loglike = lambda params: np.zeros(np.asarray(params).shape[0])
            self.kwargs = None

        def run(self, **kwargs):
            self.kwargs = kwargs
            return {"status": "ok"}

    sampler = FakeSampler()
    run_reactive_sampler(sampler, verbose=False)

    assert sampler.ndraw_min == (
        72
        * ultranest_utils.HIGH_AUTO_POINTS_PER_WORKER
        * ultranest_utils.AUTO_POINTS_PER_WORKER_MULTIPLIER
    )
    assert sampler.ndraw_max == 65536


def test_run_reactive_sampler_auto_uses_smaller_chunks_when_ram_per_cpu_is_low(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("EXOTIC_ULTRANEST_WORKERS", "72")
    monkeypatch.setenv("EXOTIC_ULTRANEST_WORKER_BACKEND", "thread")
    monkeypatch.setattr(ultranest_utils, "_available_cpu_count", lambda: 72)
    monkeypatch.setattr(ultranest_utils, "_system_total_memory_bytes", lambda: 16 * 1024 ** 3)

    class FakeSampler:
        def __init__(self):
            self.ndraw_min = 128
            self.ndraw_max = 65536
            self.draw_multiple = True
            self.x_dim = 6
            self.num_params = 6
            self.loglike = lambda params: np.zeros(np.asarray(params).shape[0])
            self.kwargs = None

        def run(self, **kwargs):
            self.kwargs = kwargs
            return {"status": "ok"}

    sampler = FakeSampler()
    run_reactive_sampler(sampler, verbose=False)

    assert sampler.ndraw_min == (
        72
        * ultranest_utils.MIN_AUTO_POINTS_PER_WORKER
        * ultranest_utils.AUTO_POINTS_PER_WORKER_MULTIPLIER
    )


def test_configured_ultranest_workers_defaults_to_available_cpu_count(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setattr(ultranest_utils.os, "process_cpu_count", lambda: 12, raising=False)

    assert ultranest_utils._configured_ultranest_workers() == 12


def test_configured_ultranest_workers_accepts_auto_override(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("EXOTIC_ULTRANEST_WORKERS", "auto")
    monkeypatch.setattr(ultranest_utils.os, "process_cpu_count", lambda: 10, raising=False)

    assert ultranest_utils._configured_ultranest_workers() == 10


def test_configured_ultranest_workers_preserves_numeric_override(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("EXOTIC_ULTRANEST_WORKERS", "3")
    monkeypatch.setattr(ultranest_utils.os, "process_cpu_count", lambda: 12, raising=False)

    assert ultranest_utils._configured_ultranest_workers() == 3


def test_run_reactive_sampler_parallelizes_vectorized_loglike_batches(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("EXOTIC_ULTRANEST_WORKERS", "3")
    monkeypatch.setenv("EXOTIC_ULTRANEST_WORKER_BACKEND", "thread")

    class FakeSampler:
        def __init__(self):
            self.kwargs = None
            self.chunk_sizes = []

            def loglike(points):
                self.chunk_sizes.append(len(points))
                return points[:, 0]

            self.loglike = loglike

        def run(self, **kwargs):
            self.kwargs = kwargs
            values = self.loglike(np.arange(12, dtype=float).reshape(6, 2))
            return {"values": values}

    sampler = FakeSampler()
    result = run_reactive_sampler(sampler, verbose=False)

    assert result["values"].tolist() == [0, 2, 4, 6, 8, 10]
    assert sorted(sampler.chunk_sizes) == [2, 2, 2]


def test_run_reactive_sampler_auto_workers_uses_available_cpu_count(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("EXOTIC_ULTRANEST_WORKER_BACKEND", "thread")
    monkeypatch.setattr(ultranest_utils.os, "process_cpu_count", lambda: 4, raising=False)

    class FakeSampler:
        def __init__(self):
            self.chunk_sizes = []

            def loglike(points):
                self.chunk_sizes.append(len(points))
                return points[:, 0]

            self.loglike = loglike

        def run(self, **kwargs):
            return {"values": self.loglike(np.arange(16, dtype=float).reshape(8, 2))}

    sampler = FakeSampler()
    result = run_reactive_sampler(sampler, verbose=False)

    assert result["values"].tolist() == [0, 2, 4, 6, 8, 10, 12, 14]
    assert sorted(sampler.chunk_sizes) == [2, 2, 2, 2]


def test_process_backend_disables_parent_gc_while_pool_is_active(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("EXOTIC_ULTRANEST_WORKERS", "2")
    monkeypatch.setenv("EXOTIC_ULTRANEST_WORKER_BACKEND", "process")
    monkeypatch.setattr(ultranest_utils.sys, "platform", "linux")
    monkeypatch.setattr(ultranest_utils, "_is_colab_runtime", lambda: False)
    monkeypatch.setattr(
        ultranest_utils,
        "get_mpi_status",
        lambda: {"available": False, "size": 1, "rank": 0, "source": "test", "error": None},
    )

    cleanup_calls = []
    monkeypatch.setattr(
        ultranest_utils,
        "suppress_inherited_tk_cleanup_in_worker",
        lambda: cleanup_calls.append("suppress"),
    )

    pool_events = []

    class FakePool:
        def __init__(self, processes, initializer):
            pool_events.append(("init", processes, initializer))

        def map(self, func, chunks):
            pool_events.append(("map", ultranest_utils.gc.isenabled()))
            return [func(chunk) for chunk in chunks]

        def close(self):
            pool_events.append(("close", None))

        def join(self):
            pool_events.append(("join", None))

    class FakeContext:
        Pool = FakePool

    monkeypatch.setattr(ultranest_utils.multiprocessing, "get_context", lambda _method: FakeContext())

    class FakeSampler:
        def __init__(self):
            self.chunk_sizes = []

            def loglike(points):
                self.chunk_sizes.append(len(points))
                return points[:, 0]

            self.loglike = loglike

        def run(self, **kwargs):
            assert ultranest_utils.gc.isenabled() is False
            return {"values": self.loglike(np.arange(8, dtype=float).reshape(4, 2))}

    gc_was_enabled = ultranest_utils.gc.isenabled()
    ultranest_utils.gc.enable()
    try:
        sampler = FakeSampler()
        result = run_reactive_sampler(sampler, verbose=False)
        assert ultranest_utils.gc.isenabled() is True
    finally:
        if not gc_was_enabled:
            ultranest_utils.gc.disable()

    assert result["values"].tolist() == [0, 2, 4, 6]
    assert sorted(sampler.chunk_sizes) == [2, 2]
    assert cleanup_calls == ["suppress"]
    assert pool_events[0] == ("init", 2, ultranest_utils.suppress_inherited_tk_cleanup_in_worker)
    assert ("map", False) in pool_events
    assert pool_events[-2:] == [("close", None), ("join", None)]


def test_process_worker_initializer_suppresses_inherited_tk_destructors(monkeypatch):
    class FakeImage:
        def __del__(self):
            raise RuntimeError("main thread is not in main loop")

    class FakeVariable:
        def __del__(self):
            raise RuntimeError("main thread is not in main loop")

    original_image_del = FakeImage.__del__
    original_variable_del = FakeVariable.__del__
    fake_tkinter = types.SimpleNamespace(Image=FakeImage, Variable=FakeVariable)
    monkeypatch.setitem(sys.modules, "tkinter", fake_tkinter)

    assert ultranest_utils._suppress_inherited_tk_cleanup() is True

    assert FakeImage._exotic_worker_original_del is original_image_del
    assert FakeVariable._exotic_worker_original_del is original_variable_del
    assert FakeImage().__del__() is None
    assert FakeVariable().__del__() is None
    assert ultranest_utils._suppress_inherited_tk_cleanup() is False


def test_run_reactive_sampler_preserves_explicit_live_point_override(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("EXOTIC_ULTRANEST_MIN_NUM_LIVE_POINTS", "320")

    class FakeSampler:
        def __init__(self):
            self.kwargs = None

        def run(self, **kwargs):
            self.kwargs = kwargs
            return {"status": "ok"}

    sampler = FakeSampler()
    run_reactive_sampler(
        sampler,
        run_kwargs={"min_num_live_points": 450},
        verbose=False,
    )

    assert sampler.kwargs["min_num_live_points"] == 450


def test_run_reactive_sampler_silences_mpi_worker_rank(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setattr(
        ultranest_utils,
        "get_mpi_status",
        lambda: {
            "available": True,
            "size": 4,
            "rank": 2,
            "source": "mpi4py",
            "error": None,
        },
    )
    stream = _FakeStream(tty=True)

    class FakeSampler:
        def __init__(self):
            self.kwargs = None

        def run(self, **kwargs):
            self.kwargs = kwargs
            return {"status": "ok"}

    sampler = FakeSampler()
    result = run_reactive_sampler(
        sampler,
        verbose=True,
        stream=stream,
    )

    assert result == {"status": "ok"}
    assert sampler.kwargs["show_status"] is False
    assert sampler.kwargs["viz_callback"] is False
    assert stream.getvalue() == ""


def test_run_reactive_sampler_tty_defaults_to_simple_status(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    stream = _FakeStream(tty=True)

    class FakeSampler:
        def __init__(self):
            self.kwargs = None

        def run(self, **kwargs):
            self.kwargs = kwargs
            return {"status": "ok"}

    sampler = FakeSampler()
    result = run_reactive_sampler(
        sampler,
        run_kwargs={"max_ncalls": 1000},
        verbose=True,
        stream=stream,
    )

    assert result == {"status": "ok"}
    assert sampler.kwargs["show_status"] is False
    assert callable(sampler.kwargs["viz_callback"])
    output = stream.getvalue()
    assert "Using simple progress updates" in output
    assert "10s heartbeat" in output


def test_run_reactive_sampler_plain_override_uses_simple_status(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("EXOTIC_ULTRANEST_PLAIN_PROGRESS", "1")
    stream = _FakeStream(tty=True)

    class FakeSampler:
        def __init__(self):
            self.kwargs = None

        def run(self, **kwargs):
            self.kwargs = kwargs
            return {"status": "ok"}

    sampler = FakeSampler()
    result = run_reactive_sampler(
        sampler,
        run_kwargs={"max_ncalls": 1000},
        verbose=True,
        stream=stream,
    )

    assert result == {"status": "ok"}
    assert sampler.kwargs["show_status"] is False
    assert callable(sampler.kwargs["viz_callback"])
    output = stream.getvalue()
    assert "Using simple progress updates" in output
    assert "10s heartbeat" in output


def test_run_reactive_sampler_rich_override_uses_native_status(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("EXOTIC_ULTRANEST_RICH_PROGRESS", "1")
    monkeypatch.setenv("TERM", "xterm-256color")
    stream = _FakeStream(tty=True)

    class FakeSampler:
        def __init__(self):
            self.kwargs = None

        def run(self, **kwargs):
            self.kwargs = kwargs
            return {"status": "ok"}

    sampler = FakeSampler()
    result = run_reactive_sampler(
        sampler,
        run_kwargs={"max_ncalls": 1000},
        verbose=True,
        stream=stream,
    )

    assert result == {"status": "ok"}
    assert sampler.kwargs["show_status"] is True
    assert "viz_callback" not in sampler.kwargs
    assert stream.getvalue() == ""


def test_run_reactive_sampler_silent_mode_mutes_ultranest_logger(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    progress_stream = _FakeStream(tty=False)
    log_stream = io.StringIO()

    class FakeSampler:
        def __init__(self):
            self.kwargs = None
            self.logger = logging.getLogger("ultranest.tests.silent_mode")
            self.logger.handlers = []
            self.logger.propagate = False
            handler = logging.StreamHandler(log_stream)
            handler.setLevel(logging.INFO)
            self.logger.addHandler(handler)
            self.logger.setLevel(logging.INFO)

        def run(self, **kwargs):
            self.kwargs = kwargs
            self.logger.info("native ultranest noise")
            return {"status": "ok"}

    sampler = FakeSampler()
    result = run_reactive_sampler(
        sampler,
        run_kwargs={"max_ncalls": 1000},
        verbose=False,
        stream=progress_stream,
    )

    assert result == {"status": "ok"}
    assert sampler.kwargs["show_status"] is False
    assert sampler.kwargs["viz_callback"] is False
    assert progress_stream.getvalue() == ""
    assert log_stream.getvalue() == ""


def test_run_reactive_sampler_plain_mode_mutes_ultranest_logger(monkeypatch):
    _reset_ultranest_env(monkeypatch)
    monkeypatch.setenv("EXOTIC_ULTRANEST_PLAIN_PROGRESS", "1")
    progress_stream = _FakeStream(tty=False)
    log_stream = io.StringIO()

    class FakeSampler:
        def __init__(self):
            self.kwargs = None
            self.logger = logging.getLogger("ultranest.tests.plain_mode")
            self.logger.handlers = []
            self.logger.propagate = False
            handler = logging.StreamHandler(log_stream)
            handler.setLevel(logging.INFO)
            self.logger.addHandler(handler)
            self.logger.setLevel(logging.INFO)

        def run(self, **kwargs):
            self.kwargs = kwargs
            callback = kwargs["viz_callback"]
            self.logger.info("native ultranest noise")
            callback(None, {"it": 25, "ncall": 40, "logz": -10.0, "logz_remain": -9.5})
            return {"status": "ok"}

    sampler = FakeSampler()
    result = run_reactive_sampler(
        sampler,
        run_kwargs={"max_ncalls": 1000},
        verbose=True,
        stream=progress_stream,
        interval_seconds=0.0,
    )

    assert result == {"status": "ok"}
    assert sampler.kwargs["show_status"] is False
    assert callable(sampler.kwargs["viz_callback"])
    assert "Using simple progress updates" in progress_stream.getvalue()
    assert "native ultranest noise" not in progress_stream.getvalue()
    assert log_stream.getvalue() == ""
