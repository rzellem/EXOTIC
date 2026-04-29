import logging
import math
import multiprocessing
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager

import numpy as np


_TRUTHY = {"1", "true", "yes", "on"}
_FALSEY = {"0", "false", "no", "off", "n"}
DEFAULT_PROGRESS_INTERVAL_SECONDS = 10.0
DEFAULT_MIN_NUM_LIVE_POINTS = 200
DEFAULT_RUN_KWARGS = {
    "min_num_live_points": DEFAULT_MIN_NUM_LIVE_POINTS,
    "min_ess": 200,
    "dlogz": 1.0,
    "dKL": 1.0,
    "frac_remain": 0.05,
    "max_num_improvement_loops": 1,
}
MIN_LIVE_POINTS_ENV_KEYS = (
    "EXOTIC_ULTRANEST_MIN_NUM_LIVE_POINTS",
    "EXOTIC_ULTRANEST_MIN_LIVE_POINTS",
)
MPI_SIZE_ENV_KEYS = (
    "OMPI_COMM_WORLD_SIZE",
    "PMI_SIZE",
    "PMIX_SIZE",
    "MV2_COMM_WORLD_SIZE",
)
MPI_RANK_ENV_KEYS = (
    "OMPI_COMM_WORLD_RANK",
    "PMI_RANK",
    "PMIX_RANK",
    "MV2_COMM_WORLD_RANK",
)
ULTRANEST_WORKER_ENV_KEYS = (
    "EXOTIC_ULTRANEST_WORKERS",
    "NEXTASTRO_EXOTIC_ULTRANEST_WORKERS",
)
ULTRANEST_WORKER_BACKEND_ENV = "EXOTIC_ULTRANEST_WORKER_BACKEND"
_PROCESS_LOGLIKE = None


def _is_enabled(value):
    return str(value).strip().lower() in _TRUTHY


def _is_disabled(value):
    return str(value).strip().lower() in _FALSEY


def _coerce_positive_int(value, default=None):
    try:
        parsed = int(float(str(value).strip()))
    except (TypeError, ValueError):
        return default

    if parsed <= 0:
        return default
    return parsed


def _coerce_int(value, default=None):
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return default


def _configured_mpi_int(env_keys, default=None):
    for env_key in env_keys:
        value = _coerce_int(os.environ.get(env_key), default=None)
        if value is not None:
            return value
    return default


def get_mpi_status():
    """Return basic MPI status without requiring MPI to be installed."""
    env_size = _configured_mpi_int(MPI_SIZE_ENV_KEYS, default=1)
    env_rank = _configured_mpi_int(MPI_RANK_ENV_KEYS, default=0)
    try:
        from mpi4py import MPI

        comm = MPI.COMM_WORLD
        return {
            "available": True,
            "size": int(comm.Get_size()),
            "rank": int(comm.Get_rank()),
            "source": "mpi4py",
            "error": None,
        }
    except Exception as exc:
        return {
            "available": False,
            "size": max(int(env_size or 1), 1),
            "rank": max(int(env_rank or 0), 0),
            "source": "environment",
            "error": str(exc),
        }


def is_mpi_worker_process():
    status = get_mpi_status()
    return status["size"] > 1 and status["rank"] > 0


def _is_colab_runtime():
    return bool(
        os.environ.get("COLAB_RELEASE_TAG")
        or os.environ.get("GOOGLE_COLAB")
        or "google.colab" in sys.modules
    )


def _configured_ultranest_workers():
    for env_key in ULTRANEST_WORKER_ENV_KEYS:
        value = os.environ.get(env_key)
        if value in (None, ""):
            continue
        if _is_disabled(value):
            return 1
        parsed = _coerce_positive_int(value, default=None)
        if parsed is not None:
            return parsed

    return 1


def _configured_ultranest_worker_backend():
    backend = str(os.environ.get(ULTRANEST_WORKER_BACKEND_ENV, "")).strip().lower()
    if backend in {"process", "processes", "multiprocessing"}:
        return "process"
    if backend in {"thread", "threads", "threading"}:
        return "thread"
    return "none"


def _process_loglike_chunk(chunk):
    if _PROCESS_LOGLIKE is None:
        raise RuntimeError("UltraNest process worker was not initialized.")
    return _PROCESS_LOGLIKE(chunk)


@contextmanager
def _parallel_vectorized_loglike(sampler, workers=None):
    worker_count = max(int(workers or _configured_ultranest_workers()), 1)
    backend = _configured_ultranest_worker_backend()
    status = get_mpi_status()
    if int(status.get("size") or 1) > 1:
        worker_count = 1

    original_loglike = getattr(sampler, "loglike", None)
    if worker_count <= 1 or backend == "none" or not callable(original_loglike):
        yield 1, "single"
        return

    pool = None
    executor = None
    if backend == "process":
        if not sys.platform.startswith("linux") or _is_colab_runtime():
            yield 1, "single"
            return
        global _PROCESS_LOGLIKE
        _PROCESS_LOGLIKE = original_loglike
        ctx = multiprocessing.get_context("fork")
        pool = ctx.Pool(processes=worker_count)
    elif backend == "thread":
        executor = ThreadPoolExecutor(max_workers=worker_count, thread_name_prefix="exotic-ultranest")
    else:
        yield 1, "single"
        return

    def parallel_loglike(params):
        params_array = np.asarray(params)
        if params_array.ndim != 2 or params_array.shape[0] < 2:
            return original_loglike(params)

        chunk_count = min(worker_count, params_array.shape[0])
        chunks = [chunk for chunk in np.array_split(params_array, chunk_count) if len(chunk)]
        if pool is not None:
            results = pool.map(_process_loglike_chunk, chunks)
        else:
            results = list(executor.map(original_loglike, chunks))
        return np.concatenate([np.atleast_1d(result) for result in results])

    sampler.loglike = parallel_loglike
    try:
        yield worker_count, backend
    finally:
        sampler.loglike = original_loglike
        if pool is not None:
            pool.close()
            pool.join()
            _PROCESS_LOGLIKE = None
        if executor is not None:
            executor.shutdown(wait=True)


def _configured_min_num_live_points(default=DEFAULT_MIN_NUM_LIVE_POINTS):
    for env_key in MIN_LIVE_POINTS_ENV_KEYS:
        value = _coerce_positive_int(os.environ.get(env_key), default=None)
        if value is not None:
            return value
    return default


def _apply_default_run_kwargs(run_kwargs):
    kwargs = dict(DEFAULT_RUN_KWARGS)
    kwargs["min_num_live_points"] = _configured_min_num_live_points()
    kwargs.update({} if run_kwargs is None else dict(run_kwargs))
    return kwargs


def supports_ultranest_live_status(stream=None):
    """Return True only when rich UltraNest status is explicitly enabled."""
    if _is_enabled(os.environ.get("EXOTIC_ULTRANEST_PLAIN_PROGRESS", "")):
        return False
    if not _is_enabled(os.environ.get("EXOTIC_ULTRANEST_RICH_PROGRESS", "")):
        return False

    if stream is None:
        stream = sys.stdout
    if stream is None:
        return False

    isatty = getattr(stream, "isatty", None)
    if not callable(isatty) or not isatty():
        return False

    if _is_enabled(os.environ.get("CI", "")):
        return False

    term = str(os.environ.get("TERM", "")).strip().lower()
    if term == "dumb":
        return False

    return True


def supports_ultranest_simple_status():
    """Return True only when simple UltraNest status is explicitly enabled."""
    return _is_enabled(os.environ.get("EXOTIC_ULTRANEST_PLAIN_PROGRESS", ""))


def _progress_mode(verbose, stream=None):
    if not verbose:
        return "silent"
    if supports_ultranest_live_status(stream=stream):
        return "rich"
    if supports_ultranest_simple_status():
        return "simple"
    return "simple"


@contextmanager
def _mute_ultranest_logging(sampler):
    muted = []
    seen = set()

    sampler_logger = getattr(sampler, "logger", None)
    if isinstance(sampler_logger, logging.Logger):
        muted.append(sampler_logger)
        seen.add(id(sampler_logger))

    for logger in logging.root.manager.loggerDict.values():
        if not isinstance(logger, logging.Logger):
            continue
        if not logger.name.startswith("ultranest"):
            continue
        if id(logger) in seen:
            continue
        muted.append(logger)
        seen.add(id(logger))

    states = [(logger, logger.disabled) for logger in muted]
    try:
        for logger in muted:
            logger.disabled = True
        yield
    finally:
        for logger, disabled in states:
            logger.disabled = disabled


def _read_float(mapping, *keys):
    for key in keys:
        if key not in mapping:
            continue
        try:
            return float(mapping[key])
        except (TypeError, ValueError):
            continue
    return None


def _read_int(mapping, *keys):
    for key in keys:
        if key not in mapping:
            continue
        try:
            return int(mapping[key])
        except (TypeError, ValueError):
            continue
    return None


def _extract_info(args, kwargs):
    info = kwargs.get("info")
    if isinstance(info, dict):
        return info
    for item in reversed(args):
        if isinstance(item, dict):
            return item
    return {}


class _UltraNestSimpleProgress:
    def __init__(self, stream=None, interval_seconds=DEFAULT_PROGRESS_INTERVAL_SECONDS, bar_width=28):
        self.stream = stream if stream is not None else sys.stdout
        self.interval_seconds = max(float(interval_seconds), 0.0)
        self.bar_width = max(int(bar_width), 8)
        self.start_time = time.monotonic()
        self.last_emit = self.start_time - self.interval_seconds
        self.iteration = None
        self.evaluations = None
        self.progress = 0.0

    def _write(self, message):
        if self.stream is None:
            return
        self.stream.write(message + "\n")
        self.stream.flush()

    @staticmethod
    def _progress_from_info(info):
        progress = _read_float(info, "progress", "fraction_done")
        if progress is not None:
            if progress > 1.0:
                progress /= 100.0
            return min(max(progress, 0.0), 1.0)

        remainder = _read_float(info, "remainder_fraction", "remaining_fraction", "frac_remain")
        if remainder is not None:
            if remainder > 1.0:
                remainder /= 100.0
            return min(max(1.0 - remainder, 0.0), 1.0)

        logz = _read_float(info, "logz")
        logz_remain = _read_float(info, "logz_remain", "logzremain")
        if logz is None or logz_remain is None or not math.isfinite(logz_remain):
            return None
        if not math.isfinite(logz):
            return 0.0

        delta = logz_remain - logz
        if delta >= 50:
            return 0.0
        if delta <= -50:
            return 1.0
        return 1.0 / (1.0 + math.exp(delta))

    def _elapsed(self):
        total = max(int(time.monotonic() - self.start_time), 0)
        hours, rem = divmod(total, 3600)
        minutes, seconds = divmod(rem, 60)
        if hours:
            return f"{hours:d}:{minutes:02d}:{seconds:02d}"
        return f"{minutes:02d}:{seconds:02d}"

    def _line(self, done=False):
        fraction = 1.0 if done else min(max(self.progress, 0.0), 1.0)
        filled = int(round(fraction * self.bar_width))
        bar = "#" * filled + "-" * (self.bar_width - filled)
        pct = f"{100.0 * fraction:6.2f}%"
        iteration = "?" if self.iteration is None else str(self.iteration)
        evaluations = "?" if self.evaluations is None else str(self.evaluations)
        state = "done" if done else "running"
        return (
            f"[ultranest] {state} {pct} [{bar}] "
            f"it={iteration} evals={evaluations} elapsed={self._elapsed()}"
        )

    def start(self):
        heartbeat = f"{self.interval_seconds:g}s"
        self._write(f"[ultranest] Using simple progress updates ({heartbeat} heartbeat).")

    def update(self, *args, **kwargs):
        info = _extract_info(args, kwargs)
        if not info:
            return

        progress = self._progress_from_info(info)
        if progress is not None:
            self.progress = progress

        iteration = _read_int(info, "it", "iteration")
        if iteration is not None:
            self.iteration = iteration

        evaluations = _read_int(info, "ncall", "ncalls", "evals")
        if evaluations is not None:
            self.evaluations = evaluations

    def maybe_emit(self, force=False):
        now = time.monotonic()
        if not force and (now - self.last_emit) < self.interval_seconds:
            return
        self.last_emit = now
        self._write(self._line(done=False))

    def finish(self):
        self._write(self._line(done=True))


def run_reactive_sampler(
    sampler,
    run_kwargs=None,
    verbose=True,
    stream=None,
    interval_seconds=DEFAULT_PROGRESS_INTERVAL_SECONDS,
):
    """
    Run an UltraNest sampler with simple text progress by default.

    Set EXOTIC_ULTRANEST_RICH_PROGRESS=1 for UltraNest's native status output.
    Use verbose=False to silence all progress updates.
    """
    kwargs = _apply_default_run_kwargs(run_kwargs)
    mode = "silent" if is_mpi_worker_process() else _progress_mode(verbose=verbose, stream=stream)
    output_stream = stream if stream is not None else sys.stdout

    with _parallel_vectorized_loglike(sampler) as (worker_count, worker_backend):
        if worker_count > 1 and mode != "silent":
            worker_label = "processes" if worker_backend == "process" else "threads"
            print(
                f"[ultranest] Using {worker_count} worker {worker_label} for vectorized likelihood batches.",
                file=output_stream,
                flush=True,
            )

        if mode == "silent":
            kwargs["show_status"] = False
            kwargs["viz_callback"] = False
            with _mute_ultranest_logging(sampler):
                return sampler.run(**kwargs)

        if mode == "rich":
            kwargs.setdefault("show_status", True)
            return sampler.run(**kwargs)

        progress = _UltraNestSimpleProgress(stream=stream, interval_seconds=interval_seconds)
        upstream_callback = kwargs.get("viz_callback")

        def callback(*args, **callback_kwargs):
            progress.update(*args, **callback_kwargs)
            progress.maybe_emit(force=False)
            if callable(upstream_callback):
                upstream_callback(*args, **callback_kwargs)

        kwargs["show_status"] = False
        kwargs["viz_callback"] = callback

        progress.start()
        try:
            with _mute_ultranest_logging(sampler):
                result = sampler.run(**kwargs)
        finally:
            progress.finish()
        return result
