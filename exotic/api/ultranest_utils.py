import gc
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
_AUTO_WORKER_VALUES = {"auto", "all", "available", "cpu", "cpus", "core", "cores"}
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
BYTES_PER_GIB = 1024 ** 3
MIN_AUTO_POINTS_PER_WORKER = 8
MEDIUM_AUTO_POINTS_PER_WORKER = 16
HIGH_AUTO_POINTS_PER_WORKER = 24
MAX_AUTO_POINTS_PER_WORKER = 32
AUTO_POINTS_PER_WORKER_MULTIPLIER = 2
AUTO_DRAW_RAM_FRACTION = 0.005
MIN_AUTO_DRAW_RAM_BUDGET_BYTES = 64 * 1024 ** 2
_PROCESS_LOGLIKE = None
_TK_CLEANUP_CLASSES = ("Image", "Variable")


def _is_enabled(value):
    return str(value).strip().lower() in _TRUTHY


def _is_disabled(value):
    return str(value).strip().lower() in _FALSEY


def _is_auto_worker_count(value):
    return str(value).strip().lower() in _AUTO_WORKER_VALUES


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


def _available_cpu_count():
    process_cpu_count = getattr(os, "process_cpu_count", None)
    if callable(process_cpu_count):
        count = process_cpu_count()
    else:
        count = os.cpu_count()
    return max(_coerce_positive_int(count, default=1), 1)


def _configured_ultranest_workers():
    for env_key in ULTRANEST_WORKER_ENV_KEYS:
        value = os.environ.get(env_key)
        if value in (None, ""):
            continue
        if _is_disabled(value):
            return 1
        if _is_auto_worker_count(value):
            return _available_cpu_count()
        parsed = _coerce_positive_int(value, default=None)
        if parsed is not None:
            return parsed

    return _available_cpu_count()


def _configured_ultranest_worker_backend():
    backend = str(os.environ.get(ULTRANEST_WORKER_BACKEND_ENV, "")).strip().lower()
    if backend in {"process", "processes", "multiprocessing"}:
        return "process"
    if backend in {"thread", "threads", "threading"}:
        return "thread"
    return "none"


def _system_total_memory_bytes():
    if hasattr(os, "sysconf"):
        try:
            page_size = int(os.sysconf("SC_PAGE_SIZE"))
            page_count = int(os.sysconf("SC_PHYS_PAGES"))
            if page_size > 0 and page_count > 0:
                return page_size * page_count
        except (AttributeError, OSError, TypeError, ValueError):
            pass

    if sys.platform.startswith("win"):
        try:
            import ctypes

            class MEMORYSTATUSEX(ctypes.Structure):
                _fields_ = [
                    ("dwLength", ctypes.c_ulong),
                    ("dwMemoryLoad", ctypes.c_ulong),
                    ("ullTotalPhys", ctypes.c_ulonglong),
                    ("ullAvailPhys", ctypes.c_ulonglong),
                    ("ullTotalPageFile", ctypes.c_ulonglong),
                    ("ullAvailPageFile", ctypes.c_ulonglong),
                    ("ullTotalVirtual", ctypes.c_ulonglong),
                    ("ullAvailVirtual", ctypes.c_ulonglong),
                    ("sullAvailExtendedVirtual", ctypes.c_ulonglong),
                ]

            memory_status = MEMORYSTATUSEX()
            memory_status.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(memory_status)):
                return int(memory_status.ullTotalPhys)
        except Exception:
            pass

    return None


def _auto_points_per_worker(cpu_count, total_memory_bytes):
    if total_memory_bytes is None:
        return MEDIUM_AUTO_POINTS_PER_WORKER * AUTO_POINTS_PER_WORKER_MULTIPLIER

    ram_per_cpu_gib = total_memory_bytes / max(int(cpu_count or 1), 1) / BYTES_PER_GIB
    if ram_per_cpu_gib >= 2.0:
        return MAX_AUTO_POINTS_PER_WORKER * AUTO_POINTS_PER_WORKER_MULTIPLIER
    if ram_per_cpu_gib >= 1.0:
        return HIGH_AUTO_POINTS_PER_WORKER * AUTO_POINTS_PER_WORKER_MULTIPLIER
    if ram_per_cpu_gib >= 0.5:
        return MEDIUM_AUTO_POINTS_PER_WORKER * AUTO_POINTS_PER_WORKER_MULTIPLIER
    return MIN_AUTO_POINTS_PER_WORKER * AUTO_POINTS_PER_WORKER_MULTIPLIER


def _estimate_ultranest_draw_point_bytes(sampler):
    x_dim = _coerce_positive_int(getattr(sampler, "x_dim", None), default=None)
    num_params = _coerce_positive_int(getattr(sampler, "num_params", None), default=None)
    coordinate_count = 16
    if x_dim is not None and num_params is not None:
        coordinate_count = max(16, 3 + x_dim + num_params)

    return max(1024, coordinate_count * np.dtype(float).itemsize * 16)


def _apply_auto_ultranest_draw_sizes(sampler, worker_count):
    worker_count = max(int(worker_count or 1), 1)
    if worker_count <= 1 or getattr(sampler, "draw_multiple", True) is False:
        return None

    current_min = _coerce_positive_int(getattr(sampler, "ndraw_min", None), default=None)
    current_max = _coerce_positive_int(getattr(sampler, "ndraw_max", None), default=None)
    if current_min is None and current_max is None:
        return None
    if current_min is None:
        current_min = 128
    if current_max is None:
        current_max = max(current_min, 65536)
    if current_max < current_min:
        current_max = current_min

    total_memory_bytes = _system_total_memory_bytes()
    cpu_count = _available_cpu_count()
    points_per_worker = _auto_points_per_worker(cpu_count, total_memory_bytes)
    target_min = max(current_min, worker_count * points_per_worker)
    memory_limited_max = current_max

    if total_memory_bytes is not None:
        memory_budget = max(
            int(total_memory_bytes * AUTO_DRAW_RAM_FRACTION),
            MIN_AUTO_DRAW_RAM_BUDGET_BYTES,
        )
        point_bytes = _estimate_ultranest_draw_point_bytes(sampler)
        memory_limited_max = max(current_min, min(current_max, memory_budget // point_bytes))
        target_min = min(target_min, memory_limited_max)

    target_min = int(max(current_min, min(current_max, target_min)))
    target_max = int(max(target_min, min(current_max, memory_limited_max)))
    if target_min == current_min and target_max == current_max:
        return None

    sampler.ndraw_min = target_min
    sampler.ndraw_max = target_max

    ram_gib = None
    if total_memory_bytes is not None:
        ram_gib = total_memory_bytes / BYTES_PER_GIB

    return {
        "ndraw_min": target_min,
        "ndraw_max": target_max,
        "workers": worker_count,
        "cpu_count": cpu_count,
        "ram_gib": ram_gib,
        "points_per_worker": points_per_worker,
    }


def _process_loglike_chunk(chunk):
    if _PROCESS_LOGLIKE is None:
        raise RuntimeError("UltraNest process worker was not initialized.")
    return _PROCESS_LOGLIKE(chunk)


def _noop_tk_destructor(_instance):
    return None


def _suppress_inherited_tk_cleanup():
    """Avoid noisy Tk destructor calls in forked worker processes."""
    tkinter_module = sys.modules.get("tkinter")
    if tkinter_module is None:
        return False

    patched = False
    for class_name in _TK_CLEANUP_CLASSES:
        tk_class = getattr(tkinter_module, class_name, None)
        if tk_class is None or getattr(tk_class, "_exotic_worker_tk_cleanup_suppressed", False):
            continue

        try:
            original_del = getattr(tk_class, "__del__", None)
            if original_del is None:
                continue
            setattr(tk_class, "_exotic_worker_original_del", original_del)
            setattr(tk_class, "__del__", _noop_tk_destructor)
            setattr(tk_class, "_exotic_worker_tk_cleanup_suppressed", True)
            patched = True
        except Exception:
            continue

    return patched


def suppress_inherited_tk_cleanup_in_worker():
    return _suppress_inherited_tk_cleanup()


@contextmanager
def suppress_tk_cleanup_during_process_pool():
    suppress_inherited_tk_cleanup_in_worker()
    restore_gc_after_pool = gc.isenabled()
    if restore_gc_after_pool:
        gc.disable()
    try:
        yield
    finally:
        if restore_gc_after_pool:
            gc.enable()


@contextmanager
def _parallel_vectorized_loglike(sampler, workers=None):
    worker_count = max(int(workers or _configured_ultranest_workers()), 1)
    backend = _configured_ultranest_worker_backend()
    status = get_mpi_status()
    if int(status.get("size") or 1) > 1:
        worker_count = 1

    original_loglike = getattr(sampler, "loglike", None)
    if worker_count <= 1 or backend == "none" or not callable(original_loglike):
        yield 1, "single", None
        return

    pool = None
    executor = None
    process_pool_guard = None
    if backend == "process":
        if not sys.platform.startswith("linux") or _is_colab_runtime():
            yield 1, "single", None
            return
        global _PROCESS_LOGLIKE
        process_pool_guard = suppress_tk_cleanup_during_process_pool()
        process_pool_guard.__enter__()
        _PROCESS_LOGLIKE = original_loglike
        try:
            ctx = multiprocessing.get_context("fork")
            pool = ctx.Pool(processes=worker_count, initializer=suppress_inherited_tk_cleanup_in_worker)
        except Exception:
            _PROCESS_LOGLIKE = None
            process_pool_guard.__exit__(*sys.exc_info())
            raise
    elif backend == "thread":
        executor = ThreadPoolExecutor(max_workers=worker_count, thread_name_prefix="exotic-ultranest")
    else:
        yield 1, "single", None
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
    draw_sizes = _apply_auto_ultranest_draw_sizes(sampler, worker_count)
    try:
        yield worker_count, backend, draw_sizes
    finally:
        sampler.loglike = original_loglike
        if pool is not None:
            try:
                pool.close()
                pool.join()
            finally:
                _PROCESS_LOGLIKE = None
                if process_pool_guard is not None:
                    process_pool_guard.__exit__(None, None, None)
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


def _traceback_mentions_ultranest_mlfriends(exc):
    traceback = exc.__traceback__
    while traceback is not None:
        filename = str(traceback.tb_frame.f_code.co_filename).replace("\\", "/")
        if "ultranest/mlfriends" in filename:
            return True
        traceback = traceback.tb_next
    return False


def _is_ultranest_degenerate_region_error(exc):
    if not isinstance(exc, ValueError):
        return False

    message = str(exc)
    if "Buffer has wrong number of dimensions" not in message:
        return False
    if "expected 2" not in message or "got 0" not in message:
        return False
    return _traceback_mentions_ultranest_mlfriends(exc)


def _run_sampler_with_degenerate_region_guard(sampler, kwargs):
    try:
        return sampler.run(**kwargs)
    except ValueError as exc:
        if not _is_ultranest_degenerate_region_error(exc):
            raise
        raise np.linalg.LinAlgError(
            "UltraNest failed while building a degenerate sampling region."
        ) from exc


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

    with _parallel_vectorized_loglike(sampler) as (worker_count, worker_backend, draw_sizes):
        if draw_sizes is not None and mode != "silent":
            memory_label = "unknown RAM"
            if draw_sizes["ram_gib"] is not None:
                memory_label = f"{draw_sizes['ram_gib']:.1f} GiB RAM"
            print(
                "[ultranest] Auto proposal draw sizes: "
                f"ndraw_min={draw_sizes['ndraw_min']}, ndraw_max={draw_sizes['ndraw_max']} "
                f"({draw_sizes['workers']} workers, {draw_sizes['points_per_worker']} points/worker, "
                f"{memory_label}).",
                file=output_stream,
                flush=True,
            )
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
                return _run_sampler_with_degenerate_region_guard(sampler, kwargs)

        if mode == "rich":
            kwargs.setdefault("show_status", True)
            return _run_sampler_with_degenerate_region_guard(sampler, kwargs)

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
                result = _run_sampler_with_degenerate_region_guard(sampler, kwargs)
        finally:
            progress.finish()
        return result
