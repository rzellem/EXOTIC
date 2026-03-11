import logging
import math
import os
import sys
import time
from contextlib import contextmanager


_TRUTHY = {"1", "true", "yes", "on"}
DEFAULT_PROGRESS_INTERVAL_SECONDS = 10.0


def _is_enabled(value):
    return str(value).strip().lower() in _TRUTHY


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
        interval_label = f"{self.interval_seconds:g}s"
        self._write(
            "[ultranest] Using simple progress updates "
            f"({interval_label} heartbeat, set EXOTIC_ULTRANEST_RICH_PROGRESS=1 for native status)."
        )

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
    kwargs = {} if run_kwargs is None else dict(run_kwargs)
    mode = _progress_mode(verbose=verbose, stream=stream)

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
