from datetime import datetime

import pytest

from exotic import exotic as exotic_module


def test_format_clock_log_message_uses_hour_and_minute_and_preserves_leading_newlines():
    message = exotic_module.format_clock_log_message(
        "\n\nStarting reduction",
        clock_time=datetime(2026, 7, 16, 7, 5, 42),
    )

    assert message == "\n\n[07:05] Starting reduction"


def test_reduction_stage_timer_logs_step_and_total_elapsed_seconds(monkeypatch):
    clock_values = iter([100.0, 102.5, 109.0])
    logged = []
    monkeypatch.setattr(exotic_module, "log_info", logged.append)

    timer = exotic_module.ReductionStageTimer(time_source=lambda: next(clock_values))

    assert timer.checkpoint("first stage") == pytest.approx(2.5)
    assert timer.checkpoint("second stage") == pytest.approx(6.5)
    assert logged == [
        "STEP TIMING | first stage | elapsed_s=2.50 | total_s=2.50",
        "STEP TIMING | second stage | elapsed_s=6.50 | total_s=9.00",
    ]
