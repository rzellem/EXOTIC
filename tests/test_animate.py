import threading

from exotic import animate


def test_stopping_spinner_wakes_display_wait_and_supports_restart(monkeypatch, capsys):
    entered_wait = threading.Event()

    class DisplayEvent(threading.Event):
        def wait(self, timeout=None):
            entered_wait.set()
            return super().wait(timeout)

    stop = DisplayEvent()
    monkeypatch.setattr(animate, '_animation_stopped', stop)
    monkeypatch.setattr(animate, 'th_animate_exotic', None)
    monkeypatch.setattr(animate, 'done_flag_animate_exotic', None)
    try:
        for _ in range(3):
            entered_wait.clear()
            assert animate.animate_toggle(True) is True
            worker = animate.th_animate_exotic
            assert entered_wait.wait(2)
            assert animate.animate_toggle(False) is False
            assert not worker.is_alive()
            assert animate.th_animate_exotic is None
            assert stop.is_set()
    finally:
        animate.animate_toggle(False)
    assert capsys.readouterr().out.count('Thinking ... DONE!') == 3
