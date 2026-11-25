import threading
import sys


done_flag_animate_exotic = None
th_animate_exotic = None
_animation_stopped = threading.Event()


def animate():
    import itertools
    for c in itertools.cycle(['|', '/', '-', '\\']):
        if done_flag_animate_exotic:
            sys.stdout.write('\rThinking ... DONE!')
            sys.stdout.write('\n')
            break
        sys.stdout.write('\rThinking ' + c + ' ... \r')
        sys.stdout.flush()
        # Pace the display without delaying shutdown when the work finishes.
        _animation_stopped.wait(0.11)
    sys.stdout.flush()


def animate_toggle(enabled=False):
    """
    Console feature intended to be run in a single synchronous block.
    """
    global done_flag_animate_exotic
    global th_animate_exotic
    if enabled:
        done_flag_animate_exotic = False
        _animation_stopped.clear()
        th_animate_exotic = threading.Thread(target=animate, daemon=True)
        th_animate_exotic.start()
    else:
        done_flag_animate_exotic = True
        _animation_stopped.set()
        if th_animate_exotic is not None:  # kill animate thread before proceeding
            th_animate_exotic.join()
        th_animate_exotic = None
    return not done_flag_animate_exotic
