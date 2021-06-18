import threading
import time
import sys


done_flag_animate_exotic = None
th_animate_exotic = None


def animate():
    import itertools
    for c in itertools.cycle(['|', '/', '-', '\\']):
        if done_flag_animate_exotic:
            sys.stdout.write('\rThinking ... DONE!')
            sys.stdout.write('\n')
            break
        sys.stdout.write('\rThinking ' + c + ' ... \r')
        sys.stdout.flush()
        time.sleep(0.11)
    sys.stdout.flush()


def animate_toggle(enabled=False):
    """
    Console feature intended to be run in a single synchronous block.
    """
    global done_flag_animate_exotic
    global th_animate_exotic
    counter_wait = 0
    if enabled:
        done_flag_animate_exotic = False
        th_animate_exotic = threading.Thread(target=animate, daemon=True)
        th_animate_exotic.start()
    else:
        done_flag_animate_exotic = True
        if th_animate_exotic is not None:  # kill animate thread before proceeding
            while th_animate_exotic.is_alive() and counter_wait <= 30:
                time.sleep(.37)
                counter_wait += 1
            th_animate_exotic.join()  # lightest solution
        th_animate_exotic = None
    return not done_flag_animate_exotic
