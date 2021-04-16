import logging
import requests
from numpy import floor, log10
from tenacity import retry, retry_if_exception_type, retry_if_result, \
    stop_after_attempt, wait_exponential

try:
    from api.plate_solution import is_false, result_if_max_retry_count
except ImportError:
    from .api.plate_solution import is_false, result_if_max_retry_count


log = logging.getLogger(__name__)


def user_input(prompt, type_, val1=None, val2=None, val3=None):
    while True:
        try:
            result = type_(input(prompt))
            log.debug(f"{prompt}{result}")
        except ValueError:
            log.info("Sorry, not a valid datatype.")
            continue
        if type_ == str and val1 and val2 and val3:
            result = result.lower().strip()
            if result not in (val1, val2, val3):
                log.info("Sorry, your response was not valid.")
            else:
                return result
        elif type_ == int and val1 and val2 and val3:
            if result not in (val1, val2, val3):
                log.info("Sorry, your response was not valid.")
            else:
                return result
        elif type_ == int and val1 and val2 and val3:
            if result not in (val1, val2, val3):
                log.info("Sorry, your response was not valid.")
            else:
                return result
        else:
            return result


def init_params(comp, dict1, dict2):
    for key, value in comp.items():
        try:
            dict1[key] = dict2[value]
        except KeyError:
            pass
    return dict1


def typecast_check(type_, val):
    try:
        return type_(val)
    except (ValueError, TypeError):
        return False


def round_to_2(*args):
    x = args[0]
    if len(args) == 1:
        y = args[0]
    else:
        y = args[1]
    if floor(y) >= 1. or y == 0.0:
        roundval = 2
    else:
        roundval = -int(floor(log10(abs(y)))) + 1
    return round(x, roundval)


def dms_to_dd(dms_in):
    """
    Quick helper method to convert long/lat values in degree-minute-second (dms) form
    (using ':' separators) to decimal (dd) form
    :param dms_in: DMS long/lat value, colon separated
    :return float: Properly signed long/lat value in decimal float form
    """
    if dms_in is None or isinstance(dms_in, str) is False or str(dms_in).count(":") != 2:
        raise ValueError("Invalid DMS input provided for calculations. ...")
    # clean string of errant leading/trailing/internal spaces
    dms = str(dms_in).strip().replace(" ", "")
    degrees, minutes, seconds = dms.split(":")
    dec = abs(float(degrees)) + float(minutes) / 60. + float(seconds) / 3600.
    if float(degrees) < 0.:
        dec = dec * -1.
    return dec


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
       retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
       retry_error_callback=result_if_max_retry_count)
def open_elevation(lat, long):
    query = f"https://api.open-elevation.com/api/v1/lookup?locations={lat},{long}"
    try:
        r = requests.get(query).json()
        return r['results'][0]['elevation']
    except requests.exceptions.RequestException:
        return False
