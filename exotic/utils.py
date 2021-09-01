import logging
import re
import requests
from numpy import floor, log10
from tenacity import retry, retry_if_exception_type, retry_if_result, \
    stop_after_attempt, wait_exponential

try:
    from api.plate_solution import is_false, result_if_max_retry_count
except ImportError:
    from .api.plate_solution import is_false, result_if_max_retry_count

log = logging.getLogger(__name__)


def user_input(prompt, type_, values=None):
    while True:
        try:
            result = type_(input(prompt))
            log.debug(f"{prompt}{result}")
        except ValueError:
            print("Sorry, not a valid datatype.")
            continue
        if type_ == str and values is not None:
            result = result.lower().strip()
            if result not in values:
                print("Sorry, your response was not valid.")
            else:
                return result
        elif type_ == int and values is not None:
            if result not in values:
                print("Sorry, your response was not valid.")
            else:
                return result
        else:
            return result


def init_params(comp, dict1, dict2):
    for key, value in comp.items():
        try:
            if not isinstance(value, tuple):
                dict1[key] = dict2[value]
            else:
                for val in value:
                    try:
                        dict1[key] = dict2[val]
                    except KeyError:
                        pass
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


# Credit: Kalee Tock
def get_val(hdr, ks):
    for key in ks:
        if key in hdr.keys():
            return hdr[key]
        if key.lower() in hdr.keys():
            return hdr[key.lower()]
        new_key = key[0] + key[1:len(key)].lower()  # first letter capitalized
        if new_key in hdr.keys():
            return hdr[new_key]
    return None


# Credit: Kalee Tock
def add_sign(var):
    str_var = str(var)
    m = re.search(r"^[+\-]", str_var)

    if m:
        return str_var
    if float(var) >= 0:
        return f"+{float(var):.6f}"
    else:
        return f"-{float(var):.6f}"


# Credit: Kalee Tock
def process_lat_long(val, key):
    m = re.search(r"\'?([+-]?\d+)[\s:](\d+)[\s:](\d+\.?\d*)", val) or \
        re.search(r"\'?([+-]?\d+)[\s:](\d+\.\d*)", val)
    if m:
        try:
            deg, min, sec = float(m.group(1)), float(m.group(2)), float(m.group(3))
        except IndexError:
            deg, min, sec = float(m.group(1)), float(m.group(2)), 0
        if deg < 0:
            v = deg - (((60 * min) + sec) / 3600)
        else:
            v = deg + (((60 * min) + sec) / 3600)
        return add_sign(v)

    m = re.search("^\'?([+-]?\d+\.\d+)", val)

    if m:
        v = float(m.group(1))
        return add_sign(v)
    else:
        print(f"Cannot match value {val}, which is meant to be {key}.")


# Credit: Kalee Tock
def find(hdr, ks, obs=None):
    # Special stuff for MObs and Boyce-Astro Observatories
    boyce = {"LATITUDE": "+32.6135", "LONGITUD": "-116.3334", "HEIGHT": 1405}
    mobs = {"LATITUDE": "+37.04", "LONGITUD": "-110.73", "HEIGHT": 2606}

    if "OBSERVAT" in hdr.keys() and hdr["OBSERVAT"] == 'Whipple Observatory':
        obs = "MObs"

    #  if "USERID" in hdr.keys() and hdr["USERID"] == 'PatBoyce':
    #    obs = "Boyce"

    if obs == "Boyce":
        boyce_val = get_val(boyce, ks)
        if boyce_val:
            return boyce_val
    if obs == "MObs":
        mobs_val = get_val(mobs, ks)
        if mobs_val:
            return mobs_val

    val = get_val(hdr, ks)

    if ks[0] == "LATITUDE" and val:
        return process_lat_long(str(val), "latitude")
    if ks[0] == "LONGITUD" and val:
        return process_lat_long(str(val), "longitude")

    return val


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
