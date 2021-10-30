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


def user_input(prompt, type_, values=None, max_tries=1000):
    """
    Captures user_input and casts it to the expected type


    Parameters
    ----------
    prompt : str
        A message shown to the user to get a desired answer in the right type
    type_ : type
        The type expected to be captured from the user. The user's response is
        attempted to be cast to this type.
    values : list[type_]
        Acceptable values to receive from the user. If the response from the user
        is valid after the type check BUT the response is not in this list then
        the user will be prompted to try again.
    max_tries : int
        The maximum number of times the user should be prompted to provide valid
        input. Defaults to 1000. Inserted to the function's signature to aid in
        simplicity of tests.

    Returns
    -------
    any
        The user's response cast to the type provided by the `type_` argument to
        the function.
    """

    tries_count = 0

    while True:
        if tries_count >= max_tries:
            print("You have exceeded the maximum number of retries")
            return None

        try:
            result = type_(input(prompt))
            log.debug(f"{prompt}{result}")
        except ValueError:
            tries_count = tries_count + 1
            print("Sorry, not a valid datatype.")
            continue

        if type_ == str and values is not None:
            result = result.lower().strip()
            if result not in values:
                tries_count = tries_count + 1
                print("Sorry, your response was not valid.")
            else:
                return result
        elif type_ == int and values is not None:
            if result not in values:
                tries_count = tries_count + 1
                print("Sorry, your response was not valid.")
            else:
                return result
        else:
            return result


def init_params(comp, dict1, dict2):
    """
    Populates dict1 to be used by the reduction program code


    Uses comp as a source of acceptable keys to populate. Iterates over the keys
    in comp and populates dict1 with values from dict2. The values for each of
    comp's keys can be a string or a tuple. If a comp key has a value that is a
    string, then dict1 is populated by looking up the value of the key for dict2
    directly using comps key's value. If a comp key is a tuple then the tuple
    values are iterated over and dict1 is populated by looking for values in dict2.
    If both values in comp's tuple are found in dict2 then the last value in the
    tuple is populated.

    Examples:

    dict1["foo"] is set to 123 when comp = {"foo": "bar"} and dict2 = {"bar": 123}

    dict1["foo"] is set to 123 when comp = {"foo": ("bar", "baz")} and dict2 has
    a value of 123 where the key is _either_ "bar" or "baz"

    dict1["foo"] is set to 123 when comp = {"foo": ("bar", "baz")} and dict2 has
    this structure: {"bar": 345, "baz": 123}

    Parameters
    ----------
    comp : dict
        Used to map dictionaries used in the reduction program code to human
        readable and sensical input provided by humans.
    dict1 : dict
        Dictionary to be populated and used by the reduction program
    dict2 : dict
        Dictionary provided by other sources like an init file. The keys are more
        sensical for planetary scientists to provide expected values. In
        practice, these values are provided predominantly by an init file ?? or
        an API call for planet_dict ?? FIXME: needs fact checking

    Returns
    -------
    dict
      Populated dict1 with values from dict2
    """

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
    """
    Casts `val` into `type_`

    Parameters
    ----------
    type_ : type
        type to cast val. ex: float
    val : any

    Returns
    -------
    any
        value casted to type_. ex 4.0. Returns False if val cannot be casted.
    """

    try:
        return type_(val)
    except (ValueError, TypeError):
        return False


def round_to_2(*args):
    """
    Rounds a number to the first two non-zero figures after the decimal point


    If is a number is more than or equal to one or less than or equal to negative
    1, the number is rounded to the hundredths place.  If the number is between
    1 and -1 (exlusive) then the number is rounded such that the zeros after the
    decimal and next two non-zero numbers in the decimal are returned.

    Parameters
    ----------
    args : float
        An arbitrary number of numeric args. Expects one or two args. When
        one argument is passed in it rounds according to the docs above. When
        two arguments are passed in, the first number is rounded to either the
        second number's two significant figures' decimal place. Arguments beyond
        two are discarded.

    Returns
    -------
    float
        the original number rounded to two non-zero numbers after the decimal place
    """

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
    """
    Pluck the value for a certain key from myriad possible known keys

    See pull request #882 for good details provided by Kalee Tock. Astronomers
    refer to various pieces of data in non-standard ways. For example, we need
    to use the latitude of the observation to build a reference frame to fit a
    light curve.

    Astronomers use different values to refer to latitude. This function gets the
    desired value by searching through a list of known keys.

    This function can be used to look up the latitude of an observation by
    passing in the headers of the FITS file as the hdr argument for this function
    and passing in ["LATITUDE", "LAT", "SITELAT"] as a list of known values via
    the ks argument.

    Parameters
    ----------
    hdr : dict
        a dictionary of details about the observatory originally embedded in the
        header of the FITS image header.
    ks : list[str]
        a list of known values that astronomers use for a piece of information.

    Returns
    -------
    str
        _first_ match found from the hdr dictionary from the ks list
    """

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
    """
    Adds a + or - to the coordinate if one isn't there already

    Parameters
    ----------
    var : str
        Coordinate, in degrees, of a latitude or longitude

    Returns
    -------
    str
        var as a string if +/- already present. Otherwise it adds a +/- depending
        on the value of var. Returns precision of six digits after the decimal point
        if +/- not already present in `var`
    """
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
    """
    Converts a longitude or latitude into standardized a value

    Parameters
    ----------
    val : str
        either a longitude or latitude coordinate, with a preceding + or -,
        expressed in _either_ HH:MM:SS or degree values. ex: +152.51 or +37:2:24.
    key : str
        expects "longitude" or "latitude"

    Returns
    -------
    str
        longitude or latitude expressed in degree coordinates with a preceding
        + or -. Six digits of precision after the decimal. ex: +152.510000
    """
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
    """
    finds stuff

    Parameters
    ----------
    hdr : dict
        a dictionary of details about the observatory originally embedded in the
        header of the FITS image header.
    ks : list[str]
        a list of known values that astronomers use for a piece of information.
    obs : string
        A specific observatory. Should be one of 'Boyce' or 'MObs' (no quotes).
        Other values are ignored.

    Returns
    -------
    any
        Most often returns a string but can return anything. Designed to return
        the latitude or longitude of an observation as a string.
    """
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
