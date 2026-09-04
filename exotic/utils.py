import logging
from math import isfinite
from pathlib import Path
import re
import requests
from numpy import floor, log10
from tenacity import retry, retry_if_exception_type, retry_if_result, \
    stop_after_attempt, wait_exponential

try:
    from api.plate_solution import is_false
except ImportError:
    from .api.plate_solution import is_false

log = logging.getLogger(__name__)


_WINDOWS_RESERVED_FILENAME_STEMS = {
    'CON',
    'PRN',
    'AUX',
    'NUL',
    *(f'COM{i}' for i in range(1, 10)),
    *(f'LPT{i}' for i in range(1, 10)),
}
_WINDOWS_ILLEGAL_FILENAME_CHARS_RE = re.compile(r'[<>:"/\\|?*\x00-\x1f\x7f]')
_FILENAME_WHITESPACE_RE = re.compile(r'\s+')
_COMPACT_EXOPLANET_SUFFIX_RE = re.compile(r'(?<=[0-9A-Z])([b-z])$')
MAX_APPARENT_MAGNITUDE = 30.0
MAGNITUDE_DECIMAL_PLACES = 4
MINIMUM_MAGNITUDE_ERROR = 0.001
AAVSO_OUTPUT_FOLDER_NAME = 'AAVSO_Files'
BOOLEAN_CONFIG_TRUE_STRINGS = frozenset(('y', 'yes', 'true', '1', 'on'))
BOOLEAN_CONFIG_FALSE_STRINGS = frozenset(('n', 'no', 'false', '0', 'off', ''))
OPEN_ELEVATION_URL = 'https://api.open-elevation.com/api/v1/lookup'
OPEN_ELEVATION_TIMEOUT = 30


def coerce_boolean_config_value(value):
    """Return a configured boolean, or ``None`` when the value is not boolean-like.

    JSON booleans and numeric 1/0 are accepted directly. String values are
    case-insensitive and accept y/n, yes/no, true/false, 1/0, and on/off.
    """

    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        if not isfinite(value):
            return None
        if value == 1:
            return True
        if value == 0:
            return False
        return None
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in BOOLEAN_CONFIG_TRUE_STRINGS:
            return True
        if normalized in BOOLEAN_CONFIG_FALSE_STRINGS:
            return False
    return None


def aavso_output_directory(root):
    """Return the dedicated AAVSO output directory, creating it when needed."""

    output_directory = Path(root) / AAVSO_OUTPUT_FOLDER_NAME
    output_directory.mkdir(parents=True, exist_ok=True)
    return output_directory


def format_aavso_exoplanet_name(value):
    """Separate a compact trailing planet letter for the AAVSO header."""

    name = str(value or '').strip()
    return _COMPACT_EXOPLANET_SUFFIX_RE.sub(r' \1', name)


def _clean_filename_text(value):
    cleaned = _WINDOWS_ILLEGAL_FILENAME_CHARS_RE.sub('-', str(value or ''))
    cleaned = _FILENAME_WHITESPACE_RE.sub('', cleaned)
    return cleaned.rstrip(' .')


def sanitize_filename_component(value, fallback='output'):
    """Return one filename component that is safe on Windows, macOS, and Linux."""

    cleaned = _clean_filename_text(value)
    if cleaned in {'', '.', '..'}:
        cleaned = _clean_filename_text(fallback)
        if cleaned in {'', '.', '..'}:
            cleaned = 'output'
    device_stem = cleaned.split('.', 1)[0].upper()
    if device_stem in _WINDOWS_RESERVED_FILENAME_STEMS:
        cleaned = f'_{cleaned}'
    return cleaned


def filename_date_token(value):
    """Return YYYY-MM-DD when a filename date includes a time component."""

    text = str(value or '').strip()
    match = re.match(r'(\d{4})[-/]?(\d{2})[-/]?(\d{2})', text)
    if match:
        return f'{match.group(1)}-{match.group(2)}-{match.group(3)}'
    return text


def safe_output_filename(prefix, *parts, extension):
    """Build a filename from EXOTIC output labels without illegal path characters."""

    stem_parts = [str(prefix), *(str(part) for part in parts)]
    safe_stem = sanitize_filename_component('_'.join(stem_parts), fallback=str(prefix or 'output'))
    ext = _FILENAME_WHITESPACE_RE.sub('', str(extension or ''))
    if ext and not ext.startswith('.'):
        ext = f'.{ext}'
    return f'{safe_stem}{ext}'


def parse_finite_float(value, default=None):
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return default
    return parsed if isfinite(parsed) else default


def is_usable_apparent_magnitude(value, max_magnitude=MAX_APPARENT_MAGNITUDE):
    parsed = parse_finite_float(value)
    return parsed is not None and parsed <= max_magnitude


def format_magnitude(value, default="na", digits=MAGNITUDE_DECIMAL_PLACES,
                     max_magnitude=MAX_APPARENT_MAGNITUDE):
    parsed = parse_finite_float(value)
    if parsed is None or parsed > max_magnitude:
        return default
    return f"{parsed:.{digits}f}"


def rounded_magnitude_value(value, default=None, digits=MAGNITUDE_DECIMAL_PLACES,
                            max_magnitude=MAX_APPARENT_MAGNITUDE):
    parsed = parse_finite_float(value)
    if parsed is None or parsed > max_magnitude:
        return default
    return round(parsed, digits)


def normalized_magnitude_error(value, default=None, minimum=MINIMUM_MAGNITUDE_ERROR,
                               max_magnitude=MAX_APPARENT_MAGNITUDE):
    parsed = parse_finite_float(value)
    if parsed is None:
        return default
    parsed = abs(parsed)
    if parsed > max_magnitude:
        return default
    return max(parsed, minimum)


def format_magnitude_error(value, default="na", digits=MAGNITUDE_DECIMAL_PLACES,
                           minimum=MINIMUM_MAGNITUDE_ERROR,
                           max_magnitude=MAX_APPARENT_MAGNITUDE):
    parsed = normalized_magnitude_error(
        value,
        default=None,
        minimum=minimum,
        max_magnitude=max_magnitude,
    )
    if parsed is None:
        return default
    return f"{parsed:.{digits}f}"


def rounded_magnitude_error(value, default=None, digits=MAGNITUDE_DECIMAL_PLACES,
                            minimum=MINIMUM_MAGNITUDE_ERROR,
                            max_magnitude=MAX_APPARENT_MAGNITUDE):
    parsed = normalized_magnitude_error(
        value,
        default=None,
        minimum=minimum,
        max_magnitude=max_magnitude,
    )
    if parsed is None:
        return default
    return round(parsed, digits)


def magnitude_text(band, magnitude, magnitude_error=None):
    formatted_mag = format_magnitude(magnitude, default=None)
    if formatted_mag is None:
        return None

    formatted_error = format_magnitude_error(magnitude_error, default=None)
    if formatted_error is None:
        return f"{band}={formatted_mag}"
    return f"{band}={formatted_mag} +/- {formatted_error}"


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


def _two_significant_figure_decimal_places(uncertainty):
    """Return the decimal place needed to show an uncertainty with two sig figs."""

    uncertainty = float(uncertainty)
    if not isfinite(uncertainty) or uncertainty < 0:
        raise ValueError("uncertainty must be a finite, non-negative number")
    if uncertainty == 0:
        return 2

    exponent = int(floor(log10(abs(uncertainty))))
    decimal_places = 1 - exponent

    # A carry can change the exponent (for example, 0.00999 -> 0.010).
    rounded_uncertainty = round(uncertainty, decimal_places)
    if rounded_uncertainty:
        rounded_exponent = int(floor(log10(abs(rounded_uncertainty))))
        decimal_places = 1 - rounded_exponent
    return decimal_places


def _format_at_decimal_place(value, decimal_places):
    """Format a number at a decimal place, including insignificant zeroes."""

    value = float(value)
    if not isfinite(value):
        raise ValueError("value must be a finite number")
    if decimal_places >= 0:
        return f"{value:.{decimal_places}f}"
    return f"{round(value, decimal_places):.0f}"


def format_value_and_uncertainty(value, uncertainty):
    """Return value/error text with a two-significant-figure uncertainty.

    Both strings end at the same decimal place. Unlike ``round_to_2``, this is
    a reporting helper: it deliberately retains trailing zeroes which carry
    precision information.
    """

    decimal_places = _two_significant_figure_decimal_places(uncertainty)
    return (
        _format_at_decimal_place(value, decimal_places),
        format_uncertainty(uncertainty),
    )


def format_uncertainty(uncertainty):
    """Format an uncertainty with exactly two significant figures."""

    decimal_places = _two_significant_figure_decimal_places(uncertainty)
    if decimal_places < 0:
        return f"{float(uncertainty):.1e}"
    return _format_at_decimal_place(uncertainty, decimal_places)


def format_value_with_uncertainty(value, uncertainty):
    """Return ``value +/- uncertainty`` using matched two-sig-fig precision."""

    value_text, uncertainty_text = format_value_and_uncertainty(value, uncertainty)
    return f"{value_text} +/- {uncertainty_text}"


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
        Either a longitude or latitude coordinate expressed in HH:MM:SS or
        decimal degrees. It may use a leading + or - or a FITS-style N/S/E/W
        hemisphere letter. Examples: +152.51, +37:2:24, or 16 30 39.7 W.
    key : str
        expects "longitude" or "latitude"

    Returns
    -------
    str
        longitude or latitude expressed in degree coordinates with a preceding
        + or -. Six digits of precision after the decimal. ex: +152.510000
    """
    text = str(val).strip()
    coordinate_type = str(key).strip().lower()
    valid_hemispheres = {
        "latitude": {"N", "S"},
        "longitude": {"E", "W"},
    }.get(coordinate_type)
    hemisphere = None

    # FITS writers commonly append a hemisphere letter to an otherwise
    # unsigned decimal or sexagesimal coordinate.  A hemisphere overrides a
    # redundant leading sign so that ``-16 30 W`` is not double-negated.
    trailing_hemisphere = re.search(r"([NSEW])\s*$", text)
    leading_hemisphere = re.match(r"\s*([NSEW])(?=\s|[+-]?\d)", text)
    hemisphere_match = trailing_hemisphere or leading_hemisphere
    if hemisphere_match:
        hemisphere = hemisphere_match.group(1)
        if valid_hemispheres is not None and hemisphere not in valid_hemispheres:
            print(f"Cannot match value {val}, which is meant to be {key}.")
            return None
        start, end = hemisphere_match.span(1)
        text = f"{text[:start]}{text[end:]}".strip()

    number_tokens = re.findall(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)", text)
    if not 1 <= len(number_tokens) <= 3:
        print(f"Cannot match value {val}, which is meant to be {key}.")
        return None
    if len(number_tokens) == 1 and hemisphere is None \
            and "." not in number_tokens[0] and number_tokens[0][0] not in "+-":
        # Preserve the historical rejection of an unsigned integer while
        # accepting one when a hemisphere supplies the otherwise missing sign.
        print(f"Cannot match value {val}, which is meant to be {key}.")
        return None

    degrees = float(number_tokens[0])
    minutes = abs(float(number_tokens[1])) if len(number_tokens) >= 2 else 0.0
    seconds = abs(float(number_tokens[2])) if len(number_tokens) >= 3 else 0.0
    magnitude = abs(degrees) + minutes / 60.0 + seconds / 3600.0
    if hemisphere:
        sign = -1.0 if hemisphere in {"S", "W"} else 1.0
    else:
        sign = -1.0 if number_tokens[0].startswith("-") else 1.0
    return add_sign(sign * magnitude)


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
    # MicroObservatory telescopes sit at the Whipple Observatory base camp
    # (Amado, AZ), not on the Mount Hopkins summit (2606 m). See PR #1382.
    mobs = {"LATITUDE": "+31.675467", "LONGITUD": "-110.951376", "HEIGHT": 1268}

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


def _return_false_after_retries(retry_state):
    return False


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10),
       retry=(retry_if_result(is_false) | retry_if_exception_type(requests.exceptions.RequestException)),
       retry_error_callback=_return_false_after_retries)
def open_elevation(lat, long):
    try:
        latitude = float(lat)
        longitude = float(long)
        if not isfinite(latitude) or not isfinite(longitude):
            return False
        if not -90.0 <= latitude <= 90.0 or not -180.0 <= longitude <= 180.0:
            return False

        response = requests.get(
            OPEN_ELEVATION_URL,
            params={'locations': f'{latitude},{longitude}'},
            timeout=OPEN_ELEVATION_TIMEOUT,
        )
        response.raise_for_status()
        result = response.json()['results'][0]['elevation']
        result = float(result)
        return result if isfinite(result) else False
    except (requests.exceptions.RequestException, KeyError, IndexError, TypeError, ValueError) as exc:
        log.debug("Open-Elevation lookup failed: %s", exc)
        return False
