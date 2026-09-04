from exotic.utils import *
from unittest.mock import patch

import pytest


def test_coerce_boolean_config_value_accepts_all_supported_forms():
    for value in (True, 1, "1", "y", "Y", "yes", "TRUE", "on"):
        assert coerce_boolean_config_value(value) is True

    for value in (False, 0, "0", "n", "N", "no", "FALSE", "off"):
        assert coerce_boolean_config_value(value) is False


def test_coerce_boolean_config_value_rejects_non_boolean_values():
    for value in (None, 2, -1, "sometimes", [], {}):
        assert coerce_boolean_config_value(value) is None


def test_filename_date_token_uses_date_only_for_iso_timestamp():
    assert filename_date_token("2026-05-06T19:51:13.964-0700") == "2026-05-06"
    assert filename_date_token("20260506T195113") == "2026-05-06"
    assert filename_date_token("2026/05/06 19:51:13") == "2026-05-06"


def test_safe_output_filename_sanitizes_filename_chars():
    filename = safe_output_filename(
        "BestFit",
        "XO-1/b ",
        filename_date_token("2026-05-06T19:51:13.964-0700"),
        extension=" png",
    )

    assert filename == "BestFit_XO-1-b_2026-05-06.png"


def test_safe_output_filename_removes_spaces_from_planet_names():
    filename = safe_output_filename(
        "FinalLightCurve",
        "Kepler-12 b",
        "03-JUN-2026",
        extension="png",
    )

    assert filename == "FinalLightCurve_Kepler-12b_03-JUN-2026.png"
    assert " " not in filename


@pytest.mark.parametrize(
    ("planet_name", "expected"),
    (
        ("TOI-4010b", "TOI-4010 b"),
        ("Kepler-11c", "Kepler-11 c"),
        ("HD 41004Ag", "HD 41004A g"),
        ("TOI-4010 b", "TOI-4010 b"),
        ("Candidate", "Candidate"),
    ),
)
def test_format_aavso_exoplanet_name_separates_planet_suffix(planet_name, expected):
    assert format_aavso_exoplanet_name(planet_name) == expected


def test_sanitize_filename_component_cleans_fallback():
    filename = sanitize_filename_component("   ", fallback="bad fallback")

    assert filename == "badfallback"


class TestUserInput:
    """tests the `user_input()` function"""

    @patch("builtins.print", autospec=True)
    @patch("builtins.input", autospec=True)
    def test_max_retries_exceeded(self, mock_input, mock_print):
        """
        added in order to get the `while True` to expire and make the function
        testable
        """

        # NOTE: foo is not in the accepted `values` arg
        user_provided_input = "foo"
        max_retry_count = 4
        mock_input.return_value = user_provided_input

        result = user_input("Enter a y or n",
                            type_=str,
                            values=["y", "n"],
                            max_tries=max_retry_count)

        # NOTE: n+1. n = number or retries allows, +1 to message the user
        # that the max retries have expired
        assert mock_input.call_count == max_retry_count
        assert mock_print.call_count == (max_retry_count + 1)
        assert result is None

    @patch("builtins.print", autospec=True)
    @patch("builtins.input", autospec=True)
    def test_yes_no_responses(self, mock_input, mock_print):
        """NOTE: used in various places in the code. Not arbitrary"""

        user_provided_input = "y"
        mock_input.return_value = user_provided_input
        assert user_provided_input == user_input("Enter a y or n",
                                                 type_=str,
                                                 values=["y", "n"],
                                                 max_tries=1)

        # NOTE: foo is not in the accepted `values` arg
        user_provided_input = "foo"
        mock_input.reset_mock()
        mock_input.return_value = user_provided_input
        result = user_input("Enter a y or n",
                            type_=str,
                            values=["y", "n"],
                            max_tries=1)

        assert result is None

        acceptable_values = ["foo", "bar"]
        mock_input.reset_mock()
        mock_input.return_value = user_provided_input
        result = user_input("More abstractly, provide an acceptable value",
                            type_=str,
                            values=acceptable_values,
                            max_tries=1)

        assert result == user_provided_input

    @patch("builtins.print", autospec=True)
    @patch("builtins.input", autospec=True)
    def test_when_floats_are_expected(self, mock_input, mock_print):
        """
        Test _not_ ints and _not_ strs (special cases in the function under test)
        Commonly used to ask for floats
        """

        # golden path case:
        user_provided_input = 3.14
        mock_input.return_value = user_provided_input
        result = user_input("More floaty now, provide an acceptable value",
                            type_=float,
                            max_tries=1)

        mock_print.assert_not_called()
        assert result == user_provided_input

        user_provided_input = 3  # int
        expected_result = 3.0  # float
        mock_input.reset_mock()
        mock_print.reset_mock()
        mock_input.return_value = user_provided_input
        result = user_input("More maddening, provide something that can be cast",
                            type_=float,
                            max_tries=1)

        mock_print.assert_not_called()
        assert 3 == 3.0  # this passes. That's kind of annoying b/c lhs is an int
        assert type(result) == float
        assert result == expected_result

        # NOTE: raises a value error
        user_provided_input = "foo"
        mock_input.reset_mock()
        mock_print.reset_mock()
        mock_input.return_value = user_provided_input
        result = user_input("More maddening, provide something that can be cast",
                            type_=float,
                            max_tries=1)

        assert 2 == mock_print.call_count  # called once for invalid + max try expiry
        assert result is None

    @patch("builtins.print", autospec=True)
    @patch("builtins.input", autospec=True)
    def test_when_int_is_expected(self, mock_input, mock_print):
        user_provided_input = 123
        mock_input.return_value = user_provided_input
        result = user_input("Give me an int, any int",
                            type_=int,
                            max_tries=1)

        mock_print.assert_not_called()
        assert result == user_provided_input

        mock_input.reset_mock()
        mock_print.reset_mock()
        allowed_values = [123, 234]
        user_provided_input = allowed_values[0]
        mock_input.return_value = user_provided_input

        result = user_input("Give me an int, any int",
                            type_=int,
                            values=allowed_values,
                            max_tries=1)

        mock_print.assert_not_called()
        assert result == user_provided_input

        mock_input.reset_mock()
        mock_print.reset_mock()
        user_provided_input = 456  # not allowed
        mock_input.return_value = user_provided_input
        result = user_input("Give me an int, any int",
                            type_=int,
                            values=allowed_values,
                            max_tries=1)

        assert 2 == mock_print.call_count  # called once for invalid + max try expiry
        assert result is None

        mock_input.reset_mock()
        mock_print.reset_mock()
        user_provided_input = "foo"  # can't be cast to an int
        mock_input.return_value = user_provided_input
        result = user_input("Give me an int, any int",
                            type_=int,
                            max_tries=1)

        assert 2 == mock_print.call_count # called once for invalid + max try expiry
        assert result is None

    @patch("builtins.print", autospec=True)
    @patch("builtins.input", autospec=True)
    def test_when_str_is_expected(self, mock_input, mock_print):

        allowed_values = ["foo", "bar"]

        user_provided_input = allowed_values[0]
        mock_input.return_value = user_provided_input
        result = user_input("Give me a str, any str",
                            type_=str,
                            values=allowed_values,
                            max_tries=1)

        mock_print.assert_not_called()
        assert result == user_provided_input

        mock_input.reset_mock()
        mock_print.reset_mock()
        user_provided_input = "not allowed!"
        mock_input.return_value = user_provided_input
        result = user_input("Give me a str, any str",
                            type_=str,
                            values=allowed_values,
                            max_tries=1)

        assert 2 == mock_print.call_count
        assert result is None

        # # with spaces and weird capitalization
        mock_input.reset_mock()
        mock_print.reset_mock()
        user_provided_input = " FoO     "
        mock_input.return_value = user_provided_input
        result = user_input("Give me a str, any str",
                            type_=str,
                            values=allowed_values,
                            max_tries=1)

        # this is used in the function. Pretty brittle test
        mock_print.assert_not_called()
        assert result == user_provided_input.lower().strip()

        mock_input.reset_mock()
        mock_print.reset_mock()
        user_provided_input = "@llowed!"
        mock_input.return_value = user_provided_input
        result = user_input("Give me a str with non alpha-nums, any str",
                            type_=str,
                            max_tries=1)

        mock_print.assert_not_called()
        assert result == user_provided_input


class TestInitParams:
    """tests the init_params() function"""

    def test_populate_key(self):
        comp = {"foo": "This is used to make the init file make sense"}
        dict1 = {"foo": None}
        dict2 = {"This is used to make the init file make sense": 123}

        result = init_params(comp, dict1, dict2)
        assert type(result) == dict
        assert result.get("foo") == 123

    def test_key_error(self):
        comp = {"foo": "bar"}
        dict1 = {"herp": None}
        dict2 = {"derp": 123}

        result = init_params(comp, dict1, dict2)
        assert result.get("foo") is None
        assert result == dict1

    def test_when_val_in_comp_is_tuple(self):
        # NOTE: Accepts bar or baz as keys in dict2
        comp = {"foo": ("bar", "baz")}
        dict1 = {"foo": None}
        dict2 = {"bar": 123}

        result = init_params(comp, dict1, dict2)
        assert result.get("foo") == 123

        comp = {"foo": ("bar", "baz")}
        dict1 = {"foo": None}
        dict2 = {"baz": 123}

        result = init_params(comp, dict1, dict2)
        assert result.get("foo") == 123

        comp = {"foo": ("bar", "baz")}
        dict1 = {"foo": None}
        dict2 = {"bar": 234, "baz": 123}

        result = init_params(comp, dict1, dict2)
        assert result.get("foo") == 123

        comp = {"foo": ("bar", "baz")}
        dict1 = {"foo": None}
        dict2 = {"herp": 123}

        result = init_params(comp, dict1, dict2)
        assert result.get("foo") is None


class TestTypecastCheck:
    """tests the `typecase_check()` function"""

    @staticmethod
    def _returns_four_point_oh(val_to_check):
        assert 4.0 == typecast_check(float, val_to_check)

    def test_checking_for_floats(self):
        # NOTE: there are two usages (as of 2021-09-20) of the `typecast_check`
        # function and both check for floats

        # floats return floats
        self._returns_four_point_oh(4.0)

        # strings that look like floats return floats
        self._returns_four_point_oh("4.0")

        # ints can be converted to floats
        self._returns_four_point_oh(4)

        # strings that look like ints can be converted to floats
        self._returns_four_point_oh("4")

        # really nutty things like 4x10^0 are okay too
        self._returns_four_point_oh(4e0)

    def test_uncastable_value(self):
        assert typecast_check(float, "foo") is False


class TestRoundToTwo:
    """tests the round_to_2() function"""

    _ARBITRARY_NUMBER = 4.0

    def test_with_one_arg(self):
        # One arg may have been passed in
        result = round_to_2(self._ARBITRARY_NUMBER)
        assert self._ARBITRARY_NUMBER == result

    def test_second_arg_special_zero_case(self):

        # when second argument passed in is 0.0.
        result = round_to_2(self._ARBITRARY_NUMBER, 0.0)
        assert self._ARBITRARY_NUMBER == result

    def test_round_to_two_decimal_places(self):

        float_with_long_fractional_part = 2.34567
        result = round_to_2(float_with_long_fractional_part, 2)
        assert 2.35 == result

        result = round_to_2(2.000123456)
        assert 2.0 == result

        # NOTE: it's kind of weird that the second argument would
        # just be discarded here. Fix it later once I understand
        # the code better
        result = round_to_2(float_with_long_fractional_part, 3)
        assert 2.35 == result

    def test_small_numbers(self):
        # the meat of the function. This gets into testing the
        # -int(floor(log10(abs(y)))) expression

        # for numbers where log10(n) * -1 is 1, round to two places
        result = round_to_2(0.1234567890)
        assert 0.12 == result

        # for numbers that are very small, round to two sig figs
        result = round_to_2(0.000123456789)
        assert 0.00012 == result

        # for numbers that are very small, and negative, round to two sig figs and keep the negativity
        result = round_to_2(-0.000123456789)
        assert -0.00012 == result

        result = round_to_2(2.123, 0.000123456789)
        assert 2.123 == result

        # not sure if this is an acceptable edge case? 0.00195 may be desired.
        result = round_to_2(0.0001951234)
        assert 0.0002 == result


class TestFormatValueAndUncertainty:
    def test_preserves_two_significant_figures_and_matches_value_precision(self):
        assert format_value_and_uncertainty(0.073, 0.01) == ("0.073", "0.010")
        assert format_value_and_uncertainty(1.0, 0.00023) == ("1.00000", "0.00023")
        assert format_value_and_uncertainty(0.0, 0.0031) == ("0.0000", "0.0031")

    def test_formats_uncertainties_above_one_to_two_significant_figures(self):
        assert format_value_and_uncertainty(89.3511, 2.16) == ("89.4", "2.2")
        assert format_value_and_uncertainty(1234, 100) == ("1230", "1.0e+02")

    def test_recomputes_precision_when_rounding_crosses_a_decade(self):
        assert format_value_and_uncertainty(0.0732, 0.00999) == ("0.073", "0.010")

    def test_full_report_text_uses_the_same_precision(self):
        assert format_value_with_uncertainty(12.0, 0.4) == "12.00 +/- 0.40"


class TestGetVal:
    """tests the get_val() function

    NOTE: this could be changed to a private method. The callers are all
    internal to this module
    """

    def test_key_not_lowered(self):
        ks = ["LONGITUD", "LONG", "LONGITUDE", "SITELONG"]
        _expected_value = "a hat"
        hdr = {"LONGITUD": _expected_value}

        assert _expected_value == get_val(hdr, ks)

    def test_lower_key_before_find(self):
        ks = ["LONGITUD", "LONG", "LONGITUDE", "SITELONG"]
        _expected_value = "a hat"
        hdr = {"longitud": _expected_value}

        assert _expected_value == get_val(hdr, ks)

    def test_capitalized_key(self):
        ks = ["LONGITUD", "LONG", "LONGITUDE", "SITELONG"]

        _expected_value = "a hat"
        hdr = {"Longitud": _expected_value}

        assert _expected_value == get_val(hdr, ks)

    def test_key_not_found_at_all(self):
        ks = ["LONGITUD", "LONG", "LONGITUDE", "SITELONG"]
        _expected_value = "a hat"
        hdr = {"foo": _expected_value}

        assert get_val(hdr, ks) is None

    # NOTE: an edge case, but maybe we should make the return
    # more explicit
    def test_key_in_dict_more_than_once(self):
        ks = ["LONGITUD", "LONG", "LONGITUDE", "SITELONG"]

        _expected_value = "a hat"
        hdr = {"LONG": _expected_value,
               "LONG": "foo" }

        assert _expected_value != get_val(hdr, ks)


class TestAddSign:
    """tests the `add_sign()` function

    NOTE: this could be changed to a private method. The callers are all
    internal to this module
    """

    def test_plus_minus_already_present(self):
        input = "+120"
        output = "+120"
        assert output == add_sign(input)

        input = "-120"
        output = "-120"
        assert output == add_sign(input)

    def test_adding_plus_to_coordinate(self):

        assert "+120.000000" == add_sign(120)

        # NOTE: may be a bug here where we want to raise awareness that
        # the coordinate is beyond the coordinate system of planets
        assert "+820.000000" == add_sign(820)

    def test_adding_minus_to_coordinate(self):

        # NOTE: I don't think the else statement which returns a negative
        # coordinate with 6 decimal place precision is ever reached
        assert "-120.000000" != add_sign(-120)
        assert "-120" == add_sign(-120)


class TestProcessLatLong:
    """tests the process_lat_long() function"""

    _ARBITRARY_LONGITUDE = "+152.51"
    _EXPECTED_LONGITUDE_RESULT = "+152.510000"

    _ARBITRARY_LATITUDE = "+37.04"
    _EXPECTED_LATITUDE_RESULT = "+37.040000"

    def test_process_lat_long_degree_inputs(self):
        assert self._EXPECTED_LONGITUDE_RESULT == process_lat_long(self._ARBITRARY_LONGITUDE, "longitude")
        assert self._EXPECTED_LATITUDE_RESULT == process_lat_long(self._ARBITRARY_LATITUDE, "latitude")

    def test_process_lat_long_dms_inputs(self):
        assert self._EXPECTED_LONGITUDE_RESULT == process_lat_long("+152:30:36", "longitude")
        assert self._EXPECTED_LATITUDE_RESULT == process_lat_long("+37:2:24", "latitude")

    @pytest.mark.parametrize(
        ("value", "coordinate_type", "expected"),
        (
            ("28 17 58.8 N", "latitude", 28.2996666667),
            ("28 17 58.8 S", "latitude", -28.2996666667),
            ("16 30 39.7 E", "longitude", 16.5110277778),
            ("16 30 39.7 W", "longitude", -16.5110277778),
            ("-16 30 39.7 W", "longitude", -16.5110277778),
            ("S28:17:58.8", "latitude", -28.2996666667),
        ),
    )
    def test_process_lat_long_hemisphere_inputs(self, value, coordinate_type, expected):
        assert float(process_lat_long(value, coordinate_type)) == pytest.approx(expected)

    def test_process_lat_long_rejects_wrong_hemisphere_for_axis(self):
        assert process_lat_long("28 17 58.8 W", "latitude") is None

    @patch("builtins.print")
    def test_bad_inputs(self, mock_print):
        result = process_lat_long("foo", "longitude")
        self._assert_prints_output_and_returns_none(mock_print, result)

    # NOTE: The following two tests might be bugs. Might want to tighten this up a bit
    def test_when_key_is_not_long_or_lat(self):
        assert self._EXPECTED_LONGITUDE_RESULT == process_lat_long(self._ARBITRARY_LONGITUDE, "HERP")
        assert self._EXPECTED_LATITUDE_RESULT == process_lat_long(self._ARBITRARY_LATITUDE, "DERP")

    @patch("builtins.print")
    def test_process_out_of_range(self, mock_print):

        # When the long and lat are way outside the acceptable values
        assert "+999.000000" == process_lat_long("+999.0", "longitude")
        assert "+999.000000" == process_lat_long("+999.0", "latitude")

        # When a plus or minus sign are missing
        assert "+999.000000" == process_lat_long("999.0", "longitude")

        # When a number without a sign or decimal is passed in
        mock_print.reset_mock()
        result = process_lat_long("999", "longitude")
        self._assert_prints_output_and_returns_none(mock_print, result)

    @staticmethod
    def _assert_prints_output_and_returns_none(mock_print, result):
        mock_print.assert_called()
        assert result is None


class OpenElevationResponse:
    def __init__(self, payload):
        self.payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self.payload


@patch("exotic.utils.requests.get")
def test_open_elevation_uses_encoded_parameters_timeout_and_numeric_result(mock_get):
    mock_get.return_value = OpenElevationResponse({
        "results": [{"latitude": 32.5, "longitude": 151.2, "elevation": 87.0}],
    })

    assert open_elevation("+32.5", "+151.2") == pytest.approx(87.0)
    mock_get.assert_called_once_with(
        OPEN_ELEVATION_URL,
        params={"locations": "32.5,151.2"},
        timeout=OPEN_ELEVATION_TIMEOUT,
    )


@patch("exotic.utils.requests.get")
def test_open_elevation_treats_invalid_response_as_failed_lookup(mock_get):
    mock_get.return_value = OpenElevationResponse({"results": []})

    lookup_without_wait = open_elevation.retry_with(wait=lambda retry_state: 0)

    assert lookup_without_wait(-33.86, 151.21) is False
    assert mock_get.call_count == 3


@patch("exotic.utils.requests.get")
def test_open_elevation_treats_nonfinite_elevation_as_failed_lookup(mock_get):
    mock_get.return_value = OpenElevationResponse({"results": [{"elevation": "nan"}]})

    assert open_elevation.__wrapped__(-33.86, 151.21) is False


class TestFind:
    """tests the find() function"""

    def test_whipple_special_case(self):

        hdr = {"OBSERVAT": "Whipple Observatory",
               "LONG": 4,
               "LAT": 3}

        # these search keys are copied from the implementation code
        search_keys = ['LONGITUD', 'LONG', 'LONGITUDE', 'SITELONG']
        whipple_observatory_longitude = "-110.73"
        result = find(hdr, search_keys)

        assert result == whipple_observatory_longitude

        # these search keys are copied from the implementation code
        search_keys = ['LATITUDE', 'LAT', 'SITELAT']
        whipple_observatory_latitude = "+37.04"
        result = find(hdr, search_keys)

        assert result == whipple_observatory_latitude

        # these search keys are copied from the implementation code
        search_keys = ['HEIGHT', 'ELEVATION', 'ELE', 'EL', 'OBSGEO-H', 'ALT-OBS', 'SITEELEV']
        whipple_observatory_height = 2606
        result = find(hdr, search_keys)

        assert result == whipple_observatory_height

    def test_boyce_observatory(self):
        """This does not appear to used in the implementation code"""

        hdr = {"OBSERVAT": "NOT Whipple Observatory",
               "LONG": "-123.45",
               "LAT": "+34.56"}

        search_keys = ['LONGITUD', 'LONG', 'LONGITUDE', 'SITELONG']
        result = find(hdr, search_keys, obs="Boyce")

        assert result == "-116.3334"  # this value is hard coded in the function

        search_keys = ['LATITUDE', 'LAT', 'SITELAT']
        result = find(hdr, search_keys, obs="Boyce")

        assert result == "+32.6135"  # this value is hard coded in the function

    def test_mobs_observatory(self):
        """This does not appear to used in the implementation code"""

        hdr = {"OBSERVAT": "NOT Whipple Observatory",
               "LONG": "-123.45",
               "LAT": "+34.56"}

        search_keys = ['LONGITUD', 'LONG', 'LONGITUDE', 'SITELONG']
        result = find(hdr, search_keys, obs="MObs")

        assert result == "-110.73"  # this value is hard coded in the function

        search_keys = ['LATITUDE', 'LAT', 'SITELAT']
        result = find(hdr, search_keys, obs="MObs")

        assert result == "+37.04"  # this value is hard coded in the function

    @patch("exotic.utils.process_lat_long")
    def test_generic_hdr(self, mock_pll):
        """This mimics calls by the implementation code"""
        # NOTE: this test is coupled to the implementation of
        # exotic.utils.process_lat_long as the result of that function is used
        # in the return value of this function. That's fine for now but a future
        # improvement could be made to decouple the two functions

        hdr = {"OBSERVAT": "NOT Whipple Observatory",
               "LONG": "-123.45",
               "LAT": "+34.56"}

        # NOTE: the order of these keys matters!
        search_keys = ['LONGITUD', 'LONG', 'LONGITUDE', 'SITELONG']
        mock_pll.return_value = hdr["LONG"]
        result = find(hdr, search_keys)

        mock_pll.assert_called_once()
        assert result == hdr["LONG"]

        mock_pll.reset_mock()
        # NOTE: the order of these keys matters!
        search_keys = ['LATITUDE', 'LAT', 'SITELAT']
        mock_pll.return_value = hdr["LAT"]
        result = find(hdr, search_keys)
        mock_pll.assert_called_once()
        assert result == hdr["LAT"]
        # NOTE: actual return value is "+34.560000" but I mocked this call

    def test_generic_hdr_interprets_coordinate_hemispheres(self):
        hdr = {
            "SITELAT": "28 17 58.8 S",
            "SITELONG": "16 30 39.7 W",
        }

        latitude_result = find(hdr, ['LATITUDE', 'LAT', 'SITELAT'])
        longitude_result = find(hdr, ['LONGITUD', 'LONG', 'LONGITUDE', 'SITELONG'])

        assert float(latitude_result) == pytest.approx(-28.2996666667)
        assert float(longitude_result) == pytest.approx(-16.5110277778)

    @patch("exotic.utils.get_val")
    def test_ks_zero_not_expected(self, mock_get_val):
        # NOTE: returns whatever is returned in `val = get_val()`

        hdr = {"OBSERVAT": "NOT Whipple Observatory",
               "LONG": "-123.45",
               "LAT": "+34.56"}
        # NOTE: changed search key order
        search_keys = ['FOO', 'LONGITUDE', 'SITELONG']
        get_val_returns = 3
        mock_get_val.return_value = get_val_returns

        result = find(hdr, search_keys)

        mock_get_val.assert_called_once()
        assert result == get_val_returns
        assert type(result) == int
