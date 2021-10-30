from exotic.utils import *
from unittest.mock import patch


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
