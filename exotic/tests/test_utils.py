from exotic.utils import *
from unittest.mock import patch


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
