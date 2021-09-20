from exotic.utils import round_to_2, typecast_check

class TestTypecastCheck:
    """tests the `typecase_check()` function"""

    @staticmethod
    def _returns_four_point_oh(val_to_check):
        assert 4.0 == typecast_check(float, val_to_check)

    def test_checking_for_floats(self):
        # NOTE: there are two usages (as of 2021-09-20) of the `typecast_check`
        # function and both check for floats

        # returns_four_point_oh = lambda val_to_check: assert 4.0 == typecast_check(float, val_to_check)

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

        # not sure if this is a bug?
        result = round_to_2(2.123, 0.000123456789)
        assert 2.123 == result

        # not sure if this is an acceptable edge case? 0.00195 may be desired.
        result = round_to_2(0.0001951234)
        assert 0.0002 == result

