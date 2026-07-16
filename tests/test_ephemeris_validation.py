import numpy as np
import pytest

from exotic.exotic import resolve_required_transit_ephemeris


def test_required_ephemeris_keeps_valid_initialization_values_without_archive_lookup():
    lookup_called = False

    def archive_lookup():
        nonlocal lookup_called
        lookup_called = True
        return {'pPer': 9.0, 'midT': 2469999.0}

    result = resolve_required_transit_ephemeris(
        {'pName': 'Example b', 'pPer': 2.5, 'midT': 2460000.25},
        archive_lookup=archive_lookup,
    )

    assert result['pPer'] == 2.5
    assert result['midT'] == 2460000.25
    assert lookup_called is False


@pytest.mark.parametrize(
    ('initial_values', 'expected_period', 'expected_tmid'),
    [
        ({'pPer': None, 'midT': 2460000.25}, 2.5, 2460000.25),
        ({'pPer': 0.0, 'midT': 2460000.25}, 2.5, 2460000.25),
        ({'pPer': np.nan, 'midT': 2460000.25}, 2.5, 2460000.25),
        ({'pPer': True, 'midT': 2460000.25}, 2.5, 2460000.25),
        ({'pPer': 2.5, 'midT': None}, 2.5, 2460000.25),
        ({'pPer': 2.5, 'midT': 0.0}, 2.5, 2460000.25),
        ({'pPer': 2.5, 'midT': np.nan}, 2.5, 2460000.25),
    ],
)
def test_required_ephemeris_fills_only_invalid_values_from_archive(
        initial_values, expected_period, expected_tmid):
    result = resolve_required_transit_ephemeris(
        {'pName': 'Example b', **initial_values},
        archive_planet_dict={
            'pPer': 2.5,
            'pPerUnc': 0.001,
            'midT': 2460000.25,
            'midTUnc': 0.002,
        },
    )

    assert result['pPer'] == expected_period
    assert result['midT'] == expected_tmid


def test_required_ephemeris_copies_archive_uncertainty_with_fallback_value():
    result = resolve_required_transit_ephemeris(
        {
            'pName': 'Example b',
            'pPer': None,
            'pPerUnc': None,
            'midT': 2460000.25,
            'midTUnc': 0.005,
        },
        archive_planet_dict={
            'pPer': 2.5,
            'pPerUnc': 0.001,
            'midT': 2461111.0,
            'midTUnc': 0.002,
        },
    )

    assert result['pPer'] == 2.5
    assert result['pPerUnc'] == 0.001
    assert result['midT'] == 2460000.25
    assert result['midTUnc'] == 0.005


def test_required_ephemeris_fails_before_reduction_when_archive_values_are_unusable():
    with pytest.raises(ValueError, match=r"Cannot start EXOTIC reduction.*pPer.*midT"):
        resolve_required_transit_ephemeris(
            {'pName': 'Example b', 'pPer': None, 'midT': 0.0},
            archive_planet_dict={'pPer': np.nan, 'midT': None},
        )


def test_required_ephemeris_reports_archive_lookup_failure():
    def archive_lookup():
        raise RuntimeError('archive unavailable')

    with pytest.raises(ValueError, match=r"NASA Exoplanet Archive fallback failed.*archive unavailable"):
        resolve_required_transit_ephemeris(
            {'pName': 'Example b', 'pPer': None, 'midT': 2460000.25},
            archive_lookup=archive_lookup,
        )
