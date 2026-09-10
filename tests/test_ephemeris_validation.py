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


def test_elca_transit_duration_matches_model_for_eccentric_orbits():
    """Regression for the fourth site of issue #1383: elca.transit_duration
    retained the inverted-eccentricity formula after PR #1384 fixed the three
    copies in exotic.py. The analytic duration must match the transit model's
    own first-to-fourth-contact duration for eccentric geometries (it was off
    by 1.95x at e=0.3/omega=90 and 3.48x at e=0.5/omega=90)."""
    import numpy as np
    from exotic.api.elca import transit, transit_duration

    for ecc, omega in [(0.0, 90.0), (0.3, 90.0), (0.3, 270.0), (0.5, 90.0)]:
        values = dict(u0=0.5, u1=0.1, u2=0.3, u3=-0.1, rprs=0.1, per=3.0,
                      ars=8.0, tmid=1.0, ecc=ecc, omega=omega, inc=88.0)
        times = 1.0 + np.linspace(-0.35, 0.35, 300001)
        flux = transit(times, values)
        in_transit = np.flatnonzero(flux < 1 - 1e-9)
        measured = times[in_transit[-1]] - times[in_transit[0]]
        analytic = transit_duration(values)
        assert analytic == pytest.approx(measured, rel=0.01), (
            f"e={ecc} omega={omega}: analytic {analytic*24:.4f} h vs model {measured*24:.4f} h")
