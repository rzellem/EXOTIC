"""Age affects refresh and diagnostics, never admission of transit data."""
import logging

import numpy as np
import pytest
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers

from exotic import earth_orientation as orientation
from exotic.timing import reference_to_bjd_tdb


@pytest.fixture
def stale_table(monkeypatch):
    table = iers.IERS_Auto.read()
    prediction = table.meta['predictive_mjd']
    now = Time(prediction + 20.62, format='mjd', scale='utc')
    monkeypatch.setattr(Time, 'now', classmethod(lambda cls: now))
    orientation.reset_earth_orientation_refresh()
    with iers.earth_orientation_table.set(table), \
            iers.conf.set_temp('auto_max_age', 15), \
            iers.conf.set_temp('auto_download', True):
        yield table
    orientation.reset_earth_orientation_refresh()


def epochs(table):
    return Time(table.meta['predictive_mjd'] + np.array([8., 19.]), format='mjd', scale='utc')


def test_mirrors_then_nextastro_then_original_publisher():
    with iers.conf.set_temp('iers_auto_url', 'https://publisher.example/iers'), \
            iers.conf.set_temp('iers_auto_url_mirror', 'https://mirror.example/iers'):
        assert orientation._sources() == (
            'https://mirror.example/iers',
            *orientation.IERS_ADDITIONAL_SOURCES,
            'https://downloads.nextastro.org/iers/finals2000A.all',
            orientation.IERS_STANDARD_URL,
            'https://publisher.example/iers',
        )


def test_duplicate_primary_is_still_tried_last():
    with iers.conf.set_temp('iers_auto_url', orientation.IERS_NEXTASTRO_URL), \
            iers.conf.set_temp('iers_auto_url_mirror', orientation.IERS_NEXTASTRO_URL):
        sources = orientation._sources()
    assert sources[-1] == orientation.IERS_NEXTASTRO_URL
    assert len(sources) == len(set(sources))


def test_all_sources_failing_still_converts_and_does_not_retry_per_frame(stale_table, monkeypatch, caplog):
    attempts = []

    def unavailable(url, **kwargs):
        attempts.append(url)
        raise OSError('HTTP 400 from publisher')

    monkeypatch.setattr(orientation, 'download_file', unavailable)
    original = iers.earth_orientation_table.get()
    # Reproduce the reported exception before exercising the repair.
    with monkeypatch.context() as reproduction:
        # Older development Astropy versions refresh even with auto_download
        # disabled. Isolate that unrelated network access from the age check.
        reproduction.setattr(stale_table, '_refresh_table_as_needed', lambda mjd: None)
        with pytest.raises(ValueError, match='predictive values'):
            _ = epochs(stale_table).delta_ut1_utc
    with caplog.at_level(logging.INFO):
        for instant in epochs(stale_table):
            with orientation.earth_orientation_context(instant):
                assert iers.conf.auto_max_age == 15
                assert iers.conf.auto_download is True
                assert np.isfinite(instant.delta_ut1_utc)
    assert len(attempts) >= 3
    assert len(attempts) == len(set(attempts))
    assert 'HTTP 400 from publisher' in caplog.text
    assert caplog.text.count('Continuing with available IERS predictions') == 1
    assert '20.62 days' in caplog.text
    assert iers.conf.auto_max_age == 15
    assert iers.conf.auto_download is True
    assert iers.earth_orientation_table.get() is original


def test_unavailable_and_malformed_sources_do_not_prevent_using_nextastro(
        stale_table, monkeypatch, caplog):
    monkeypatch.setattr(orientation, 'IERS_ADDITIONAL_SOURCES', ('https://independent.example/iers',))
    fresh = stale_table.copy()
    fresh.meta['predictive_mjd'] += 17
    fresh['UT1_UTC'] += .01 * u.s
    attempts = []

    def download(url, **kwargs):
        attempts.append(url)
        if len(attempts) == 1:
            raise OSError('connection reset')
        return url

    def read(cls, file):
        if file == 'https://independent.example/iers':
            raise ValueError('HTML instead of IERS data')
        return fresh

    monkeypatch.setattr(orientation, 'download_file', download)
    monkeypatch.setattr(iers.IERS_Auto, 'read', classmethod(read))
    with caplog.at_level(logging.INFO), orientation.earth_orientation_context(epochs(stale_table)) as selected:
        actual = epochs(stale_table).delta_ut1_utc
        expected = fresh.ut1_utc(epochs(stale_table), return_status=True)[0].to_value(u.s)
        np.testing.assert_array_equal(actual, expected)
        assert selected is fresh
        assert selected.meta['data_url'] == orientation.IERS_NEXTASTRO_URL
    assert len(attempts) == 3
    assert 'HTML instead of IERS data' in caplog.text
    assert '3.62 days' in caplog.text
    assert 'Continuing with available IERS predictions' not in caplog.text


def test_readable_stale_mirrors_are_used_without_replacing_newer_data_with_older(
        stale_table, monkeypatch, caplog):
    attempts = []

    def download(url, **kwargs):
        attempts.append(url)
        return url

    newer = stale_table.copy()
    newer.meta['predictive_mjd'] += 2
    older = stale_table.copy()
    older.meta['predictive_mjd'] -= 2

    def read(cls, file):
        return newer if file == iers.conf.iers_auto_url_mirror else older.copy()

    monkeypatch.setattr(orientation, 'download_file', download)
    monkeypatch.setattr(iers.IERS_Auto, 'read', classmethod(read))
    with caplog.at_level(logging.INFO), orientation.earth_orientation_context(epochs(stale_table)) as selected:
        assert selected is newer
        assert np.all(np.isfinite(epochs(stale_table).delta_ut1_utc))
    assert len(attempts) >= 3  # A readable but old first response still tries other hosts.
    count = len(attempts)
    for dates in (epochs(stale_table), Time('2019-07-16', scale='utc')):
        with orientation.earth_orientation_context(dates) as selected:
            assert selected is newer
            assert np.all(np.isfinite(dates.delta_ut1_utc))
    assert len(attempts) == count
    assert '18.62 days' in caplog.text
    assert 'Continuing with available IERS predictions' in caplog.text


@pytest.mark.parametrize('nextastro_stale', [False, True])
def test_original_publisher_is_used_when_mirrors_fail_or_are_stale(
        stale_table, monkeypatch, nextastro_stale):
    attempts = []
    fresh = stale_table.copy()
    fresh.meta['predictive_mjd'] += 17

    def download(url, **kwargs):
        attempts.append(url)
        if url == orientation.IERS_NEXTASTRO_URL and nextastro_stale:
            return url
        if url != iers.conf.iers_auto_url:
            raise OSError('mirror unavailable')
        return url

    def read(cls, file):
        return stale_table.copy() if file == orientation.IERS_NEXTASTRO_URL else fresh

    monkeypatch.setattr(orientation, 'download_file', download)
    monkeypatch.setattr(iers.IERS_Auto, 'read', classmethod(read))
    with orientation.earth_orientation_context(epochs(stale_table)) as selected:
        assert selected is fresh
        assert selected.meta['data_url'] == iers.conf.iers_auto_url
        assert np.all(np.isfinite(epochs(stale_table).delta_ut1_utc))
    assert attempts == list(orientation._sources())
    assert attempts.index(orientation.IERS_NEXTASTRO_URL) < len(attempts) - 1


@pytest.mark.parametrize('case', ['offline', 'historical', 'fresh'])
def test_no_network_when_offline_or_refresh_is_unnecessary(stale_table, monkeypatch, case):
    monkeypatch.setattr(orientation, 'download_file', lambda *a, **k: pytest.fail('unnecessary download'))
    dates = epochs(stale_table)
    if case == 'offline':
        monkeypatch.setattr(iers.conf, 'auto_download', False)
    elif case == 'historical':
        dates = Time('2019-07-16', scale='utc')
    else:
        stale_table.meta['predictive_mjd'] += 17
    with orientation.earth_orientation_context(dates):
        assert np.all(np.isfinite(dates.delta_ut1_utc))
    assert iers.conf.auto_download is (case != 'offline')
    assert iers.conf.auto_max_age == 15


def test_real_processing_errors_propagate_and_configuration_is_restored(stale_table, monkeypatch):
    monkeypatch.setattr(iers.conf, 'auto_download', False)
    with pytest.raises(ValueError, match='actual calculation failure'):
        with orientation.earth_orientation_context(epochs(stale_table)):
            raise ValueError('actual calculation failure')
    assert iers.conf.auto_max_age == 15
    assert iers.conf.auto_download is False


def test_predictions_far_beyond_the_refresh_age_still_remain_permitted(stale_table, monkeypatch):
    later = Time(stale_table.meta['predictive_mjd'] + 90, format='mjd')
    monkeypatch.setattr(Time, 'now', classmethod(lambda cls: later))
    monkeypatch.setattr(iers.conf, 'auto_download', False)
    with orientation.earth_orientation_context(epochs(stale_table)):
        assert np.all(np.isfinite(epochs(stale_table).delta_ut1_utc))


def test_primary_barycentric_converter_keeps_running_without_switching_calculators(stale_table, monkeypatch):
    from exotic import exotic as reduction

    monkeypatch.setattr(iers.conf, 'auto_download', False)
    calls = []

    def barycentric(dates, **kwargs):
        times = Time(dates, format='jd', scale='utc')
        calls.append(times.delta_ut1_utc)
        return times.tdb.jd, [], 0

    monkeypatch.setattr(reduction, 'JDUTC_to_BJDTDB', barycentric)
    monkeypatch.setattr(reduction, 'SkyCoord', lambda *a, **k: pytest.fail('calculator fallback on age'))
    dates = epochs(stale_table).jd
    actual = reduction.convert_jd_to_bjd(
        dates, {'ra': 291.77, 'dec': 1.38}, {'lat': 31.68, 'long': -110.88, 'elev': 1268.})
    np.testing.assert_array_equal(actual, Time(dates, format='jd', scale='utc').tdb.jd)
    assert len(calls) == 1
    assert np.all(np.isfinite(calls[0]))


def test_calculated_airmass_uses_available_predictions(stale_table, monkeypatch):
    from exotic import exotic as reduction

    monkeypatch.setattr(iers.conf, 'auto_download', False)
    result = reduction.air_mass({}, 291.77, 1.38, 31.68, -110.88, 1268., epochs(stale_table)[0].jd)
    assert np.isfinite(result)


def test_header_airmass_still_needs_no_time_conversion(monkeypatch):
    from exotic import exotic as reduction

    monkeypatch.setattr(reduction, 'earth_orientation_context',
                        lambda *a, **k: pytest.fail('header airmass requires no IERS data'))
    assert reduction.air_mass({'AIRMASS': 1.4}, None, None, None, None, None, None) == 1.4


@pytest.mark.parametrize('standard', ['JD_UTC', 'HJD_UTC'])
def test_reference_light_time_conversion_remains_permitted_with_stale_data(stale_table, monkeypatch, standard):
    monkeypatch.setattr(iers.conf, 'auto_download', False)
    result = reference_to_bjd_tdb(epochs(stale_table)[-1].jd, standard, 214.68, -20.27)
    assert np.isfinite(result)


def test_new_reduction_can_retry_previously_failed_sources(stale_table, monkeypatch):
    attempts = []

    def unavailable(url, **kwargs):
        attempts.append(url)
        raise OSError('offline')

    monkeypatch.setattr(orientation, 'download_file', unavailable)
    with orientation.earth_orientation_context(epochs(stale_table)):
        pass
    count = len(attempts)
    orientation.reset_earth_orientation_refresh()
    with orientation.earth_orientation_context(epochs(stale_table)):
        pass
    assert len(attempts) == 2 * count


def test_one_second_ut1_error_changes_geometric_light_time_by_less_than_two_microseconds(stale_table):
    site = EarthLocation.from_geodetic(lon=0 * u.deg, lat=0 * u.deg, height=0 * u.m)
    target = SkyCoord(ra=np.arange(0, 360, 15) * u.deg, dec=np.zeros(24) * u.deg)
    first = Time('2019-07-16', scale='utc', location=site)
    second = Time(first)
    first.delta_ut1_utc = 0.
    second.delta_ut1_utc = 1.
    differences = (second.light_travel_time(target) - first.light_travel_time(target)).to_value(u.s)
    maximum = np.max(np.abs(differences))
    assert 1.4e-6 < maximum < 1.6e-6
