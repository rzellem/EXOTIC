"""Refresh Earth-orientation data without rejecting transit data on age alone.

The existing Astropy age setting triggers a refresh, not a science-data gate.
Use the newest readable IERS table when sources are unavailable or stale, and
report that choice. Astropy's global download and age settings are preserved.
"""
from contextlib import contextmanager
import logging
from threading import RLock

import numpy as np
from astropy.time import Time
from astropy.utils import iers
from astropy.utils.data import download_file, import_file_to_cache, is_url_in_cache


# Same IAU2000A product: IERS and USNO are the publishers; the last source is
# Astropy's maintained package-data snapshot, hosted independently on GitHub.
IERS_ADDITIONAL_SOURCES = (
    'https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all',
    'https://raw.githubusercontent.com/astropy/astropy-iers-data/main/'
    'astropy_iers_data/data/finals2000A.all',
)
_REFRESH_RESULTS = {}
_REPORTED_STALE = set()
_PREPARED_TABLES = {}
_STATE_LOCK = RLock()


def load_barycentric_converter():
    """Bootstrap barycorrpy from local reference data, without a startup download.

    barycorrpy 0.4.4 downloads IERS-B at import even in Astropy's offline mode.
    For an empty cache, supply the official IERS-B snapshot already installed
    with astropy-iers-data. Keep any existing cached table. Calculation-time
    IERS-A refreshes remain the responsibility of earth_orientation_context.
    """
    if not is_url_in_cache(iers.IERS_B_URL):
        # Astropy hashes only the canonical URL into a safe cache directory name,
        # as required by barycorrpy's lookup; no file contents are hashed.
        import_file_to_cache(iers.IERS_B_URL, iers.IERS_B_FILE, replace=False)
        logging.getLogger('exotic.exotic').info(
            'Initializing barycorrpy with bundled IERS-B reference data: %s', iers.IERS_B_FILE)
    # barycorrpy sets auto_max_age=15 on import; preserve the user's setting.
    with iers.conf.set_temp('auto_max_age', iers.conf.auto_max_age):
        from barycorrpy.utc_tdb import JDUTC_to_BJDTDB
    return JDUTC_to_BJDTDB


class _PreparedIERS(iers.IERS_Auto):
    """Interpolate an already selected table without repeating freshness work."""

    def _refresh_table_as_needed(self, mjd):
        # EXOTIC already tried the sources once for this reduction.
        pass

    def _check_interpolate_indices(self, indices_orig, indices_clipped, max_input_mjd):
        # IERS_Auto's check raises only for prediction age. Its actual time
        # interpolation and the coordinate transformations remain unchanged.
        pass


def reset_earth_orientation_refresh():
    """Allow another refresh on the next reduction in the same notebook."""
    with _STATE_LOCK:
        _REFRESH_RESULTS.clear()
        _REPORTED_STALE.clear()
        _PREPARED_TABLES.clear()


def _revision(table):
    return (str(table.meta.get('data_path')), table.meta.get('data_url'),
            table.meta.get('predictive_mjd'), float(table['MJD'][-1].value))


def _needs_refresh(table, latest_mjd, now_mjd, max_age):
    prediction = table.meta.get('predictive_mjd')
    return (prediction is not None and max_age is not None
            and latest_mjd > prediction and now_mjd - prediction > max_age)


def _sources():
    return tuple(dict.fromkeys((iers.conf.iers_auto_url,
                               iers.conf.iers_auto_url_mirror,
                               *IERS_ADDITIONAL_SOURCES)))


def _newest_table(table, latest_mjd, now_mjd, max_age, logger):
    sources = _sources()
    key = (_revision(table), sources)
    if key in _REFRESH_RESULTS:
        return _REFRESH_RESULTS[key]

    best = table
    for url in sources:
        try:
            # Refresh only this small reference table. Astropy's cache maps
            # source URLs to safe filenames; it does not identify science
            # products by content. No science products are scanned or reopened.
            path = download_file(url, cache='update', timeout=iers.conf.remote_timeout)
            candidate = iers.IERS_Auto.read(file=path)
        except Exception as exc:
            logger.warning('IERS refresh failed from %s: %s: %s', url, type(exc).__name__, exc)
            continue
        candidate.meta['data_url'] = url
        prediction = candidate.meta['predictive_mjd']
        logger.info('IERS source %s: predictions start %s; prediction age %.2f days.',
                    url, Time(prediction, format='mjd').utc.isot, now_mjd - prediction)
        if (prediction, float(candidate['MJD'][-1].value)) > (
                best.meta['predictive_mjd'], float(best['MJD'][-1].value)):
            best = candidate
        if not _needs_refresh(best, latest_mjd, now_mjd, max_age):
            break

    # Share the result across frame timestamps, airmass, and reference-epoch
    # conversions. A failed source is not retried for every individual frame.
    _REFRESH_RESULTS[key] = best
    return best


def _select_table(times, logger):
    epochs = times if isinstance(times, Time) else Time(times, format='jd', scale='utc')
    latest_mjd = float(np.max(epochs.utc.mjd, initial=-np.inf))
    now_mjd = float(Time.now().mjd)
    max_age = iers.conf.auto_max_age
    table = iers.earth_orientation_table.get()
    key = (_revision(table), _sources())
    if key in _REFRESH_RESULTS:
        table = _REFRESH_RESULTS[key]
    elif _needs_refresh(table, latest_mjd, now_mjd, max_age) and iers.conf.auto_download:
        table = _newest_table(table, latest_mjd, now_mjd, max_age, logger)
    if _needs_refresh(table, latest_mjd, now_mjd, max_age):
        key = _revision(table)
        if key not in _REPORTED_STALE:
            _REPORTED_STALE.add(key)
            prediction = table.meta['predictive_mjd']
            logger.warning(
                'Continuing with available IERS predictions: start %s, age %.2f days, '
                'refresh age %.2f days, source %s. Table age alone will not abort '
                'the reduction; Earth-orientation prediction uncertainty remains.',
                Time(prediction, format='mjd').utc.isot, now_mjd - prediction, max_age,
                table.meta.get('data_url') or table.meta.get('data_path') or 'bundled table',
            )
    return table


@contextmanager
def earth_orientation_context(times, logger=None):
    """Try all reference sources as needed, then calculate with available data.

    Age alone never aborts a transit reduction. Download and age settings,
    including leap-second update behavior, are left untouched. Explicit offline
    settings are respected. The selected table is scoped to the calculation.
    """
    logger = logger if logger is not None else logging.getLogger('exotic.exotic')
    # Astropy's table is process-wide. Keep overlapping EXOTIC lookups from
    # restoring each other's temporary table prematurely.
    with _STATE_LOCK:
        table = _select_table(times, logger)
        prepared = table
        if isinstance(table, iers.IERS_Auto):
            key = _revision(table)
            if key not in _PREPARED_TABLES:
                _PREPARED_TABLES[key] = _PreparedIERS(table, copy=False)
            prepared = _PREPARED_TABLES[key]
        with iers.earth_orientation_table.set(prepared):
            yield table
