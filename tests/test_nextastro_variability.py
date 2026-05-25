import importlib.util
import sys
import types
import gzip
import json

import numpy as np
import pytest
from tenacity import Future, RetryError

fake_barycorrpy = types.ModuleType('barycorrpy')
fake_utc_tdb = types.ModuleType('barycorrpy.utc_tdb')
fake_utc_tdb.JDUTC_to_BJDTDB = lambda *args, **kwargs: None
fake_barycorrpy.utc_tdb = fake_utc_tdb


def _module_available(name: str) -> bool:
    try:
        return importlib.util.find_spec(name) is not None
    except (ModuleNotFoundError, ValueError):
        return False


def _set_stub_if_missing(name: str, module: types.ModuleType) -> None:
    if not _module_available(name):
        sys.modules.setdefault(name, module)


fake_astroalign = types.ModuleType("astroalign")
fake_astroalign.PIXEL_TOL = 1
fake_astroquery = types.ModuleType("astroquery")
fake_astroquery_simbad = types.ModuleType("astroquery.simbad")
fake_astroquery_simbad.Simbad = type("Simbad", (), {})
fake_astroquery_gaia = types.ModuleType("astroquery.gaia")
fake_astroquery_gaia.Gaia = type("Gaia", (), {})
fake_imreg_dft = types.ModuleType("imreg_dft")
fake_colour_demosaicing = types.ModuleType("colour_demosaicing")
fake_colour_demosaicing.demosaicing_CFA_Bayer_bilinear = lambda *args, **kwargs: None
fake_photutils = types.ModuleType("photutils")
fake_photutils_aperture = types.ModuleType("photutils.aperture")
fake_photutils_aperture.CircularAperture = type("CircularAperture", (), {})
fake_photutils_aperture.CircularAnnulus = type("CircularAnnulus", (), {})
fake_photutils_detection = types.ModuleType("photutils.detection")
fake_photutils_detection.DAOStarFinder = type("DAOStarFinder", (), {})
fake_ldtk = types.ModuleType("ldtk")
fake_ldtk.LDPSet = type("LDPSet", (), {})
fake_ldtk.ldtk = types.SimpleNamespace(LDPSet=fake_ldtk.LDPSet)
fake_ldtk_ldmodel = types.ModuleType("ldtk.ldmodel")
fake_ldtk_ldmodel.LinearModel = type("LinearModel", (), {})
fake_ldtk_ldmodel.QuadraticModel = type("QuadraticModel", (), {})
fake_ldtk_ldmodel.NonlinearModel = type("NonlinearModel", (), {})
fake_lmfit = types.ModuleType("lmfit")
fake_pylightcurve = types.ModuleType("pylightcurve")
fake_pylightcurve_models = types.ModuleType("pylightcurve.models")
fake_pylightcurve_exoplanet = types.ModuleType("pylightcurve.models.exoplanet_lc")
fake_pylightcurve_exoplanet.transit = lambda *args, **kwargs: None
fake_pyvo = types.ModuleType("pyvo")
fake_ultranest = types.ModuleType("ultranest")
fake_ultranest.ReactiveNestedSampler = type("ReactiveNestedSampler", (), {})
fake_elca = types.ModuleType("exotic.api.elca")
fake_elca.lc_fitter = lambda *args, **kwargs: None
fake_elca.binner = lambda *args, **kwargs: None
fake_elca.transit = lambda *args, **kwargs: None
fake_elca.get_phase = lambda *args, **kwargs: None
fake_ld = types.ModuleType("exotic.api.ld")
fake_ld.LimbDarkening = type("LimbDarkening", (), {})
fake_ld.ld_re_punct_p = lambda *args, **kwargs: None

_set_stub_if_missing("astroalign", fake_astroalign)
_set_stub_if_missing("astroquery", fake_astroquery)
_set_stub_if_missing("astroquery.simbad", fake_astroquery_simbad)
_set_stub_if_missing("astroquery.gaia", fake_astroquery_gaia)
_set_stub_if_missing("imreg_dft", fake_imreg_dft)
_set_stub_if_missing("colour_demosaicing", fake_colour_demosaicing)
_set_stub_if_missing("photutils", fake_photutils)
_set_stub_if_missing("photutils.aperture", fake_photutils_aperture)
_set_stub_if_missing("photutils.detection", fake_photutils_detection)
_set_stub_if_missing("ldtk", fake_ldtk)
_set_stub_if_missing("ldtk.ldmodel", fake_ldtk_ldmodel)
_set_stub_if_missing("lmfit", fake_lmfit)
_set_stub_if_missing("pylightcurve", fake_pylightcurve)
_set_stub_if_missing("pylightcurve.models", fake_pylightcurve_models)
_set_stub_if_missing("pylightcurve.models.exoplanet_lc", fake_pylightcurve_exoplanet)
_set_stub_if_missing("pyvo", fake_pyvo)
_set_stub_if_missing("ultranest", fake_ultranest)
sys.modules.setdefault('barycorrpy', fake_barycorrpy)
sys.modules.setdefault('barycorrpy.utc_tdb', fake_utc_tdb)
sys.modules.setdefault("exotic.api.elca", fake_elca)
sys.modules.setdefault("exotic.api.ld", fake_ld)

from exotic import exotic as exotic_module


class DummyResponse:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")


def _decode_request_body(body, headers):
    encoding = headers["Content-Encoding"]
    if encoding == "gzip":
        return json.loads(gzip.decompress(body).decode("utf-8"))
    if encoding == "zstd":
        zstandard = pytest.importorskip("zstandard")
        return json.loads(zstandard.ZstdDecompressor().decompress(body).decode("utf-8"))
    raise AssertionError(f"Unexpected content encoding: {encoding}")


def test_nextastro_variability_logs_json_request_and_response(monkeypatch):
    captured = {}
    logged = []

    def fake_post(url, data, headers, timeout):
        captured['url'] = url
        captured['json'] = _decode_request_body(data, headers)
        captured['headers'] = headers
        captured['timeout'] = timeout
        return DummyResponse([
            {'is_in_vsx': 0},
            {'is_in_vsx': 1},
        ])

    monkeypatch.setattr(exotic_module.requests, 'post', fake_post)
    monkeypatch.setattr(exotic_module, 'log_info', lambda message, warn=False, error=False: logged.append(message))

    variability_flags = exotic_module.nextastro_variability_test([(10.1, -11.2), (22.3, -33.4)])

    assert captured['url'].endswith('/variability_test')
    assert captured['timeout'] == 30
    assert captured['headers']['Content-Type'] == 'application/json'
    assert captured['headers']['Content-Encoding'] in {'gzip', 'zstd'}
    assert captured['json'] == [{'ra': 10.1, 'dec': -11.2}, {'ra': 22.3, 'dec': -33.4}]
    assert variability_flags == [False, True]
    assert any('NextAstro variability request JSON:' in message for message in logged)
    assert any('NextAstro variability request compression:' in message for message in logged)
    assert any('NextAstro variability response JSON:' in message for message in logged)


def test_nextastro_variability_retries_zstd_415_once_with_gzip(monkeypatch):
    logged = []
    encodings = []

    def fake_build_compressed_json_request(payload, content_encoding=None):
        encoding = content_encoding or 'zstd'
        body = json.dumps(payload).encode('utf-8')
        headers = {
            'Content-Type': 'application/json',
            'Content-Encoding': encoding,
        }
        return body, headers, encoding, len(body), len(body)

    def fake_post(url, data, headers, timeout):
        encodings.append(headers['Content-Encoding'])
        if headers['Content-Encoding'] == 'zstd':
            return DummyResponse(None, status_code=415)
        return DummyResponse([{'is_in_vsx': 1}])

    monkeypatch.setattr(exotic_module, 'build_compressed_json_request', fake_build_compressed_json_request)
    monkeypatch.setattr(exotic_module.requests, 'post', fake_post)
    monkeypatch.setattr(exotic_module, 'log_info', lambda message, warn=False, error=False: logged.append(message))

    variability_flags = exotic_module.nextastro_variability_test([(10.1, -11.2)])

    assert variability_flags == [True]
    assert encodings == ['zstd', 'gzip']
    assert any('rejected zstd-compressed request (HTTP 415)' in message for message in logged)


def test_nextastro_variability_caps_retry_attempts_at_five(monkeypatch):
    attempts = []

    def fake_build_compressed_json_request(payload, content_encoding=None):
        encoding = content_encoding or 'gzip'
        body = b'{}'
        headers = {
            'Content-Type': 'application/json',
            'Content-Encoding': encoding,
        }
        return body, headers, encoding, len(body), len(body)

    def fake_post(url, data, headers, timeout):
        attempts.append(headers['Content-Encoding'])
        return DummyResponse(None, status_code=502)

    monkeypatch.setattr(exotic_module, 'build_compressed_json_request', fake_build_compressed_json_request)
    monkeypatch.setattr(exotic_module.requests, 'post', fake_post)
    monkeypatch.setattr(exotic_module.nextastro_variability_test.retry, 'sleep', lambda _: None)

    with pytest.raises(RetryError) as excinfo:
        exotic_module.nextastro_variability_test([(10.1, -11.2)])

    assert len(attempts) == 5
    assert excinfo.value.last_attempt.attempt_number == 5


def test_nextastro_photometry_catalog_match_prefers_requested_filter():
    catalog = {
        'columns': ['id', 'source_id', 'ra', 'dec', 'Vmag', 'err_Vmag', 'g', 'dg'],
        'count': 2,
        'row_format': 'objects',
        'rows': [
            {
                'id': 1,
                'source_id': 111,
                'ra': 10.0001,
                'dec': 20.0001,
                'Vmag': None,
                'err_Vmag': None,
                'g': 12.1,
                'dg': 0.02,
            },
            {
                'id': 2,
                'source_id': 222,
                'ra': 10.0002,
                'dec': 20.0002,
                'Vmag': 12.3,
                'err_Vmag': 0.04,
                'g': 12.0,
                'dg': 0.02,
            },
        ],
    }

    match = exotic_module.nextastro_photometry_catalog_match(catalog, 10.0, 20.0, 'CV')

    assert match['source_id'] == 222
    assert match['mag'] == pytest.approx(12.3)
    assert match['error'] == pytest.approx(0.04)
    assert match['mag_band'] == 'V'
    assert match['separation_arcsec'] > 0


def test_nextastro_photometry_catalog_match_ignores_over_30_magnitudes():
    catalog = {
        'columns': ['id', 'source_id', 'ra', 'dec', 'Vmag', 'err_Vmag'],
        'count': 1,
        'row_format': 'objects',
        'rows': [
            {
                'id': 1,
                'source_id': 111,
                'ra': 10.0001,
                'dec': 20.0001,
                'Vmag': 99.99,
                'err_Vmag': 99.99,
            },
        ],
    }

    match = exotic_module.nextastro_photometry_catalog_match(catalog, 10.0, 20.0, 'CV')

    assert match is None


def test_nextastro_photometry_catalog_match_rejects_separations_over_two_arcsec():
    catalog = {
        'columns': ['id', 'source_id', 'ra', 'dec', 'Vmag', 'err_Vmag'],
        'count': 1,
        'row_format': 'objects',
        'rows': [
            {
                'id': 1,
                'source_id': 111,
                'ra': 10.001,
                'dec': 20.0,
                'Vmag': 12.3,
                'err_Vmag': 0.02,
            },
        ],
    }

    match = exotic_module.nextastro_photometry_catalog_match(
        catalog,
        10.0,
        20.0,
        'CV',
        max_separation_arcsec=30.0,
    )

    assert match is None


def test_merge_nextastro_calibration_stars_adds_non_vsp_metadata():
    catalog = {
        'columns': ['id', 'source_id', 'ra', 'dec', 'Vmag', 'err_Vmag'],
        'count': 1,
        'row_format': 'objects',
        'rows': [
            {
                'id': 9,
                'source_id': 12345,
                'ra': 10.00001,
                'dec': -20.00001,
                'Vmag': 11.2,
                'err_Vmag': 0.03,
            }
        ],
    }

    calibration_stars = exotic_module.merge_nextastro_calibration_stars(
        comp_stars=[[100, 200]],
        comp_ra_dec=[(10.0, -20.0)],
        obs_filter='V',
        existing_comp_stars={},
        field_catalog=catalog,
    )

    assert list(calibration_stars) == ['NextAstro-12345']
    calibration = calibration_stars['NextAstro-12345']
    assert calibration['is_aavso_vsp'] is False
    assert calibration['catalog_source'] == 'NextAstro photometry catalog'
    assert calibration['ra'] == pytest.approx(10.0)
    assert calibration['dec'] == pytest.approx(-20.0)
    assert calibration['mag'] == pytest.approx(11.2)
    assert calibration['error'] == pytest.approx(0.03)
    assert calibration['observed_filter'] == 'V'


def test_build_stellar_variability_params_records_nextastro_reference(monkeypatch, tmp_path):
    captured = {}

    class DummyFit:
        data = np.array([1.0, 1.02, 0.98], dtype=float)
        airmass_model = np.ones(3, dtype=float)
        airmass = np.array([1.1, 1.2, 1.3], dtype=float)
        jd_times = np.array([2450000.1, 2450000.2, 2450000.3], dtype=float)
        transit = np.ones(3, dtype=float)

    def fake_plot(params, save, s_name, label):
        captured['params'] = params
        captured['label'] = label

    monkeypatch.setattr(exotic_module, 'plot_stellar_variability', fake_plot)

    calibration_star = {
        'mag': 12.0,
        'error': 0.05,
        'ra': 10.1,
        'dec': -20.2,
        'catalog_ra': 10.10001,
        'catalog_dec': -20.20001,
        'catalog_source': 'NextAstro photometry catalog',
        'is_aavso_vsp': False,
        'mag_band': 'V',
        'observed_filter': 'V',
        'source_id': 123,
        'separation_arcsec': 0.2,
    }

    params = exotic_module.build_stellar_variability_params_from_fit(
        DummyFit(),
        calibration_star,
        [100, 200],
        'NextAstro-123',
        tmp_path,
        'Host Star',
        observed_filter='CV',
    )

    assert captured['label'] == 'RA=10.1000000 Dec=-20.2000000'
    assert len(params) == 3
    assert params[0]['catalog_source'] == 'NextAstro photometry catalog'
    assert params[0]['is_aavso_vsp'] is False
    assert params[0]['comp_ra'] == pytest.approx(10.1)
    assert params[0]['comp_dec'] == pytest.approx(-20.2)
    assert params[0]['cmag'] == pytest.approx(12.0)
    assert params[0]['cmag_err'] == pytest.approx(0.05)
    assert params[0]['observed_filter'] == 'CV'


def test_check_for_variable_stars_uses_nextastro_flags_to_filter(monkeypatch):
    logged = []

    ra_wcs = np.array([[100.1, 100.2], [100.3, 100.4]])
    dec_wcs = np.array([[-10.1, -10.2], [-10.3, -10.4]])
    comp_stars = [[0, 0], [1, 1]]

    monkeypatch.setattr(exotic_module, 'nextastro_variability_test', lambda payload: [False, True])
    monkeypatch.setattr(exotic_module, 'log_info', lambda message, warn=False, error=False: logged.append(message))

    exotic_module.check_for_variable_stars(
        ra_wcs, dec_wcs, comp_stars, use_nextastro_variability_server=True
    )

    assert comp_stars == [[0, 0]]
    assert any('NextAstro flagged variable: False' in message for message in logged)
    assert any('NextAstro flagged variable: True' in message for message in logged)


def test_check_for_variable_stars_logs_underlying_nextastro_retry_error(monkeypatch):
    logged = []

    ra_wcs = np.array([[100.1]])
    dec_wcs = np.array([[-10.1]])
    comp_stars = [[0, 0]]

    last_attempt = Future(5)
    last_attempt.set_exception(RuntimeError('HTTP 502'))

    def raise_retry_error(payload):
        raise RetryError(last_attempt)

    monkeypatch.setattr(exotic_module, 'nextastro_variability_test', raise_retry_error)
    monkeypatch.setattr(exotic_module, 'query_variable_star_apis', lambda ra, dec: False)
    monkeypatch.setattr(exotic_module, 'log_info', lambda message, warn=False, error=False: logged.append(message))

    exotic_module.check_for_variable_stars(
        ra_wcs, dec_wcs, comp_stars, use_nextastro_variability_server=True
    )

    assert any('RetryError after 5 attempts (RuntimeError: HTTP 502)' in message for message in logged)


def test_get_wcs_falls_back_to_nextastro_when_nova_fails(monkeypatch):
    service_calls = []

    class DummyPlateSolution:
        def __init__(self, **kwargs):
            pass

        def plate_solution(self):
            service_calls.append('nova')
            return False

    class DummyNextAstroSolution:
        def __init__(self, **kwargs):
            pass

        def plate_solution(self):
            service_calls.append('nextastro')
            return 'nextastro-wcs'

    monkeypatch.setattr(exotic_module, 'PlateSolution', DummyPlateSolution)
    monkeypatch.setattr(exotic_module, 'NextAstroPlateSolution', DummyNextAstroSolution)
    monkeypatch.setattr(exotic_module, 'animate_toggle', lambda *args, **kwargs: None)

    solved_wcs = exotic_module.get_wcs('frame.fits', directory='.')

    assert solved_wcs == 'nextastro-wcs'
    assert service_calls == ['nova', 'nextastro']


def test_get_wcs_logs_nextastro_bad_gateway_before_nova_fallback(monkeypatch):
    logged = []
    service_calls = []

    class DummyPlateSolution:
        def __init__(self, **kwargs):
            self.last_error_type = None

        def plate_solution(self):
            service_calls.append('nova')
            return 'nova-wcs'

    class DummyNextAstroSolution:
        def __init__(self, **kwargs):
            self.last_http_status = 502
            self.last_error_type = 'NextAstro solve submission'
            self.api_url = 'https://astrometry.nextastro.org'

        def plate_solution(self):
            service_calls.append('nextastro')
            return False

    monkeypatch.setattr(exotic_module, 'PlateSolution', DummyPlateSolution)
    monkeypatch.setattr(exotic_module, 'NextAstroPlateSolution', DummyNextAstroSolution)
    monkeypatch.setattr(exotic_module, 'animate_toggle', lambda *args, **kwargs: None)
    monkeypatch.setattr(exotic_module, 'log_info', lambda message, warn=False, error=False: logged.append(message))

    solved_wcs = exotic_module.get_wcs('frame.fits', directory='.', use_nextastro_astrometry=True)

    assert solved_wcs == 'nova-wcs'
    assert service_calls == ['nextastro', 'nova']
    assert any(message == 'NextAstro Server not responding. Will try nova.astrometry.net' for message in logged)


def test_get_wcs_logs_nextastro_bad_gateway_after_both_methods_fail(monkeypatch):
    logged = []
    service_calls = []

    class DummyPlateSolution:
        def __init__(self, **kwargs):
            self.last_error_type = 'Upload'

        def plate_solution(self):
            service_calls.append('nova')
            return False

    class DummyNextAstroSolution:
        def __init__(self, **kwargs):
            self.last_http_status = 502
            self.last_error_type = 'NextAstro solve submission'
            self.api_url = 'https://astrometry.nextastro.org'

        def plate_solution(self):
            service_calls.append('nextastro')
            return False

    monkeypatch.setattr(exotic_module, 'PlateSolution', DummyPlateSolution)
    monkeypatch.setattr(exotic_module, 'NextAstroPlateSolution', DummyNextAstroSolution)
    monkeypatch.setattr(exotic_module, 'animate_toggle', lambda *args, **kwargs: None)
    monkeypatch.setattr(exotic_module, 'log_info', lambda message, warn=False, error=False: logged.append(message))

    solved_wcs = exotic_module.get_wcs('frame.fits', directory='.')

    assert solved_wcs is False
    assert service_calls == ['nova', 'nextastro']
    assert any(message == 'NextAstro Server not responding. Both astrometry methods trialed, pushing forward without astrometry solution'
               for message in logged)


def test_vsx_variable_falls_back_to_nextastro(monkeypatch):
    class DummyFailedResponse:
        def raise_for_status(self):
            raise RuntimeError('VSX unavailable')

    monkeypatch.setattr(exotic_module.requests, 'get', lambda *args, **kwargs: DummyFailedResponse())
    monkeypatch.setattr(exotic_module, 'nextastro_variability_test', lambda payload: [True])

    is_variable = exotic_module.vsx_variable(ra=10.0, dec=20.0)

    assert is_variable is True


def test_vsx_variable_handles_empty_default_response_without_fallback(monkeypatch):
    nextastro_called = {'value': False}

    monkeypatch.setattr(
        exotic_module.requests,
        'get',
        lambda *args, **kwargs: DummyResponse({'VSXObjects': []}),
    )
    monkeypatch.setattr(
        exotic_module,
        'nextastro_variability_test',
        lambda payload: nextastro_called.__setitem__('value', True),
    )

    is_variable = exotic_module.vsx_variable(ra=10.0, dec=20.0)

    assert is_variable is False
    assert nextastro_called['value'] is False


def test_vsx_variable_parses_default_vsx_object_list(monkeypatch):
    monkeypatch.setattr(
        exotic_module.requests,
        'get',
        lambda *args, **kwargs: DummyResponse({
            'VSXObjects': {
                'VSXObject': [
                    {
                        'Name': 'alf Ori',
                        'RA2000': '88.79292',
                        'Declination2000': '7.40706',
                        'Category': 'Variable',
                    }
                ]
            }
        }),
    )
    monkeypatch.setattr(exotic_module, 'nextastro_variability_test', lambda payload: [False])

    is_variable = exotic_module.vsx_variable(ra=88.79292, dec=7.40706)

    assert is_variable is True
