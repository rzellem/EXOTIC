import importlib.util
import sys
import types
import gzip
import json

import numpy as np
import pytest
from tenacity import Future, RetryError

from exotic.plate_status import PlateStatus

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


def test_nextastro_vsx_query_boxes_split_ra_wrap():
    boxes = exotic_module.nextastro_vsx_query_boxes(359.9, 0.0, 0.2)

    np.testing.assert_allclose(boxes, [
        (359.7, 360.0, -0.2, 0.2),
        (0.0, 0.1, -0.2, 0.2),
    ])


def test_nextastro_vsx_field_query_normalizes_rows(monkeypatch):
    captured = {}

    def fake_post(url, json, timeout):
        captured['url'] = url
        captured['json'] = json
        captured['timeout'] = timeout
        return DummyResponse({
            'columns': ['oid', 'name', 'ra_deg', 'dec_deg', 'var_type', 'mag1', 'mag1_band'],
            'count': 1,
            'row_format': 'objects',
            'rows': [{
                'oid': 123,
                'name': 'Cached Variable',
                'ra_deg': 10.1,
                'dec_deg': -20.2,
                'var_type': 'EA',
                'mag1': 12.3,
                'mag1_band': 'V',
            }],
        })

    monkeypatch.setattr(exotic_module.requests, 'post', fake_post)

    rows = exotic_module.nextastro_vsx_field_query(10.0, -20.0, 0.25)

    assert captured['url'].endswith('/vsx_query')
    assert captured['timeout'] == 30
    assert captured['json']['compact'] is False
    assert rows[0]['Name'] == 'Cached Variable'
    assert rows[0]['OID'] == 123
    assert rows[0]['RA2000'] == pytest.approx(10.1)
    assert rows[0]['Declination2000'] == pytest.approx(-20.2)
    assert rows[0]['VariabilityType'] == 'EA'
    assert rows[0]['MaxMag'] == '12.3 V'
    assert rows[0]['_vsx_source'] == 'nextastro_cache'
    assert rows[0]['_vsx_has_full_metadata'] is False


def test_vsx_field_query_cache_first_falls_back_when_cache_empty(monkeypatch):
    calls = []
    monkeypatch.setattr(exotic_module, 'nextastro_vsx_field_query', lambda *args: [])
    monkeypatch.setattr(
        exotic_module,
        'vsx_field_query',
        lambda *args, **kwargs: calls.append((args, kwargs)) or [{'Name': 'AAVSO Variable'}],
    )

    rows = exotic_module.vsx_field_query_with_preference(
        10.0,
        -20.0,
        0.25,
        use_nextastro_vsx_cache_first=True,
    )

    assert rows == [{'Name': 'AAVSO Variable'}]
    assert len(calls) == 1


def test_vsx_field_query_cache_first_enriches_period_and_amplitude(monkeypatch):
    monkeypatch.setattr(
        exotic_module,
        'nextastro_vsx_field_query',
        lambda *args: [{
            'OID': 123,
            'Name': 'Cached Name',
            'RA2000': 10.1,
            'Declination2000': -20.2,
            'VariabilityType': 'EA',
            '_vsx_source': 'nextastro_cache',
        }],
    )
    monkeypatch.setattr(
        exotic_module,
        'vsx_field_query',
        lambda *args, **kwargs: [{
            'OID': '123',
            'Name': 'AAVSO Name',
            'RA2000': '10.1000',
            'Declination2000': '-20.2000',
            'Period': '2.5',
            'MaxMag': '12.0 V',
            'MinMag': '12.4 V',
        }],
    )

    rows = exotic_module.vsx_field_query_with_preference(
        10.0,
        -20.0,
        0.25,
        use_nextastro_vsx_cache_first=True,
    )

    assert rows[0]['Name'] == 'AAVSO Name'
    assert rows[0]['Period'] == '2.5'
    assert exotic_module.vsx_object_amplitude_mag(rows[0]) == pytest.approx(0.4)
    assert rows[0]['_vsx_source'] == 'nextastro_cache+aavso_metadata'


def test_vsx_field_query_cache_first_uses_full_nextastro_metadata_without_aavso(monkeypatch):
    cached_row = exotic_module.normalize_nextastro_vsx_row({
        'oid': 123,
        'name': 'Cached Full Variable',
        'ra_deg': 10.1,
        'dec_deg': -20.2,
        'var_type': 'EA',
        'period_days': 2.5,
        'amplitude_mag': 0.4,
        'max_mag': 12.0,
        'max_passband': 'V',
        'min_mag': 12.4,
        'min_passband': 'V',
    })
    monkeypatch.setattr(
        exotic_module,
        'nextastro_vsx_field_query',
        lambda *args: [cached_row],
    )

    def unexpected_aavso_call(*args, **kwargs):
        raise AssertionError('AAVSO should not be called for the full NextAstro schema')

    monkeypatch.setattr(exotic_module, 'vsx_field_query', unexpected_aavso_call)

    rows = exotic_module.vsx_field_query_with_preference(
        10.0,
        -20.0,
        0.25,
        use_nextastro_vsx_cache_first=True,
    )

    assert rows[0]['_vsx_has_full_metadata'] is True
    assert rows[0]['Period'] == 2.5
    assert rows[0]['Amplitude'] == 0.4
    assert rows[0]['MaxMag'] == '12.0 V'
    assert rows[0]['MinMag'] == '12.4 V'


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


def test_nextastro_photometry_catalog_match_floors_zero_magnitude_error():
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
                'Vmag': 12.3,
                'err_Vmag': 0.0,
            },
        ],
    }

    match = exotic_module.nextastro_photometry_catalog_match(catalog, 10.0, 20.0, 'CV')

    assert match['error'] == pytest.approx(0.001)


def test_nextastro_photometry_catalog_match_accepts_relaxed_v_magnitude_error():
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
                'Vmag': 12.3,
                'err_Vmag': 0.051,
            },
        ],
    }

    match = exotic_module.nextastro_photometry_catalog_match(catalog, 10.0, 20.0, 'CV')

    assert match['mag_band'] == 'V'
    assert match['error'] == pytest.approx(0.051)
    assert match['uses_relaxed_bv_error_limit'] is True


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


def test_nextastro_photometry_catalog_match_honors_scale_aware_radius():
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

    default_match = exotic_module.nextastro_photometry_catalog_match(
        catalog,
        10.0,
        20.0,
        'CV',
    )
    scale_aware_match = exotic_module.nextastro_photometry_catalog_match(
        catalog,
        10.0,
        20.0,
        'CV',
        max_separation_arcsec=5.2,
    )

    assert default_match is None
    assert scale_aware_match['source_id'] == 111
    assert 2.0 < scale_aware_match['separation_arcsec'] < 5.2


def test_nextastro_catalog_match_radius_uses_one_pixel_with_two_arcsec_floor():
    assert exotic_module.nextastro_catalog_match_radius_arcsec(None) == pytest.approx(2.0)
    assert exotic_module.nextastro_catalog_match_radius_arcsec(1.4) == pytest.approx(2.0)
    assert exotic_module.nextastro_catalog_match_radius_arcsec(5.153485) == pytest.approx(5.153485)


def test_direct_selected_catalog_candidate_uses_targeted_scale_aware_lookup(monkeypatch):
    calls = []

    def fake_lookup(ra, dec, observed_filter, radius_arcsec=2.0):
        calls.append((ra, dec, observed_filter, radius_arcsec))
        return {
            'id': 185212647,
            'source_id': None,
            'mag': 10.038,
            'error': 0.027,
            'mag_band': 'V',
            'catalog_ra': 294.68751,
            'catalog_dec': 31.360037,
            'separation_arcsec': 2.643594,
        }

    monkeypatch.setattr(exotic_module, 'nextastro_photometry_for_coordinate', fake_lookup)

    candidate = exotic_module.build_direct_selected_catalog_candidate(
        [[232.0, 348.0]],
        [(294.68834158880503, 31.36022406806484)],
        {'row_format': 'objects', 'rows': []},
        0,
        observed_filter='CV',
        match_radius_arcsec=5.153485,
    )

    assert calls == [(294.68834158880503, 31.36022406806484, 'CV', 5.153485)]
    assert candidate['source'] == 'direct_catalog'
    assert candidate['star']['mag'] == pytest.approx(10.038)
    assert candidate['star']['mag_band'] == 'V'
    assert candidate['star']['pos'] == [232.0, 348.0]


def test_reported_stellar_variability_band_distinguishes_clearv_from_catalog_v():
    for observed_filter in (
        'MObs CV',
        'CV',
        'Clear',
        'Luminance',
        'Photographic G',
        'Gaia G',
        'G',
        'G1',
        'G2',
    ):
        assert exotic_module.reported_stellar_variability_band(
            observed_filter, 'V'
        ) == 'ClearV'
    assert exotic_module.reported_stellar_variability_band('bv', 'V') == 'V'
    for observed_filter in ('R', 'SR', 'rp'):
        assert exotic_module.reported_stellar_variability_band(
            observed_filter, 'r'
        ) == 'rp'


def test_selected_comparison_direct_catalog_match_precedes_derived_fallback():
    selected = exotic_module.choose_selected_comp_catalog_reference_candidate(
        [
            {
                'source': 'field_derived',
                'error': 0.010,
                'star': {'mag_band': 'V', 'error': 0.010},
            },
            {
                'source': 'direct_catalog',
                'error': 0.027,
                'star': {'mag_band': 'V', 'error': 0.027},
            },
        ],
        observed_filter='CV',
    )

    assert selected['source'] == 'direct_catalog'


def test_aavso_vsp_band_for_filter_uses_observed_filter_aliases():
    assert exotic_module.aavso_vsp_band_for_filter('MObs CV') == 'V'
    assert exotic_module.aavso_vsp_band_for_filter('Clear (unfiltered) reduced to V sequence') == 'V'
    assert exotic_module.aavso_vsp_band_for_filter('Photographic G') == 'V'
    assert exotic_module.aavso_vsp_band_for_filter('Gaia G') == 'V'
    assert exotic_module.aavso_vsp_band_for_filter('Cousins R') == 'Rc'
    assert exotic_module.aavso_vsp_band_for_filter('Sloan g') == 'SG'


def test_nextastro_catalog_band_tokens_keep_bessell_sloan_and_clearv_distinct():
    assert exotic_module.nextastro_photometry_band_candidates('bv') == [
        ('Vmag', 'err_Vmag', 'V')
    ]
    assert exotic_module.nextastro_photometry_band_candidates('bb') == [
        ('Bmag', 'err_Bmag', 'B')
    ]
    assert exotic_module.nextastro_photometry_band_candidates('Sloan g') == [('g', 'dg', 'g')]
    for observed_filter in ('G', 'Photographic G', 'Gaia G', 'G1', 'G2'):
        assert exotic_module.nextastro_photometry_band_candidates(observed_filter) == [
            ('Vmag', 'err_Vmag', 'V')
        ]
    assert exotic_module.catalog_band_priority('V', 'bv') == 0
    assert exotic_module.catalog_band_priority('B', 'bb') == 0
    assert exotic_module.catalog_band_priority('g', 'Sloan g') == 0
    assert exotic_module.catalog_band_priority('g', 'G') == 1
    assert exotic_module.catalog_band_priority('g', 'Gaia G') == 1
    assert exotic_module.catalog_band_priority('V', 'G') == 0
    assert exotic_module.catalog_band_priority('V', 'Gaia G') == 0


@pytest.mark.parametrize(
    ('observed_filter', 'expected_labels'),
    [
        ('u', ['u-g', 'B-V', 'BP-RP']),
        ('Johnson U', ['u-g', 'B-V', 'BP-RP']),
        ('B', ['B-V', 'BP-RP']),
        ('Photographic B', ['B-V', 'BP-RP']),
        ('V', ['B-V', 'BP-RP']),
        ('Sloan g', ['g-r', 'B-V', 'BP-RP']),
        ('Sloan r', ['r-i', 'B-V', 'BP-RP']),
        ('Cousins R', ['r-i', 'B-V', 'BP-RP']),
        ('Sloan i', ['r-i', 'B-V', 'BP-RP']),
        ('Cousins I', ['r-i', 'B-V', 'BP-RP']),
        ('Sloan z', ['i-z', 'B-V', 'BP-RP']),
        ('CV', ['B-V', 'BP-RP']),
        ('Clear', ['B-V', 'BP-RP']),
        ('Luminance', ['B-V', 'BP-RP']),
        ('Photographic G', ['B-V', 'BP-RP']),
        ('Gaia G', ['B-V', 'BP-RP']),
        ('G', ['B-V', 'BP-RP']),
        ('G1', ['B-V', 'BP-RP']),
        ('G2', ['B-V', 'BP-RP']),
        ('Unknown', ['B-V', 'BP-RP']),
    ],
)
def test_all_filters_use_filter_specific_color_then_universal_fallbacks(
        observed_filter, expected_labels):
    pairs = exotic_module.nextastro_color_candidate_pairs(observed_filter)
    assert [pair[2] for pair in pairs] == expected_labels


@pytest.mark.parametrize(
    ('observed_filter', 'primary_label', 'primary_columns'),
    [
        ('u', 'u-g', ('umag',)),
        ('B', 'B-V', ()),
        ('V', 'B-V', ()),
        ('Sloan g', 'g-r', ('g',)),
        ('Sloan r', 'r-i', ('r',)),
        ('Sloan i', 'r-i', ('r',)),
        ('Sloan z', 'i-z', ('i',)),
        ('Clear', 'B-V', ()),
        ('Photographic G', 'B-V', ()),
        ('Gaia G', 'B-V', ()),
    ],
)
def test_all_filters_fall_back_from_specific_color_to_bv_then_bp_rp(
        observed_filter, primary_label, primary_columns):
    row = {
        'Bmag': 13.0,
        'Vmag': 12.5,
        'umag': 13.8,
        'g': 12.8,
        'r': 12.2,
        'i': 12.0,
        'z': 11.8,
        'bp_rp': 1.1,
    }
    assert exotic_module.nextastro_catalog_color(row, observed_filter)['label'] == primary_label

    for column in primary_columns:
        row[column] = None
    assert exotic_module.nextastro_catalog_color(row, observed_filter)['label'] == 'B-V'

    row['Bmag'] = None
    assert exotic_module.nextastro_catalog_color(row, observed_filter) == {
        'color': pytest.approx(1.1),
        'label': 'BP-RP',
        'first_column': 'bp_rp',
        'second_column': None,
    }


def test_clearv_catalog_color_derives_bp_rp_from_gaia_magnitudes():
    color = exotic_module.nextastro_catalog_color(
        {
            'Vmag': 12.5,
            'phot_bp_mean_mag': 13.4,
            'phot_rp_mean_mag': 12.1,
        },
        'Gaia G',
    )

    assert color == {
        'color': pytest.approx(1.3),
        'label': 'BP-RP',
        'first_column': 'phot_bp_mean_mag',
        'second_column': 'phot_rp_mean_mag',
    }


def test_nextastro_gaia_bp_rp_lookup_is_cached(monkeypatch):
    captured = []

    def fake_get(url, params, timeout):
        captured.append((url, params, timeout))
        return DummyResponse({
            'gaia': {
                'source_id': 123456,
                'separation_arcsec': 0.2,
                'phot_bp_mean_mag': 13.4,
                'phot_rp_mean_mag': 12.1,
                'bp_rp': 1.3,
            },
        })

    exotic_module._cached_nextastro_gaia_bp_rp.cache_clear()
    monkeypatch.setattr(exotic_module.requests, 'get', fake_get)

    first = exotic_module.nextastro_gaia_bp_rp_for_coordinate(10.12345678, -20.25)
    second = exotic_module.nextastro_gaia_bp_rp_for_coordinate(10.12345678, -20.25)

    assert first == second
    assert first['color'] == pytest.approx(1.3)
    assert first['label'] == 'BP-RP'
    assert first['catalog_source'] == 'NextAstro Gaia DR3'
    assert captured == [(
        exotic_module.NEXTASTRO_GAIA_DISTPM_ENDPOINT,
        {'ra': 10.1234568, 'dec': -20.25},
        exotic_module.NEXTASTRO_GAIA_COLOR_LOOKUP_TIMEOUT_SECONDS,
    )]


def test_local_catalog_color_does_not_query_gaia(monkeypatch):
    monkeypatch.setattr(
        exotic_module,
        'nextastro_gaia_bp_rp_for_coordinate',
        lambda *args, **kwargs: pytest.fail('Gaia should not be queried when B-V is available.'),
    )
    state = {'remaining': 3, 'attempted': 0, 'matched': 0}

    color = exotic_module.nextastro_catalog_color_with_gaia_fallback(
        {'ra': 10.0, 'dec': 20.0, 'Bmag': 13.0, 'Vmag': 12.5},
        'V',
        lookup_state=state,
    )

    assert color['label'] == 'B-V'
    assert color['color'] == pytest.approx(0.5)
    assert state == {'remaining': 3, 'attempted': 0, 'matched': 0}


def test_nearest_catalog_color_row_uses_gaia_bp_rp_last_resort(monkeypatch):
    calls = []

    def fake_gaia_lookup(ra, dec, max_separation_arcsec):
        calls.append((ra, dec, max_separation_arcsec))
        return {'color': 1.25, 'label': 'BP-RP'}

    monkeypatch.setattr(exotic_module, 'nextastro_gaia_bp_rp_for_coordinate', fake_gaia_lookup)
    state = {'remaining': 3, 'attempted': 0, 'matched': 0}
    catalog = {
        'rows': [
            {'source_id': 42, 'ra': 10.00001, 'dec': 20.0, 'Vmag': 12.5},
        ],
    }

    match = exotic_module.nextastro_catalog_nearest_color_row(
        catalog,
        10.0,
        20.0,
        'V',
        gaia_lookup_state=state,
        gaia_match_radius_arcsec=1.5,
    )

    assert match['source_id'] == 42
    assert match['color'] == {'color': 1.25, 'label': 'BP-RP'}
    assert calls == [(10.00001, 20.0, 1.5)]
    assert state == {'remaining': 2, 'attempted': 1, 'matched': 1}


@pytest.mark.parametrize(
    ('observed_filter', 'magnitude_column', 'error_column'),
    [
        ('u', 'umag', 'err_umag'),
        ('B', 'Bmag', 'err_Bmag'),
        ('V', 'Vmag', 'err_Vmag'),
        ('Sloan g', 'g', 'dg'),
        ('Sloan r', 'r', 'dr'),
        ('Sloan i', 'i', 'di'),
        ('Sloan z', 'z', 'dz'),
    ],
)
def test_nextastro_catalog_match_never_falls_back_to_another_band(
        observed_filter, magnitude_column, error_column):
    row = {
        'id': 1,
        'ra': 10.0,
        'dec': 20.0,
        'Bmag': 12.1,
        'err_Bmag': 0.01,
        'Vmag': 12.2,
        'err_Vmag': 0.01,
        'umag': 12.3,
        'err_umag': 0.01,
        'g': 12.4,
        'dg': 0.01,
        'r': 12.5,
        'dr': 0.01,
        'i': 12.6,
        'di': 0.01,
        'z': 12.7,
        'dz': 0.01,
    }
    row[magnitude_column] = None
    row[error_column] = None

    match = exotic_module.nextastro_photometry_catalog_match(
        {'row_format': 'objects', 'rows': [row]},
        10.0,
        20.0,
        observed_filter,
    )

    assert match is None


def test_nextastro_catalog_match_uses_relaxed_error_only_as_bv_fallback():
    catalog = {
        'row_format': 'objects',
        'rows': [
            {
                'id': 1, 'ra': 10.00001, 'dec': 20.0,
                'Vmag': 12.1, 'err_Vmag': 0.07,
            },
            {
                'id': 2, 'ra': 10.00010, 'dec': 20.0,
                'Vmag': 12.2, 'err_Vmag': 0.03,
            },
        ],
    }

    preferred = exotic_module.nextastro_photometry_catalog_match(catalog, 10.0, 20.0, 'V')
    assert preferred['id'] == 2
    assert preferred['uses_relaxed_bv_error_limit'] is False

    relaxed_v = exotic_module.nextastro_photometry_catalog_match(
        {'row_format': 'objects', 'rows': [catalog['rows'][0]]},
        10.0,
        20.0,
        'bv',
    )
    assert relaxed_v['id'] == 1
    assert relaxed_v['error'] == pytest.approx(0.07)
    assert relaxed_v['uses_relaxed_bv_error_limit'] is True

    relaxed_b = exotic_module.nextastro_photometry_catalog_match(
        {
            'row_format': 'objects',
            'rows': [{'id': 3, 'ra': 10.0, 'dec': 20.0, 'Bmag': 13.0, 'err_Bmag': 0.10}],
        },
        10.0,
        20.0,
        'bb',
    )
    assert relaxed_b['id'] == 3
    assert relaxed_b['uses_relaxed_bv_error_limit'] is True

    rejected_g = exotic_module.nextastro_photometry_catalog_match(
        {
            'row_format': 'objects',
            'rows': [{'id': 4, 'ra': 10.0, 'dec': 20.0, 'g': 13.0, 'dg': 0.051}],
        },
        10.0,
        20.0,
        'Sloan g',
    )
    assert rejected_g is None

    rejected_v = exotic_module.nextastro_photometry_catalog_match(
        {
            'row_format': 'objects',
            'rows': [{'id': 5, 'ra': 10.0, 'dec': 20.0, 'Vmag': 13.0, 'err_Vmag': 0.101}],
        },
        10.0,
        20.0,
        'V',
    )
    assert rejected_v is None


def test_nextastro_photometry_for_coordinate_uses_single_object_endpoint(monkeypatch):
    captured = {}

    def fake_post(url, json, timeout):
        captured['url'] = url
        captured['json'] = json
        captured['timeout'] = timeout
        return DummyResponse({
            'columns': list(exotic_module.NEXTASTRO_PHOTOMETRY_COLUMNS),
            'match': {
                'id': 9,
                'source_id': 12345,
                'ra': 10.00001,
                'dec': -20.00001,
                'Vmag': 11.2,
                'err_Vmag': 0.03,
            },
            'separation_arcsec': 0.05,
        })

    monkeypatch.setattr(exotic_module.requests, 'post', fake_post)
    monkeypatch.setattr(exotic_module, 'log_info', lambda *args, **kwargs: None)

    match = exotic_module.nextastro_photometry_for_coordinate(10.0, -20.0, 'V')

    assert captured['url'] == 'https://photometry.nextastro.org/single_object'
    assert captured['json']['ra'] == pytest.approx(10.0)
    assert captured['json']['dec'] == pytest.approx(-20.0)
    assert captured['json']['radius_arcsec'] == pytest.approx(2.0)
    assert captured['json']['columns'] == [
        'id', 'source_id', 'ra', 'dec', 'Vmag', 'err_Vmag',
    ]
    assert captured['json']['required_columns'] == ['Vmag', 'err_Vmag']
    assert captured['timeout'] == 30
    assert match['source_id'] == 12345
    assert match['mag'] == pytest.approx(11.2)
    assert match['mag_band'] == 'V'


def test_nextastro_photometry_for_coordinates_uses_objects_query_endpoint(monkeypatch):
    captured = {}

    def fake_post(url, json, timeout):
        captured['url'] = url
        captured['json'] = json
        captured['timeout'] = timeout
        return DummyResponse({
            'columns': list(exotic_module.NEXTASTRO_PHOTOMETRY_COLUMNS),
            'count': 1,
            'results': [
                {
                    'key': '0',
                    'ra': 10.0,
                    'dec': -20.0,
                    'match': {
                        'id': 9,
                        'source_id': 12345,
                        'ra': 10.00001,
                        'dec': -20.00001,
                        'Vmag': 11.2,
                        'err_Vmag': 0.03,
                    },
                    'separation_arcsec': 0.05,
                },
                {
                    'key': '1',
                    'ra': 11.0,
                    'dec': -21.0,
                    'match': None,
                    'separation_arcsec': None,
                },
            ],
        })

    monkeypatch.setattr(exotic_module.requests, 'post', fake_post)
    monkeypatch.setattr(exotic_module, 'log_info', lambda *args, **kwargs: None)

    matches = exotic_module.nextastro_photometry_for_coordinates(
        [(10.0, -20.0), (11.0, -21.0)],
        'V',
    )

    assert captured['url'] == 'https://photometry.nextastro.org/objects_query'
    assert captured['json']['objects'] == [
        {'key': '0', 'ra': 10.0, 'dec': -20.0},
        {'key': '1', 'ra': 11.0, 'dec': -21.0},
    ]
    assert captured['json']['radius_arcsec'] == pytest.approx(2.0)
    assert captured['json']['columns'] == [
        'id', 'source_id', 'ra', 'dec', 'Vmag', 'err_Vmag',
    ]
    assert captured['json']['required_columns'] == ['Vmag', 'err_Vmag']
    assert captured['timeout'] == 30
    assert matches[0]['source_id'] == 12345
    assert matches[0]['mag_band'] == 'V'
    assert matches[1] is None


def test_nextastro_photometry_for_coordinates_routes_single_target_to_single_object(monkeypatch):
    calls = []
    expected_match = {'source_id': 12345, 'mag': 11.2, 'error': 0.03, 'mag_band': 'V'}

    def fake_single(ra, dec, obs_filter, radius_arcsec=2.0):
        calls.append((ra, dec, obs_filter, radius_arcsec))
        return expected_match

    monkeypatch.setattr(exotic_module, 'nextastro_photometry_for_coordinate', fake_single)
    monkeypatch.setattr(
        exotic_module,
        'nextastro_photometry_objects_query',
        lambda *args, **kwargs: pytest.fail('one target should not use /objects_query'),
    )

    matches = exotic_module.nextastro_photometry_for_coordinates([(10.0, -20.0)], 'V')

    assert calls == [(10.0, -20.0, 'V', 2.0)]
    assert matches == [expected_match]


def test_merge_nextastro_calibration_stars_batches_missing_field_matches(monkeypatch):
    calls = []

    def fake_matches(coordinates, obs_filter, radius_arcsec=2.0):
        calls.append((coordinates, obs_filter, radius_arcsec))
        return [
            {
                'source_id': 101,
                'id': 1,
                'mag': 11.2,
                'error': 0.03,
                'mag_band': 'V',
                'catalog_ra': 10.0,
                'catalog_dec': -20.0,
                'separation_arcsec': 0.1,
                'catalog_row': {'source_id': 101, 'ra': 10.0, 'dec': -20.0},
            },
            {
                'source_id': 202,
                'id': 2,
                'mag': 12.1,
                'error': 0.04,
                'mag_band': 'V',
                'catalog_ra': 11.0,
                'catalog_dec': -21.0,
                'separation_arcsec': 0.2,
                'catalog_row': {'source_id': 202, 'ra': 11.0, 'dec': -21.0},
            },
        ]

    monkeypatch.setattr(
        exotic_module,
        'nextastro_photometry_for_coordinates',
        fake_matches,
    )
    monkeypatch.setattr(exotic_module, 'log_info', lambda *args, **kwargs: None)

    calibration_stars = exotic_module.merge_nextastro_calibration_stars(
        comp_stars=[[100, 200], [130, 230]],
        comp_ra_dec=[(10.0, -20.0), (11.0, -21.0)],
        obs_filter='V',
        existing_comp_stars={},
        field_catalog=None,
    )

    assert calls == [([(10.0, -20.0), (11.0, -21.0)], 'V', 2.0)]
    assert list(calibration_stars) == ['NextAstro-101', 'NextAstro-202']


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


def test_merge_nextastro_calibration_stars_uses_image_scale_match_radius():
    catalog = {
        'row_format': 'objects',
        'rows': [{
            'id': 9,
            'source_id': 12345,
            'ra': 10.001,
            'dec': 20.0,
            'Vmag': 10.038,
            'err_Vmag': 0.027,
        }],
    }

    calibration_stars = exotic_module.merge_nextastro_calibration_stars(
        comp_stars=[[232, 348]],
        comp_ra_dec=[(10.0, 20.0)],
        obs_filter='MObs CV',
        field_catalog=catalog,
        match_radius_arcsec=5.153485,
    )

    assert list(calibration_stars) == ['NextAstro-12345']
    assert calibration_stars['NextAstro-12345']['mag'] == pytest.approx(10.038)
    assert calibration_stars['NextAstro-12345']['mag_band'] == 'V'
    assert calibration_stars['NextAstro-12345']['separation_arcsec'] > 2.0


def test_merge_nextastro_calibration_stars_deduplicates_catalog_source_ids():
    catalog = {
        'columns': ['id', 'source_id', 'ra', 'dec', 'Vmag', 'err_Vmag'],
        'count': 1,
        'row_format': 'objects',
        'rows': [{
            'id': 9,
            'source_id': 12345,
            'ra': 10.0,
            'dec': -20.0,
            'Vmag': 11.2,
            'err_Vmag': 0.03,
        }],
    }

    calibration_stars = exotic_module.merge_nextastro_calibration_stars(
        comp_stars=[[100, 200], [130, 230]],
        comp_ra_dec=[(10.0, -20.0), (10.0001, -20.0001)],
        obs_filter='V',
        existing_comp_stars={},
        field_catalog=catalog,
    )

    assert list(calibration_stars) == ['NextAstro-12345']


def test_fetch_aavso_vsp_chart_retries_malformed_json_five_times_then_succeeds(monkeypatch):
    payload = {'chartid': 'X-RETRY', 'photometry': []}
    responses = [None] * exotic_module.AAVSO_VSP_MAX_RETRIES + [DummyResponse(payload)]
    request_timeouts = []
    sleep_delays = []
    log_messages = []

    class InvalidJSONResponse:
        def raise_for_status(self):
            return None

        def json(self):
            return json.loads('')

    def fake_get(url, timeout):
        request_timeouts.append(timeout)
        response = responses.pop(0)
        return InvalidJSONResponse() if response is None else response

    monkeypatch.setattr(exotic_module.requests, 'get', fake_get)
    monkeypatch.setattr(exotic_module, 'sleep', sleep_delays.append)
    monkeypatch.setattr(
        exotic_module,
        'log_info',
        lambda message, **kwargs: log_messages.append((message, kwargs)),
    )

    assert exotic_module.fetch_aavso_vsp_chart('https://example.invalid/vsp') == payload
    assert request_timeouts == [
        exotic_module.AAVSO_VSP_REQUEST_TIMEOUT_SECONDS
    ] * (exotic_module.AAVSO_VSP_MAX_RETRIES + 1)
    assert sleep_delays == [
        exotic_module.AAVSO_VSP_RETRY_DELAY_SECONDS
    ] * exotic_module.AAVSO_VSP_MAX_RETRIES
    assert 'attempt 1/6' in log_messages[0][0]
    assert 'attempt 5/6' in log_messages[-1][0]
    assert all(kwargs.get('warn') is True for _, kwargs in log_messages)


def test_fetch_aavso_vsp_chart_raises_after_five_failed_retries(monkeypatch):
    request_count = 0
    sleep_delays = []

    class InvalidJSONResponse:
        def raise_for_status(self):
            return None

        def json(self):
            return json.loads('')

    def fake_get(url, timeout):
        nonlocal request_count
        request_count += 1
        assert timeout == exotic_module.AAVSO_VSP_REQUEST_TIMEOUT_SECONDS
        return InvalidJSONResponse()

    monkeypatch.setattr(exotic_module.requests, 'get', fake_get)
    monkeypatch.setattr(exotic_module, 'sleep', sleep_delays.append)
    monkeypatch.setattr(exotic_module, 'log_info', lambda *args, **kwargs: None)

    with pytest.raises(exotic_module.AAVSOVSPUnavailableError) as exc_info:
        exotic_module.fetch_aavso_vsp_chart('https://example.invalid/vsp')

    assert request_count == exotic_module.AAVSO_VSP_MAX_RETRIES + 1
    assert sleep_delays == [
        exotic_module.AAVSO_VSP_RETRY_DELAY_SECONDS
    ] * exotic_module.AAVSO_VSP_MAX_RETRIES
    assert 'after 6 attempts (5 retries)' in str(exc_info.value)
    assert 'JSONDecodeError' in str(exc_info.value)


def test_vsp_query_rejects_band_errors_over_limit(monkeypatch):
    class DummyWCS:
        def pixel_to_world_values(self, x_pixel, y_pixel):
            return 10.0, 20.0

        def world_to_pixel_values(self, ra_deg, dec_deg):
            return np.array([40.0]), np.array([50.0])

    payload = {
        'chartid': 'X123',
        'photometry': [
            {
                'auid': 'HIGH',
                'ra': '00:00:00.0',
                'dec': '+00:00:00.0',
                'bands': [{'band': 'V', 'mag': 12.0, 'error': 0.051}],
            },
            {
                'auid': 'LOW',
                'ra': '00:00:00.0',
                'dec': '+00:00:00.0',
                'bands': [{'band': 'V', 'mag': 12.1, 'error': 0.05}],
            },
        ],
    }
    user_comp_stars = []

    monkeypatch.setattr(exotic_module, 'search_wcs', lambda file: DummyWCS())
    monkeypatch.setattr(exotic_module, 'radec_hours_to_degree', lambda ra, dec: (10.0, 20.0))
    monkeypatch.setattr(
        exotic_module.requests,
        'get',
        lambda url, timeout: DummyResponse(payload),
    )
    monkeypatch.setattr(exotic_module, 'log_info', lambda *args, **kwargs: None)

    vsp_comp_stars, chart_id = exotic_module.vsp_query(
        'frame.fits',
        [100, 100],
        'MObs CV',
        1.0,
        user_comp_stars=user_comp_stars,
        user_targ_star=[10, 10],
    )

    assert chart_id == 'X123'
    assert list(vsp_comp_stars) == ['LOW']
    assert vsp_comp_stars['LOW']['error'] == pytest.approx(0.05)
    assert user_comp_stars == [[40, 50]]


def test_vsp_query_keeps_late_supplied_matches_after_new_star_limit(monkeypatch):
    class DummyWCS:
        def pixel_to_world_values(self, x_pixel, y_pixel):
            return 10.0, 20.0

        def world_to_pixel_values(self, ra_deg, dec_deg):
            return np.array([ra_deg]), np.array([dec_deg])

    payload = {
        'chartid': 'X-LIMIT',
        'photometry': [
            {
                'auid': 'NEW-1',
                'ra': '20',
                'dec': '20',
                'bands': [{'band': 'V', 'mag': 11.0, 'error': 0.01}],
            },
            {
                'auid': 'NEW-2',
                'ra': '40',
                'dec': '40',
                'bands': [{'band': 'V', 'mag': 12.0, 'error': 0.01}],
            },
            {
                'auid': 'SUPPLIED',
                'ra': '80',
                'dec': '80',
                'bands': [{'band': 'V', 'mag': 13.0, 'error': 0.02}],
            },
        ],
    }
    user_comp_stars = [[80, 80]]

    monkeypatch.setattr(exotic_module, 'search_wcs', lambda file: DummyWCS())
    monkeypatch.setattr(
        exotic_module,
        'radec_hours_to_degree',
        lambda ra, dec: (float(ra), float(dec)),
    )
    monkeypatch.setattr(
        exotic_module.requests,
        'get',
        lambda url, timeout: DummyResponse(payload),
    )
    monkeypatch.setattr(exotic_module, 'log_info', lambda *args, **kwargs: None)

    vsp_comp_stars, chart_id = exotic_module.vsp_query(
        'frame.fits',
        [100, 100],
        'Clear',
        1.0,
        user_comp_stars=user_comp_stars,
        max_new_comp_stars=1,
    )

    assert chart_id == 'X-LIMIT'
    assert list(vsp_comp_stars) == ['NEW-1', 'SUPPLIED']
    assert vsp_comp_stars['SUPPLIED']['pos'] == [80, 80]
    assert user_comp_stars == [[80, 80], [20, 20]]


def test_vsp_query_assigns_only_nearest_catalog_source_to_supplied_coordinate(monkeypatch):
    class DummyWCS:
        def pixel_to_world_values(self, x_pixel, y_pixel):
            return 10.0, 20.0

        def world_to_pixel_values(self, ra_deg, dec_deg):
            return np.array([ra_deg]), np.array([dec_deg])

    payload = {
        'chartid': 'X-NEAREST',
        'photometry': [
            {
                'auid': 'FARTHER',
                'ra': '47',
                'dec': '50',
                'bands': [{'band': 'V', 'mag': 11.0, 'error': 0.01}],
            },
            {
                'auid': 'NEAREST',
                'ra': '51',
                'dec': '50',
                'bands': [{'band': 'V', 'mag': 12.0, 'error': 0.02}],
            },
        ],
    }
    user_comp_stars = [[50, 50]]

    monkeypatch.setattr(exotic_module, 'search_wcs', lambda file: DummyWCS())
    monkeypatch.setattr(
        exotic_module,
        'radec_hours_to_degree',
        lambda ra, dec: (float(ra), float(dec)),
    )
    monkeypatch.setattr(
        exotic_module.requests,
        'get',
        lambda url, timeout: DummyResponse(payload),
    )
    monkeypatch.setattr(exotic_module, 'log_info', lambda *args, **kwargs: None)

    vsp_comp_stars, _ = exotic_module.vsp_query(
        'frame.fits',
        [100, 100],
        'Clear',
        1.0,
        user_comp_stars=user_comp_stars,
        max_new_comp_stars=0,
    )

    assert list(vsp_comp_stars) == ['NEAREST']
    assert vsp_comp_stars['NEAREST']['pos'] == [50, 50]
    assert user_comp_stars == [[50, 50]]


def test_tracked_comparison_position_keeps_full_field_anchor_index_after_science_reset():
    science_comp_stars = [[217.0, 210.0], [408.0, 261.0]]
    tracked_calibration_stars = [
        *science_comp_stars,
        [415.0, 203.0],
        [449.0, 267.0],
    ]

    # The target fit restores the shorter science list, while catalog-anchor
    # indices retain the full tracking-list index space.
    restored_science_comp_stars = list(science_comp_stars)
    assert len(restored_science_comp_stars) == 2
    assert exotic_module.tracked_comparison_position(
        tracked_calibration_stars,
        2,
    ) == [415.0, 203.0]
    assert exotic_module.tracked_comparison_position(
        tracked_calibration_stars,
        3,
    ) == [449.0, 267.0]


def test_selected_comparison_finder_entries_include_every_ensemble_member():
    entries = exotic_module.selected_comparison_finder_entries(
        [[100.0, 200.0], [300.0, 400.0], [500.0, 600.0]],
        ensemble_member_keys=['comp3', 'comp1'],
    )

    assert entries == [
        {'key': 'comp3', 'label': 'Comp 3', 'position': [500.0, 600.0]},
        {'key': 'comp1', 'label': 'Comp 1', 'position': [100.0, 200.0]},
    ]


def test_selected_comparison_finder_entries_do_not_depend_on_aavso_metadata():
    entries = exotic_module.selected_comparison_finder_entries(
        [[217.0, 210.0]],
        comp_index=0,
    )

    assert entries == [
        {'key': 'comp1', 'label': 'Comp 1', 'position': [217.0, 210.0]},
    ]


def test_clear_v_calibration_fallback_merges_aavso_with_existing_pool(monkeypatch):
    calls = []
    supplied_positions = [[100, 200]]

    def fake_vsp_query(file, axis, obs_filter, img_scale, **kwargs):
        calls.append((file, axis, obs_filter, img_scale, kwargs))
        kwargs['user_comp_stars'].append([300, 400])
        return {
            '000-BPW-929': {
                'pos': [100, 200],
                'mag': 11.85,
                'error': 0.046,
                'mag_band': 'V',
                'catalog_source': 'AAVSO VSP',
                'is_aavso_vsp': True,
            },
            '000-BMX-191': {
                'pos': [300, 400],
                'mag': 12.121,
                'error': 0.005,
                'mag_band': 'V',
                'catalog_source': 'AAVSO VSP',
                'is_aavso_vsp': True,
            },
        }, 'X42753ZU'

    monkeypatch.setattr(exotic_module, 'vsp_query', fake_vsp_query)
    monkeypatch.setattr(exotic_module, 'log_info', lambda *args, **kwargs: None)

    combined, fallback_stars, chart_id, queried = (
        exotic_module.merge_aavso_vsp_v_calibration_fallback(
            'frame.fits',
            [512, 512],
            'Clear',
            1.2,
            {
                'NextAstro-g-only': {
                    'pos': [100, 200],
                    'mag': 11.7,
                    'error': 0.01,
                    'mag_band': 'g',
                    'catalog_source': 'NextAstro photometry catalog',
                },
            },
            supplied_positions,
            user_targ_star=[250, 250],
        )
    )

    assert queried is True
    assert chart_id == 'X42753ZU'
    assert len(calls) == 1
    assert calls[0][2] == 'Clear'
    assert calls[0][4]['max_new_comp_stars'] == 5
    assert supplied_positions == [[100, 200], [300, 400]]
    assert set(fallback_stars) == {'000-BPW-929', '000-BMX-191'}
    assert set(combined) == {'NextAstro-g-only', '000-BPW-929', '000-BMX-191'}


def test_clear_v_calibration_fallback_skips_vsp_when_nextastro_has_usable_v(monkeypatch):
    def unexpected_vsp_query(*args, **kwargs):
        raise AssertionError('VSP must not be queried when NextAstro supplied usable V')

    monkeypatch.setattr(exotic_module, 'vsp_query', unexpected_vsp_query)

    existing = {
        'NextAstro-123': {
            'pos': [100, 200],
            'mag': 11.7,
            'error': 0.02,
            'mag_band': 'V',
            'catalog_source': 'NextAstro photometry catalog',
        },
    }
    combined, fallback_stars, chart_id, queried = (
        exotic_module.merge_aavso_vsp_v_calibration_fallback(
            'frame.fits',
            [512, 512],
            'Clear',
            1.2,
            existing,
            [[100, 200]],
        )
    )

    assert combined == existing
    assert fallback_stars == {}
    assert chart_id is None
    assert queried is False


def test_clear_v_calibration_fallback_does_not_repeat_exhausted_vsp_query(monkeypatch):
    def unexpected_vsp_query(*args, **kwargs):
        raise AssertionError('An exhausted VSP request must not start another retry cycle')

    log_messages = []
    monkeypatch.setattr(exotic_module, 'vsp_query', unexpected_vsp_query)
    monkeypatch.setattr(
        exotic_module,
        'log_info',
        lambda message, **kwargs: log_messages.append((message, kwargs)),
    )

    combined, fallback_stars, chart_id, queried = (
        exotic_module.merge_aavso_vsp_v_calibration_fallback(
            'frame.fits',
            [512, 512],
            'Clear',
            1.2,
            {},
            [[100, 200]],
            vsp_query_available=False,
        )
    )

    assert combined == {}
    assert fallback_stars == {}
    assert chart_id is None
    assert queried is False
    assert 'already exhausted all retries' in log_messages[-1][0]
    assert log_messages[-1][1].get('warn') is True


@pytest.mark.parametrize(
    'observed_filter',
    ['CV', 'Clear', 'Luminance', 'Photographic G', 'Gaia G'],
)
def test_build_stellar_variability_params_records_nextastro_reference(
        monkeypatch, tmp_path, observed_filter):
    captured = {}

    class DummyFit:
        time = np.array([2450000.105, 2450000.205, 2450000.305], dtype=float)
        data = np.array([1.0, 1.02, 0.98], dtype=float)
        dataerr = np.full(3, 0.01, dtype=float)
        airmass_model = np.ones(3, dtype=float)
        airmass = np.array([1.1, 1.2, 1.3], dtype=float)
        jd_times = np.array([2450000.1, 2450000.2, 2450000.3], dtype=float)
        transit = np.ones(3, dtype=float)
        stellar_variability_target_flux = np.array([1000.0, 1020.0, 980.0], dtype=float)
        stellar_variability_comp_flux = np.full(3, 1000.0, dtype=float)
        stellar_variability_target_flux_error = np.full(3, 2.0, dtype=float)
        stellar_variability_comp_flux_error = np.full(3, 2.0, dtype=float)

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
        observed_filter=observed_filter,
    )

    assert captured['label'] == 'RA=10.1000000 Dec=-20.2000000'
    assert len(params) == 3
    assert params[0]['catalog_source'] == 'NextAstro photometry catalog'
    assert params[0]['is_aavso_vsp'] is False
    assert params[0]['comp_ra'] == pytest.approx(10.1)
    assert params[0]['comp_dec'] == pytest.approx(-20.2)
    assert params[0]['cmag'] == pytest.approx(12.0)
    assert params[0]['cmag_err'] == pytest.approx(0.05)
    assert params[0]['observed_filter'] == observed_filter
    assert params[0]['mag_band'] == 'ClearV'
    assert params[0]['catalog_mag_band'] == 'V'
    assert [row['time'] for row in params] == pytest.approx(DummyFit.time)
    assert [row['jd_time'] for row in params] == pytest.approx(DummyFit.jd_times)


def test_build_stellar_variability_params_rejects_cross_band_calibration(tmp_path):
    class DummyFit:
        data = np.ones(3, dtype=float)
        dataerr = np.full(3, 0.01, dtype=float)
        airmass_model = np.ones(3, dtype=float)
        airmass = np.array([1.1, 1.2, 1.3], dtype=float)
        jd_times = np.array([2450000.1, 2450000.2, 2450000.3], dtype=float)
        transit = np.ones(3, dtype=float)
        stellar_variability_target_flux = np.full(3, 1000.0, dtype=float)
        stellar_variability_comp_flux = np.full(3, 1000.0, dtype=float)
        stellar_variability_target_flux_error = np.full(3, 2.0, dtype=float)
        stellar_variability_comp_flux_error = np.full(3, 2.0, dtype=float)

    with pytest.raises(RuntimeError, match='cross-band absolute calibration is not permitted'):
        exotic_module.build_stellar_variability_params_from_fit(
            DummyFit(),
            {
                'mag': 11.615,
                'error': 0.001,
                'mag_band': 'g',
                'observed_filter': 'V',
            },
            [100, 200],
            'NextAstro-invalid-g-reference',
            tmp_path,
            'Host Star',
            observed_filter='V',
            observation_date='2026-08-02',
        )

    differential_csv = next(
        tmp_path.glob('StellarVariabilityDifferentialMagnitude_HostStar_2026-08-02.csv')
    )
    assert '# AIRMASS_CORRECTION=NO' in differential_csv.read_text(encoding='utf-8')
    assert (
        tmp_path / 'StellarVariabilityDifferentialMagnitude_HostStar_2026-08-02.png'
    ).exists()


def test_build_stellar_variability_params_uses_raw_ratio_and_per_exposure_errors(monkeypatch, tmp_path):
    comp_mag = 9.751
    comp_mag_error = 0.018
    target_mag = 13.1
    flux_ratio = 10 ** ((comp_mag - target_mag) / 2.5)
    detrended = flux_ratio * np.array([0.94, 1.0, 1.06], dtype=float)

    class DummyFit:
        # The fitted series is intentionally normalized: the absolute target
        # magnitude must come from the retained raw target/comparison fluxes.
        data = detrended / np.nanmedian(detrended)
        dataerr = np.full(3, 0.01, dtype=float)
        airmass_model = np.array([0.94, 1.0, 1.06], dtype=float)
        airmass = np.array([1.1, 1.2, 1.3], dtype=float)
        jd_times = np.array([2450000.1, 2450000.2, 2450000.3], dtype=float)
        transit = np.ones(3, dtype=float)
        stellar_variability_comp_flux = np.full(3, 100000.0, dtype=float)
        stellar_variability_target_flux = stellar_variability_comp_flux * detrended
        stellar_variability_target_flux_error = np.array([20.0, 21.0, 22.0], dtype=float)
        stellar_variability_comp_flux_error = np.array([30.0, 31.0, 32.0], dtype=float)

    monkeypatch.setattr(exotic_module, 'plot_stellar_variability', lambda *args, **kwargs: None)

    calibration_star = {
        'mag': comp_mag,
        'error': comp_mag_error,
        'catalog_source': 'AAVSO VSP',
        'is_aavso_vsp': True,
        'mag_band': 'V',
        'observed_filter': 'V',
    }

    params = exotic_module.build_stellar_variability_params_from_fit(
        DummyFit(),
        calibration_star,
        [100, 200],
        '000-BJX-718',
        tmp_path,
        'HAT-P-37',
        observed_filter='CV',
    )

    expected_mag_error = np.hypot(
        comp_mag_error,
        (2.5 / np.log(10.0)) * np.hypot(
            DummyFit.stellar_variability_target_flux_error[1]
            / DummyFit.stellar_variability_target_flux[1],
            DummyFit.stellar_variability_comp_flux_error[1]
            / DummyFit.stellar_variability_comp_flux[1],
        ),
    )

    expected_raw_magnitudes = comp_mag - (2.5 * np.log10(detrended))
    np.testing.assert_allclose(
        [row['mag'] for row in params],
        expected_raw_magnitudes,
        atol=1.0e-10,
    )
    np.testing.assert_allclose(
        [row['differential_mag'] for row in params],
        -2.5 * np.log10(detrended),
        atol=1.0e-10,
    )
    assert params[0]['mag'] != pytest.approx(target_mag)
    assert params[1]['mag_err'] == pytest.approx(expected_mag_error)
    assert params[1]['differential_mag_err'] == pytest.approx(
        np.sqrt(expected_mag_error ** 2 - comp_mag_error ** 2)
    )
    assert params[1]['mag_err'] < 0.08


def test_annotate_stellar_variability_raw_photometry_restores_final_selected_fit_fluxes():
    class DummyFit:
        data = np.array([0.99, 1.0, 1.01], dtype=float)

    target_flux = np.array([9900.0, 10000.0, 10100.0], dtype=float)
    comp_flux = np.full(3, 20000.0, dtype=float)
    target_error = np.array([10.0, 11.0, 12.0], dtype=float)
    comp_error = np.array([20.0, 21.0, 22.0], dtype=float)
    fit = DummyFit()

    exotic_module.annotate_stellar_variability_raw_photometry(
        fit,
        target_flux,
        comp_flux,
        target_flux_error=target_error,
        comp_flux_error=comp_error,
    )

    retained = exotic_module.stellar_variability_raw_photometry(fit)
    for actual, expected in zip(
        retained,
        (target_flux, comp_flux, target_error, comp_error),
    ):
        np.testing.assert_array_equal(actual, expected)


def test_build_stellar_variability_params_rejects_normalized_only_absolute_calibration(
        monkeypatch, tmp_path):
    class DummyFit:
        data = np.array([0.99, 1.0, 1.01], dtype=float)
        dataerr = np.full(3, 0.01, dtype=float)
        airmass = np.ones(3, dtype=float)
        jd_times = np.array([2450000.1, 2450000.2, 2450000.3], dtype=float)
        transit = np.ones(3, dtype=float)

    monkeypatch.setattr(exotic_module, 'plot_stellar_variability', lambda *args, **kwargs: None)

    with pytest.raises(RuntimeError, match='cannot be recovered from a normalized light curve'):
        exotic_module.build_stellar_variability_params_from_fit(
            DummyFit(),
            {'mag': 12.0, 'error': 0.02, 'mag_band': 'V'},
            [100, 200],
            'COMP',
            tmp_path,
            'Host Star',
            observed_filter='V',
        )


def test_stellar_variability_requires_selected_transit_comparison(monkeypatch, tmp_path):
    logged = []

    class DummyFit:
        data = np.array([1.0, 1.01, 0.99], dtype=float)
        airmass_model = np.ones(3, dtype=float)
        airmass = np.ones(3, dtype=float)
        jd_times = np.array([2450000.1, 2450000.2, 2450000.3], dtype=float)
        transit = np.ones(3, dtype=float)

    monkeypatch.setattr(exotic_module, 'log_info', lambda message, warn=False, error=False: logged.append(message))

    params = exotic_module.stellar_variability(
        {0: {'myfit': DummyFit(), 'pos': [100, 200]}},
        DummyFit(),
        [[100, 200]],
        {'REF': {'pos': [100, 200], 'mag': 12.0, 'error': 0.02}},
        [0],
        None,
        tmp_path,
        'Host Star',
    )

    assert params == []
    assert any('no transit-fit comparison star' in message for message in logged)


def test_stellar_variability_derives_selected_comparison_catalog_magnitude(monkeypatch, tmp_path):
    logged = []
    captured = {}

    class DummyFit:
        def __init__(self, data):
            self.data = np.array(data, dtype=float)
            self.dataerr = np.full(3, 0.01, dtype=float)
            self.airmass_model = np.ones(3, dtype=float)
            self.airmass = np.ones(3, dtype=float)
            self.jd_times = np.array([2450000.1, 2450000.2, 2450000.3], dtype=float)
            self.transit = np.ones(3, dtype=float)
            self.stellar_variability_target_flux = self.data * 1000.0
            self.stellar_variability_comp_flux = np.full(3, 1000.0, dtype=float)
            self.stellar_variability_target_flux_error = np.full(3, 2.0, dtype=float)
            self.stellar_variability_comp_flux_error = np.full(3, 2.0, dtype=float)

    monkeypatch.setattr(exotic_module, 'log_info', lambda message, warn=False, error=False: logged.append(message))
    monkeypatch.setattr(
        exotic_module,
        'plot_stellar_variability',
        lambda params, save, s_name, label: captured.update(params=params, label=label),
    )

    selected_mag = 11.0
    anchor_mag = 12.0
    selected_to_anchor_flux_ratio = 10 ** ((anchor_mag - selected_mag) / 2.5)
    selected_fit = DummyFit([1.0, 1.01, 0.99])
    anchor_fit = DummyFit(selected_to_anchor_flux_ratio * np.array([1.0, 1.01, 0.99]))

    params = exotic_module.stellar_variability(
        {
            0: {'myfit': selected_fit, 'pos': [100, 200]},
            1: {'myfit': anchor_fit, 'pos': [300, 400]},
        },
        DummyFit([1.0, 1.01, 0.99]),
        [[100, 200], [300, 400]],
        {'REF': {'pos': [300, 400], 'mag': anchor_mag, 'error': 0.02}},
        [1],
        0,
        tmp_path,
        'Host Star',
        comp_ra_dec=[(10.0, -20.0), (11.0, -21.0)],
    )

    assert len(params) == 3
    assert params[0]['cmag'] == pytest.approx(selected_mag)
    assert params[0]['cmag_err'] == pytest.approx(0.02)
    assert params[0]['comp_ra'] == pytest.approx(10.0)
    assert params[0]['comp_dec'] == pytest.approx(-20.0)
    assert params[0]['derived_catalog_reference'] is True
    assert params[0]['derived_reference_anchor_count'] == 1
    assert params[0]['derived_reference_anchor_labels'] == ['REF']
    assert captured['label'] == 'RA=10.0000000 Dec=-20.0000000'
    assert any('derived catalog magnitude' in message for message in logged)


def test_stellar_variability_uses_direct_catalog_without_shared_fit_oot_points(monkeypatch, tmp_path):
    logged = []

    class DummyFit:
        def __init__(self, data):
            self.data = np.array(data, dtype=float)
            self.dataerr = np.full(4, 0.01, dtype=float)
            self.airmass_model = np.ones(4, dtype=float)
            self.airmass = np.ones(4, dtype=float)
            self.jd_times = np.array([2450000.1, 2450000.2, 2450000.3, 2450000.4], dtype=float)
            self.transit = np.ones(4, dtype=float)
            self.stellar_variability_target_flux = self.data * 1000.0
            self.stellar_variability_comp_flux = np.full(4, 1000.0, dtype=float)
            self.stellar_variability_target_flux_error = np.full(4, 2.0, dtype=float)
            self.stellar_variability_comp_flux_error = np.full(4, 2.0, dtype=float)

    monkeypatch.setattr(exotic_module, 'log_info', lambda message, warn=False, error=False: logged.append(message))
    monkeypatch.setattr(exotic_module, 'plot_stellar_variability', lambda *args, **kwargs: None)

    selected_fit = DummyFit([1.0, 1.0, 1.0, 1.0])
    noisy_anchor_fit = DummyFit([2.0, 0.8, 2.2, 0.7])
    best_fit = DummyFit([1.0, 1.0, 1.0, 1.0])
    best_fit.transit = np.zeros(4, dtype=float)
    params = exotic_module.stellar_variability(
        {
            0: {'myfit': selected_fit, 'pos': [100, 200]},
            1: {'myfit': noisy_anchor_fit, 'pos': [300, 400]},
        },
        best_fit,
        [[100, 200], [300, 400]],
        {'ANCHOR': {'pos': [300, 400], 'mag': 12.0, 'error': 0.02, 'mag_band': 'g'}},
        [1],
        0,
        tmp_path,
        'Host Star',
        observed_filter='g',
        comp_ra_dec=[(10.0, -20.0), (11.0, -21.0)],
        field_catalog={
            'rows': [{
                'ra': 10.0,
                'dec': -20.0,
                'g': 11.5,
                'dg': 0.08,
                'source_id': 12345,
            }]
        },
    )

    assert len(params) == 4
    assert params[0]['cmag'] == pytest.approx(11.5)
    assert params[0]['cmag_err'] == pytest.approx(0.08)
    assert params[0]['derived_catalog_reference'] is False
    assert params[0]['allow_high_error_catalog_reference'] is True
    assert any('direct selected-comparison catalog magnitude' in message for message in logged)


def test_stellar_variability_derives_catalog_magnitude_from_full_field(monkeypatch, tmp_path):
    class DummyFit:
        data = np.array([1.0, 1.01, 0.99], dtype=float)
        dataerr = np.full(3, 0.01, dtype=float)
        airmass_model = np.ones(3, dtype=float)
        airmass = np.ones(3, dtype=float)
        jd_times = np.array([2450000.1, 2450000.2, 2450000.3], dtype=float)
        transit = np.ones(3, dtype=float)
        stellar_variability_target_flux = data * 1000.0
        stellar_variability_comp_flux = np.full(3, 1000.0, dtype=float)
        stellar_variability_target_flux_error = np.full(3, 2.0, dtype=float)
        stellar_variability_comp_flux_error = np.full(3, 2.0, dtype=float)

    class DummyWcs:
        def world_to_pixel_values(self, ra, dec):
            return float(ra), float(dec)

        def pixel_to_world_values(self, x, y):
            return 123.4, -45.6

    image = np.full((60, 60), 10.0, dtype=float)
    image[20, 20] = 50.0
    image[40, 40] = 110.0

    monkeypatch.setattr(exotic_module, 'plot_stellar_variability', lambda *args, **kwargs: None)
    monkeypatch.setattr(exotic_module, 'search_wcs', lambda _path: DummyWcs())

    params = exotic_module.stellar_variability(
        {0: {'myfit': DummyFit(), 'pos': [20, 20]}},
        DummyFit(),
        [[20, 20]],
        {},
        [],
        0,
        tmp_path,
        'Host Star',
        observed_filter='g',
        field_catalog={
            'rows': [{
                'ra': 40.0,
                'dec': 40.0,
                'g': 12.0,
                'dg': 0.03,
                'source_id': 67890,
            }]
        },
        reference_image=image,
        wcs_file='dummy.wcs',
    )

    expected_mag = 12.0 - 2.5 * np.log10(40.0 / 100.0)
    assert len(params) == 3
    assert params[0]['cmag'] == pytest.approx(expected_mag)
    assert params[0]['cmag_err'] == pytest.approx(0.03)
    assert params[0]['comp_ra'] == pytest.approx(123.4)
    assert params[0]['comp_dec'] == pytest.approx(-45.6)
    assert params[0]['derived_catalog_reference'] is True
    assert params[0]['derived_reference_anchor_count'] == 1
    assert params[0]['derived_reference_anchor_labels'] == ['NextAstro-67890']


def test_derived_catalog_reference_skips_missing_anchor_fit():
    class DummyFit:
        data = np.array([1.0, 1.01, 0.99], dtype=float)

    label, star = exotic_module.derived_catalog_reference_for_selected_comp(
        {
            0: {'myfit': DummyFit(), 'pos': [10, 10]},
            1: None,
        },
        [[10, 10], [20, 20]],
        {
            'Anchor': {
                'pos': [20, 20],
                'mag': 12.0,
                'error': 0.03,
                'mag_band': 'V',
            },
        },
        [1],
        0,
        observed_filter='V',
    )

    assert label is None
    assert star is None


def test_derived_catalog_reference_ensembles_multiple_aavso_v_anchors():
    class DummyFit:
        def __init__(self, reference_curve):
            reference_curve = np.asarray(reference_curve, dtype=float)
            self.data = reference_curve
            self.time = np.array([1.0, 2.0, 3.0], dtype=float)
            self.transit = np.ones(3, dtype=float)
            self.stellar_variability_target_flux = reference_curve * 1000.0
            self.stellar_variability_comp_flux = np.full(3, 1000.0, dtype=float)
            self.stellar_variability_target_flux_error = np.full(3, 2.0, dtype=float)
            self.stellar_variability_comp_flux_error = np.full(3, 2.0, dtype=float)

    selected_mag = 11.0
    first_anchor_mag = 12.0
    second_anchor_mag = 13.0
    first_ratio = 10.0 ** ((first_anchor_mag - selected_mag) / 2.5)
    second_ratio = 10.0 ** ((second_anchor_mag - selected_mag) / 2.5)

    label, star = exotic_module.derived_catalog_reference_for_selected_comp(
        {
            0: {'myfit': DummyFit(np.ones(3)), 'pos': [10, 10]},
            1: {'myfit': DummyFit(np.full(3, first_ratio)), 'pos': [20, 20]},
            2: {'myfit': DummyFit(np.full(3, second_ratio)), 'pos': [30, 30]},
        },
        [[10, 10], [20, 20], [30, 30]],
        {
            'AAVSO-1': {
                'pos': [20, 20],
                'mag': first_anchor_mag,
                'error': 0.02,
                'mag_band': 'V',
                'catalog_source': 'AAVSO VSP',
                'is_aavso_vsp': True,
            },
            'AAVSO-2': {
                'pos': [30, 30],
                'mag': second_anchor_mag,
                'error': 0.04,
                'mag_band': 'V',
                'catalog_source': 'AAVSO VSP',
                'is_aavso_vsp': True,
            },
        },
        [1, 2],
        0,
        observed_filter='Clear',
    )

    assert label == 'Derived Comp 1'
    assert star['mag'] == pytest.approx(selected_mag)
    assert star['error'] == pytest.approx((1.0 / (1.0 / 0.02 ** 2 + 1.0 / 0.04 ** 2)) ** 0.5)
    assert star['mag_band'] == 'V'
    assert star['derived_catalog_reference'] is True
    assert star['derived_reference_anchor_count'] == 2
    assert star['derived_reference_anchor_labels'] == ['AAVSO-1', 'AAVSO-2']


def test_stellar_variability_rejects_g_catalog_anchor_for_clearv(monkeypatch, tmp_path):
    logged = []

    class DummyFit:
        data = np.array([1.0, 1.01, 0.99], dtype=float)
        dataerr = np.full(3, 0.01, dtype=float)
        airmass_model = np.ones(3, dtype=float)
        airmass = np.ones(3, dtype=float)
        jd_times = np.array([2450000.1, 2450000.2, 2450000.3], dtype=float)
        transit = np.ones(3, dtype=float)
        stellar_variability_target_flux = data * 1000.0
        stellar_variability_comp_flux = np.full(3, 1000.0, dtype=float)
        stellar_variability_target_flux_error = np.full(3, 2.0, dtype=float)
        stellar_variability_comp_flux_error = np.full(3, 2.0, dtype=float)

    class DummyWcs:
        def world_to_pixel_values(self, ra, dec):
            return float(ra), float(dec)

        def pixel_to_world_values(self, x, y):
            return 123.4, -45.6

    image = np.full((60, 60), 10.0, dtype=float)
    image[20, 20] = 50.0
    image[40, 40] = 110.0

    monkeypatch.setattr(exotic_module, 'plot_stellar_variability', lambda *args, **kwargs: None)
    monkeypatch.setattr(exotic_module, 'search_wcs', lambda _path: DummyWcs())
    monkeypatch.setattr(exotic_module, 'log_info', lambda message, warn=False, error=False: logged.append(message))

    params = exotic_module.stellar_variability(
        {0: {'myfit': DummyFit(), 'pos': [20, 20]}},
        DummyFit(),
        [[20, 20]],
        {},
        [],
        0,
        tmp_path,
        'Host Star',
        observed_filter='CV',
        field_catalog={
            'rows': [{
                'ra': 40.0,
                'dec': 40.0,
                'g': 12.0,
                'dg': 0.03,
                'source_id': 67890,
            }]
        },
        reference_image=image,
        wcs_file='dummy.wcs',
    )

    assert params == []
    assert any('no derived magnitude could be inferred' in message for message in logged)


def test_stellar_variability_uses_v_catalog_anchor_for_clearv(monkeypatch, tmp_path):
    class DummyFit:
        data = np.array([1.0, 1.01, 0.99], dtype=float)
        dataerr = np.full(3, 0.01, dtype=float)
        airmass_model = np.ones(3, dtype=float)
        airmass = np.ones(3, dtype=float)
        jd_times = np.array([2450000.1, 2450000.2, 2450000.3], dtype=float)
        transit = np.ones(3, dtype=float)
        stellar_variability_target_flux = data * 1000.0
        stellar_variability_comp_flux = np.full(3, 1000.0, dtype=float)
        stellar_variability_target_flux_error = np.full(3, 2.0, dtype=float)
        stellar_variability_comp_flux_error = np.full(3, 2.0, dtype=float)

    class DummyWcs:
        def world_to_pixel_values(self, ra, dec):
            return float(ra), float(dec)

        def pixel_to_world_values(self, x, y):
            return 123.4, -45.6

    image = np.full((60, 60), 10.0, dtype=float)
    image[20, 20] = 50.0
    image[40, 40] = 110.0

    monkeypatch.setattr(exotic_module, 'plot_stellar_variability', lambda *args, **kwargs: None)
    monkeypatch.setattr(exotic_module, 'search_wcs', lambda _path: DummyWcs())

    params = exotic_module.stellar_variability(
        {0: {'myfit': DummyFit(), 'pos': [20, 20]}},
        DummyFit(),
        [[20, 20]],
        {},
        [],
        0,
        tmp_path,
        'Host Star',
        observed_filter='CV',
        field_catalog={
            'rows': [{
                'ra': 40.0,
                'dec': 40.0,
                'Vmag': 12.0,
                'err_Vmag': 0.03,
                'g': 11.7,
                'dg': 0.01,
                'source_id': 67890,
            }]
        },
        reference_image=image,
        wcs_file='dummy.wcs',
    )

    assert len(params) == 3
    assert params[0]['mag_band'] == 'ClearV'
    assert params[0]['catalog_mag_band'] == 'V'
    assert params[0]['cmag_err'] == pytest.approx(0.03)


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


def _stellar_variability_only_planet_dict():
    return {
        'pName': 'Synthetic b',
        'sName': 'Synthetic',
        'pPer': 1.0,
        'pPerUnc': 0.0001,
        'midT': 10.0,
        'midTUnc': 0.0001,
        'rprs': 0.1,
        'rprsUnc': 0.001,
        'aRs': 12.0,
        'aRsUnc': 0.2,
        'inc': 89.0,
        'incUnc': 0.1,
        'ecc': 0.0,
        'omega': 90.0,
    }


def test_stellar_variability_out_of_transit_mask_excludes_predicted_transit_window():
    p_dict = _stellar_variability_only_planet_dict()
    duration = exotic_module.estimate_transit_duration_from_prior_geometry(
        exotic_module.stellar_variability_transit_prior_from_planet_dict(p_dict)
    )
    times = np.array([
        p_dict['midT'] - duration,
        p_dict['midT'],
        p_dict['midT'] + 0.49 * duration,
        p_dict['midT'] + duration,
    ])

    keep_mask, summary = exotic_module.stellar_variability_out_of_transit_mask(times, p_dict)

    assert keep_mask.tolist() == [True, False, False, True]
    assert summary['applied'] is True
    assert summary['rejected_point_count'] == 2
    assert summary['duration_days'] == pytest.approx(duration)


def test_build_stellar_variability_only_lightcurve_discards_transit_points(monkeypatch):
    p_dict = _stellar_variability_only_planet_dict()
    duration = exotic_module.estimate_transit_duration_from_prior_geometry(
        exotic_module.stellar_variability_transit_prior_from_planet_dict(p_dict)
    )
    offsets = np.array([-3.0, -2.2, -1.4, -0.7, -0.1, 0.0, 0.1, 0.7, 1.4, 2.2, 3.0]) * duration
    times = p_dict['midT'] + offsets
    target_flux = np.full(times.shape, 10000.0)
    comp_flux = np.full(times.shape, 10000.0)
    flux_err = np.full(times.shape, 20.0)
    airmass = np.ones(times.shape)

    monkeypatch.setattr(
        exotic_module,
        'get_phase',
        lambda t, per, tmid: ((np.asarray(t, dtype=float) - tmid) / per + 0.5) % 1.0 - 0.5,
    )

    fit, prepared = exotic_module.build_stellar_variability_only_lightcurve_from_fluxes(
        times,
        target_flux,
        comp_flux,
        airmass,
        p_dict,
        jd_times=times,
        target_flux_error=flux_err,
        comp_flux_error=flux_err,
        exposure_times_seconds=np.full(times.shape, 60.0),
        gain_e_per_adu=1.0,
        comp_index=0,
        comp_label="Comp 1",
        comp_position=[1, 2],
        method_label="PSF photometry",
    )

    assert prepared['applied'] is True
    assert fit is not None
    assert fit.stellar_variability_only is True
    assert np.all(fit.transit == 1.0)
    assert fit.stellar_variability_transit_exclusion['rejected_point_count'] == 3
    assert not np.any(np.isclose(fit.time, p_dict['midT']))
    assert len(fit.time) == times.size - 3


def test_stellar_variability_ensemble_masks_saturated_frames_and_error_clips_members():
    frame_count = 8
    ranked_summaries = [
        {
            'key': 'comp1', 'comp_index': 0, 'label': 'Comp 1', 'position': [10, 20],
            'overexposure_rejected_count': 0,
        },
        {
            'key': 'comp2', 'comp_index': 1, 'label': 'Comp 2', 'position': [30, 40],
            'overexposure_rejected_count': 0,
        },
        {
            'key': 'comp3', 'comp_index': 2, 'label': 'Comp 3', 'position': [50, 60],
            'overexposure_rejected_count': 1,
        },
        {
            'key': 'comp4', 'comp_index': 3, 'label': 'Comp 4', 'position': [70, 80],
            'overexposure_rejected_count': 0,
        },
    ]
    comp_flux_map = {
        'comp1': np.full(frame_count, 1000.0),
        'comp2': np.full(frame_count, 2000.0),
        'comp3': np.full(frame_count, 3000.0),
        'comp4': np.full(frame_count, 1500.0),
    }
    comp_flux_map['comp3'][0] = np.nan
    calibration_stars = {
        'C1': {'pos': [10, 20], 'mag': 12.0, 'error': 0.010, 'mag_band': 'V'},
        'C2': {'pos': [30, 40], 'mag': 12.5, 'error': 0.011, 'mag_band': 'V'},
        'C3': {'pos': [50, 60], 'mag': 12.2, 'error': 0.010, 'mag_band': 'V'},
        'C4': {'pos': [70, 80], 'mag': 12.3, 'error': 0.200, 'mag_band': 'V'},
    }

    selection = exotic_module.select_stellar_variability_ensemble_members(
        ranked_summaries,
        calibration_stars,
        comp_flux_map,
        observed_filter='V',
    )

    assert [member['key'] for member in selection['members']] == ['comp3', 'comp2', 'comp1']
    assert selection['members'][0]['summary']['overexposure_rejected_count'] == 1
    rejected_reasons = {item['key']: item['reason'] for item in selection['rejected']}
    assert 'sigma-clip' in rejected_reasons['comp4']
    assert selection['calibration_error_clip']['high_threshold'] < 0.2
    assert selection['calibration_error_clip']['high_threshold'] >= 0.01
    assert selection['calibration_error_clip']['minimum_high_threshold'] == pytest.approx(0.01)


def test_stellar_variability_ensemble_skips_cross_band_catalog_reference():
    frame_count = 8
    ranked_summaries = [
        {
            'key': 'comp1', 'comp_index': 0, 'label': 'Comp 1', 'position': [10, 20],
            'overexposure_rejected_count': 0,
        },
        {
            'key': 'comp2', 'comp_index': 1, 'label': 'Comp 2', 'position': [30, 40],
            'overexposure_rejected_count': 0,
        },
    ]
    calibration_stars = {
        'C1': {'pos': [10, 20], 'mag': 11.0, 'error': 0.001, 'mag_band': 'g'},
        'C2': {'pos': [30, 40], 'mag': 12.5, 'error': 0.011, 'mag_band': 'V'},
    }
    selection = exotic_module.select_stellar_variability_ensemble_members(
        ranked_summaries,
        calibration_stars,
        {
            'comp1': np.full(frame_count, 2000.0),
            'comp2': np.full(frame_count, 1000.0),
        },
        observed_filter='V',
        min_members=1,
        max_members=1,
    )

    assert [member['key'] for member in selection['members']] == ['comp2']
    assert selection['members'][0]['star']['mag_band'] == 'V'
    rejected = {item['key']: item['reason'] for item in selection['rejected']}
    assert rejected['comp1'] == 'no usable catalog calibration'


def test_stellar_variability_ensemble_combines_nextastro_and_aavso_v_members():
    frame_count = 8
    ranked_summaries = [
        {
            'key': 'comp1', 'comp_index': 0, 'label': 'Comp 1', 'position': [10, 20],
            'overexposure_rejected_count': 0,
        },
        {
            'key': 'comp2', 'comp_index': 1, 'label': 'Comp 2', 'position': [30, 40],
            'overexposure_rejected_count': 0,
        },
    ]
    calibration_stars = {
        'NextAstro-123': {
            'pos': [10, 20],
            'mag': 11.8,
            'error': 0.02,
            'mag_band': 'V',
            'catalog_source': 'NextAstro photometry catalog',
            'source_id': 123,
        },
        '000-BMX-191': {
            'pos': [30, 40],
            'mag': 12.121,
            'error': 0.005,
            'mag_band': 'V',
            'catalog_source': 'AAVSO VSP',
            'is_aavso_vsp': True,
            'catalog_ra': 18.0,
            'catalog_dec': 35.0,
        },
    }

    selection = exotic_module.select_stellar_variability_ensemble_members(
        ranked_summaries,
        calibration_stars,
        {
            'comp1': np.full(frame_count, 2000.0),
            'comp2': np.full(frame_count, 1000.0),
        },
        observed_filter='Clear',
    )

    assert [member['key'] for member in selection['members']] == ['comp1', 'comp2']
    assert [member['label'] for member in selection['members']] == [
        'NextAstro-123',
        '000-BMX-191',
    ]
    assert [member['star']['catalog_source'] for member in selection['members']] == [
        'NextAstro photometry catalog',
        'AAVSO VSP',
    ]


def test_stellar_variability_single_mode_can_select_aavso_v_fallback_member():
    frame_count = 8
    selection = exotic_module.select_stellar_variability_ensemble_members(
        [
            {
                'key': 'comp1', 'comp_index': 0, 'label': 'Comp 1', 'position': [10, 20],
                'overexposure_rejected_count': 0,
            },
            {
                'key': 'comp2', 'comp_index': 1, 'label': 'Comp 2', 'position': [30, 40],
                'overexposure_rejected_count': 0,
            },
        ],
        {
            '000-BPW-929': {
                'pos': [10, 20],
                'mag': 11.85,
                'error': 0.046,
                'mag_band': 'V',
                'catalog_source': 'AAVSO VSP',
                'is_aavso_vsp': True,
                'catalog_ra': 18.0,
                'catalog_dec': 35.0,
            },
            '000-BMX-191': {
                'pos': [30, 40],
                'mag': 12.121,
                'error': 0.005,
                'mag_band': 'V',
                'catalog_source': 'AAVSO VSP',
                'is_aavso_vsp': True,
                'catalog_ra': 18.1,
                'catalog_dec': 35.1,
            },
        },
        {
            'comp1': np.full(frame_count, 2000.0),
            'comp2': np.full(frame_count, 1000.0),
        },
        observed_filter='Clear',
        min_members=1,
        max_members=1,
    )

    assert [member['key'] for member in selection['members']] == ['comp1']
    assert selection['members'][0]['label'] == '000-BPW-929'
    assert selection['members'][0]['star']['catalog_source'] == 'AAVSO VSP'


def test_stellar_variability_rejects_comparison_that_steps_across_acquisition_gap():
    frame_count = 180
    cadence_days = 6.0 / 86400.0
    times = 2460000.0 + (np.arange(frame_count, dtype=float) * cadence_days)
    times[90:] += 120.0 / 86400.0
    ranked_summaries = [
        {
            'key': f'comp{index}',
            'comp_index': index - 1,
            'label': f'Comp {index}',
            'position': [10 * index, 20 * index],
            'overexposure_rejected_count': 0,
        }
        for index in range(1, 5)
    ]
    calibration_stars = {
            f'C{index}': {
                'pos': [10 * index, 20 * index],
                'mag': 12.0 + (0.1 * index),
                'error': 0.01,
                'mag_band': 'V',
                'catalog_source': 'Synthetic catalog',
            }
        for index in range(1, 5)
    }
    phase = np.linspace(0, 4 * np.pi, frame_count)
    common = 1000.0 * (1.0 + (0.01 * np.sin(phase)))
    comp_flux_map = {
        f'comp{index}': common * (1.0 + (0.001 * index * np.cos(phase)))
        for index in range(1, 5)
    }
    # A +0.06 mag instrumental discontinuity in comparison 1.
    comp_flux_map['comp1'] = comp_flux_map['comp1'].copy()
    comp_flux_map['comp1'][90:] *= 10.0 ** (-0.4 * 0.06)

    selection = exotic_module.select_stellar_variability_ensemble_members(
        ranked_summaries,
        calibration_stars,
        comp_flux_map,
        observed_filter='V',
        min_members=1,
        max_members=1,
        times=times,
    )

    assert selection['gap_stability']['applied'] is True
    assert selection['gap_stability']['boundaries'][0]['source_index'] == 90
    assert selection['gap_stability']['candidates']['comp1']['rejected'] is True
    assert selection['gap_stability']['candidates']['comp1'][
        'maximum_absolute_step_magnitude'
    ] == pytest.approx(0.06, abs=0.003)
    assert selection['gap_stability']['candidates']['comp1']['label'] == 'C1'
    assert selection['gap_stability']['candidates']['comp1']['position'] == [10, 20]
    assert selection['gap_stability']['candidates']['comp1']['catalog_magnitude_band'] == 'V'
    assert selection['gap_stability']['candidates']['comp1']['catalog_source'] == 'Synthetic catalog'
    assert selection['members'][0]['key'] != 'comp1'
    rejected = {item['key']: item for item in selection['rejected']}
    assert rejected['comp1']['label'] == 'C1'
    assert rejected['comp1']['position'] == [10, 20]
    assert 'changed discontinuously' in rejected['comp1']['reason']


def test_fortuitous_output_error_limit_relaxes_only_for_flagged_bv_catalog_reference():
    assert exotic_module.fortuitous_output_magnitude_error_limit([
        {'star': {'mag_band': 'V', 'uses_relaxed_bv_error_limit': True}}
    ]) == pytest.approx(0.10)
    assert exotic_module.fortuitous_output_magnitude_error_limit([
        {'star': {'mag_band': 'B', 'uses_relaxed_bv_error_limit': True}}
    ]) == pytest.approx(0.10)
    assert exotic_module.fortuitous_output_magnitude_error_limit([
        {'star': {'mag_band': 'g', 'uses_relaxed_bv_error_limit': True}}
    ]) == pytest.approx(0.05)
    assert exotic_module.fortuitous_output_magnitude_error_limit([
        {'star': {'mag_band': 'V', 'uses_relaxed_bv_error_limit': False}}
    ]) == pytest.approx(0.05)

    metadata = exotic_module.fortuitous_variable_target_metadata({
        'name': 'Synthetic',
        'reference_mode': 'single_comparison',
        'output_magnitude_error_limit': 0.10,
    })
    assert metadata['detection_magnitude_error_limit'] == pytest.approx(0.05)
    assert metadata['output_magnitude_error_limit'] == pytest.approx(0.10)


def test_stellar_variability_ensemble_error_clip_does_not_reject_below_point_zero_one_mag():
    candidates = [
        {'magnitude_error': error}
        for error in (0.0010, 0.0011, 0.0012, 0.0090, 0.0110)
    ]

    keep, summary = exotic_module.stellar_variability_ensemble_calibration_error_clip(candidates)

    assert keep.tolist() == [True, True, True, True, False]
    assert summary['high_threshold'] == pytest.approx(0.01)
    assert summary['minimum_high_threshold'] == pytest.approx(0.01)


def test_automatic_comparison_merge_deduplicates_only_added_sources():
    merged, messages = exotic_module.merge_automatic_comparison_star_coords(
        [[10.0, 20.0], [11.0, 20.0]],
        [[10.4, 20.3], [50.0, 60.0], [50.5, 60.2]],
        duplicate_radius_pixels=2.0,
    )

    # Nearby primary/user selections remain intentional; automatic additions
    # cannot repeat either a primary source or an earlier automatic source.
    assert merged == [[10.0, 20.0], [11.0, 20.0], [50.0, 60.0]]
    assert len(messages) == 2


def test_full_field_vsx_variables_are_removed_from_science_comparisons():
    comparison_stars = [
        [4969.0, 1695.0],
        [5221.0, 2714.0],
        [3923.0, 1362.0],
    ]
    fortuitous_variables = [
        {'name': 'DI Her', 'pos': [4971.55, 1695.42]},
        {
            'name': 'ASASSN-V J185327.35+241158.6',
            'x': 3926.38,
            'y': 1360.36,
        },
    ]

    retained, rejected = exotic_module.filter_comparison_stars_against_fortuitous_variables(
        comparison_stars,
        fortuitous_variables,
        duplicate_radius_pixels=10.0,
    )

    assert retained == [[5221.0, 2714.0]]
    assert [item['comparison_index'] for item in rejected] == [0, 2]
    assert [item['variable_name'] for item in rejected] == [
        'DI Her',
        'ASASSN-V J185327.35+241158.6',
    ]
    assert rejected[0]['distance_pixels'] == pytest.approx(2.5840, abs=1.0e-3)
    assert rejected[1]['distance_pixels'] == pytest.approx(3.7563, abs=1.0e-3)


def test_tracked_vsx_overexposure_warning_does_not_call_variable_a_comparison_star():
    messages = []
    status = PlateStatus(lambda message, **kwargs: messages.append(message))
    status.setCurrentFilename('frame.fits')

    status.overexposedWarning(
        12,
        4973.2,
        1676.6,
        58981.5,
        starLabel='Tracked VSX variable DI Her',
    )

    assert messages[0] == (
        'Tracked VSX variable DI Her is overexposed in file frame.fits; '
        'aperture pixels near [4973.2, 1676.6] exceeded 58981.5.'
    )
    assert 'repeated frame-level star warnings are aggregated' in messages[1]
    assert 'Comparison star' not in messages[0]


def test_stellar_variability_ensemble_caps_at_five_by_target_color_and_magnitude():
    frame_count = 8
    target_match = {
        'mag': 12.0,
        'error': 0.01,
        'mag_band': 'V',
        'catalog_row': {'Bmag': 12.5, 'Vmag': 12.0},
    }
    candidate_values = [
        (10.0, -0.5),
        (11.0, 0.0),
        (12.1, 0.55),
        (12.2, 0.60),
        (11.9, 0.45),
        (12.3, 0.40),
        (12.0, 0.52),
    ]
    ranked_summaries = []
    calibration_stars = {}
    comp_flux_map = {}
    for index, (magnitude, color) in enumerate(candidate_values, start=1):
        key = f'comp{index}'
        position = [index * 10, index * 10 + 1]
        ranked_summaries.append({
            'key': key,
            'comp_index': index - 1,
            'label': f'Comp {index}',
            'position': position,
            'overexposure_rejected_count': 0,
        })
        calibration_stars[f'C{index}'] = {
            'pos': position,
            'mag': magnitude,
            'error': 0.01,
            'mag_band': 'V',
            'catalog_row': {'Bmag': magnitude + color, 'Vmag': magnitude},
        }
        comp_flux_map[key] = np.full(frame_count, 10000.0 - index * 100.0)

    selection = exotic_module.select_stellar_variability_ensemble_members(
        ranked_summaries,
        calibration_stars,
        comp_flux_map,
        observed_filter='V',
        target_catalog_match=target_match,
    )

    assert selection['prelimit_member_count'] == 7
    assert selection['member_limit'] == 5
    assert [member['key'] for member in selection['members']] == [
        'comp7', 'comp3', 'comp5', 'comp4', 'comp6',
    ]
    assert all(member['color_delta'] is not None for member in selection['members'])
    assert all(member['magnitude_delta'] is not None for member in selection['members'])
    limited_keys = {
        rejected['key']
        for rejected in selection['rejected']
        if 'closest to the target' in rejected['reason']
    }
    assert limited_keys == {'comp1', 'comp2'}

    expanded_selection = exotic_module.select_stellar_variability_ensemble_members(
        ranked_summaries,
        calibration_stars,
        comp_flux_map,
        observed_filter='V',
        target_catalog_match=target_match,
        max_members=7,
    )

    assert expanded_selection['member_limit'] == 7
    assert len(expanded_selection['members']) == 7
    assert not any(
        'closest to the target' in rejected['reason']
        for rejected in expanded_selection['rejected']
    )


def test_stellar_variability_ensemble_uses_gaia_bp_rp_when_local_colors_are_missing(
        monkeypatch):
    calls = []

    def fake_gaia_lookup(ra, dec, max_separation_arcsec):
        calls.append((ra, dec, max_separation_arcsec))
        return {
            'color': 1.2 if ra == 10.0 else 1.25,
            'label': 'BP-RP',
        }

    monkeypatch.setattr(exotic_module, 'nextastro_gaia_bp_rp_for_coordinate', fake_gaia_lookup)
    selection = exotic_module.select_stellar_variability_ensemble_members(
        [{
            'key': 'comp1',
            'comp_index': 0,
            'label': 'Comp 1',
            'position': [10, 20],
            'overexposure_rejected_count': 0,
        }],
        {
            'NextAstro-1': {
                'pos': [10, 20],
                'mag': 12.1,
                'error': 0.01,
                'mag_band': 'V',
                'catalog_row': {'ra': 11.0, 'dec': 20.0, 'Vmag': 12.1},
            },
        },
        {'comp1': np.full(8, 1000.0)},
        observed_filter='V',
        target_catalog_match={
            'mag': 12.0,
            'error': 0.01,
            'mag_band': 'V',
            'catalog_row': {'ra': 10.0, 'dec': 20.0, 'Vmag': 12.0},
        },
        min_members=1,
        max_members=1,
    )

    assert selection['target_catalog_profile']['color_label'] == 'BP-RP'
    assert selection['target_catalog_profile']['color'] == pytest.approx(1.2)
    assert selection['members'][0]['color_label'] == 'BP-RP'
    assert selection['members'][0]['color_delta'] == pytest.approx(0.05)
    assert calls == [(10.0, 20.0, 2.0), (11.0, 20.0, 2.0)]


def test_stellar_variability_ensemble_rejects_duplicate_catalog_sources():
    frame_count = 8
    ranked_summaries = [
        {
            'key': 'comp1', 'comp_index': 0, 'label': 'Comp 1', 'position': [10, 20],
            'overexposure_rejected_count': 0,
        },
        {
            'key': 'comp2', 'comp_index': 1, 'label': 'Comp 2', 'position': [30, 40],
            'overexposure_rejected_count': 0,
        },
        {
            'key': 'comp3', 'comp_index': 2, 'label': 'Comp 3', 'position': [50, 60],
            'overexposure_rejected_count': 0,
        },
    ]
    calibration_stars = {
        'NextAstro-111': {
            'pos': [10, 20], 'mag': 12.0, 'error': 0.01, 'mag_band': 'V', 'source_id': 111,
        },
        'NextAstro-111-2': {
            'pos': [30, 40], 'mag': 12.0, 'error': 0.01, 'mag_band': 'V', 'source_id': 111,
        },
        'NextAstro-222': {
            'pos': [50, 60], 'mag': 12.5, 'error': 0.01, 'mag_band': 'V', 'source_id': 222,
        },
    }
    comp_flux_map = {
        'comp1': np.full(frame_count, 3000.0),
        'comp2': np.full(frame_count, 2000.0),
        'comp3': np.full(frame_count, 1000.0),
    }

    selection = exotic_module.select_stellar_variability_ensemble_members(
        ranked_summaries,
        calibration_stars,
        comp_flux_map,
        observed_filter='V',
    )

    assert [member['key'] for member in selection['members']] == ['comp1', 'comp3']
    duplicate = next(item for item in selection['rejected'] if item['key'] == 'comp2')
    assert 'duplicate catalog source' in duplicate['reason']


def test_discover_fortuitous_vsx_variables_filters_on_count_rate_error_and_classifies(monkeypatch):
    class FakeWcs:
        def pixel_to_world_values(self, x_value, y_value):
            return x_value, y_value

        def world_to_pixel_values(self, ra_value, dec_value):
            return ra_value, dec_value

    reference_image = np.zeros((80, 80), dtype=float)
    reference_image[20, 20] = 10000.0
    reference_image[40, 40] = 100.0
    monkeypatch.setattr(exotic_module, 'search_wcs', lambda _path: FakeWcs())
    monkeypatch.setattr(
        exotic_module,
        'vsx_field_query',
        lambda *args, **kwargs: [
            {
                'Name': 'Bright VSX', 'AUID': '000-AAA-001',
                'RA2000': 20.0, 'Declination2000': 20.0,
                'Period': '5.0', 'MaxMag': '12.0 V', 'MinMag': '12.5 V',
            },
            {
                'Name': 'Faint VSX', 'AUID': '000-AAA-002',
                'RA2000': 40.0, 'Declination2000': 40.0,
                'Period': '20.0', 'MaxMag': '15.0 V', 'MinMag': '15.2 V',
            },
        ],
    )

    variables = exotic_module.discover_fortuitous_vsx_variables(
        'synthetic.wcs',
        reference_image.shape,
        1.0,
        reference_image,
        'V',
        target_pixel=[60, 60],
        exposure_seconds=60.0,
    )

    assert len(variables) == 1
    assert variables[0]['name'] == 'Bright VSX'
    assert variables[0]['estimated_magnitude_error'] < 0.05
    assert variables[0]['category'] == 'optimal_variables'
    assert variables[0]['period_days'] == pytest.approx(5.0)
    assert variables[0]['amplitude_mag'] == pytest.approx(0.5)
    assert exotic_module.fortuitous_variable_category(20.0, 0.2) == 'normal'


def test_fortuitous_reference_error_estimate_includes_sky_noise(monkeypatch):
    reference_image = np.full((80, 80), 1000.0, dtype=float)
    reference_image[40, 40] += 5000.0
    monkeypatch.setattr(
        exotic_module,
        'skybg_phot',
        lambda *args, **kwargs: (1000.0, 100.0, 500.0),
    )

    estimate = exotic_module.estimated_magnitude_error_from_reference_count_rate(
        reference_image,
        40.0,
        40.0,
        exposure_seconds=60.0,
        gain_e_per_adu=1.0,
    )

    source_only_error = (
        (2.5 / np.log(10.0))
        * exotic_module.source_flux_uncertainty_from_counts(5000.0, gain_e_per_adu=1.0)
        / 5000.0
    )
    assert source_only_error < 0.05
    assert estimate['estimated_magnitude_error'] > 0.05
    assert estimate['reference_noise_components_adu']['sky_aperture'] > 0
    assert estimate['reference_noise_components_adu']['sky_estimate'] > 0


def test_calibrated_stellar_variability_ensemble_combines_catalog_zero_points():
    frame_count = 7
    target_flux = np.full(frame_count, 1000.0)
    target_error = np.full(frame_count, 1.0)
    comp_flux_map = {
        'comp1': np.full(frame_count, 500.0),
        'comp2': np.full(frame_count, 250.0),
    }
    comp_error_map = {
        'comp1': np.full(frame_count, 1.0),
        'comp2': np.full(frame_count, 1.0),
    }
    members = [
        {
            'key': 'comp1',
            'magnitude': 12.0,
            'magnitude_error': 0.01,
            'summary': {'ensemble_frame_keep_mask': np.ones(frame_count, dtype=bool)},
        },
        {
            'key': 'comp2',
            'magnitude': 12.0 + 2.5 * np.log10(2.0),
            'magnitude_error': 0.01,
            'summary': {'ensemble_frame_keep_mask': np.ones(frame_count, dtype=bool)},
        },
    ]

    result = exotic_module.build_stellar_variability_calibrated_ensemble_series(
        target_flux,
        target_error,
        comp_flux_map,
        comp_error_map,
        members,
    )

    expected_target_magnitude = 12.0 - 2.5 * np.log10(2.0)
    assert result['applied'] is True
    np.testing.assert_allclose(result['magnitude'], expected_target_magnitude, atol=1.0e-10)
    np.testing.assert_allclose(result['relative_flux'], 1.0, atol=1.0e-10)
    np.testing.assert_allclose(result['raw_reference_flux'], 375.0, atol=1.0e-10)
    assert np.all(np.isfinite(result['raw_reference_flux_error']))
    np.testing.assert_array_equal(result['valid_member_count'], np.full(frame_count, 2))
    assert np.all(result['magnitude_error'] > 0)


def test_calibrated_stellar_variability_ensemble_rejects_frame_missing_any_member():
    frame_count = 7
    target_flux = np.full(frame_count, 1000.0)
    target_error = np.full(frame_count, 1.0)
    comp_flux_map = {
        'comp1': np.full(frame_count, 500.0),
        'comp2': np.full(frame_count, 250.0),
        'comp3': np.full(frame_count, 400.0),
    }
    comp_error_map = {
        key: np.full(frame_count, 1.0)
        for key in comp_flux_map
    }
    members = [
        {
            'key': 'comp1',
            'magnitude': 12.0,
            'magnitude_error': 0.01,
            'summary': {'ensemble_frame_keep_mask': np.ones(frame_count, dtype=bool)},
        },
        {
            'key': 'comp2',
            'magnitude': 12.0 + 2.5 * np.log10(2.0),
            'magnitude_error': 0.01,
            'summary': {'ensemble_frame_keep_mask': np.ones(frame_count, dtype=bool)},
        },
        {
            'key': 'comp3',
            'magnitude': 12.0 + 2.5 * np.log10(1.25),
            'magnitude_error': 0.01,
            'summary': {
                'ensemble_frame_keep_mask': np.array(
                    [True, True, True, False, True, True, True],
                    dtype=bool,
                ),
            },
        },
    ]

    result = exotic_module.build_stellar_variability_calibrated_ensemble_series(
        target_flux,
        target_error,
        comp_flux_map,
        comp_error_map,
        members,
        minimum_members=2,
    )

    assert result['applied'] is True
    assert result['valid_member_count'][3] == 2
    assert np.isnan(result['magnitude'][3])
    assert np.isnan(result['relative_flux'][3])
    finite_indices = np.flatnonzero(np.isfinite(result['magnitude']))
    np.testing.assert_array_equal(finite_indices, [0, 1, 2, 4, 5, 6])


def test_build_stellar_variability_ensemble_params_preserves_member_metadata(monkeypatch, tmp_path):
    captured = {}
    monkeypatch.setattr(
        exotic_module,
        'plot_stellar_variability',
        lambda params, save, target, label: captured.update(
            params=params,
            save=save,
            target=target,
            label=label,
        ),
    )
    fit = types.SimpleNamespace(
        time=np.array([2460000.105, 2460000.205]),
        jd_times=np.array([2460000.1, 2460000.2]),
        data=np.array([1.0, 1.1]),
        dataerr=np.array([0.01, 0.011]),
        airmass=np.array([1.1, 1.2]),
        airmass_model=np.ones(2),
        transit=np.ones(2),
        stellar_variability_only=True,
        stellar_variability_target_flux=np.array([500.0, 550.0]),
        stellar_variability_comp_flux=np.full(2, 1000.0),
        stellar_variability_target_flux_error=np.full(2, 2.0),
        stellar_variability_comp_flux_error=np.full(2, 3.0),
        stellar_variability_ensemble_magnitudes=np.array([12.30, 12.31]),
        stellar_variability_ensemble_magnitude_errors=np.array([0.01, 0.011]),
        stellar_variability_ensemble_members=[
            {
                'label': 'C1', 'position': [10, 20], 'magnitude': 12.0,
                'magnitude_error': 0.01,
                'star': {'catalog_source': 'Catalog A', 'ra': 10.1, 'dec': -20.1},
            },
            {
                'label': 'C2', 'position': [30, 40], 'magnitude': 12.5,
                'magnitude_error': 0.011,
                'star': {'catalog_source': 'Catalog B', 'ra': 10.2, 'dec': -20.2},
            },
        ],
    )

    params = exotic_module.build_stellar_variability_ensemble_params_from_fit(
        fit,
        tmp_path,
        'Target Star',
        observed_filter='MObs CV',
    )

    assert len(params) == 2
    assert [row['time'] for row in params] == pytest.approx([2460000.105, 2460000.205])
    assert [row['jd_time'] for row in params] == pytest.approx([2460000.1, 2460000.2])
    assert params[0]['cname'] == 'ENSEMBLE (2 stars)'
    assert params[0]['cmag'] is None
    assert params[0]['mag_band'] == 'ClearV'
    assert params[0]['catalog_mag_band'] == 'V'
    assert params[0]['differential_mag'] == pytest.approx(-2.5 * np.log10(0.5))
    assert params[1]['differential_mag'] == pytest.approx(-2.5 * np.log10(0.55))
    assert params[0]['differential_mag_err'] > 0
    assert params[0]['ensemble_member_labels'] == ['C1', 'C2']
    assert params[0]['ensemble_member_catalog_errors'] == [0.01, 0.011]
    assert params[0]['ensemble_member_ra_degs'] == [10.1, 10.2]
    assert params[0]['ensemble_member_dec_degs'] == [-20.1, -20.2]
    assert params[0]['ensemble_members'][0]['ra_deg'] == pytest.approx(10.1)
    assert params[0]['ensemble_members'][1]['dec_deg'] == pytest.approx(-20.2)
    assert fit.stellar_variability_params == params
    assert captured['label'] == 'ENSEMBLE (2 stars)'

    r_params = exotic_module.build_stellar_variability_ensemble_params_from_fit(
        fit,
        tmp_path,
        'Target Star',
        observed_filter='R',
    )

    assert r_params[0]['mag_band'] == 'rp'
    assert r_params[0]['catalog_mag_band'] == 'r'


def test_stellar_variability_ensemble_selection_json_lists_color_and_magnitude(monkeypatch, tmp_path):
    monkeypatch.setattr(exotic_module, 'plot_stellar_variability', lambda *args, **kwargs: None)
    member = {
        'selection_rank': 1,
        'key': 'comp1',
        'label': 'C1',
        'position': [10, 20],
        'magnitude': 12.1,
        'magnitude_error': 0.01,
        'color': 0.55,
        'color_label': 'B-V',
        'target_color': 0.50,
        'target_color_label': 'B-V',
        'color_delta': 0.05,
        'target_magnitude': 12.0,
        'magnitude_delta': 0.1,
        'color_magnitude_similarity_score': np.hypot(0.05, 0.1),
        'median_flux': 5000.0,
        'star': {
            'ra': 10.1,
            'dec': -20.2,
            'mag_band': 'V',
            'catalog_source': 'Synthetic catalog',
        },
    }
    fit = types.SimpleNamespace(
        jd_times=np.array([2460000.1]),
        airmass=np.array([1.1]),
        stellar_variability_ensemble_magnitudes=np.array([12.3]),
        stellar_variability_ensemble_magnitude_errors=np.array([0.02]),
        stellar_variability_ensemble_members=[member],
        stellar_variability_ensemble_selection={
            'members': [member],
            'rejected': [{'key': 'comp2', 'reason': 'not among the 5 closest'}],
            'member_limit': 5,
            'prelimit_member_count': 6,
            'calibration_error_clip': {'high_threshold': 0.03},
            'target_catalog_profile': {
                'magnitude': 12.0,
                'magnitude_band': 'V',
                'color': 0.5,
                'color_label': 'B-V',
            },
        },
    )

    exotic_module.build_stellar_variability_ensemble_params_from_fit(
        fit,
        tmp_path,
        'Target Star',
        observed_filter='V',
        observation_date='2024-01-02',
    )

    output_path = next(tmp_path.glob('EnsembleSelection_TargetStar_2024-01-02.json'))
    payload = json.loads(output_path.read_text(encoding='utf-8'))
    assert payload['ensemble']['maximum_members'] == 5
    assert payload['ensemble']['member_count_before_five_star_limit'] == 6
    assert payload['ensemble']['per_frame_required_members'] == 1
    assert 'every selected ensemble member is valid' in payload['ensemble']['selection_rule']
    assert payload['target']['catalog_profile']['color'] == pytest.approx(0.5)
    assert payload['ensemble']['members'][0]['color_delta'] == pytest.approx(0.05)
    assert payload['ensemble']['members'][0]['magnitude_delta'] == pytest.approx(0.1)


def test_process_fortuitous_variables_write_independent_and_combined_aid_products(monkeypatch, tmp_path):
    monkeypatch.setattr(exotic_module, 'plot_stellar_variability', lambda *args, **kwargs: None)
    monkeypatch.setattr(
        exotic_module,
        'psf_quality_mask_for_key',
        lambda psf_data, key, frame_count, psf_flux_data=None: np.ones(frame_count, dtype=bool),
    )
    frame_count = 12
    times = np.linspace(2460000.0, 2460000.1, frame_count)
    quality_mask = np.ones(frame_count, dtype=bool)
    comparison_calibration = {
        'method': 'aperture',
        'method_label': 'Aperture photometry',
        'a': 0,
        'an': 0,
        'field_image_keep_mask': quality_mask,
        'comp_summaries': [
            {
                'key': 'comp1', 'comp_index': 0, 'label': 'Comp 1', 'position': [10, 20],
                'aggregate_score': 0.001, 'coverage_rejected': False,
                'suitability_outlier_rejected': False, 'overexposure_rejected_count': 0,
                'psf_quality_keep_mask': quality_mask,
                'ensemble_frame_keep_mask': quality_mask,
            },
            {
                'key': 'comp2', 'comp_index': 1, 'label': 'Comp 2', 'position': [30, 40],
                'aggregate_score': 0.002, 'coverage_rejected': False,
                'suitability_outlier_rejected': False, 'overexposure_rejected_count': 0,
                'psf_quality_keep_mask': quality_mask,
                'ensemble_frame_keep_mask': quality_mask,
            },
        ],
    }
    psf_data = {
        # The exoplanet target is unusable in every frame. Fortuitous-variable
        # processing must remain independent of that target-specific mask.
        'target': np.full((frame_count, 7), np.nan),
        'comp1': np.ones((frame_count, 7)),
        'comp2': np.ones((frame_count, 7)),
        'comp3': np.ones((frame_count, 7)),
        'comp4': np.ones((frame_count, 7)),
    }
    aper_data = {
        'target': np.full((frame_count, 1, 1), 1000.0),
        'comp1': np.full((frame_count, 1, 1), 500.0),
        'comp1_unc': np.full((frame_count, 1, 1), 1.0),
        'comp2': np.full((frame_count, 1, 1), 250.0),
        'comp2_unc': np.full((frame_count, 1, 1), 1.0),
        'comp3': (
            800.0 * (1.0 + 0.02 * np.sin(np.linspace(0, 2 * np.pi, frame_count)))
        )[:, None, None],
        'comp3_unc': np.full((frame_count, 1, 1), 1.0),
        'comp4': (
            700.0 * (1.0 + 0.01 * np.cos(np.linspace(0, 2 * np.pi, frame_count)))
        )[:, None, None],
        'comp4_unc': np.full((frame_count, 1, 1), 1.0),
    }
    # One otherwise valid frame has an internal target error above 0.05 mag.
    aper_data['comp3_unc'][2, 0, 0] = 80.0
    calibrations = {
        'C1': {
            'pos': [10, 20], 'mag': 12.0, 'error': 0.01, 'mag_band': 'V',
            'ra': 10.1, 'dec': -20.1,
            'catalog_source': 'Synthetic catalog',
            'catalog_row': {'Bmag': 12.5, 'Vmag': 12.0},
        },
        'C2': {
            'pos': [30, 40], 'mag': 12.75, 'error': 0.011, 'mag_band': 'V',
            'ra': 10.2, 'dec': -20.2,
            'catalog_source': 'Synthetic catalog',
            'catalog_row': {'Bmag': 13.35, 'Vmag': 12.75},
        },
    }
    variable = {
        'name': 'Synthetic VSX',
        'auid': '000-AAA-001',
        'variable_type': 'EA',
        'period_days': 5.0,
        'amplitude_mag': 0.5,
        'category': 'optimal_variables',
        'ra': 10.0,
        'dec': -20.0,
        'pos': [50, 60],
        'tracking_key': 'comp3',
        'aperture_flux_adu': 800.0,
        'count_rate_adu_per_second': 13.3,
        'estimated_magnitude_error': 0.04,
        'catalog_match': {
            'mag': 12.4,
            'error': 0.02,
            'mag_band': 'V',
            'catalog_row': {'Bmag': 12.95, 'Vmag': 12.4},
        },
    }
    second_variable = {
        **variable,
        'name': 'Synthetic VSX 2',
        'auid': '000-AAA-002',
        'pos': [70, 80],
        'tracking_key': 'comp4',
        'aperture_flux_adu': 700.0,
        'count_rate_adu_per_second': 11.7,
    }
    info_dict = {
        'save': str(tmp_path),
        'date': '2024-01-02',
        'aavso_num': 'RTZ',
        'camera': 'CCD',
        'filter': 'V',
        'lat': '+32.4',
        'long': '-110.7',
        'elev': 2600,
    }

    variable_overexposed = np.zeros(frame_count, dtype=bool)
    variable_overexposed[:2] = True
    results = exotic_module.process_fortuitous_variables(
        [variable, second_variable],
        comparison_calibration,
        calibrations,
        times,
        times,
        np.linspace(1.1, 1.3, frame_count),
        psf_data,
        aper_data,
        info_dict,
        comp_overexposed_masks={'comp3': variable_overexposed},
        exposure_times_seconds=np.full(frame_count, 60.0),
        observed_filter='V',
        use_single_comparison=False,
        maximum_number_of_ensemble_comparisons_for_stellar_variability=2,
    )

    assert results[0]['status'] == 'completed'
    assert results[0]['input_frame_count'] == frame_count
    assert results[0]['target_overexposure_rejected_frame_count'] == 2
    assert results[0]['output_magnitude_error_rejected_frame_count'] == 1
    assert results[0]['point_count'] == frame_count - 3
    variable_dir = tmp_path / 'variables' / 'optimal_variables' / 'SyntheticVSX'
    assert next(
        (variable_dir / 'AAVSO_Files').glob('AID_AAVSO_SyntheticVSX_2024-01-02.txt')
    ).is_file()
    second_variable_dir = (
        tmp_path / 'variables' / 'optimal_variables' / 'SyntheticVSX2'
    )
    assert next(
        (second_variable_dir / 'AAVSO_Files').glob('AID_AAVSO_SyntheticVSX2_2024-01-02.txt')
    ).is_file()
    assert next(variable_dir.glob('EnsembleSelection_SyntheticVSX_2024-01-02.json')).is_file()
    assert next(variable_dir.glob('StellarVariability_SyntheticVSX_2024-01-02.csv')).is_file()
    combined_aid_path = (
        tmp_path / 'variables' / 'AAVSO_Files' / 'AID_AAVSO_FortuitousVariables_2024-01-02.txt'
    )
    combined_aid_text = combined_aid_path.read_text(encoding='utf-8')
    combined_aid_rows = [
        line for line in combined_aid_text.splitlines()
        if line and not line.startswith('#')
    ]
    assert combined_aid_text.count('#TYPE=EXTENDED') == 1
    assert '#ENSEMBLE-COMPARISONS-XC=' not in combined_aid_text
    assert len(combined_aid_rows) == sum(result['point_count'] for result in results)
    assert {row.split(',', 1)[0] for row in combined_aid_rows} == {
        '000-AAA-001',
        '000-AAA-002',
    }
    manifest = json.loads(
        next((tmp_path / 'variables').glob('FortuitousVariables_2024-01-02.json')).read_text(
            encoding='utf-8'
        )
    )
    assert manifest['combined_aid'] == str(combined_aid_path)
    assert manifest['variables'][0]['ensemble_member_count'] == 2
    assert manifest['variables'][0]['comparison_gap_stability']['applied'] is False
    assert manifest['variables'][0]['comparison_gap_rejected_candidates'] == []
    assert manifest['variables'][0]['output_magnitude_error_rejected_frame_count'] == 1
    assert manifest['variables'][0]['output_magnitude_error_max'] > 0.05
    selection = json.loads(
        next(variable_dir.glob('EnsembleSelection_SyntheticVSX_2024-01-02.json')).read_text(
            encoding='utf-8'
        )
    )
    assert selection['target']['input_frame_count'] == frame_count
    assert selection['target']['target_overexposure_rejected_frame_count'] == 2
    assert selection['target']['output_magnitude_error_rejected_frame_count'] == 1
    assert selection['target']['valid_output_frame_count'] == frame_count - 3
    assert 'exoplanet target overexposure mask is not applied' in (
        selection['target']['saturation_rejection_scope']
    )
    aid_text = next(
        (variable_dir / 'AAVSO_Files').glob('AID_AAVSO_SyntheticVSX_2024-01-02.txt')
    ).read_text(encoding='utf-8')
    assert (
        '#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS,'
        'GROUP,CHART,NOTES'
    ) in aid_text
    assert '|DIFFMAG=' in aid_text
    assert '|DIFFERR=' in aid_text
    aid_data_row = next(line for line in aid_text.splitlines() if not line.startswith('#'))
    assert len(aid_data_row.split(',')) == 15
    assert aid_data_row.split(',')[-1].startswith('|DIFFMAG=')
    assert '|DIFFERR=' in aid_data_row.split(',')[-1]
    ensemble_header = next(
        line for line in aid_text.splitlines()
        if line.startswith('#ENSEMBLE-COMPARISONS-XC=')
    )
    ensemble_metadata = json.loads(ensemble_header.split('=', 1)[1])
    assert ensemble_metadata['member_count'] == 2
    assert ensemble_metadata['members'][0]['ra_deg'] == pytest.approx(10.1)
    assert ensemble_metadata['members'][0]['dec_deg'] == pytest.approx(-20.1)
    assert ensemble_metadata['members'][1]['ra_deg'] == pytest.approx(10.2)
    assert ensemble_metadata['members'][1]['dec_deg'] == pytest.approx(-20.2)
    csv_path = next(variable_dir.glob('StellarVariability_SyntheticVSX_2024-01-02.csv'))
    csv_text = csv_path.read_text(encoding='utf-8')
    assert 'Apparent Magnitude' in csv_text.splitlines()[0]
    assert 'Raw Differential Magnitude' in csv_text.splitlines()[0]
    exported_errors = [
        float(row.split(',')[3])
        for row in csv_text.splitlines()[1:]
        if row.strip()
    ]
    assert exported_errors
    assert max(exported_errors) < 0.05

    single_root = tmp_path / 'single'
    single_info = {**info_dict, 'save': str(single_root)}
    single_results = exotic_module.process_fortuitous_variables(
        [variable],
        comparison_calibration,
        calibrations,
        times,
        times - 0.005,
        np.linspace(1.1, 1.3, frame_count),
        psf_data,
        aper_data,
        single_info,
        comp_overexposed_masks={'comp3': variable_overexposed},
        exposure_times_seconds=np.full(frame_count, 60.0),
        observed_filter='V',
    )

    assert single_results[0]['status'] == 'completed'
    assert single_results[0]['reference_mode'] == 'single_comparison'
    assert single_results[0]['comparison_member_count'] == 1
    assert single_results[0]['comparison_label'] == 'C2'
    single_dir = (
        single_root / 'variables' / 'optimal_variables' / 'SyntheticVSX'
    )
    assert not list(single_dir.glob('EnsembleSelection_*.json'))
    single_csv = next(single_dir.glob('StellarVariability_SyntheticVSX_2024-01-02.csv'))
    single_csv_rows = [
        line for line in single_csv.read_text(encoding='utf-8').splitlines()[1:]
        if line.strip()
    ]
    assert float(single_csv_rows[0].split(',')[0]) == pytest.approx(times[3])
    single_aid = next(
        (single_dir / 'AAVSO_Files').glob('AID_AAVSO_SyntheticVSX_2024-01-02.txt')
    )
    single_aid_text = single_aid.read_text(encoding='utf-8')
    assert '#ENSEMBLE-COMPARISONS-XC=' not in single_aid_text
    single_aid_row = next(
        line for line in single_aid_text.splitlines()
        if line and not line.startswith('#')
    )
    assert single_aid_row.split(',')[7] == 'C2'

    failed_root = tmp_path / 'failed'
    failed_results = exotic_module.process_fortuitous_variables(
        [variable],
        comparison_calibration,
        {},
        times,
        times,
        np.linspace(1.1, 1.3, frame_count),
        psf_data,
        aper_data,
        {**info_dict, 'save': str(failed_root)},
        comp_overexposed_masks={'comp3': variable_overexposed},
        exposure_times_seconds=np.full(frame_count, 60.0),
        observed_filter='V',
    )

    assert failed_results[0]['status'] == 'completed'
    assert failed_results[0]['apparent_magnitude_point_count'] == 0
    assert failed_results[0]['apparent_magnitude_error']
    differential_only_dir = (
        failed_root / 'variables' / 'optimal_variables' / 'SyntheticVSX'
    )
    assert differential_only_dir.exists()
    assert next(differential_only_dir.glob('DifferentialMagnitude_*.csv')).is_file()
    assert not list(differential_only_dir.glob('StellarVariability_*.csv'))
    assert not list((differential_only_dir / 'AAVSO_Files').glob('AID_AAVSO_*.txt'))
    failed_manifest = json.loads(
        next((failed_root / 'variables').glob('FortuitousVariables_2024-01-02.json')).read_text(
            encoding='utf-8'
        )
    )
    assert failed_manifest['variables'][0]['status'] == 'completed'
    assert failed_manifest['variables'][0]['apparent_magnitude_point_count'] == 0
    assert failed_manifest['variables'][0]['differential_magnitude_csv']


def test_stellar_variability_selector_uses_calibrated_ensemble_by_default(monkeypatch, tmp_path):
    monkeypatch.setattr(exotic_module, 'plot_stellar_variability', lambda *args, **kwargs: None)
    logged = []
    monkeypatch.setattr(exotic_module, 'log_info', lambda message, **kwargs: logged.append(message))
    frame_count = 12
    times = np.linspace(10.2, 10.3, frame_count)
    target_flux = 1000.0 * (1.0 + np.linspace(-0.002, 0.002, frame_count))
    comp1_flux = np.full(frame_count, 500.0)
    comp2_flux = np.full(frame_count, 250.0)
    quality_mask = np.ones(frame_count, dtype=bool)
    comp_summaries = [
        {
            'key': 'comp1', 'comp_index': 0, 'label': 'Comp 1', 'position': [10, 20],
            'aggregate_score': 0.001, 'coverage_rejected': False,
            'suitability_outlier_rejected': False, 'overexposure_rejected_count': 0,
            'psf_quality_keep_mask': quality_mask,
            'ensemble_frame_keep_mask': quality_mask,
        },
        {
            'key': 'comp2', 'comp_index': 1, 'label': 'Comp 2', 'position': [30, 40],
            'aggregate_score': 0.002, 'coverage_rejected': False,
            'suitability_outlier_rejected': False, 'overexposure_rejected_count': 0,
            'psf_quality_keep_mask': quality_mask,
            'ensemble_frame_keep_mask': quality_mask,
        },
    ]
    comparison_calibration = {
        'method': 'aperture',
        'method_label': 'Aperture photometry (aper=5px, annulus=10px)',
        'a': 0,
        'an': 0,
        'aper': 5.0,
        'annulus': 10.0,
        'field_score': 0.0015,
        'field_image_keep_mask': quality_mask,
        'comp_summaries': comp_summaries,
    }
    psf_data = {
        'target': np.ones((frame_count, 7), dtype=float),
        'comp1': np.ones((frame_count, 7), dtype=float),
        'comp2': np.ones((frame_count, 7), dtype=float),
    }
    aper_data = {
        'target': target_flux[:, None, None],
        'target_unc': np.full((frame_count, 1, 1), 1.0),
        'comp1': comp1_flux[:, None, None],
        'comp1_unc': np.full((frame_count, 1, 1), 1.0),
        'comp2': comp2_flux[:, None, None],
        'comp2_unc': np.full((frame_count, 1, 1), 1.0),
    }
    calibration_stars = {
        'C1': {
            'pos': [10, 20], 'mag': 12.0, 'error': 0.01, 'mag_band': 'V',
            'catalog_source': 'Synthetic catalog',
        },
        'C2': {
            'pos': [30, 40], 'mag': 12.0 + 2.5 * np.log10(2.0),
            'error': 0.011, 'mag_band': 'V', 'catalog_source': 'Synthetic catalog',
        },
    }

    result = exotic_module.select_stellar_variability_only_photometry(
        times,
        times,
        np.ones(frame_count),
        _stellar_variability_only_planet_dict(),
        comparison_calibration,
        psf_data,
        aper_data,
        target_flux,
        use_ensemble_photometry=True,
        maximum_number_of_ensemble_comparisons_for_stellar_variability=2,
        calibration_stars=calibration_stars,
        observed_filter='V',
    )

    selected = result['selected_result']
    assert result['selection_metric'] == 'stellar_variability_ensemble'
    assert selected['comp_index'] is None
    assert selected['ensemble_member_keys'] == ['comp1', 'comp2']
    assert selected['fit'].stellar_variability_ensemble_members
    assert any(
        'Using a 2-star comparison ensemble for stellar-variability products only:' in message
        for message in logged
    )
    assert len(selected['fit'].stellar_variability_ensemble_magnitudes) == len(selected['fit'].time)
    np.testing.assert_allclose(
        selected['fit'].differential_magnitude_reference_flux,
        375.0,
        atol=1.0e-10,
    )
    differential_series = exotic_module.differential_magnitude_series_from_fit(
        selected['fit'],
        apply_airmass_correction=False,
    )
    expected_differential = -2.5 * np.log10(target_flux / 375.0)
    np.testing.assert_allclose(
        differential_series['magnitude'],
        expected_differential,
        atol=1.0e-10,
    )
    assert abs(float(np.nanmedian(differential_series['magnitude']))) > 0.5

    # Final output preparation may refresh the normalized fitting photometry.
    # That must never overwrite the separately retained raw ensemble reference.
    exotic_module.annotate_stellar_variability_raw_photometry(
        selected['fit'],
        target_flux,
        target_flux,
        target_flux_error=np.ones(frame_count),
        comp_flux_error=np.ones(frame_count),
    )
    differential_after_fit_refresh = exotic_module.differential_magnitude_series_from_fit(
        selected['fit'],
        apply_airmass_correction=False,
    )
    np.testing.assert_allclose(
        differential_after_fit_refresh['magnitude'],
        expected_differential,
        atol=1.0e-10,
    )

    vsp_params = exotic_module.build_stellar_variability_params_from_photometry_selection(
        result,
        calibration_stars,
        tmp_path,
        'Synthetic',
        observed_filter='V',
        observation_date='2024-01-02',
    )
    aid_path = exotic_module.AIDOutputFiles(
        selected['fit'],
        {'sName': 'Synthetic', 'pName': 'Synthetic b'},
        {
            'save': tmp_path,
            'date': '2024-01-02',
            'aavso_num': 'TEST',
            'camera': 'CCD',
            'lat': 0.0,
            'long': 0.0,
            'elev': 0.0,
            'filter': 'V',
        },
        'AUID-TEST',
        None,
        vsp_params,
    ).aavso()

    assert len(vsp_params) == frame_count
    np.testing.assert_allclose(
        [row['mag'] for row in vsp_params],
        selected['fit'].stellar_variability_ensemble_magnitudes,
        atol=1.0e-10,
    )
    np.testing.assert_allclose(
        [row['differential_mag'] for row in vsp_params],
        expected_differential,
        atol=1.0e-10,
    )
    assert aid_path.is_file()
    aid_text = aid_path.read_text(encoding='utf-8')
    assert '#ENSEMBLE-COMPARISONS-XC=' in aid_text
    assert 'AUID-TEST,' in aid_text


def test_stellar_variability_exact_comparisons_use_every_supplied_member_without_catalogue():
    frame_count = 8
    times = np.linspace(10.2, 10.3, frame_count)
    quality_mask = np.ones(frame_count, dtype=bool)
    target_flux = np.linspace(990.0, 1010.0, frame_count)
    comparison_calibration = {
        'method': 'aperture',
        'method_label': 'Aperture photometry',
        'a': 0,
        'an': 0,
        'field_score': np.inf,
        'field_image_keep_mask': quality_mask,
        'comp_summaries': [
            {
                'key': 'comp1', 'comp_index': 0, 'label': 'Comp 1', 'position': [10, 20],
                'aggregate_score': np.inf, 'coverage_rejected': True,
                'suitability_outlier_rejected': True,
                'psf_quality_keep_mask': quality_mask,
                'ensemble_frame_keep_mask': quality_mask,
            },
            {
                'key': 'comp2', 'comp_index': 1, 'label': 'Comp 2', 'position': [30, 40],
                'aggregate_score': np.inf, 'coverage_rejected': True,
                'suitability_outlier_rejected': True,
                'psf_quality_keep_mask': quality_mask,
                'ensemble_frame_keep_mask': quality_mask,
            },
        ],
    }
    psf_data = {
        'target': np.ones((frame_count, 7), dtype=float),
        'comp1': np.ones((frame_count, 7), dtype=float),
        'comp2': np.ones((frame_count, 7), dtype=float),
    }
    aper_data = {
        'target': target_flux[:, None, None],
        'target_unc': np.ones((frame_count, 1, 1)),
        'comp1': np.full((frame_count, 1, 1), 500.0),
        'comp1_unc': np.ones((frame_count, 1, 1)),
        'comp2': np.full((frame_count, 1, 1), 250.0),
        'comp2_unc': np.ones((frame_count, 1, 1)),
    }

    result = exotic_module.select_stellar_variability_only_photometry(
        times,
        times,
        np.linspace(1.0, 1.5, frame_count),
        _stellar_variability_only_planet_dict(),
        comparison_calibration,
        psf_data,
        aper_data,
        target_flux,
        use_ensemble_photometry=True,
        calibration_stars={},
        observed_filter='V',
        require_apparent_magnitudes=False,
        use_exactly_the_comps_provided=True,
    )

    selected = result['selected_result']
    assert result['selection_metric'] == 'exact_stellar_variability_ensemble'
    assert selected['ensemble_member_keys'] == ['comp1', 'comp2']
    assert [
        member['key'] for member in selected['fit'].stellar_variability_ensemble_members
    ] == ['comp1', 'comp2']
    assert np.all(np.isnan(selected['fit'].stellar_variability_ensemble_magnitudes))


def test_stellar_variability_selector_opt_out_restores_single_comp_selection():
    frame_count = 24
    times = np.linspace(10.2, 10.3, frame_count)
    quality_mask = np.ones(frame_count, dtype=bool)
    comparison_calibration = {
        'method': 'aperture',
        'method_label': 'Aperture photometry',
        'a': 0,
        'an': 0,
        'aper': 5.0,
        'annulus': 10.0,
        'field_score': 0.001,
        'field_image_keep_mask': quality_mask,
        'comp_summaries': [{
            'key': 'comp1', 'comp_index': 0, 'label': 'Comp 1', 'position': [10, 20],
            'aggregate_score': 0.001, 'coverage_rejected': False,
            'suitability_outlier_rejected': False, 'overexposure_rejected_count': 0,
            'psf_quality_keep_mask': quality_mask,
            'ensemble_frame_keep_mask': quality_mask,
        }],
    }
    target_flux = 1000.0 * (1.0 + np.linspace(-0.001, 0.001, frame_count))
    comp_flux = np.full(frame_count, 500.0)
    psf_data = {
        'target': np.ones((frame_count, 7), dtype=float),
        'comp1': np.ones((frame_count, 7), dtype=float),
    }
    aper_data = {
        'target': target_flux[:, None, None],
        'target_unc': np.full((frame_count, 1, 1), 1.0),
        'comp1': comp_flux[:, None, None],
        'comp1_unc': np.full((frame_count, 1, 1), 1.0),
    }

    result = exotic_module.select_stellar_variability_only_photometry(
        times,
        times,
        np.ones(frame_count),
        _stellar_variability_only_planet_dict(),
        comparison_calibration,
        psf_data,
        aper_data,
        target_flux,
        use_ensemble_photometry=False,
    )

    assert result['selection_metric'] == 'stellar_variability_scatter'
    assert result['selected_result']['comp_index'] == 0
