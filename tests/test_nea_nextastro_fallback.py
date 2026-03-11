import requests

from exotic.api.nea import NASAExoplanetArchive


class DummyResponse:
    def __init__(self, payload):
        self._payload = payload

    def raise_for_status(self):
        return None

    def json(self):
        return self._payload


def test_planet_info_uses_nextastro_fallback_when_nasa_archive_unavailable(monkeypatch):
    nea = NASAExoplanetArchive('WASP-12 b')

    def raise_direct_failure(*args, **kwargs):
        raise requests.exceptions.RequestException('ipac unavailable')

    monkeypatch.setattr(nea, '_new_scrape', raise_direct_failure)

    payload = {
        'params': {
            'name': 'WASP-12 b',
            'hostStarName': 'WASP-12',
            'raDeg': 180.0,
            'decDeg': 29.0,
            'orbitalPeriodDays': {'value': 1.09, 'errPlus': 0.001, 'errMinus': 0.001},
            'midTransitTimeDays': {'value': 2450000.5, 'errPlus': 0.0002, 'errMinus': 0.0003},
            'rpOverRs': {'value': 0.12, 'errPlus': 0.004, 'errMinus': 0.003},
            'aOverRs': {'value': 3.0, 'errPlus': 0.2, 'errMinus': 0.1},
            'inclinationDeg': {'value': 83.5, 'errPlus': 0.8, 'errMinus': 0.7},
            'eccentricity': 0.0,
            'argPeriastronDeg': 0.0,
            'starTeffK': {'value': 6300.0, 'errPlus': 50.0, 'errMinus': 40.0},
            'starFeh': {'value': 0.2, 'errPlus': 0.03, 'errMinus': 0.02},
            'starLogg': {'value': 4.1, 'errPlus': 0.05, 'errMinus': 0.04},
            'source': {'table': 'pscomppars', 'localCache': True},
        }
    }

    called = {}

    def fake_get(url, params, timeout):
        called['url'] = url
        called['params'] = params
        called['timeout'] = timeout
        return DummyResponse(payload)

    monkeypatch.setattr(requests, 'get', fake_get)

    planet_name, candidate, pl_dict = nea.planet_info()

    assert planet_name == 'WASP-12 b'
    assert candidate is False
    assert called['url'] == 'https://archive.nextastro.org/api/exoplanet_params'
    assert called['params'] == {'name': 'WASP-12 b'}
    assert pl_dict['pName'] == 'WASP-12 b'
    assert pl_dict['sName'] == 'WASP-12'
    assert pl_dict['ra'] == 180.0
    assert pl_dict['dec'] == 29.0
    assert pl_dict['pPer'] == 1.09
    assert pl_dict['midT'] == 2450000.5
    assert pl_dict['rprs'] == 0.12
    assert pl_dict['aRs'] == 3.0
