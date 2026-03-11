import sys
import types

import numpy as np

fake_barycorrpy = types.ModuleType('barycorrpy')
fake_utc_tdb = types.ModuleType('barycorrpy.utc_tdb')
fake_utc_tdb.JDUTC_to_BJDTDB = lambda *args, **kwargs: None
fake_barycorrpy.utc_tdb = fake_utc_tdb
sys.modules.setdefault('barycorrpy', fake_barycorrpy)
sys.modules.setdefault('barycorrpy.utc_tdb', fake_utc_tdb)

from exotic import exotic as exotic_module


class DummyResponse:
    def __init__(self, payload, status_code=200):
        self._payload = payload
        self.status_code = status_code

    def json(self):
        return self._payload


def test_nextastro_variability_logs_json_request_and_response(monkeypatch):
    captured = {}
    logged = []

    def fake_post(url, json, timeout):
        captured['url'] = url
        captured['json'] = json
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
    assert captured['json'] == [{'ra': 10.1, 'dec': -11.2}, {'ra': 22.3, 'dec': -33.4}]
    assert variability_flags == [False, True]
    assert any('NextAstro variability request JSON:' in message for message in logged)
    assert any('NextAstro variability response JSON:' in message for message in logged)


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


def test_vsx_variable_falls_back_to_nextastro(monkeypatch):
    class DummyFailedResponse:
        def raise_for_status(self):
            raise RuntimeError('VSX unavailable')

    monkeypatch.setattr(exotic_module.requests, 'get', lambda *args, **kwargs: DummyFailedResponse())
    monkeypatch.setattr(exotic_module, 'nextastro_variability_test', lambda payload: [True])

    is_variable = exotic_module.vsx_variable(ra=10.0, dec=20.0)

    assert is_variable is True
