import pandas
import pytest
import requests

from exotic.api.nea import NASAExoplanetArchive, planet_name_lookup_candidates


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


@pytest.mark.parametrize(
    ("planet_name", "reason"),
    [
        ("TOI-3889.01", "the name ends with a decimal suffix"),
        ("TIC 123456789", "the name starts with 'TIC'"),
    ],
)
def test_new_scrape_auto_marks_candidate_like_names_without_prompt(monkeypatch, tmp_path, capsys,
                                                                   planet_name, reason):
    nea = NASAExoplanetArchive(planet_name)

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(nea, 'planet_names', lambda filename="pl_names.json": None)
    monkeypatch.setattr(nea, '_tap_query', lambda *args, **kwargs: pandas.DataFrame())
    monkeypatch.setattr(
        'builtins.input',
        lambda prompt: pytest.fail("interactive prompt should not run for candidate-like targets"),
    )

    resolved_name, candidate = nea._new_scrape()

    assert resolved_name == planet_name
    assert candidate is True

    output = capsys.readouterr().out
    assert f"Cannot find target ({planet_name}) in NASA Exoplanet Archive." in output
    assert f"Assuming {planet_name} is a planet candidate because {reason}." in output


def test_new_scrape_non_interactive_unknown_target_aborts_without_prompt(monkeypatch, tmp_path):
    planet_name = 'Definitely Not A Planet b'
    nea = NASAExoplanetArchive(planet_name, non_interactive=True)

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(nea, 'planet_names', lambda filename="pl_names.json": None)
    monkeypatch.setattr(nea, '_tap_query', lambda *args, **kwargs: pandas.DataFrame())
    monkeypatch.setattr(
        'builtins.input',
        lambda prompt: pytest.fail("non-interactive NASA lookup must not prompt"),
    )

    with pytest.raises(
        RuntimeError,
        match=r"Non-interactive run cancelled: target \(Definitely Not A Planet b\) was not found",
    ):
        nea._new_scrape()


def test_planet_name_lookup_candidates_strip_phase_and_preserve_planet_letter():
    candidates = planet_name_lookup_candidates('field WASP-164 b ingress')

    assert 'field WASP-164 b' in candidates
    assert 'field WASP-164b' in candidates
    assert 'WASP-164b' in candidates
    assert 'ingress' not in candidates


@pytest.mark.parametrize('phase', ['ingress', 'EGRESS'])
def test_new_scrape_resolves_phase_labeled_name_from_planet_cache(
    monkeypatch,
    tmp_path,
    phase,
):
    nea = NASAExoplanetArchive(f'WASP-164 b {phase}', non_interactive=True)
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'pl_names.json').write_text(
        '{"wasp164b": "WASP-164 b"}',
        encoding='utf-8',
    )

    class LookupResolved(Exception):
        pass

    def stop_after_name_resolution(*args, **kwargs):
        assert nea.planet == 'WASP-164 b'
        raise LookupResolved

    monkeypatch.setattr(nea, '_tap_query', stop_after_name_resolution)

    with pytest.raises(LookupResolved):
        nea._new_scrape()
