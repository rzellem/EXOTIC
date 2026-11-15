import json
from types import SimpleNamespace

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


def test_tap_query_preserves_empty_text_columns_and_numeric_epoch(monkeypatch):
    response = SimpleNamespace(
        status_code=200,
        headers={'content-type': 'text/csv'},
        text=('pl_name,hostname,pl_tsystemref,pl_refname,pl_pubdate,pl_tranmid\n'
              'WASP-50 b,,,,,2458411.841111111\n'),
    )
    monkeypatch.setattr('exotic.api.nea.requests.get', lambda *args, **kwargs: response)
    nea = NASAExoplanetArchive('WASP-50 b')

    table = nea._tap_query('https://example.org/?query=', {'select': '*', 'from': 'ps'})

    for column in ('hostname', 'pl_tsystemref', 'pl_refname', 'pl_pubdate'):
        assert table[column].dtype == object
        assert pandas.isna(table[column].iloc[0])
    assert pandas.api.types.is_numeric_dtype(table['pl_tranmid'])
    assert table['pl_tranmid'].iloc[0] == float('2458411.841111111')


@pytest.mark.parametrize('label_source', ['newest', 'older', 'default_and_newest', 'missing'])
def test_new_scrape_merges_time_standard_text_into_empty_default_column(
    monkeypatch, tmp_path, label_source,
):
    nea = NASAExoplanetArchive('WASP-50 b')
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'pl_names.json').write_text('{"wasp50b": "WASP-50 b"}', encoding='utf-8')
    epoch = float('2458411.841111111')
    default = {
        'pl_name': 'WASP-50 b', 'hostname': 'WASP-50', 'pl_pubdate': '2019-07',
        'pl_refname': 'Default paper', 'pl_tsystemref': None, 'pl_tranmid': None,
        'pl_tranmiderr1': .0001, 'pl_tranmiderr2': -.0001, 'pl_orbper': 1.9550962,
        'pl_ratror': .12, 'pl_ratrorerr1': .001, 'pl_ratrorerr2': -.001, 'pl_ratdor': 8.0,
    }
    newest = dict(default, pl_pubdate='2024-10', pl_refname='Newest paper')
    older = dict(default, pl_pubdate='2024-05', pl_refname='Older paper', pl_tranmid=epoch)
    if label_source in ('newest', 'default_and_newest'):
        newest.update(pl_tsystemref='BJD-TDB', pl_tranmid=epoch)
    elif label_source == 'older':
        older['pl_tsystemref'] = 'BJD-TDB'
    if label_source == 'default_and_newest':
        default['pl_tsystemref'] = 'BJD-TDB'
    csv_tables = iter([
        pandas.DataFrame([default]).to_csv(index=False, float_format='%.17g'),
        pandas.DataFrame([newest, older]).to_csv(index=False, float_format='%.17g'),
    ])

    def get_response(*args, **kwargs):
        return SimpleNamespace(status_code=200, headers={'content-type': 'text/csv'},
                               text=next(csv_tables))

    monkeypatch.setattr('exotic.api.nea.requests.get', get_response)

    planet, candidate, parameters = nea.planet_info()

    assert planet == 'WASP-50 b'
    assert candidate is False
    assert parameters['midT'] == epoch
    expected_standard = None if label_source == 'missing' else 'BJD-TDB'
    assert parameters['midTStandard'] == expected_standard
    stored = json.loads((tmp_path / 'eaConf.json').read_text())[0]
    assert stored['pl_tsystemref'] == expected_standard
    exported = json.loads(nea.planet_info(fancy=True))
    assert exported['Published Mid-Transit Time Standard'] == (expected_standard or 'UNKNOWN')


def test_legacy_epoch_standard_lookup_uses_matching_publication_and_preserves_precision(monkeypatch):
    nea = NASAExoplanetArchive('WASP-16 b')
    captured = {}
    def tap(base, query, dataframe=True):
        captured.update(query)
        assert dataframe is False
        return ('pl_tranmid,pl_orbper,pl_tsystemref,pl_refname\n'
                '2456355.79653,3.11860349,BJD-TDB,Kokori 2023\n'
                '2454584.42878,3.1186009,HJD,Lister 2009\n')
    monkeypatch.setattr(nea,'_tap_query',tap)
    result=nea.lookup_ephemeris_time_standard('WASP-16 b',2456355.79653,3.11860349)
    assert result=={'midT':2456355.79653,'midTStandard':'BJD-TDB','midTSource':'Kokori 2023'}
    assert 'pl_tsystemref' in captured['select']
    assert nea.lookup_ephemeris_time_standard('WASP-16 b',2456355.7,3.11860349) is None


def test_loaded_archive_parameters_resolve_and_preserve_epoch_metadata(monkeypatch):
    nea=NASAExoplanetArchive('WASP-16 b')
    data={'pl_name':'WASP-16 b','hostname':'WASP-16','ra':214.68301,'dec':-20.27544,
          'pl_ratror':.12,'pl_ratrorerr1':.001,'pl_ratrorerr2':-.001,
          'pl_ratdor':8.21,'pl_tranmid':2456355.79653,'pl_tsystemref':'BJD-TDB',
          'pl_refname':'Kokori 2023'}
    nea._get_params(data)
    monkeypatch.setattr(
        nea, '_tap_query',
        lambda *args, **kwargs: pytest.fail('loaded planetary parameters should resolve without another NASA query'),
    )
    assert nea.resolve_name() is True
    result=json.loads(nea.planet_info(fancy=True))
    assert result['Published Mid-Transit Time Standard']=='BJD-TDB'
    assert result['Published Mid-Transit Time Reference']=='Kokori 2023'
    assert result['Published Mid-Transit Time']==2456355.79653
    assert 'Published Mid-Transit Time (BJD-UTC)' not in result


def test_colab_lookup_resolves_nextastro_parameters_when_nasa_archive_unavailable(monkeypatch, capsys):
    nea = NASAExoplanetArchive('WASP-12 b')

    def raise_direct_failure(*args, **kwargs):
        raise requests.exceptions.RequestException('ipac unavailable')

    monkeypatch.setattr(nea, '_new_scrape', raise_direct_failure)
    monkeypatch.setattr(nea, '_tap_query', raise_direct_failure)

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
        output = capsys.readouterr().out
        assert "Falling back to the NextAstro NEA cache server." in output
        assert "Successfully loaded" not in output
        called['url'] = url
        called['params'] = params
        called['timeout'] = timeout
        return DummyResponse(payload)

    monkeypatch.setattr(requests, 'get', fake_get)

    planet_name, candidate, pl_dict = nea.planet_info()

    assert (
        "Successfully loaded cached NASA Exoplanet Archive parameters for WASP-12 b "
        "from the NextAstro NEA cache server!"
    ) in capsys.readouterr().out
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

    # Standard Colab validates the name and exports parameters after planet_info().
    assert nea.resolve_name() is True
    exported = json.loads(nea.planet_info(fancy=True))
    assert exported['Planet Name'] == 'WASP-12 b'
    assert exported['Orbital Period (days)'] == 1.09


@pytest.mark.parametrize('loaded_params', [None, {}, {'pName': 'WASP-16 b'}])
@pytest.mark.parametrize('found', [False, True])
def test_resolve_name_queries_nasa_without_matching_loaded_parameters(monkeypatch, loaded_params, found):
    nea = NASAExoplanetArchive('WASP-12 b')
    nea.pl_dict = loaded_params
    queries = []

    def tap(base_url, query):
        queries.append(query)
        return pandas.DataFrame({'pl_name': ['WASP-12 b'] if found else []})

    monkeypatch.setattr(nea, '_tap_query', tap)

    assert nea.resolve_name() is found
    assert len(queries) == 1
    assert queries[0]['where'] == "tran_flag = 1 and pl_name = 'WASP-12 b'"


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
