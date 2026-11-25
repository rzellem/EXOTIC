import ast
import json
from pathlib import Path

import pytest
import requests

from exotic import inputs as inputs_module
from exotic.inputs import Inputs, populate_missing_candidate_parameters


def response(parameters, **extra):
    class Response:
        def raise_for_status(self):
            pass
        def json(self):
            return {'parameters': parameters, 'sources': {}, 'missing': [], 'warnings': [], **extra}
    return Response()


@pytest.mark.parametrize('configured', [False, 'n', 'off', 'false', 0])
def test_unattended_runs_respect_explicit_disable_without_prompt(monkeypatch, configured):
    monkeypatch.setattr(inputs_module.requests, 'post', lambda *a, **k: pytest.fail('explicitly disabled'))
    monkeypatch.setattr(inputs_module, 'user_input', lambda *a, **k: pytest.fail('must not prompt'))
    assert populate_missing_candidate_parameters({}, {'use_nextastro_candidate_estimator': configured},
                                                non_interactive_run=True) == set()


@pytest.mark.parametrize('info', [{}, {'use_nextastro_candidate_estimator': None},
                                {'use_nextastro_candidate_estimator': True},
                                {'use_nextastro_candidate_estimator': 'y'}])
def test_unattended_candidate_auto_fills_without_prompt(monkeypatch, info):
    monkeypatch.setattr(inputs_module, 'user_input', lambda *a, **k: pytest.fail('unattended must never prompt'))
    calls = []
    monkeypatch.setattr(inputs_module.requests, 'post', lambda *a, **k: calls.append(k) or response({'aRs': 12.3}))
    parameters = {'pPer': 10}
    assert populate_missing_candidate_parameters(parameters, info, non_interactive_run=True) == {'aRs'}
    assert parameters['aRs'] == 12.3
    assert len(calls) == 1


def test_interactive_decline_makes_no_request(monkeypatch):
    monkeypatch.setattr(inputs_module, 'user_input', lambda *a, **k: 'n')
    monkeypatch.setattr(inputs_module.requests, 'post', lambda *a, **k: pytest.fail('declined'))
    assert populate_missing_candidate_parameters({}, {}) == set()


def test_partial_coordinate_is_sent_and_nonfinite_values_are_missing(monkeypatch):
    calls = []
    monkeypatch.setattr(inputs_module.requests, 'post', lambda *a, **k: calls.append(k) or response(
        {'teff': 5100, 'dec': -20}))
    parameters = {'ra': '08:00:00', 'dec': None, 'teff': float('nan')}
    filled = populate_missing_candidate_parameters(parameters, {'use_nextastro_candidate_estimator': True})
    assert calls[0]['json']['parameters']['ra'] == pytest.approx(120)
    assert calls[0]['json']['parameters']['dec'] is None
    assert calls[0]['json']['parameters']['teff'] is None
    assert parameters['ra'] == '08:00:00'
    assert parameters['teff'] == 5100
    assert filled == {'teff', 'dec'}


def test_interactive_accept_fills_missing_and_saves_sources(monkeypatch, tmp_path):
    monkeypatch.setattr(inputs_module, 'user_input', lambda *a, **k: 'y')
    calls = []
    monkeypatch.setattr(inputs_module.requests, 'post', lambda *a, **k: calls.append((a, k)) or response(
        {'teff': 5100, 'rprs': 0.12, 'met': -0.3, 'ecc': 0.2, 'midTStandard': 'BJD_TDB'},
        sources={'teff': {'kind': 'gaia', 'reference': 'Gaia DR3 42'}},
    ))
    parameters = {'pName': 'candidate', 'ra': '08:00:00', 'dec': '-20:00:00',
                  'rprs': 0.1, 'met': 0, 'ecc': 0, 'teff': None, 'midTStandard': 'JD_UTC'}
    info = {'save': str(tmp_path)}
    assert populate_missing_candidate_parameters(parameters, info) == {'teff'}
    assert parameters['teff'] == 5100
    assert parameters['rprs'] == 0.1 and parameters['met'] == 0 and parameters['ecc'] == 0
    assert parameters['midTStandard'] == 'JD_UTC'
    assert calls[0][0][0] == 'https://archive.nextastro.org/api/candidate_parameters'
    assert calls[0][1]['json']['parameters']['ra'] == pytest.approx(120)
    report = tmp_path / 'working_artifacts' / 'NextAstroCandidateParameters.json'
    assert json.loads(report.read_text())['sources']['teff']['kind'] == 'gaia'


@pytest.mark.parametrize('standard', ['BJD_TDB', None])
def test_new_epoch_replaces_old_timing_metadata(monkeypatch, standard):
    monkeypatch.setattr(inputs_module.requests, 'post', lambda *a, **k: response(
        {'midT': 2460000, 'midTStandard': standard, 'midTSource': 'TOI catalogue'}))
    parameters = {'midT': None, 'midTStandard': 'JD_UTC', 'midTSource': 'old epoch'}
    populate_missing_candidate_parameters(parameters, {'use_nextastro_candidate_estimator': True},
                                          non_interactive_run=True)
    assert parameters['midT'] == 2460000
    assert parameters['midTStandard'] == standard
    assert parameters['midTSource'] == 'TOI catalogue'


@pytest.mark.parametrize('error', [requests.HTTPError('503'), requests.Timeout('offline'), ValueError('bad JSON')])
def test_failed_request_is_reported_and_does_not_change_inputs(monkeypatch, error, capsys):
    def fail(*a, **k):
        raise error
    monkeypatch.setattr(inputs_module.requests, 'post', fail)
    parameters = {'pPer': 4, 'teff': None}
    assert populate_missing_candidate_parameters(parameters, {'use_nextastro_candidate_estimator': True}) == set()
    assert parameters == {'pPer': 4, 'teff': None}
    assert 'auto-fill failed' in capsys.readouterr().out


def test_invalid_response_does_not_apply_partial_values(monkeypatch):
    monkeypatch.setattr(inputs_module.requests, 'post', lambda *a, **k: response({'teff': 5100, 'logg': 'bad'}))
    parameters = {}
    assert populate_missing_candidate_parameters(parameters, {'use_nextastro_candidate_estimator': True}) == set()
    assert parameters == {}


def test_complete_inputs_make_no_request_or_prompt(monkeypatch):
    monkeypatch.setattr(inputs_module, 'user_input', lambda *a, **k: pytest.fail('already complete'))
    monkeypatch.setattr(inputs_module.requests, 'post', lambda *a, **k: pytest.fail('already complete'))
    parameters = dict.fromkeys(inputs_module.CANDIDATE_NUMERIC_KEYS, 1)
    parameters.update(ra=120, dec=-20, sName='star')
    assert populate_missing_candidate_parameters(parameters, {}) == set()


@pytest.mark.parametrize('enabled', [True, False, 'y', 'n'])
def test_inits_option_is_loaded(tmp_path, enabled):
    path = tmp_path / 'inits.json'
    path.write_text(json.dumps({'user_info': {}, 'planetary_parameters': {},
                               'optional_info': {'use_nextastro_candidate_estimator': enabled}}))
    inputs = Inputs(init_opt='y')
    inputs.comp_params(path, {})
    assert inputs.info_dict['use_nextastro_candidate_estimator'] == enabled


def test_supplied_depth_is_loaded_and_sent_without_adding_planetary_prompt_fields(tmp_path, monkeypatch):
    path = tmp_path / 'inits.json'
    path.write_text(json.dumps({'user_info': {}, 'planetary_parameters': {}, 'optional_info': {
        'use_nextastro_candidate_estimator': True, 'transit_depth': 10000,
        'transit_depth_unc': 100, 'transit_depth_unit': 'ppm'}}))
    inputs = Inputs(init_opt='y')
    parameters = inputs.comp_params(path, {})
    calls = []
    monkeypatch.setattr(inputs_module.requests, 'post', lambda *a, **k: calls.append(k) or response(
        {'rprs': 0.1, 'rprsUnc': 0.0005}, sources={'rprs': {'kind': 'derived', 'reference': 'sqrt(depth)'}}))
    populate_missing_candidate_parameters(parameters, inputs.info_dict, non_interactive_run=True)
    assert calls[0]['json']['parameters']['transit_depth'] == 10000
    assert calls[0]['json']['parameters']['transit_depth_unc'] == 100
    assert calls[0]['json']['parameters']['transit_depth_unit'] == 'ppm'
    assert parameters['rprs'] == 0.1
    assert 'transit_depth' not in parameters
    assert inputs.info_dict['candidate_parameter_estimates']['sources']['rprs']['kind'] == 'derived'


def test_inferred_mass_and_orbit_derivation_are_saved_for_review(monkeypatch, tmp_path):
    derivation = {'mass_solar': 0.8, 'mass_method': 'surface_gravity_radius',
                  'logg_dex': 4.4, 'radius_solar': 0.9,
                  'missing_uncertainties': ['logg_dex']}
    monkeypatch.setattr(inputs_module.requests, 'post', lambda *a, **k: response(
        {'aRs': 12.3, 'aRsUnc': 0.5}, derivations={'aRs': derivation},
        sources={'aRs': {'kind': 'derived', 'reference': 'Gaia gravity/radius + period'}}))
    parameters = {'pPer': 10}
    info = {'use_nextastro_candidate_estimator': True, 'save': str(tmp_path)}
    assert populate_missing_candidate_parameters(parameters, info, non_interactive_run=True) == {'aRs', 'aRsUnc'}
    assert parameters['aRs'] == 12.3
    report = json.loads((tmp_path / 'working_artifacts' / 'NextAstroCandidateParameters.json').read_text())
    assert report['derivations']['aRs'] == derivation
    assert 'mass_solar' not in parameters


def test_duration_is_loaded_and_sent_and_derived_inclination_is_applied(tmp_path, monkeypatch):
    path = tmp_path / 'inits.json'
    path.write_text(json.dumps({'user_info': {}, 'planetary_parameters': {}, 'optional_info': {
        'use_nextastro_candidate_estimator': True, 'transit_duration': 180,
        'transit_duration_unc': 5, 'transit_duration_unit': 'minutes', 'ecc_unc': 0.02, 'omega_unc': 3}}))
    inputs = Inputs(init_opt='y')
    parameters = inputs.comp_params(path, {})
    calls = []
    monkeypatch.setattr(inputs_module.requests, 'post', lambda *a, **k: calls.append(k) or response(
        {'inc': 87.3, 'incUnc': 0.4}, sources={'inc': {'kind': 'derived', 'reference': 'circular T14 geometry'}},
        derivations={'inc': {'model': 'circular_T14', 'duration_days': 0.125}}))
    populate_missing_candidate_parameters(parameters, inputs.info_dict, non_interactive_run=True)
    assert calls[0]['json']['parameters']['transit_duration'] == 180
    assert calls[0]['json']['parameters']['transit_duration_unc'] == 5
    assert calls[0]['json']['parameters']['transit_duration_unit'] == 'minutes'
    assert calls[0]['json']['parameters']['eccUnc'] == 0.02
    assert calls[0]['json']['parameters']['omegaUnc'] == 3
    assert parameters['inc'] == 87.3 and parameters['incUnc'] == 0.4
    assert 'transit_duration' not in parameters
    assert inputs.info_dict['candidate_parameter_estimates']['derivations']['inc']['model'] == 'circular_T14'


def test_candidate_confirmation_does_not_reprompt_for_accepted_estimates(monkeypatch):
    from exotic import exotic as reduction
    parameters = dict.fromkeys(inputs_module.CANDIDATE_NUMERIC_KEYS, 1)
    parameters.update(ra=120, dec=-20, pName='candidate', sName='star')
    monkeypatch.setattr(reduction, 'user_input', lambda *a, **k: pytest.fail('already accepted'))
    assert reduction.get_planetary_parameters(True, parameters, autofilled_keys=set(parameters)) == parameters


@pytest.mark.parametrize('candidate,override,option,unattended,expected', [
    (True, False, None, True, True), (False, True, True, True, True),
    (False, True, False, True, False), (False, False, True, True, False),
    (False, True, None, True, True), (False, True, None, False, False),
    (False, True, True, False, True), (False, False, None, True, False),
])
def test_main_wires_candidate_and_unattended_override_defaults(candidate, override, option, unattended, expected):
    # Execute the actual main-flow statement without starting a reduction.
    from types import SimpleNamespace
    from exotic import exotic as reduction
    tree = ast.parse(Path(reduction.__file__).read_text(encoding='utf-8'))
    matches = [node for node in ast.walk(tree) if isinstance(node, ast.If)
               and 'CandidatePlanetBool or' in ast.unparse(node.test)]
    statement = matches[0]
    calls = []
    namespace = dict(CandidatePlanetBool=candidate, args=SimpleNamespace(override=override, non_interactive_run=unattended),
        exotic_infoDict={'use_nextastro_candidate_estimator': option}, userpDict={},
        parse_bool_config_value=reduction.parse_bool_config_value,
        populate_missing_candidate_parameters=lambda *a, **k: calls.append(1) or set())
    exec(compile(ast.Module(body=[statement], type_ignores=[]), '<main candidate flow>', 'exec'), namespace)
    assert bool(calls) is expected
