"""Exercise portable trial inputs without running nested sampling."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from scripts import expanded_prior_real_data_trial as trial


@pytest.fixture
def fitted_trials(monkeypatch):
    calls = []

    def fit(**kwargs):
        calls.append(kwargs)
        warm = kwargs.get('warmstart_source') is not None
        return object(), {
            'rprs_stdev': 0.01,
            'rprs_median': 0.20,
            'logz': 0.0,
            'logzerr': 0.1,
            'warmstart_applied': warm,
            'warmstart_expanded_keys': ['rprs'] if warm else [],
        }

    monkeypatch.setattr(trial, 'run_fit', fit)
    return calls


def write_lightcurve(path, index):
    if trial.DATASETS[index]['format'] == 'normalized':
        path.write_text(
            'BJD,Norm Flux,Norm Err,AM\n'
            '2461120.680,1.01,0.002,1.2\n'
            '2461120.681,0.99,0.003,1.3\n', encoding='utf-8',
        )
    else:
        path.write_text(
            '# BJD,phase,flux,error,model,airmass\n'
            '2461120.680,0,1.01,0.002,1,1.2\n'
            '2461120.681,0,0.99,0.003,1,1.3\n', encoding='utf-8',
        )


@pytest.mark.parametrize('index', range(len(trial.DATASETS)))
def test_cli_runs_supplied_lightcurve_with_original_preset(tmp_path, fitted_trials, index):
    lightcurve = tmp_path / 'local light curve.csv'
    output = tmp_path / 'results' / 'trial.json'
    write_lightcurve(lightcurve, index)

    status = trial.main([
        '--dataset-path', str(index), str(lightcurve),
        '--output', str(output), '--seed', '1000', '--live-points', '57',
    ])

    assert status == 0
    payload = json.loads(output.read_text())
    assert len(payload['datasets']) == 1
    result = payload['datasets'][0]
    assert result['name'] == trial.DATASETS[index]['name']
    assert result['source_path'] == str(lightcurve)
    assert result['point_count'] == 2
    assert result['initial_rprs_bounds'] == trial.DATASETS[index]['initial_rprs_bounds']
    assert result['expanded_rprs_bounds'] == trial.DATASETS[index]['expanded_rprs_bounds']
    assert result['tmid_bounds'] == trial.DATASETS[index]['tmid_bounds']
    assert len(fitted_trials) == 3
    assert [call['seed'] for call in fitted_trials] == [1000 + index * 100 + i for i in range(3)]
    assert [call['bounds'] for call in fitted_trials] == [
        trial.DATASETS[index]['initial_rprs_bounds'],
        trial.DATASETS[index]['expanded_rprs_bounds'],
        trial.DATASETS[index]['expanded_rprs_bounds'],
    ]
    for call in fitted_trials:
        assert call['live_points'] == 57
        assert all(call['prior'][key] == value for key, value in trial.DATASETS[index]['prior'].items())
        np.testing.assert_allclose(call['flux'], [1.01, 0.99])
        np.testing.assert_allclose(call['error'], [0.002, 0.003])
        np.testing.assert_allclose(call['airmass'], [1.2, 1.3])
    assert 'warmstart_source' not in fitted_trials[0]
    assert fitted_trials[1]['warmstart_source'] is not None
    assert 'warmstart_source' not in fitted_trials[2]
    assert 'path' not in trial.DATASETS[index]


def test_repeated_paths_run_in_supplied_order(tmp_path, fitted_trials):
    paths = {index: tmp_path / f'curve-{index}.csv' for index in (2, 0)}
    for index, path in paths.items():
        write_lightcurve(path, index)
    output = tmp_path / 'trial.json'

    assert trial.main([
        '--dataset-path', '2', str(paths[2]),
        '--dataset-path', '0', str(paths[0]), '--output', str(output),
    ]) == 0

    payload = json.loads(output.read_text())
    assert [row['name'] for row in payload['datasets']] == [trial.DATASETS[i]['name'] for i in (2, 0)]
    assert len(fitted_trials) == 6


def test_explicit_subset_does_not_open_unselected_path(tmp_path, fitted_trials):
    lightcurve = tmp_path / 'selected.csv'
    output = tmp_path / 'trial.json'
    write_lightcurve(lightcurve, 1)

    assert trial.main([
        '--dataset-path', '0', str(tmp_path / 'not-selected-and-absent.csv'),
        '--dataset-path', '1', str(lightcurve), '--dataset-index', '1',
        '--output', str(output),
    ]) == 0

    payload = json.loads(output.read_text())
    assert [row['name'] for row in payload['datasets']] == [trial.DATASETS[1]['name']]
    assert len(fitted_trials) == 3


def test_selected_dataset_needs_an_explicit_path(tmp_path, capsys, fitted_trials):
    output = tmp_path / 'results' / 'trial.json'
    with pytest.raises(SystemExit) as exc:
        trial.main([
            '--dataset-path', '0', 'curve.csv', '--dataset-index', '1',
            '--output', str(output),
        ])

    assert exc.value.code == 2
    assert '--dataset-path 1 PATH' in capsys.readouterr().err
    assert not output.parent.exists()
    assert not fitted_trials


def test_loading_failure_is_reported_without_substituting_another_file(tmp_path, fitted_trials):
    with pytest.raises(FileNotFoundError):
        trial.main([
            '--dataset-path', '0', str(tmp_path / 'absent.csv'),
            '--output', str(tmp_path / 'trial.json'),
        ])
    assert not fitted_trials


def test_help_works_outside_repository_without_importing_the_fitter(tmp_path):
    script = Path(trial.__file__).resolve()
    result = subprocess.run(
        [sys.executable, str(script), '--help'], cwd=tmp_path,
        capture_output=True, text=True, check=True,
    )
    assert '--dataset-path INDEX PATH' in result.stdout
    assert 'WASP-43b full transit' in result.stdout
    assert 'TrES-5b full transit' in result.stdout
    assert 'KELT-20b one-sided partial transit' in result.stdout
