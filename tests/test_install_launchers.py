"""Run the launchers with intercepted Python calls, without installing packages."""

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
LAUNCHERS = ['run_exotic_linux.sh', 'run_exotic_mac.command', 'run_exotic_windows.bat']

# sitecustomize runs before -m/-c code. Unknown invocations exit too, so these
# subprocesses cannot accidentally install packages or open the real GUI.
INTERCEPT = '''\
import json
import os
import sys

original_args = sys.orig_argv[1:]
args = original_args[1:] if original_args[:1] == ['-P'] else original_args
kind = 'unknown'
if args[:2] == ['-m', 'pip']:
    kind = 'pip_check' if args[2:] == ['--version'] else 'install'
elif args[:2] == ['-m', 'ensurepip']:
    kind = 'ensurepip'
elif args[:1] == ['-c']:
    if 'sys.version_info' in args[1]:
        kind = 'version_check'
    elif 'from exotic.exotic_gui import main; main()' == args[1]:
        kind = 'gui'

with open(os.environ['LAUNCHER_TEST_TRACE'], 'a', encoding='utf-8') as stream:
    stream.write(json.dumps({'kind': kind, 'args': original_args, 'python': sys.executable}) + '\\n')

if kind == 'version_check':
    sys.version_info = tuple(map(int, os.environ['LAUNCHER_TEST_VERSION'].split('.')))
else:
    os._exit(int(os.environ.get('LAUNCHER_TEST_' + kind.upper(), '97')))
'''


@pytest.fixture(params=LAUNCHERS)
def run_launcher(request, tmp_path):
    launcher = ROOT / request.param
    if launcher.suffix == '.bat':
        if os.name != 'nt':
            pytest.skip('Windows batch launcher requires cmd.exe')
        command = [os.environ['COMSPEC'], '/d', '/c', str(launcher)]
    else:
        bash = (Path(os.environ.get('ProgramFiles', 'C:/Program Files')) / 'Git/bin/bash.exe'
                if os.name == 'nt' else None)
        bash = str(bash) if bash and bash.exists() else shutil.which('bash')
        if not bash:
            pytest.skip('Bash is not installed')
        command = [bash, launcher.as_posix()]

    shim = tmp_path / 'python interception'
    shim.mkdir()
    (shim / 'sitecustomize.py').write_text(INTERCEPT, encoding='utf-8')
    trace = tmp_path / 'calls.jsonl'
    working = tmp_path / 'unrelated working directory'
    working.mkdir()
    env = os.environ.copy()
    env.update({
        'PYTHONPATH': str(shim),
        'PYTHONDONTWRITEBYTECODE': '1',
        'LAUNCHER_TEST_TRACE': str(trace),
        'LAUNCHER_TEST_VERSION': '3.12.0',
        'LAUNCHER_TEST_PIP_CHECK': '0',
        'LAUNCHER_TEST_INSTALL': '0',
        'LAUNCHER_TEST_ENSUREPIP': '0',
        'LAUNCHER_TEST_GUI': '0',
    })
    if os.name == 'nt':
        # Exclude Windows Store aliases and unrelated interpreters from discovery.
        git = Path(os.environ.get('ProgramFiles', 'C:/Program Files')) / 'Git'
        env['PATH'] = os.pathsep.join(map(str, [
            Path(sys.executable).parent, git / 'bin', git / 'usr/bin',
            Path(os.environ['SystemRoot']) / 'System32',
        ]))
    else:
        env['PATH'] = str(Path(sys.executable).parent) + os.pathsep + env['PATH']

    def invoke(command_args, call_env):
        return subprocess.run(
            command_args, cwd=working, env=call_env, capture_output=True, text=True,
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0,
        )

    # Verify interception before executing a launcher that could call real pip.
    preflight = invoke([sys.executable, '-c', 'raise SystemExit(98)'], env)
    assert preflight.returncode == 97, preflight.stderr
    trace.unlink()

    def run(**overrides):
        call_env = dict(env, **{'LAUNCHER_TEST_' + key.upper(): str(value)
                               for key, value in overrides.items()})
        result = invoke(command, call_env)
        calls = [json.loads(line) for line in trace.read_text(encoding='utf-8').splitlines()]
        return result, calls, launcher

    return run


@pytest.mark.parametrize('version', ['3.12.0', '3.13.0', '3.14.0'])
def test_supported_python_installs_and_launches_with_same_interpreter(run_launcher, version):
    result, calls, launcher = run_launcher(version=version)
    assert result.returncode == 0, result.stdout + result.stderr
    assert [call['kind'] for call in calls] == ['version_check', 'pip_check', 'install', 'gui']
    assert {Path(call['python']).resolve() for call in calls} == {Path(sys.executable).resolve()}
    install = calls[2]['args']
    assert install[:4] == ['-m', 'pip', 'install', '--upgrade']
    assert install[-1] == 'exotic'
    assert calls[-1]['args'][0] == '-P'
    if launcher.suffix == '.bat':
        assert Path(install[install.index('--find-links') + 1]) == ROOT / 'wheelhouse'
        assert install[install.index('--only-binary') + 1] == 'ultranest'
        assert 'ultranest>=4.5.0' in install


def test_failed_install_stops_before_gui(run_launcher):
    result, calls, _ = run_launcher(install=23)
    assert result.returncode == 23
    assert [call['kind'] for call in calls] == ['version_check', 'pip_check', 'install']


def test_missing_pip_bootstraps_selected_interpreter(run_launcher):
    result, calls, _ = run_launcher(pip_check=1)
    assert result.returncode == 0, result.stdout + result.stderr
    assert [call['kind'] for call in calls] == [
        'version_check', 'pip_check', 'ensurepip', 'install', 'gui',
    ]
    assert {Path(call['python']).resolve() for call in calls} == {Path(sys.executable).resolve()}


def test_failed_pip_bootstrap_stops_installation(run_launcher):
    result, calls, _ = run_launcher(pip_check=1, ensurepip=17)
    assert result.returncode == 17
    assert [call['kind'] for call in calls] == ['version_check', 'pip_check', 'ensurepip']


def test_gui_failure_is_reported(run_launcher):
    result, calls, _ = run_launcher(gui=19)
    assert result.returncode == 19
    assert calls[-1]['kind'] == 'gui'


def test_python_below_project_minimum_does_not_install(run_launcher):
    result, calls, _ = run_launcher(version='3.11.9')
    assert result.returncode == 65
    assert calls
    assert all(call['kind'] == 'version_check' for call in calls)
    assert 'Python 3.12 or newer is required' in result.stdout
