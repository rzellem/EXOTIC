"""Measure a full reduction and both candidate searches on identical photometry.

The modified pipeline runs normally through final products. At the candidate
boundary, save the small photometry inputs and RNG state. Then run the chosen
Git revision's candidate search against those same inputs and RNG state.
Existing Git commit IDs and the saved source/patch identify the software being
compared; file paths and modification times cannot distinguish code revisions.
Science images are processed once, without content hashes or readback QA.
"""
from __future__ import annotations

import argparse
import ast
import copy
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import pickle
import random
import subprocess
import sys
import time
import traceback
import warnings

os.environ.setdefault('MPLBACKEND', 'Agg')
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import astropy
import numpy as np
import psutil
import scipy
import ultranest

import exotic.exotic as runtime
from exotic.api import elca, ultranest_utils


def json_value(value):
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def cpu_seconds():
    current = psutil.Process().cpu_times()
    return current.user + current.system


def fit_summary(attempt):
    fit = attempt.get('fit')
    return {
        'label': attempt.get('label'), 'comp_index': attempt.get('comp_index'),
        'qc_status': attempt.get('transit_qc_status'),
        'qc_summary': attempt.get('transit_qc_summary'),
        'ktmf': attempt.get('ktmf_metric'),
        'delta_bic': attempt.get('transit_delta_bic'),
        'eebls_snr': attempt.get('eebls_snr'),
        'residual_scatter': attempt.get('residual_scatter'),
        'fit_point_count': attempt.get('fit_point_count'),
        'selection_reason': attempt.get('selection_reason'),
        'parameters': copy.deepcopy(getattr(fit, 'parameters', None)),
        'errors': copy.deepcopy(getattr(fit, 'errors', None)),
        'selected_despite_transit_qc': attempt.get('selected_despite_transit_qc', False),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inits', type=Path, required=True)
    parser.add_argument('--science-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline-ref', default='HEAD',
                        help='Git revision containing the original candidate-search function.')
    parser.add_argument('--seed', type=int, default=3853714786,
                        help="Use the seed recorded in Ivan's original run.")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[1]
    baseline_ref = subprocess.check_output(
        ['git', 'rev-parse', args.baseline_ref], cwd=repo, text=True).strip()
    current_ref = subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    baseline_source = subprocess.check_output(
        ['git', 'show', f'{baseline_ref}:exotic/exotic.py'], cwd=repo, text=True)
    (output / 'baseline_exotic.py').write_text(baseline_source, encoding='utf-8')
    (output / 'shortlist.patch').write_text(subprocess.check_output(
        ['git', 'diff', '--', 'exotic/exotic.py'], cwd=repo, text=True), encoding='utf-8')

    init = json.loads(args.inits.read_text(encoding='utf-8'))
    init['user_info']['Directory with FITS files'] = str(args.science_dir.resolve())
    init['user_info']['Directory to Save Plots'] = str(output / 'modified_full')
    init['user_info']['Directory of Darks'] = str(args.science_dir.resolve() / 'darks')
    init['optional_info']['Random Seed'] = args.seed
    init_path = output / 'inits_local.json'
    init_path.write_text(json.dumps(init, indent=2), encoding='utf-8')

    report = {
        'started_utc': datetime.now(timezone.utc).isoformat(),
        'method': 'One full modified reduction; paired candidate searches share measured photometry and initial RNG state.',
        'baseline_ref': baseline_ref,
        'current_ref': current_ref,
        'seed': args.seed,
        'environment': {
            'python': sys.version, 'numpy': np.__version__, 'scipy': scipy.__version__,
            'astropy': astropy.__version__, 'ultranest': ultranest.__version__,
            'logical_cpus': psutil.cpu_count(),
            'memory_available_bytes_at_start': psutil.virtual_memory().available,
            'ultranest_workers_configured': ultranest_utils._configured_ultranest_workers(),
            'ultranest_backend_configured': ultranest_utils._configured_ultranest_worker_backend(),
            'relevant_environment': {
                key: value for key, value in os.environ.items()
                if key.startswith(('OMP_', 'OPENBLAS_', 'MKL_')) or key in {
                    'EXOTIC_ULTRANEST_WORKERS', 'NEXTASTRO_EXOTIC_ULTRANEST_WORKERS',
                    'EXOTIC_ULTRANEST_WORKER_BACKEND', 'EXOTIC_ULTRANEST_CPU_GATE',
                    'EXOTIC_ULTRANEST_CPU_THRESHOLD',
                }
            },
        },
        'phase': 'modified_preprocessing', 'trials': {}, 'sampler_runs': [],
        'candidate_reductions': [], 'warnings': [],
    }

    def save():
        (output / 'benchmark.json').write_text(
            json.dumps(json_value(report), indent=2, allow_nan=False), encoding='utf-8')

    def event(message):
        print(f'BENCHMARK | {message}', flush=True)
        save()

    original_showwarning = warnings.showwarning

    def showwarning(message, category, filename, lineno, file=None, line=None):
        report['warnings'].append({
            'phase': report['phase'], 'candidate': report.get('candidate'),
            'category': category.__name__, 'message': str(message),
            'filename': filename, 'line': lineno,
        })
        original_showwarning(message, category, filename, lineno, file=file, line=line)

    warnings.showwarning = showwarning
    warnings.simplefilter('always', RuntimeWarning)
    original_sampler = elca.run_reactive_sampler

    def sampler_run(sampler, *sampler_args, **sampler_kwargs):
        record = {'phase': report['phase'], 'candidate': report.get('candidate'),
                  'parameters': getattr(sampler, 'paramnames', None)}
        report['sampler_runs'].append(record)
        started = time.perf_counter()
        try:
            result = original_sampler(sampler, *sampler_args, **sampler_kwargs)
            record.update({'ncall': result.get('ncall'), 'niter': result.get('niter'),
                           'logz': result.get('logz'), 'logzerr': result.get('logzerr')})
            return result
        except BaseException as exc:
            record['error'] = f'{type(exc).__name__}: {exc}'
            raise
        finally:
            record['wall_seconds'] = time.perf_counter() - started
            event(f"sampler completed: phase={record['phase']}, candidate={record['candidate']}, seconds={record['wall_seconds']:.2f}")

    elca.run_reactive_sampler = sampler_run
    original_start = runtime.log_comparison_candidate_evaluation_start

    def candidate_start(summary, *start_args, **start_kwargs):
        report['candidate'] = summary.get('label')
        event(f"candidate start: phase={report['phase']}, candidate={report['candidate']}")
        return original_start(summary, *start_args, **start_kwargs)

    runtime.log_comparison_candidate_evaluation_start = candidate_start
    original_finalize = runtime.finalize_comparison_candidate_full_reduction

    def finalize(*fit_args, **fit_kwargs):
        started = time.perf_counter()
        record = {'phase': report['phase'], 'candidate': report.get('candidate')}
        try:
            result = original_finalize(*fit_args, **fit_kwargs)
            record['applied'] = result.get('applied')
            return result
        finally:
            record['wall_seconds'] = time.perf_counter() - started
            report['candidate_reductions'].append(record)
            save()

    runtime.finalize_comparison_candidate_full_reduction = finalize
    modified_search = runtime.fit_ranked_comparison_calibration_candidates
    captured = {}

    def measure(mode, function, fit_args, fit_kwargs):
        report['phase'] = mode
        report['candidate'] = None
        before_samplers = len(report['sampler_runs'])
        started = time.perf_counter()
        cpu_start = cpu_seconds()
        load_start = psutil.cpu_percent(interval=None)
        event(f'{mode} search starting')
        try:
            result = function(*fit_args, **fit_kwargs)
        except BaseException as exc:
            report['trials'][mode] = {'error': f'{type(exc).__name__}: {exc}',
                                     'wall_seconds': time.perf_counter() - started}
            save()
            raise
        runs = report['sampler_runs'][before_samplers:]
        report['trials'][mode] = {
            'wall_seconds': time.perf_counter() - started,
            'process_cpu_seconds': cpu_seconds() - cpu_start,
            'system_cpu_percent_start': load_start,
            'system_cpu_percent_over_trial': psutil.cpu_percent(interval=None),
            'ultranest_runs': len(runs),
            'ultranest_evaluations': sum(int(run.get('ncall') or 0) for run in runs),
            'selection_metric': result.get('selection_metric'),
            'attempts': [fit_summary(attempt) for attempt in result['attempts']],
            'selected': fit_summary(result['selected_result']) if result['selected_result'] else None,
            'stopped_after_acceptable_shortlist': result.get('stopped_after_acceptable_shortlist', False),
            'stopped_after_first_qc_pass': result.get('stopped_after_first_qc_pass', False),
        }
        event(f"{mode} search complete: {report['trials'][mode]['wall_seconds']:.2f} seconds, {len(runs)} UltraNest runs")
        return result

    def capture_and_run(*fit_args, **fit_kwargs):
        captured['args'] = copy.deepcopy(fit_args)
        captured['kwargs'] = copy.deepcopy(fit_kwargs)
        captured['numpy_rng_state'] = np.random.get_state()
        captured['python_rng_state'] = random.getstate()
        # Snapshot the original function after normal runtime configuration.
        node = next(node for node in ast.parse(baseline_source).body
                    if isinstance(node, ast.FunctionDef) and node.name == modified_search.__name__)
        namespace = vars(runtime).copy()
        exec(compile(ast.Module(body=[node], type_ignores=[]),
                     f'git:{baseline_ref}:exotic/exotic.py', 'exec'), namespace)
        captured['baseline_function'] = namespace[node.name]
        checkpoint = {key: value for key, value in captured.items() if key != 'baseline_function'}
        checkpoint['baseline_ref'] = baseline_ref
        # These are the actual measured photometry inputs, not science-image copies.
        with (output / 'candidate_inputs.pkl').open('wb') as handle:
            pickle.dump(checkpoint, handle, protocol=pickle.HIGHEST_PROTOCOL)
        result = measure('shortlist', modified_search, fit_args, fit_kwargs)
        report['phase'] = 'modified_final_outputs'
        report['candidate'] = (result.get('selected_result') or {}).get('label')
        return result

    runtime.fit_ranked_comparison_calibration_candidates = capture_and_run
    command = ['exotic', '-red', str(init_path), '-ov', '--non-interactive-run',
               '--use-nextastro-astrometry']
    report['command'] = command
    sys.argv = command
    event('full modified pipeline starting')
    full_start = time.perf_counter()
    full_cpu_start = cpu_seconds()
    full_error = None
    try:
        result_code = runtime.cli()
        report['modified_full_return_code'] = result_code
    except BaseException as exc:
        full_error = exc
        report['modified_full_error'] = f'{type(exc).__name__}: {exc}'
        traceback.print_exc()
    report['modified_full_wall_seconds'] = time.perf_counter() - full_start
    report['modified_full_process_cpu_seconds'] = cpu_seconds() - full_cpu_start
    event('full modified pipeline finished')

    if captured:
        np.random.set_state(captured['numpy_rng_state'])
        random.setstate(captured['python_rng_state'])
        baseline_kwargs = copy.deepcopy(captured['kwargs'])
        baseline_output = output / 'baseline_candidates'
        baseline_output.mkdir(exist_ok=True)
        baseline_kwargs['save_dir'] = str(baseline_output)
        runtime.configure_runtime_logging(output_dir=str(baseline_output), start_new_run=True)
        try:
            measure('baseline', captured['baseline_function'],
                    copy.deepcopy(captured['args']), baseline_kwargs)
        finally:
            runtime.close_runtime_logging()
    if 'baseline' in report['trials'] and 'shortlist' in report['trials']:
        baseline = report['trials']['baseline']
        shortlist = report['trials']['shortlist']
        if 'error' not in baseline and 'error' not in shortlist:
            report['comparison'] = {
                'candidate_speedup': baseline['wall_seconds'] / shortlist['wall_seconds'],
                'candidate_seconds_saved': baseline['wall_seconds'] - shortlist['wall_seconds'],
                'selected_label_matches': (baseline['selected'] or {}).get('label') == (shortlist['selected'] or {}).get('label'),
                'parameter_differences_shortlist_minus_baseline': {
                    key: value - baseline['selected']['parameters'][key]
                    for key, value in (shortlist['selected'] or {}).get('parameters', {}).items()
                    if isinstance(value, (int, float, np.number))
                    and isinstance((baseline['selected'] or {}).get('parameters', {}).get(key), (int, float, np.number))
                },
            }
    report['phase'] = 'complete'
    report['finished_utc'] = datetime.now(timezone.utc).isoformat()
    event('benchmark complete')
    if full_error is not None:
        raise full_error


if __name__ == '__main__':
    main()
