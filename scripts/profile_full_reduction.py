"""Profile a complete EXOTIC FITS reduction without changing its science path.

Only input/output locations and a reproducible RNG seed are set in a local copy
of the initialization file. Original observations and previous outputs remain
in place. cProfile includes main-thread Python/C calls; wrapped stage times
separate nested work without counting it twice. All artifacts stay local.
"""
import argparse
import copy
import cProfile
from datetime import datetime, timezone
import functools
import json
from pathlib import Path
import pstats
import sys
import threading
from time import perf_counter
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


class StageTimes:
    def __init__(self):
        self.main_thread = threading.get_ident()
        self.stack = []
        self.groups = {}
        self.calls = {}
        self.fits = []
        self.worker_calls = {}

    def wrap(self, owner, name, category, *, fit=False):
        original = getattr(owner, name, None)
        if not callable(original):
            return
        label = f'{getattr(owner, "__name__", type(owner).__name__)}.{name}'

        @functools.wraps(original)
        def measured(*args, **kwargs):
            if threading.get_ident() != self.main_thread:
                self.worker_calls[label] = self.worker_calls.get(label, 0) + 1
                return original(*args, **kwargs)
            started = perf_counter()
            frame = dict(category=category, function=label, children=0., started=started)
            outermost = not any(f['category'] == category for f in self.stack)
            self.stack.append(frame)
            record = None
            if fit:
                obj = args[0]
                prior = getattr(obj, 'prior', {})
                record = dict(law=prior.get('_limb_darkening_law', 'claret'),
                              backend=prior.get('_hybrid4_backend'), points=len(obj.time))
            try:
                return original(*args, **kwargs)
            except BaseException as exc:
                if record is not None:
                    record['error'] = f'{type(exc).__name__}: {exc}'
                raise
            finally:
                elapsed = perf_counter()-started
                self.stack.pop()
                if self.stack:
                    self.stack[-1]['children'] += elapsed
                group = self.groups.setdefault(category, dict(calls=0, exclusive_seconds=0., inclusive_seconds=0.))
                group['calls'] += 1
                group['exclusive_seconds'] += elapsed-frame['children']
                if outermost:
                    group['inclusive_seconds'] += elapsed
                call = self.calls.setdefault(label, dict(calls=0, inclusive_seconds=0.))
                call['calls'] += 1
                call['inclusive_seconds'] += elapsed
                if record is not None:
                    record['wall_seconds'] = elapsed
                    record['ncall'] = int((getattr(obj, 'results', None) or {}).get('ncall') or 0)
                    record['likelihood_backend'] = getattr(obj, 'likelihood_backend', None)
                    self.fits.append(record)

        setattr(owner, name, measured)


class Tee:
    def __init__(self, original, stream):
        self.original, self.stream = original, stream

    def write(self, value):
        # Colorama can retain this tee for its final terminal reset at atexit,
        # after the log context has closed. The foreground console stays valid.
        if not self.stream.closed:
            self.stream.write(value)
        return self.original.write(value)

    def flush(self):
        if not self.stream.closed:
            self.stream.flush()
        self.original.flush()

    def __getattr__(self, name):
        return getattr(self.original, name)


def raw_profile_rows(profiler):
    """Preserve distinct code objects with identical source labels.

    Legacy and package imports can load elca twice. Standard pstats keys use
    only filename/line/name and overwrite such entries. Ordinals preserve the
    complete call graph without hashing code or science data.
    """
    entries = profiler.getstats()
    # lsprof recreates builtin-name strings for caller edges. Their string
    # identities can differ from the corresponding top-level entry. Python
    # code objects must still keep their separate identities across imports.
    def identity(code):
        return ('builtin',code) if isinstance(code,str) else ('python',id(code))
    identifiers = {identity(entry.code): index for index, entry in enumerate(entries)}
    # Some C-extension callback edges also have no top-level lsprof entry.
    # Retain their code identities and edge measurements without inventing a
    # top-level call count or dropping the edge.
    callee_only = []
    for entry in entries:
        for child in (entry.calls or []):
            key = identity(child.code)
            if key not in identifiers:
                identifiers[key] = len(entries) + len(callee_only)
                callee_only.append(child.code)
    rows = []
    for index, entry in enumerate(entries):
        code = entry.code
        rows.append(dict(
            code_id=index, file=getattr(code, 'co_filename', '~'),
            line=getattr(code, 'co_firstlineno', 0),
            function=getattr(code, 'co_qualname', str(code)),
            calls=entry.callcount, recursive_calls=entry.reccallcount,
            own_seconds=entry.inlinetime, cumulative_seconds=entry.totaltime,
            callees=[dict(code_id=identifiers[identity(child.code)], calls=child.callcount,
                          recursive_calls=child.reccallcount, own_seconds=child.inlinetime,
                          cumulative_seconds=child.totaltime) for child in (entry.calls or [])],
        ))
    for code in callee_only:
        rows.append(dict(code_id=identifiers[identity(code)],
                         file=getattr(code,'co_filename','~'),line=getattr(code,'co_firstlineno',0),
                         function=getattr(code,'co_qualname',str(code)),
                         entry_kind='callee_only_no_top_level_stat',calls=None,recursive_calls=None,
                         own_seconds=None,cumulative_seconds=None,callees=[]))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inits', type=Path, required=True)
    parser.add_argument('--science-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--no-cprofile', action='store_true',
                        help='Measure full wall/stage time without per-call profiler overhead')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    init = json.loads(args.inits.read_text(encoding='utf-8-sig'))
    init = copy.deepcopy(init)
    init['user_info']['Directory with FITS files'] = str(args.science_dir.resolve())
    init['user_info']['Directory to Save Plots'] = str(output/'reduction')
    init.setdefault('optional_info', {}).setdefault('Random Seed', 20261121)
    init_path = output/'inits_profile.json'
    init_path.write_text(json.dumps(init, indent=2), encoding='utf-8')
    original_stdout, original_stderr = sys.stdout, sys.stderr
    with (output/'console.log').open('w', encoding='utf-8', buffering=1) as log:
        sys.stdout, sys.stderr = Tee(original_stdout, log), Tee(original_stderr, log)
        try:
            return run(args, output, init_path)
        finally:
            sys.stdout, sys.stderr = original_stdout, original_stderr


def run(args, output, init_path):
    imported = perf_counter()
    import matplotlib
    matplotlib.use('Agg')
    import numpy as np
    import psutil
    import requests
    from matplotlib.figure import Figure
    import exotic.exotic as runtime
    from exotic.api import elca, ultranest_utils
    from exotic.api.native_build import ensure_hybrid4_native
    import_seconds = perf_counter()-imported
    prepared = perf_counter()
    ensure_hybrid4_native()
    native_setup_seconds = perf_counter()-prepared
    stages = StageTimes()
    categories = {
        'psf_centroid': ['fit_centroid'],
        'psf_flux': ['fit_psf_photometry_flux_row', 'fit_legacy_psf_photometry_flux_row'],
        'aperture': ['compute_star_aperture_grid', 'aperPhot'],
        'alignment': ['filter_pointing_outlier_frames', 'build_multiprocess_alignment_results', 'transformation'],
        'candidate_selection': ['fit_ranked_comparison_calibration_candidates',
                                'fit_lightcurve_to_every_comparison_candidate', 'run_target_driven_photometry_search'],
        'atmosphere': ['nonlinear_ld'],
        'plots': ['plot_fov', 'plot_centroids', 'plot_obs_stats', 'plot_final_lightcurve', 'plot_flux'],
        'calibration': ['process_dark_frames', 'process_bias_frames', 'process_flat_frames'],
        'image_read': ['read_reduction_image_data'],
        'image_headers': ['read_reduction_image_header'],
        'image_calibration': ['apply_cals', 'demosaic_img'],
        'bad_pixels': ['repair_bad_pixels_in_frame'],
        'frame_preparation': ['load_calibrated_reduction_frame'],
    }
    for category, functions in categories.items():
        for name in functions:
            stages.wrap(runtime, name, category)
    stages.wrap(elca.lc_fitter, 'fit_nested', 'transit_inference', fit=True)
    stages.wrap(elca.lc_fitter, 'fit_LM', 'transit_inference', fit=True)
    stages.wrap(requests.Session, 'send', 'network')
    stages.wrap(Figure, 'savefig', 'plots')
    command = ['exotic', '--reduce', str(init_path), '--override',
               '--non-interactive-run', '--use-nextastro-astrometry']
    sys.argv = command
    report = dict(status='running', started_utc=datetime.now(timezone.utc).isoformat(),
                  source_inits=str(args.inits.resolve()), science_dir=str(args.science_dir.resolve()),
                  output=str(output), command=command, import_seconds=import_seconds,
                  native_setup_seconds=native_setup_seconds, numpy=np.__version__, python=sys.version,
                  logical_cpus=psutil.cpu_count(),
                  configured_sampler_workers=ultranest_utils._configured_ultranest_workers(),
                  configured_sampler_backend=ultranest_utils._configured_ultranest_worker_backend(),
                  image_preparation_backend=runtime.bad_pixel_processing.PREPARATION_BACKEND,
                  cprofile_enabled=not args.no_cprofile,
                  profiling_scope=('Full reduction with lightweight stage clocks; cProfile disabled.'
                      if args.no_cprofile else 'Full reduction with cProfile on main thread. Profile overhead is included; worker-thread calls are counted separately. Stage exclusive times do not overlap.'))
    process = psutil.Process()
    initial_cpu = process.cpu_times()
    started = perf_counter()
    stop = threading.Event()

    def snapshot():
        data = dict(report, wall_seconds=perf_counter()-started,
                    stages=copy.deepcopy(stages.groups), functions=copy.deepcopy(stages.calls),
                    fits=copy.deepcopy(stages.fits), worker_thread_calls=dict(stages.worker_calls))
        data['active_stages'] = [dict(f) for f in stages.stack]
        (output/'profile_summary.json').write_text(json.dumps(data, indent=2), encoding='utf-8')
        return data

    def heartbeat():
        while not stop.wait(30):
            snapshot()
            active = stages.stack[-1]['function'] if stages.stack else 'reduction preprocessing / output'
            print(f'PROFILE | elapsed={perf_counter()-started:.1f}s | active={active}', flush=True)

    thread = threading.Thread(target=heartbeat, name='profile-progress')
    thread.start()
    profiler = None if args.no_cprofile else cProfile.Profile()
    snapshot()
    print('PROFILE | Starting full reduction:', args.science_dir, flush=True)
    if profiler is not None:
        profiler.enable()
    try:
        result = runtime.cli()
        report['return_code'] = result
        report['status'] = 'completed' if result in (None, 0) and any((output/'reduction').glob('FinalParams_*.json')) else 'incomplete'
    except BaseException as exc:
        report.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        traceback.print_exc()
    finally:
        if profiler is not None:
            profiler.disable()
        report['reduction_wall_seconds'] = perf_counter()-started
        stop.set()
        thread.join()
        final_cpu = process.cpu_times()
        report['process_cpu_seconds'] = final_cpu.user+final_cpu.system-initial_cpu.user-initial_cpu.system
        # Preserve the final clocks and standard profile before converting the
        # optional raw graph, so an exporter failure cannot erase measurements.
        snapshot()
        if profiler is not None:
            profiler.dump_stats(str(output/'full_reduction.prof'))
            with (output/'profile_functions.txt').open('w', encoding='utf-8') as stream:
                stats = pstats.Stats(profiler, stream=stream).strip_dirs()
                stats.sort_stats('cumulative').print_stats(100)
                stats.sort_stats('tottime').print_stats(100)
            rows = []
            for (filename, line, function), (primitive, calls, own, cumulative, callers) in pstats.Stats(profiler).stats.items():
                rows.append(dict(file=filename, line=line, function=function, calls=calls,
                                 primitive_calls=primitive, own_seconds=own, cumulative_seconds=cumulative))
            (output/'profile_functions.json').write_text(json.dumps(sorted(rows,key=lambda r:r['cumulative_seconds'],reverse=True), indent=2), encoding='utf-8')
            (output/'profile_raw_calls.json').write_text(json.dumps(raw_profile_rows(profiler), indent=2), encoding='utf-8')
        snapshot()
        print(f"PROFILE | {report['status']} | wall={report['reduction_wall_seconds']:.2f}s | output={output}", flush=True)
    return int(report['status'] != 'completed')


if __name__ == '__main__':
    raise SystemExit(main())
