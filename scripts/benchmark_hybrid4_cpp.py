"""Rerun saved Claret/hybrid4 comparison inputs with NumPy and C++ hybrid4.

Use --report for each previously generated LimbDarkeningComparison JSON and
--output for a local directory. Reuse its exact coefficients, photometry,
exposures, bounds and duration prior. No atmosphere regeneration or clipping.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from exotic.api.elca import lc_fitter
from exotic.api.hybrid_comparison import fit_summary, json_value, plot_comparison, run_hybrid4_fit
from exotic.api.limbdarkfast import hybrid_flux
from exotic.api.limbdarkfast_cpp import hybrid4_flux, library_path


SEEDS = (20261120, 20261121, 20261122)


def read_reference(report, report_directory=None):
    csv_path = Path(report['files']['photometry_and_models_csv'])
    if not csv_path.is_absolute() and report_directory is not None:
        csv_path = Path(report_directory) / csv_path
    data = np.loadtxt(csv_path, delimiter=',', skiprows=1)
    settings = report['shared_fit_settings']
    prior = copy.deepcopy(report['claret']['raw_maximum_likelihood_parameters'])
    fit = lc_fitter(
        data[:,0], data[:,1], data[:,2], data[:,3], prior,
        copy.deepcopy(settings['bounds']), mode=None,
        exposure_times_seconds=data[:,4],
        use_impactparameter_rather_than_inclination_to_fit=settings['impact_parameter_sampling'],
        duration_prior=copy.deepcopy(settings['duration_prior']),
        baseline_fit_mask=report['photometry']['baseline_fit_mask'],
        fixed_flux_baseline=settings['fixed_flux_baseline'],
        ultranest_min_num_live_points=settings['minimum_live_points'],
        exposure_smearing_supersample=settings['exposure_supersample'],
        exposure_smearing_change_tolerance=settings['exposure_change_tolerance'],
    )
    fit.max_ncalls = settings['max_ncalls']
    return fit


def uncertainty_comparison(a, b):
    result = {}
    for key in set(a['posterior']) & set(b['posterior']):
        x, y = a['posterior'][key], b['posterior'][key]
        width_a, width_b = x['upper']-x['lower'], y['upper']-y['lower']
        delta = y['median']-x['median']
        result[key] = {
            'median_difference_cpp_minus_numpy': delta,
            '68_percent_interval_width_numpy': width_a,
            '68_percent_interval_width_cpp': width_b,
            'interval_width_ratio_cpp_over_numpy': width_b/width_a if width_a else None,
            'median_difference_in_numpy_68_percent_halfwidths': delta/(0.5*width_a) if width_a else None,
        }
        if key == 'tmid':
            result[key]['median_difference_seconds'] = delta*86400.0
    return result


def microbench(reference, atmosphere):
    # Include the orbit and the existing exposure integration in a separate
    # benchmark; DLL loading is warmed before measuring either implementation.
    prior = {**reference.prior, '_limb_darkening_law': 'hybrid4',
             **dict(zip(('h0','h1','h2','h3'), atmosphere['coefficients']))}
    values = reference._values_with_exposure_smearing(prior, reference.time)
    from exotic.api.elca import transit
    z = np.linspace(0, 1.2, len(reference.time))
    timings = {}
    for backend in ('numpy','cpp'):
        p = {**values,'_hybrid4_backend':backend}
        flux = (lambda: hybrid4_flux(z, p['rprs'], atmosphere['coefficients'])) if backend == 'cpp' else (
            lambda: hybrid_flux(z, p['rprs'], atmosphere['coefficients'], 'hybrid4'))
        for name, function in [('flux_only', flux), ('orbit_and_exposure', lambda: transit(reference.time,p))]:
            function()
            durations=[]
            for _ in range(51):
                started=perf_counter(); function(); durations.append(perf_counter()-started)
            timings[f'{backend}_{name}_median_seconds']=float(np.median(durations))
    # Compare the complete integrated models at identical physical parameters,
    # including radii and geometries spanning the actual prior box.
    rng=np.random.default_rng(7931)
    worst_flux=0.0
    worst_chi2=0.0
    for index in range(17):
        p=dict(prior)
        if index:
            for key in ('rprs','ars','inc','tmid'):
                if key in reference.bounds:
                    p[key]=rng.uniform(*reference.bounds[key])
        f0=reference._transit_model(reference.time,{**p,'_hybrid4_backend':'numpy'})
        f1=reference._transit_model(reference.time,{**p,'_hybrid4_backend':'cpp'})
        worst_flux=max(worst_flux,float(np.max(np.abs(f1-f0))))
        chi0=float(np.sum(((reference.data-f0)/reference.dataerr)**2))
        chi1=float(np.sum(((reference.data-f1)/reference.dataerr)**2))
        worst_chi2=max(worst_chi2,abs(chi1-chi0))
    timings['same_parameters_max_abs_integrated_flux_difference']=worst_flux
    timings['same_parameters_max_abs_chi2_difference']=worst_chi2
    return timings


def write_pair_csv(reference, native, path):
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream)
        writer.writerow(['BJD_TDB','Shared Flux','Shared Uncertainty','NumPy Model','C++ Model',
                         'NumPy Residual','C++ Residual'])
        writer.writerows([float(x) for x in row] for row in zip(
            reference.time, reference.data, reference.dataerr, reference.model, native.model,
            reference.data-reference.model,native.data-native.model))


def run_case(path, output):
    saved=json.loads(path.read_text(encoding='utf-8'))
    reference=read_reference(saved, path.parent)
    atmosphere=saved['stellar_atmosphere']
    name=path.stem.replace('LimbDarkeningComparison_','')
    result={
        'status':'running','dataset':name,'source_report':str(path.resolve()),
        'shared_fit_settings':saved['shared_fit_settings'],
        'photometry_point_count':len(reference.time),
        'stellar_atmosphere':atmosphere,
        'validation_inputs':saved['validation_inputs'],
        'timing_scope':'Warm native library, existing CPU admission policy, full fit including posterior summaries',
        'seeds':list(SEEDS),'paired_runs':[],
    }
    result['model_benchmarks']=microbench(reference,atmosphere)
    print('\nDataset:',name,'points:',len(reference.time),flush=True)
    print('Model benchmark:',result['model_benchmarks'],flush=True)
    # A current Claret control avoids comparing a native run only with an old
    # timing captured under potentially different machine load or dependencies.
    np.random.seed(SEEDS[0])
    reference.mode='ns'
    started=perf_counter();reference.fit_nested()
    reference.fit_runtime_seconds=perf_counter()-started
    result['claret_control']=fit_summary(reference)
    for index,seed in enumerate(SEEDS):
        order=('numpy','cpp') if index%2==0 else ('cpp','numpy')
        fits={}
        row={'seed':seed,'execution_order':list(order)}
        for backend in order:
            print(name,'seed',seed,'backend',backend,flush=True)
            np.random.seed(seed)
            fits[backend]=run_hybrid4_fit(reference,atmosphere,backend=backend)
            row[backend]=fit_summary(fits[backend])
        row['uncertainty_comparison']=uncertainty_comparison(row['numpy'],row['cpp'])
        row['runtime_speedup_cpp_over_numpy']=row['numpy']['runtime_seconds']/row['cpp']['runtime_seconds']
        result['paired_runs'].append(row)
        print('Pair complete: NumPy',row['numpy']['runtime_seconds'],'C++',row['cpp']['runtime_seconds'],
              'speedup',row['runtime_speedup_cpp_over_numpy'],flush=True)
        if index==0:
            plot_path=output/f'NumPyVsCpp_{name}.png'
            csv_path=output/f'NumPyVsCpp_{name}.csv'
            plot_comparison(fits['numpy'],fits['cpp'],plot_path,f'hybrid4 NumPy and C++: {name}',
                            labels=('NumPy','C++'))
            write_pair_csv(fits['numpy'],fits['cpp'],csv_path)
            result['files']={'comparison_plot':str(plot_path),'photometry_and_models_csv':str(csv_path)}
        (output/f'CppBenchmark_{name}.json').write_text(json.dumps(json_value(result),indent=2,allow_nan=False),encoding='utf-8')
    medians={backend:float(np.median([row[backend]['runtime_seconds'] for row in result['paired_runs']]))
             for backend in ('numpy','cpp')}
    result['timing_summary']={
        'median_runtime_seconds':medians,
        'speedup_cpp_over_numpy':medians['numpy']/medians['cpp'],
        'percent_less_time_than_numpy':100*(1-medians['cpp']/medians['numpy']),
        'claret_control_runtime_seconds':result['claret_control']['runtime_seconds'],
        'speedup_cpp_over_claret_control':result['claret_control']['runtime_seconds']/medians['cpp'],
    }
    result['status']='completed'
    (output/f'CppBenchmark_{name}.json').write_text(json.dumps(json_value(result),indent=2,allow_nan=False),encoding='utf-8')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,action='append',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    hybrid4_flux([0.0],0.1,[0.3,0.2,0.1,0.05])
    build_info=library_path().parent/'build_info.json'
    metadata={'python':sys.version,'numpy':np.__version__,'platform':platform.platform(),
              'native_build':json.loads(build_info.read_text()) if build_info.is_file() else None}
    results=[]
    for path in args.report:
        try:
            result=run_case(path,args.output)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            result={'status':'failed','source_report':str(path),'error_type':type(exc).__name__,'reason':str(exc)}
        results.append(result)
        (args.output/'cpp_comparison_summary.json').write_text(
            json.dumps(json_value({'runtime':metadata,'datasets':results}),indent=2,allow_nan=False),encoding='utf-8')
    return 1 if any(row['status']!='completed' for row in results) else 0


if __name__=='__main__':
    raise SystemExit(main())
