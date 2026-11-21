"""Compare the retained Python orbit/C++ flux path with the fused C++ model.

Use saved comparison reports so atmosphere profiles, photometry, uncertainties,
exposures, bounds and sampler settings are identical. All outputs stay local.
"""
import argparse
import csv
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.benchmark_hybrid4_cpp import SEEDS, read_reference
from exotic.api.hybrid_comparison import fit_summary, json_value, plot_comparison, run_hybrid4_fit
from exotic.api.limbdarkfast_cpp import _library, library_path


def posterior_comparison(reference, compiled):
    result = {}
    for key in sorted(set(reference['posterior']) & set(compiled['posterior'])):
        a, b = reference['posterior'][key], compiled['posterior'][key]
        width_a, width_b = a['upper']-a['lower'], b['upper']-b['lower']
        delta = b['median']-a['median']
        result[key] = {
            'median_difference_compiled_minus_reference': delta,
            'reference_68_percent_interval_width': width_a,
            'compiled_68_percent_interval_width': width_b,
            'interval_width_ratio_compiled_over_reference': width_b/width_a if width_a else None,
            'median_difference_in_reference_halfwidths': delta/(.5*width_a) if width_a else None,
        }
        if key == 'tmid':
            result[key]['median_difference_seconds'] = delta*86400.
    return result


def model_comparison(reference, atmosphere):
    prior = {**reference.prior,'_limb_darkening_law':'hybrid4',
             **dict(zip(('h0','h1','h2','h3'),atmosphere['coefficients']))}
    timings = {}
    for backend in ('numpy','cpp_flux','cpp'):
        p = {**prior,'_hybrid4_backend':backend}
        reference._transit_model(reference.time,p)
        durations=[]
        for _ in range(101):
            started=perf_counter()
            reference._transit_model(reference.time,p)
            durations.append(perf_counter()-started)
        timings[backend] = float(np.median(durations))
    rng=np.random.default_rng(7931)
    worst_flux=0.; worst_chi2=0.
    for index in range(65):
        p=dict(prior)
        if index:
            for key in ('rprs','ars','inc','tmid'):
                if key in reference.bounds:
                    p[key]=rng.uniform(*reference.bounds[key])
        old=reference._transit_model(reference.time,{**p,'_hybrid4_backend':'cpp_flux'})
        new=reference._transit_model(reference.time,{**p,'_hybrid4_backend':'cpp'})
        worst_flux=max(worst_flux,float(np.max(np.abs(new-old))))
        a=float(np.sum(((reference.data-old)/reference.dataerr)**2))
        b=float(np.sum(((reference.data-new)/reference.dataerr)**2))
        worst_chi2=max(worst_chi2,abs(a-b))
    return dict(median_complete_model_seconds=timings,
                speedup_over_previous_cpp_model=timings['cpp_flux']/timings['cpp'],
                speedup_over_numpy_model=timings['numpy']/timings['cpp'],
                same_parameters_max_absolute_flux_difference=worst_flux,
                same_parameters_max_absolute_chi2_difference=worst_chi2,
                prior_draws_including_initial_parameters=65)


def run_case(path, output):
    saved=json.loads(path.read_text(encoding='utf-8'))
    reference=read_reference(saved,path.parent)
    atmosphere=saved['stellar_atmosphere']
    name=path.stem.replace('LimbDarkeningComparison_','')
    destination=output/f'OrbitCppBenchmark_{name}.json'
    result=dict(status='running',dataset=name,source_report=str(path.resolve()),
                shared_fit_settings=saved['shared_fit_settings'],stellar_atmosphere=atmosphere,
                validation_inputs=saved['validation_inputs'],photometry=saved['photometry'],
                backend_definitions={'cpp_flux':'Python orbit/exposure integration with C++ hybrid4 flux',
                                     'cpp':'C++ orbit, exposure integration and hybrid4 flux'},
                seeds=list(SEEDS),paired_runs=[],
                timing_scope='Full UltraNest fit and posterior summaries; warm library; compilation excluded; existing CPU admission unchanged')
    def save():
        destination.write_text(json.dumps(json_value(result),indent=2,allow_nan=False),encoding='utf-8')
    print('\nDataset:',name,'points:',len(reference.time),flush=True)
    result['model_comparison']=model_comparison(reference,atmosphere)
    print('Model comparison:',result['model_comparison'],flush=True)
    np.random.seed(SEEDS[0])
    reference.mode='ns'
    reference.fit_nested()
    result['claret_control']=fit_summary(reference)
    save()
    for index,seed in enumerate(SEEDS):
        order=('cpp_flux','cpp') if index%2==0 else ('cpp','cpp_flux')
        fits={}
        row=dict(seed=seed,execution_order=list(order))
        for backend in order:
            print(name,'seed',seed,'backend',backend,flush=True)
            np.random.seed(seed)
            fits[backend]=run_hybrid4_fit(reference,atmosphere,backend=backend)
            row[backend]=fit_summary(fits[backend])
        old,new=row['cpp_flux'],row['cpp']
        row['posterior_comparison']=posterior_comparison(old,new)
        row['goodness_of_fit_comparison']={
            key: new[key]-old[key] for key in ('chi2_of_reported_model','residual_rms_ppm','logz','logzerr')
            if old[key] is not None and new[key] is not None}
        row['same_seed_likelihood_call_count_difference']=new['ncall']-old['ncall']
        row['speedup_over_previous_cpp']=old['runtime_seconds']/new['runtime_seconds']
        result['paired_runs'].append(row)
        print('Pair complete: previous C++',old['runtime_seconds'],'full C++',new['runtime_seconds'],
              'speedup',row['speedup_over_previous_cpp'],flush=True)
        if index==0:
            csv_path=output/f'CppFluxVsFullCpp_{name}.csv'
            plot_path=output/f'CppFluxVsFullCpp_{name}.png'
            with csv_path.open('w',newline='',encoding='utf-8') as stream:
                writer=csv.writer(stream)
                writer.writerow(['BJD_TDB','Shared Flux','Shared Uncertainty','Previous C++ Model',
                                 'Full C++ Model','Previous C++ Residual','Full C++ Residual'])
                writer.writerows([float(x) for x in row] for row in zip(
                    reference.time,reference.data,reference.dataerr,fits['cpp_flux'].model,
                    fits['cpp'].model,reference.data-fits['cpp_flux'].model,reference.data-fits['cpp'].model))
            plot_comparison(fits['cpp_flux'],fits['cpp'],plot_path,name,
                            labels=('Previous C++','Full C++'))
            result['files']=dict(json=destination.name,comparison_plot=plot_path.name,
                                 photometry_and_models_csv=csv_path.name)
        save()
    medians={backend:float(np.median([row[backend]['runtime_seconds'] for row in result['paired_runs']]))
             for backend in ('cpp_flux','cpp')}
    result['timing_summary']=dict(
        median_runtime_seconds=medians,
        speedup_over_previous_cpp=medians['cpp_flux']/medians['cpp'],
        percent_less_time_than_previous_cpp=100.*(1.-medians['cpp']/medians['cpp_flux']),
        claret_control_runtime_seconds=result['claret_control']['runtime_seconds'],
        speedup_over_claret_control=result['claret_control']['runtime_seconds']/medians['cpp'])
    result['status']='completed'
    save()
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,action='append',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    _library()
    metadata=dict(python=sys.version,numpy=np.__version__,platform=platform.platform(),
                  native_build=json.loads((library_path().parent/'build_info.json').read_text()))
    results=[]
    for path in args.report:
        try:
            result=run_case(path,args.output)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            result=dict(status='failed',source_report=str(path),error_type=type(exc).__name__,reason=str(exc))
        results.append(result)
        (args.output/'orbit_cpp_comparison_summary.json').write_text(
            json.dumps(json_value(dict(runtime=metadata,datasets=results)),indent=2,allow_nan=False),encoding='utf-8')
    return int(any(row['status']!='completed' for row in results))


if __name__=='__main__':
    raise SystemExit(main())
