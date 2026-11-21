"""Compare native Claret fits and orbital replay on saved, identical photometry.

Full FITS reductions are measured separately with profile_full_reduction.py.
These unprofiled checks isolate model speed, paired fits and saved-draw replay.
"""
import argparse
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.benchmark_hybrid4_cpp import read_reference
from scripts.benchmark_hybrid4_orbit_cpp import posterior_comparison
from scripts.profile_full_reduction import Tee
from exotic.api.hybrid_comparison import fit_summary, json_value
from exotic.api.limbdarkfast_cpp import _library, library_path
from exotic.orbital_uncertainty import _conditional_orbit_draws, _conditional_orbit_draws_python
from exotic.posterior import fit_posterior_sample_values, summarize_posterior


def model_comparison(fit):
    timings = {}
    for backend in ('python','cpp'):
        p = {**fit.prior,'_claret_backend':backend}
        fit._transit_model(fit.time,p)
        durations=[]
        for _ in range(101):
            started=perf_counter()
            fit._transit_model(fit.time,p)
            durations.append(perf_counter()-started)
        timings[backend]=float(np.median(durations))
    rng=np.random.default_rng(7931)
    worst_flux=worst_chi2=0.
    for index in range(65):
        p=dict(fit.prior)
        if index:
            for key in ('rprs','ars','inc','tmid'):
                if key in fit.bounds:
                    p[key]=rng.uniform(*fit.bounds[key])
        old=fit._transit_model(fit.time,{**p,'_claret_backend':'python'})
        new=fit._transit_model(fit.time,{**p,'_claret_backend':'cpp'})
        worst_flux=max(worst_flux,float(np.max(np.abs(new-old))))
        a=float(np.sum(((fit.data-old)/fit.dataerr)**2))
        b=float(np.sum(((fit.data-new)/fit.dataerr)**2))
        worst_chi2=max(worst_chi2,abs(a-b))
    return dict(median_complete_model_seconds=timings,speedup=timings['python']/timings['cpp'],
                max_absolute_flux_difference=worst_flux,max_absolute_chi2_difference=worst_chi2,
                identical_prior_draws=65)


def replay_comparison(fit):
    values,weights,source=fit_posterior_sample_values(fit)
    valid=np.isfinite(weights)&(weights>0)
    for v in values.values():
        valid &= np.isfinite(v)
    values={k:v[valid] for k,v in values.items()}
    weights=weights[valid]; weights=weights/weights.sum()
    base=fit.orbital_likelihood_context['base_physical']
    original_backend=base.get('_claret_backend','cpp')
    nodes=[(float(base['ecc']),float(base['omega'])),(.025,90.),(.025,270.)]
    rows=[]
    baseline=None
    try:
        for ecc,omega in nodes:
            results={};timings={}
            for backend in ('python','scalar_cpp','batch_cpp'):
                base['_claret_backend']='python' if backend=='python' else 'cpp'
                fit._orbital_replay_backend='cpp'
                function=_conditional_orbit_draws if backend=='batch_cpp' else _conditional_orbit_draws_python
                started=perf_counter()
                results[backend]=function(fit,ecc,omega,values,weights)
                timings[backend]=perf_counter()-started
            if baseline is None:
                baseline=results['python'][1]
            row=dict(ecc=ecc,omega=omega,seconds=timings,
                     speedup_python_to_batch=timings['python']/timings['batch_cpp'],
                     speedup_scalar_cpp_to_batch=timings['scalar_cpp']/timings['batch_cpp'])
            ref=results['python']
            new=results['batch_cpp']
            finite=np.isfinite(ref[1])&np.isfinite(new[1])
            row['finite_support_disagreements']=int(np.count_nonzero(np.isfinite(ref[1])!=np.isfinite(new[1])))
            row['max_absolute_score_difference']=float(np.max(np.abs(new[1][finite]-ref[1][finite]))) if finite.any() else None
            row['max_absolute_support_expansion_difference']=float(np.max(np.abs(new[2]-ref[2])))
            if finite.any():
                summaries={};normalized={};ess={}
                for backend in ('python','batch_cpp'):
                    physical,scores,_=results[backend]
                    log_weight=np.log(weights[finite])+scores[finite]-baseline[finite]
                    w=np.exp(log_weight-np.max(log_weight)); w/=w.sum()
                    normalized[backend]=w
                    ess[backend]=float(1./np.sum(w*w))
                    summaries[backend]={k:summarize_posterior(v[finite],w) for k,v in physical.items()}
                row['max_absolute_normalized_weight_difference']=float(np.max(np.abs(normalized['python']-normalized['batch_cpp'])))
                row['effective_sample_size']=ess
                row['posterior_summaries']=summaries
            rows.append(row)
            print('Replay',ecc,omega,'draws',len(weights),'seconds',timings,flush=True)
    finally:
        base['_claret_backend']=original_backend
        del fit._orbital_replay_backend
    return dict(sample_source=source,positive_weight_draws=len(weights),nodes=rows,
                timing_scope='One unprofiled evaluation per backend and node on identical saved joint draws; no sampler')


def run_case(path,output):
    saved=json.loads(path.read_text(encoding='utf-8'))
    name=path.stem.replace('LimbDarkeningComparison_','')
    destination=output/f'ClaretCppBenchmark_{name}.json'
    reference=read_reference(saved,path.parent)
    result=dict(status='running',dataset=name,source_report=str(path.resolve()),
                shared_fit_settings=saved['shared_fit_settings'],photometry=saved['photometry'],
                stellar_atmosphere=saved['stellar_atmosphere'],paired_runs=[],
                timing_scope='Unprofiled saved-photometry fits; original budgets and live-point settings retained; warm library; compilation excluded. Fit runtimes include the existing CPU admission waits; model/replay timings have no admission gate.')
    def save():
        destination.write_text(json.dumps(json_value(result),indent=2,allow_nan=False),encoding='utf-8')
    result['model_comparison']=model_comparison(reference)
    print(name,'model',result['model_comparison'],flush=True)
    save()
    for index,seed in enumerate((20261121,20261122,20261123)):
        order=('python','cpp') if index%2==0 else ('cpp','python')
        row=dict(seed=seed,execution_order=list(order))
        for backend in order:
            fit=read_reference(saved,path.parent)
            fit.prior['_claret_backend']=backend
            fit.verbose=False
            np.random.seed(seed)
            print(name,'seed',seed,'backend',backend,flush=True)
            fit.fit_nested()
            row[backend]=fit_summary(fit)
            if index==0 and backend=='cpp':
                result['orbital_replay']=replay_comparison(fit)
        row['posterior_comparison']=posterior_comparison(row['python'],row['cpp'])
        row['goodness_of_fit_differences']={k:row['cpp'][k]-row['python'][k]
            for k in ('chi2_of_reported_model','residual_rms_ppm','logz','logzerr')
            if row['cpp'][k] is not None and row['python'][k] is not None}
        row['speedup']=row['python']['runtime_seconds']/row['cpp']['runtime_seconds']
        result['paired_runs'].append(row)
        print('Pair complete',name,'speedup',row['speedup'],flush=True)
        save()
    medians={k:float(np.median([r[k]['runtime_seconds'] for r in result['paired_runs']])) for k in ('python','cpp')}
    result['timing_summary']=dict(median_fit_seconds=medians,speedup=medians['python']/medians['cpp'])
    result['status']='completed'
    save()
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report',type=Path,action='append',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    _library()
    metadata=dict(python=sys.version,numpy=np.__version__,platform=platform.platform(),
                  native_build=json.loads((library_path().parent/'build_info.json').read_text()))
    with (args.output/'console.log').open('w',encoding='utf-8',buffering=1) as stream:
        original=sys.stdout;sys.stdout=Tee(original,stream)
        try:
            results=[run_case(path,args.output) for path in args.report]
        finally:
            sys.stdout=original
    (args.output/'summary.json').write_text(json.dumps(json_value(dict(runtime=metadata,datasets=results)),indent=2,allow_nan=False),encoding='utf-8')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
