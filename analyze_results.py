"""Summarize all outcomes; same-document paired conditional model intervals."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import argparse,itertools,json,hashlib
import numpy as np
import pandas as pd
from runner import HERE,KEY,sha,save
def collect(run,label):
    rec=[json.loads(p.read_text(encoding='utf8')) for p in sorted(run.glob('*.complete.json'))]
    assert rec and not (run/'RUNNING.lock').exists()
    for r in rec:assert sha(run/r['prediction_file'])==r['prediction_sha256']
    df=pd.DataFrame(rec);df['panel']=label
    tg=df.groupby(['panel']+KEY+['model'],as_index=False)[['mse_model','mse_mean','mae_model','mae_mean','rmse_model','rmse_mean','D','skill']].median()
    intervals=[];pairs=[]
    for k,group in df.groupby(KEY,sort=True):
        ordered=group.sort_values(['model','seed']);models=ordered.model.tolist();testids=None;loss=[]
        for r in ordered.itertuples():
            f=pd.read_csv(run/r.prediction_file)
            if testids is None:testids=f.record_uid.tolist();docs=f.document_chembl_id.to_numpy();base=f.squared_baseline_error.to_numpy()
            assert f.record_uid.tolist()==testids and np.allclose(f.squared_baseline_error,base,rtol=0,atol=1e-12)
            loss.append(f.squared_model_error.to_numpy())
        loss=np.stack(loss,axis=1);unique=np.array(sorted(set(docs)));sums=np.array([np.r_[loss[docs==d].sum(axis=0),base[docs==d].sum(),sum(docs==d)] for d in unique])
        seed=int(hashlib.sha256(('|'.join((label,)+k)).encode()).hexdigest()[:8],16)
        w=np.random.default_rng(seed).multinomial(len(unique),np.full(len(unique),1/len(unique)),size=2000);b=w@sums
        assert np.all(b[:,-2]>0) and np.all(b[:,-1]>0)
        bs=1-b[:,:len(models)]/b[:,-2,None];point=1-loss.mean(0)/base.mean();by={}
        for model in sorted(set(models)):
            ix=[i for i,n in enumerate(models) if n==model];z=np.median(bs[:,ix],axis=1);v=float(np.median(point[ix]));lo,hi=np.quantile(z,[.025,.975]);by[model]=(v,z)
            intervals.append(dict(panel=label,**dict(zip(KEY,k)),model=model,rows=len(testids),documents=len(unique),skill=v,ci_low=float(lo),ci_high=float(hi),bootstrap_seed=seed,replicates=2000,uncertainty_scope='conditional_fixed_models_and_cohort'))
        for a,c in itertools.combinations(sorted(by),2):
            delta=by[a][1]-by[c][1];lo,hi=np.quantile(delta,[.025,.975])
            pairs.append(dict(panel=label,**dict(zip(KEY,k)),model_a=a,model_b=c,skill_a_minus_b=by[a][0]-by[c][0],ci_low=float(lo),ci_high=float(hi),replicates=2000,common_document_weights=True,uncertainty_scope='conditional_fixed_models_and_cohort'))
    return df,tg,pd.DataFrame(intervals),pd.DataFrame(pairs)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--current',type=Path,required=True);ap.add_argument('--external',type=Path);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    out=args.output.resolve();out.relative_to(HERE);out.mkdir(parents=True,exist_ok=False)
    allparts=[collect(args.current,'kinase')]
    if args.external:allparts.append(collect(args.external,'nonkinase'))
    for name,i in [('run_metrics',0),('target_metrics',1),('conditional_document_intervals',2),('paired_model_intervals',3)]:pd.concat([p[i] for p in allparts],ignore_index=True).to_csv(out/(name+'.csv'),index=False)
    tg=pd.concat([p[1] for p in allparts],ignore_index=True)
    panel=tg.groupby(['panel','view','scenario','model','route']).agg(median_skill=('skill','median'),positive_targets=('skill',lambda x:int((x>0).sum())),target_count=('target','nunique'),median_rmse_model=('rmse_model','median'),median_rmse_mean=('rmse_mean','median'),median_mae_model=('mae_model','median'),median_mae_mean=('mae_mean','median')).reset_index();panel.to_csv(out/'panel_summary.csv',index=False)
    routes=tg.pivot(index=['panel','view','scenario','target','model'],columns='route',values='skill').reset_index()
    for a,b in [('temporal_exact','source_exact'),('source_scaffold','source_exact'),('temporal_scaffold','temporal_exact')]:routes[a+'_minus_'+b]=routes[a]-routes[b]
    routes.to_csv(out/'target_route_contrasts.csv',index=False)
    cols=['temporal_exact_minus_source_exact','source_scaffold_minus_source_exact','temporal_scaffold_minus_temporal_exact']
    summary=routes.groupby(['panel','view','scenario','model'])[cols].median().reset_index()
    summary['targets_with_lower_year_skill']=routes.groupby(['panel','view','scenario','model']).temporal_exact_minus_source_exact.apply(lambda x:int((x<0).sum())).values
    summary.to_csv(out/'contrast_summary.csv',index=False)
    save(out/'ANALYSIS_COMPLETE.json',dict(prediction_cases=sum(len(p[0]) for p in allparts),target_model_contexts=len(tg),bootstrap_replicates=2000,models_resampled_together=True,training_uncertainty_included=False,bootstrap_code_sha256=sha(__file__)))
    print(summary[summary.view.eq('main')].round(6).to_string(index=False))
if __name__=='__main__':main()
