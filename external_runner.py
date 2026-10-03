"""Fresh, package-relative, restartable model validation. No archived predictions loaded."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_k]='1'
from pathlib import Path
import argparse, hashlib, json, sys, time, platform
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
DEFAULT=HERE/'external/prepared_v2'
KEY=['view','scenario','target','route']
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,obj):
    p=Path(p); tmp=p.with_suffix(p.suffix+'.tmp'); tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False),encoding='utf8'); tmp.replace(p)
def sig(f): return hashlib.sha256(f[['record_uid','smiles','target']].to_csv(index=False).encode()).hexdigest()
def local(root,r):
    p=(root/'data'/str(r.view)/str(r.corrected_csv).replace('\\','/').rsplit('/',1)[-1]).resolve(); p.relative_to(root.resolve()); return p
def load_inputs(root):
    m=pd.read_csv(root/'prediction_manifest.csv'); s=pd.read_csv(root/'split_manifest.csv')
    plan=json.loads((root/'decision_record.json').read_text(encoding='utf8')); assert len(m)==len(s)==plan['contexts']*3
    assert not m.duplicated(KEY+['seed']).any() and not s.duplicated(KEY+['partition']).any()
    assert set(m.seed)=={42,53,67}
    assert set(m.target).issubset(set(plan['selected_targets']))
    assert set(m.route)=={'source_exact','source_scaffold','temporal_exact','temporal_scaffold'}
    assert set(m.view)=={'main'} and set(m.scenario)=={'revised','single_protein'}
    frames={}; hashes={}
    for r in s.itertuples():
        p=local(root,r); assert sha(p)==r.corrected_sha256; f=pd.read_csv(p)
        assert len(f)==r.corrected_rows and sig(f)==r.corrected_ordered_signature
        assert {'record_uid','smiles','target','document_chembl_id','murcko_scaffold'}.issubset(f.columns)
        assert f.record_uid.is_unique and np.isfinite(f.target).all()
        frames[tuple(getattr(r,k) for k in KEY)+ (r.partition,)]=f; hashes[str(p.relative_to(root))]=sha(p)
    for r in m.itertuples():
        k=tuple(getattr(r,x) for x in KEY); assert sig(frames[k+('train',)])==r.training_signature
        assert len(frames[k+('test',)])==r.test_rows
    for name in ['prediction_manifest.csv','split_manifest.csv','model_config.json','decision_record.json']:
        hashes[name]=sha(root/name)
    return m,s,frames,hashes
def first_nn(x,train):
    # int32 avoids uint8 overflow; np.argmax selects first original row, never labels.
    a=train.astype(np.int32); b=np.asarray(x,dtype=np.int32); inter=a@b; union=a.sum(1)+b.sum()-inter
    sim=np.divide(inter,union,out=np.ones(len(a),dtype=float),where=union!=0)
    return int(np.argmax(sim)),float(sim.max())
def metrics(y,p,mean):
    mse=float(np.mean((y-p)**2)); base=float(np.mean((y-mean)**2))
    return dict(mse_model=mse,mse_mean=base,mae_model=float(np.mean(abs(y-p))),mae_mean=float(np.mean(abs(y-mean))),rmse_model=float(np.sqrt(mse)),rmse_mean=float(np.sqrt(base)),D=base-mse,skill=1-mse/base if base>0 else None)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--input-root',type=Path,default=DEFAULT); ap.add_argument('--output',required=True,type=Path); ap.add_argument('--resume',action='store_true'); ap.add_argument('--model',choices=['all','rf','ridge','tanimoto1nn'],default='all'); ap.add_argument('--limit-cases',type=int,default=0); ap.add_argument('--include-paired',action='store_true',help='Also evaluate the preserved paired-operational sensitivity views'); args=ap.parse_args()
    root=args.input_root.resolve(); out=args.output.resolve(); out.relative_to(HERE)
    assert out!=HERE and root not in out.parents and out!=root
    if out.exists() and not args.resume: raise FileExistsError('Use a new directory or explicit --resume')
    import rdkit,sklearn,joblib
    from rdkit import Chem,DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from scipy.sparse import csr_matrix
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import Ridge
    cfg=json.loads((root/'model_config.json').read_text(encoding='utf8')); versions={k:v for k,v in [('numpy',np.__version__),('pandas',pd.__version__),('rdkit',rdkit.__version__),('sklearn',sklearn.__version__)]}
    assert all(versions[k]==cfg['environment'][k] for k in versions), versions
    m,s,frames,hashes=load_inputs(root)
    selected=m[m.scenario.isin(['revised','single_protein']) & (True if args.include_paired else m.view.eq('main'))].drop_duplicates(KEY).sort_values(KEY)
    assert len(selected)==json.loads((root/'decision_record.json').read_text(encoding='utf8'))['contexts']
    binding=dict(input_hashes=hashes,protocol_sha256=sha(HERE/'EXTERNAL_PROTOCOL.md'),runner_sha256=sha(__file__),versions=versions)
    if out.exists(): assert json.loads((out/'binding.json').read_text(encoding='utf8'))==binding
    else: out.mkdir(parents=True); save(out/'binding.json',binding)
    # Exclusive lock is intentionally retained after process kill; diagnose before resuming.
    lock=out/'RUNNING.lock'
    with lock.open('x') as f: f.write(str(os.getpid()))
    invocation=str(time.time_ns()); save(out/('command_'+invocation+'.json'),dict(argv=sys.argv,executable=sys.executable,platform=platform.platform(),versions=versions))
    gen=rdFingerprintGenerator.GetMorganGenerator(**cfg['fingerprint']); fp={}; bitvectors={}
    def matrix(frame):
        for sm in frame.smiles:
            if sm not in fp:
                mol=Chem.MolFromSmiles(sm); assert mol is not None; bitvectors[sm]=gen.GetFingerprint(mol); fp[sm]=gen.GetFingerprintAsNumPy(mol)
        return np.stack([fp[sm] for sm in frame.smiles])
    results=[]; count=0
    try:
        for r in selected.itertuples():
            k=tuple(getattr(r,z) for z in KEY); train=frames[k+('train',)]; test=frames[k+('test',)]
            for modelname in (['rf','ridge','tanimoto1nn'] if args.model=='all' else [args.model]):
                for seed in (cfg['seeds'] if modelname=='rf' else [0]):
                    if args.limit_cases and count>=args.limit_cases: break
                    modelcfg=(dict(cfg['rf'],n_jobs=2,random_state=seed) if modelname=='rf' else dict(alpha=1.0,fit_intercept=True,solver='lsqr',tol=1e-4,max_iter=10000) if modelname=='ridge' else dict(tie='first_original_train_row',zero_union_similarity=1))
                    ident=dict(model=modelname,seed=seed,training_signature=sig(train),config=modelcfg)
                    job=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest(); modelpath=out/(job+'.joblib'); ledger=out/(job+'.fit.json')
                    X=matrix(train); y=train.target.to_numpy(float)
                    if ledger.exists():
                        e=json.loads(ledger.read_text(encoding='utf8')); assert e['state']=='complete' and e['identity']==ident and sha(modelpath)==e['model_sha256']; model=joblib.load(modelpath)
                    else:
                        assert not modelpath.exists(); e=dict(state='running',identity=ident,started=time.time(),training_rows=len(train),actual_fit_executed=False); save(ledger,e)
                        model=RandomForestRegressor(**modelcfg) if modelname=='rf' else Ridge(**modelcfg) if modelname=='ridge' else ([bitvectors[sm] for sm in train.smiles],y)
                        if modelname!='tanimoto1nn': model.fit(csr_matrix(X,dtype=np.float64) if modelname=='ridge' else X,y)
                        joblib.dump(model,modelpath,compress=3); e.update(state='complete',finished=time.time(),actual_fit_executed=modelname!='tanimoto1nn',index_built=modelname=='tanimoto1nn',model_sha256=sha(modelpath)); save(ledger,e)
                    name='_'.join(k)+(f'_{modelname}_{seed}'); dest=out/(name+'.csv'); done=out/(name+'.complete.json')
                    if done.exists():
                        rec=json.loads(done.read_text(encoding='utf8')); assert sha(dest)==rec['prediction_sha256']; results.append(rec); count+=1; continue
                    assert not dest.exists(), 'Unledgered prediction: diagnose before resume'
                    T=matrix(test)
                    if modelname=='tanimoto1nn':
                        pred=np.asarray([model[1][int(np.argmax(DataStructs.BulkTanimotoSimilarity(bitvectors[sm],model[0]))) ] for sm in test.smiles])
                    else:
                        pred=model.predict(csr_matrix(T,dtype=np.float64) if modelname=='ridge' else T)
                    mean=float(y.mean()); actual=test.target.to_numpy(float); assert np.isfinite(pred).all()
                    f=test.copy(); f['prediction']=pred; f['training_mean']=mean; f['absolute_model_error']=abs(actual-pred); f['absolute_baseline_error']=abs(actual-mean); f['squared_model_error']=(actual-pred)**2; f['squared_baseline_error']=(actual-mean)**2; f.to_csv(dest,index=False)
                    rec=dict(zip(KEY,k)); rec.update(model=modelname,seed=seed,job_key=job,test_signature=sig(test),prediction_sha256=sha(dest),prediction_file=dest.name,**metrics(actual,pred,mean)); save(done,rec); results.append(rec); count+=1
                    print('COMPLETE',name,flush=True)
        # Immutable per-invocation summary; assemble all checkpointed cases on resume.
        records=[json.loads(p.read_text(encoding='utf8')) for p in sorted(out.glob('*.complete.json'))]
        df=pd.DataFrame(records); df.to_csv(out/('metrics_'+invocation+'.csv'),index=False)
        tg=df.groupby(KEY+['model'],as_index=False)[['mse_model','mse_mean','mae_model','mae_mean','rmse_model','rmse_mean','D','skill']].median(); tg.to_csv(out/('target_metrics_'+invocation+'.csv'),index=False)
        paired=tg.merge(tg,on=KEY,suffixes=('_a','_b')); paired=paired[paired.model_a<paired.model_b]; paired['skill_a_minus_b']=paired.skill_a-paired.skill_b; paired['mse_a_minus_b']=paired.mse_model_a-paired.mse_model_b; paired.to_csv(out/('paired_model_contrasts_'+invocation+'.csv'),index=False)
        routes=tg.pivot(index=['view','scenario','target','model'],columns='route',values='skill'); routes.to_csv(out/('route_skills_'+invocation+'.csv'))
        for a,b in [('temporal_exact','source_exact'),('source_scaffold','source_exact'),('temporal_scaffold','temporal_exact')]:
            if a in routes and b in routes: routes[a+'_minus_'+b]=routes[a]-routes[b]
        routes.to_csv(out/('paired_route_contrasts_'+invocation+'.csv'))
        assert load_inputs(root)[3]==hashes
        expected_cases=len(selected)*5
        save(out/('status_'+invocation+'.json'),dict(state='full_model_extension_complete' if len(records)==expected_cases else 'partial_cases_complete',completed_cases=len(records),expected_cases=expected_cases,selected_contexts=len(selected),primary_contexts=len(selected),include_paired=args.include_paired,independent_target_validation=True,public_release=False,input_hashes_unchanged=True))
    finally:
        lock.unlink()
if __name__=='__main__': main()
