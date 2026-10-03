"""Frozen common-test, equal-row-budget control. No tuning or outcome selection."""
import os
for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS'):
    os.environ[k]='1'
from pathlib import Path
import json, hashlib, datetime, argparse, sys, time, platform
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
PREP=ROOT.parent/'manuscript_provenance_closure_20261003_v1/03_representative_replay/fresh'
PUB=ROOT.parent/'algorithm_validation_extension_20261003_v1/public_release_v1'
BASE=ROOT.parent/'holdout_shift_diagnostics_20261003_v1'
TARGETS='alk braf egfr erbb2 kdr met ache dpp4 ptgs2 hdac1 pde4d'.split()
ROUTES=['source_exact','temporal_exact']
DRAWS=[101,211,307,401,503]
MS=[('rf',42),('rf',53),('rf',67),('ridge',0),('tanimoto1nn',0)]
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,obj):
    p=Path(p); tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf8'); tmp.replace(p)
def read(p): return json.loads(Path(p).read_text(encoding='utf8'))
def csv(p): return pd.read_csv(p,float_precision='round_trip')
def sig(f): return hashlib.sha256(f[['record_uid','smiles','target']].to_csv(index=False,lineterminator='\r\n').encode()).hexdigest()
def writecsv(f,p): f.to_csv(p,index=False,lineterminator='\r\n')
def stamp(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def canon(values):
    from rdkit import Chem
    out=[]
    for s in values:
        mol=Chem.MolFromSmiles(s); assert mol is not None
        out.append(Chem.MolToSmiles(mol,canonical=True,isomericSmiles=True))
    return set(out)
def environment():
    import rdkit, sklearn
    e={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'rdkit':rdkit.__version__,'sklearn':sklearn.__version__}
    assert e==read(ROOT/'ANALYSIS_PLAN.json')['environment'],e
    return dict(versions=e,executable=sys.executable,platform=platform.platform())

def prepare():
    assert not (ROOT/'PREPARED.json').exists(),'Preparation already complete; do not repeat'
    assert not (ROOT/'SOURCE_HASHES.json').exists(),'Partial preparation: inspect before restarting'
    environment()
    for name in ['inputs','predictions','models','ledgers','analysis','checkpoints']:
        (ROOT/name).mkdir(exist_ok=True)
    hashes={}; originals={}; counts=[]
    for i,t in enumerate(TARGETS):
        panel='kinase' if i<6 else 'nonkinase'; prep=PREP/('current_prepared' if i<6 else 'external_prepared')
        man=pd.read_csv(prep/'split_manifest.csv')
        for n in ['split_manifest.csv','model_config.json','decision_record.json']:
            hashes[str(prep/n)]=sha(prep/n)
        for route in ROUTES:
            cfg=read(prep/'model_config.json')
            for part in ['train','test']:
                rec=man[(man.view=='main')&(man.scenario=='revised')&(man.target==t)&(man.route==route)&(man.partition==part)]
                assert len(rec)==1; rec=rec.iloc[0]
                p=prep/'data/main'/str(rec.corrected_csv).replace('\\','/').rsplit('/',1)[-1]
                f=pd.read_csv(p); assert sha(p)==rec.corrected_sha256 and sig(f)==rec.corrected_ordered_signature
                assert len(f)==rec.corrected_rows and f.record_uid.is_unique and np.isfinite(f.target).all()
                hashes[str(p)]=sha(p); originals[(t,route,part)]=f
            for model,seed in MS:
                p=PUB/'results/predictions'/panel/f'main_revised_{t}_{route}_{model}_{seed}.csv'
                f=csv(p); te=originals[t,route,'test']; tr=originals[t,route,'train']
                assert list(f.record_uid)==list(te.record_uid) and list(f.smiles)==list(te.smiles)
                assert np.array_equal(f.target.to_numpy(),te.target.to_numpy())
                assert np.max(abs(f.training_mean-tr.target.mean()))<=1e-12
                hashes[str(p)]=sha(p); originals[t,route,model,seed]=(f,p)
    for p in (BASE/'02_editable_source').rglob('*'):
        if p.is_file(): hashes[str(p)]=sha(p)
    save(ROOT/'SOURCE_HASHES.json',hashes)
    save(ROOT/'ENVIRONMENT.json',environment())
    cfg=read(PREP/'current_prepared/model_config.json')
    assert cfg['fingerprint']==read(PREP/'external_prepared/model_config.json')['fingerprint']
    assert cfg['rf']==read(PREP/'external_prepared/model_config.json')['rf']
    save(ROOT/'MODEL_CONFIG.json',cfg)
    binding={'plan_sha256':sha(ROOT/'ANALYSIS_PLAN.json'),'source_hashes_sha256':sha(ROOT/'SOURCE_HASHES.json'),'runner_sha256':sha(__file__),'config_sha256':sha(ROOT/'MODEL_CONFIG.json')}
    save(ROOT/'BINDING.json',binding); manifest=[]; ih={}; reused={}
    for i,t in enumerate(TARGETS):
        panel='kinase' if i<6 else 'nonkinase'
        source=originals[t,'source_exact','test']; year=originals[t,'temporal_exact','test']
        common=source[source.record_uid.isin(year.record_uid)].copy().reset_index(drop=True)
        other=year.set_index('record_uid').loc[common.record_uid].reset_index()
        for c in source.columns:
            assert np.array_equal(common[c].to_numpy(),other[c].to_numpy()),(t,c)
        testfile=f'inputs/common_test_{t}.csv'; writecsv(common,ROOT/testfile); ih[testfile]=sha(ROOT/testfile)
        nsource=len(originals[t,'source_exact','train']); nyear=len(originals[t,'temporal_exact','train']); budget=min(nsource,nyear)
        counts.append(dict(panel=panel,target=t,source_full_n=nsource,year_full_n=nyear,matched_n=budget,source_original_test_n=len(source),year_original_test_n=len(year),common_test_n=len(common),common_test_molecules=common.smiles.nunique(),common_test_documents=common.document_chembl_id.nunique(),source_test_coverage=len(common)/len(source),year_test_coverage=len(common)/len(year)))
        for route in ROUTES:
            full=originals[t,route,'train']
            assert not set(common.record_uid)&set(full.record_uid)
            assert not set(common.document_chembl_id)&set(full.document_chembl_id)
            assert not canon(common.smiles)&canon(full.smiles)
            prefix='source' if route=='source_exact' else 'year'
            fullfile=f'inputs/train_{t}_{route}_full.csv'; writecsv(full,ROOT/fullfile); ih[fullfile]=sha(ROOT/fullfile)
            for draw in DRAWS:
                pos=np.sort(np.random.RandomState(draw).choice(len(full),size=budget,replace=False)) if len(full)>budget else np.arange(len(full))
                sample=full.iloc[pos].copy().reset_index(drop=True)
                samplefile=f'inputs/train_{t}_{route}_s{draw}.csv'; writecsv(sample,ROOT/samplefile); ih[samplefile]=sha(ROOT/samplefile)
                for arm,trainfile,trainframe,is_reuse in [(prefix+'_full',fullfile,full,True),(prefix+'_matched',samplefile,sample,len(full)==budget)]:
                    for model,seed in MS:
                        predfile=f'predictions/{t}_{route}_'+('full' if is_reuse else f's{draw}')+f'_{model}_{seed}.csv'
                        rec=dict(panel=panel,target=t,route=route,arm=arm,draw_seed=draw,model=model,model_seed=seed,train_file=trainfile,train_sha256=ih[trainfile],train_signature=sig(trainframe),train_n=len(trainframe),test_file=testfile,test_sha256=ih[testfile],test_signature=sig(common),test_n=len(common),prediction_file=predfile,reused_frozen_prediction=is_reuse,actual_fit_expected=not is_reuse and model!='tanimoto1nn',index_build_expected=not is_reuse and model=='tanimoto1nn')
                        if is_reuse and predfile not in reused:
                            original,p=originals[t,route,model,seed]
                            f=original.set_index('record_uid').loc[common.record_uid].reset_index()
                            writecsv(f,ROOT/predfile)
                            ledger=dict(state='complete',binding=binding,reused_frozen_prediction=True,actual_fit=False,index_built=False,source_prediction=str(p),source_prediction_sha256=sha(p),prediction_file=predfile,prediction_sha256=sha(ROOT/predfile),train_signature=sig(full),test_signature=sig(common),model=model,model_seed=seed,completed_utc=stamp())
                            save(ROOT/'ledgers'/Path(predfile).with_suffix('.complete.json').name,ledger); reused[predfile]=ledger
                        manifest.append(rec)
    writecsv(pd.DataFrame(manifest),ROOT/'manifest.csv'); writecsv(pd.DataFrame(counts),ROOT/'POOL_COUNTS.csv')
    save(ROOT/'INPUT_HASHES.json',ih)
    save(ROOT/'PREPARED.json',dict(binding=binding,manifest_sha256=sha(ROOT/'manifest.csv'),input_hashes_sha256=sha(ROOT/'INPUT_HASHES.json'),created_utc=stamp(),targets=11,logical_cases=len(manifest),unique_prediction_cases=len(set(r['prediction_file'] for r in manifest)),reused_prediction_cases=len(reused),new_prediction_cases=len(set(r['prediction_file'] for r in manifest))-len(reused),new_rf_fits=165,new_ridge_fits=55,new_nn_indices=55))
    save(ROOT/'WORK_STATUS.json',dict(state='prepared_no_new_model_executed',created_utc=stamp(),binding=binding))
    print(json.dumps(read(ROOT/'PREPARED.json'),indent=2),flush=True)

def validate_binding():
    b=read(ROOT/'BINDING.json'); p=read(ROOT/'PREPARED.json')
    assert b==p['binding']; assert b['plan_sha256']==sha(ROOT/'ANALYSIS_PLAN.json') and b['runner_sha256']==sha(__file__) and b['config_sha256']==sha(ROOT/'MODEL_CONFIG.json')
    assert sha(ROOT/'manifest.csv')==p['manifest_sha256'] and sha(ROOT/'INPUT_HASHES.json')==p['input_hashes_sha256']
    for n,h in read(ROOT/'INPUT_HASHES.json').items(): assert sha(ROOT/n)==h,n
    return b

def run():
    from rdkit import Chem,DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import Ridge
    from scipy.sparse import csr_matrix
    import joblib
    b=validate_binding(); environment(); cfg=read(ROOT/'MODEL_CONFIG.json')
    lock=ROOT/'RUNNING.lock'; lock.open('x').write(str(os.getpid()))
    save(ROOT/f'command_{time.time_ns()}.json',dict(argv=sys.argv,pid=os.getpid(),started_utc=stamp(),**environment()))
    gen=rdFingerprintGenerator.GetMorganGenerator(**cfg['fingerprint']); arrays={}; bits={}; frames={}; diagnostic=[]
    def frame(p):
        if p not in frames: frames[p]=csv(ROOT/p)
        return frames[p]
    def matrix(f):
        for s in f.smiles:
            if s not in arrays:
                mol=Chem.MolFromSmiles(s); assert mol is not None
                bits[s]=gen.GetFingerprint(mol); arrays[s]=gen.GetFingerprintAsNumPy(mol)
        return np.stack([arrays[s] for s in f.smiles])
    m=csv(ROOT/'manifest.csv').drop_duplicates('prediction_file'); done_count=0
    try:
        for r in m.itertuples():
            dest=ROOT/r.prediction_file; done=ROOT/'ledgers'/dest.with_suffix('.complete.json').name
            if done.exists():
                e=read(done); assert e['state']=='complete' and e['binding']==b and sha(dest)==e['prediction_sha256']
                done_count+=1; continue
            assert not dest.exists(),'Unledgered prediction must be diagnosed'
            tr=frame(r.train_file); te=frame(r.test_file)
            assert sig(tr)==r.train_signature and sig(te)==r.test_signature
            modelcfg=dict(cfg['rf'],n_jobs=2,random_state=int(r.model_seed)) if r.model=='rf' else dict(alpha=1.0,fit_intercept=True,solver='lsqr',tol=1e-4,max_iter=10000) if r.model=='ridge' else dict(tie='first_ordered_train_row',zero_union_similarity=1)
            ident=dict(training_signature=r.train_signature,train_sha256=r.train_sha256,model=r.model,model_seed=int(r.model_seed),config=modelcfg)
            key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()
            modelpath=ROOT/'models'/f'{key}.joblib'; fitpath=ROOT/'ledgers'/f'{key}.fit.json'
            X=matrix(tr); T=matrix(te); y=tr.target.to_numpy(float)
            if fitpath.exists():
                e=read(fitpath); assert e['state']=='complete' and e['binding']==b and e['identity']==ident and sha(modelpath)==e['model_sha256'],'Incomplete fit: do not repeat'
                model=joblib.load(modelpath)
            else:
                assert not modelpath.exists(),'Unledgered model: do not overwrite'
                e=dict(state='running',binding=b,identity=ident,started_utc=stamp(),actual_fit=False,index_built=False,pid=os.getpid())
                save(fitpath,e)
                model=RandomForestRegressor(**modelcfg) if r.model=='rf' else Ridge(**modelcfg) if r.model=='ridge' else ([bits[s] for s in tr.smiles],y)
                if r.model!='tanimoto1nn': model.fit(csr_matrix(X,dtype=np.float64) if r.model=='ridge' else X,y)
                joblib.dump(model,modelpath,compress=3)
                e.update(state='complete',completed_utc=stamp(),actual_fit=r.model!='tanimoto1nn',index_built=r.model=='tanimoto1nn',model_sha256=sha(modelpath),model_file=str(modelpath.relative_to(ROOT)))
                save(fitpath,e)
            if r.model=='tanimoto1nn':
                pred=np.array([model[1][int(np.argmax(DataStructs.BulkTanimotoSimilarity(bits[s],model[0])))] for s in te.smiles])
            else: pred=model.predict(csr_matrix(T,dtype=np.float64) if r.model=='ridge' else T)
            assert np.isfinite(pred).all()
            f=te.copy(); f['prediction']=pred; f['training_mean']=float(y.mean())
            f['squared_model_error']=(f.target-pred)**2; f['squared_baseline_error']=(f.target-y.mean())**2
            writecsv(f,dest)
            save(done,dict(state='complete',binding=b,prediction_file=r.prediction_file,prediction_sha256=sha(dest),reused_frozen_prediction=False,actual_fit=e['actual_fit'],index_built=e['index_built'],fit_key=key,fit_ledger=str(fitpath.relative_to(ROOT)),train_signature=r.train_signature,test_signature=r.test_signature,model=r.model,model_seed=int(r.model_seed),completed_utc=stamp()))
            done_count+=1; print('COMPLETE',dest.name,done_count,'/',len(m),flush=True)
            save(ROOT/'WORK_STATUS.json',dict(state='running',completed_unique_cases=done_count,expected_unique_cases=len(m),last_prediction=r.prediction_file,binding=b,updated_utc=stamp()))
        # Structural diagnostics: completed pool checkpoints avoid repeated large similarity work.
        dm=csv(ROOT/'manifest.csv'); dm=dm[dm.arm.str.endswith('matched')].drop_duplicates(['target','arm','draw_seed'])
        for r in dm.itertuples():
            cp=ROOT/'checkpoints'/f'{r.target}_{r.arm}_s{r.draw_seed}.json'
            if cp.exists():
                q=read(cp); assert q['binding']==b and q['train_sha256']==r.train_sha256 and q['test_sha256']==r.test_sha256
            else:
                tr=frame(r.train_file); te=frame(r.test_file); matrix(tr); matrix(te)
                tbits=[bits[s] for s in tr.smiles]; cache={}
                for s in te.smiles:
                    if s not in cache: cache[s]=float(max(DataStructs.BulkTanimotoSimilarity(bits[s],tbits)))
                vals=np.array([cache[s] for s in te.smiles]); years=tr.document_year.to_numpy(float)
                q=dict(binding=b,train_sha256=r.train_sha256,test_sha256=r.test_sha256,panel=r.panel,target=r.target,arm=r.arm,draw_seed=int(r.draw_seed),train_n=len(tr),test_n=len(te),train_molecules=int(tr.smiles.nunique()),train_documents=int(tr.document_chembl_id.nunique()),train_label_mean=float(tr.target.mean()),train_label_std=float(tr.target.std(ddof=0)),train_year_mean=float(years.mean()),train_year_std=float(years.std()),median_max_tanimoto=float(np.median(vals)),fraction_ge_08=float((vals>=.8).mean()))
                save(cp,q)
            diagnostic.append({k:v for k,v in q.items() if k!='binding'})
        writecsv(pd.DataFrame(diagnostic),ROOT/'analysis/structure_diagnostics.csv')
        save(ROOT/'WORK_STATUS.json',dict(state='model_execution_complete',completed_unique_cases=done_count,expected_unique_cases=len(m),binding=b,updated_utc=stamp()))
    finally:
        lock.unlink()

def summarize():
    b=validate_binding(); m=csv(ROOT/'manifest.csv'); rows=[]
    for r in m.itertuples():
        p=ROOT/r.prediction_file; e=read(ROOT/'ledgers'/p.with_suffix('.complete.json').name)
        assert e['state']=='complete' and sha(p)==e['prediction_sha256'] and e['binding']==b
        f=csv(p); tr=csv(ROOT/r.train_file); y=f.target.to_numpy(); pr=f.prediction.to_numpy(); mean=float(tr.target.mean())
        assert np.max(abs(f.training_mean-mean))<=1e-12
        mse=float(np.mean((y-pr)**2)); base=float(np.mean((y-mean)**2)); var=float(np.var(y))
        rows.append(dict(panel=r.panel,target=r.target,arm=r.arm,draw_seed=r.draw_seed,model=r.model,model_seed=r.model_seed,train_n=r.train_n,test_n=r.test_n,mse=mse,baseline_mse=base,test_variance=var,normalized_mse=mse/var if var else np.nan,training_mean_skill=1-mse/base if base else np.nan,reused_frozen_prediction=e['reused_frozen_prediction'],actual_fit=e['actual_fit'],index_built=e['index_built'],prediction_file=r.prediction_file,prediction_sha256=e['prediction_sha256']))
    runs=pd.DataFrame(rows); writecsv(runs,ROOT/'analysis/run_metrics.csv')
    keys=['panel','target','arm','draw_seed','model']
    arms=runs.groupby(keys,sort=False,as_index=False)[['train_n','test_n','mse','baseline_mse','test_variance','normalized_mse','training_mean_skill']].median()
    writecsv(arms,ROOT/'analysis/arm_metrics.csv'); contrasts=[]
    for t in TARGETS:
        for model in ['rf','ridge','tanimoto1nn']:
            for d in DRAWS:
                a=arms[(arms.target==t)&(arms.model==model)&(arms.draw_seed==d)].set_index('arm')
                sf,sm,yf,ym=[a.loc[n] for n in ['source_full','source_matched','year_full','year_matched']]
                row=dict(panel='kinase' if t in TARGETS[:6] else 'nonkinase',target=t,model=model,draw_seed=d,test_n=int(sf.test_n),matched_n=int(sm.train_n),test_variance=sf.test_variance,source_full_mse=sf.mse,source_matched_mse=sm.mse,year_full_mse=yf.mse,year_matched_mse=ym.mse,source_size=sm.mse-sf.mse,year_size=yf.mse-ym.mse,matched_pool=ym.mse-sm.mse,full_gap=yf.mse-sf.mse,source_matched_skill=sm.training_mean_skill,year_matched_skill=ym.training_mean_skill)
                for name in ['source_size','year_size','matched_pool','full_gap']:
                    row['normalized_'+name]=row[name]/sf.test_variance if sf.test_variance else np.nan
                row['identity_error']=row['full_gap']-row['source_size']-row['year_size']-row['matched_pool']
                assert abs(row['identity_error'])<=1e-12
                contrasts.append(row)
    con=pd.DataFrame(contrasts); writecsv(con,ROOT/'analysis/per_draw_contrasts.csv')
    numeric=[c for c in con.columns if c not in ['panel','target','model','draw_seed']]
    sums=[]
    for t in TARGETS:
        for model in ['rf','ridge','tanimoto1nn']:
            f=con[(con.target==t)&(con.model==model)]
            row=dict(panel=f.iloc[0].panel,target=t,model=model)
            for col in numeric:
                row[col+'_median']=float(f[col].median()); row[col+'_min']=float(f[col].min()); row[col+'_max']=float(f[col].max())
            sums.append(row)
    summary=pd.DataFrame(sums); writecsv(summary,ROOT/'analysis/target_summary.csv')
    structure=csv(ROOT/'analysis/structure_diagnostics.csv'); sc=[]
    for t in TARGETS:
        a=structure[structure.target==t].pivot(index='draw_seed',columns='arm',values='median_max_tanimoto')
        vals=a.year_matched-a.source_matched
        sc.append(dict(target=t,matched_similarity_change_median=float(vals.median()),matched_similarity_change_min=float(vals.min()),matched_similarity_change_max=float(vals.max())))
    writecsv(pd.DataFrame(sc),ROOT/'analysis/structure_summary.csv')
    results=dict(created_utc=stamp(),common_test_rows=int(csv(ROOT/'POOL_COUNTS.csv').common_test_n.sum()),common_test_n_min=int(csv(ROOT/'POOL_COUNTS.csv').common_test_n.min()),common_test_n_max=int(csv(ROOT/'POOL_COUNTS.csv').common_test_n.max()),targets=11,draws=5,models=['rf','ridge','tanimoto1nn'],unique_prediction_cases=int(m.prediction_file.nunique()),logical_prediction_references=len(m),new_rf_fits=165,new_ridge_fits=55,new_nn_indices=55,reused_unique_prediction_cases=110,matched_pool_effects={})
    for model in results['models']:
        f=summary[summary.model==model]
        results['matched_pool_effects'][model]=dict(year_higher_mse=int((f.matched_pool_median>0).sum()),year_lower_mse=int((f.matched_pool_median<0).sum()),source_subsampling_higher_mse=int((f.source_size_median>0).sum()),median_normalized_matched_pool=float(f.normalized_matched_pool_median.median()),median_normalized_source_size=float(f.normalized_source_size_median.median()),median_normalized_full_gap=float(f.normalized_full_gap_median.median()))
    save(ROOT/'RESULT_SUMMARY.json',results)
    save(ROOT/'WORK_STATUS.json',dict(state='summarized_pending_independent_verification',updated_utc=stamp(),binding=b))
    print(json.dumps(results,indent=2),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--prepare',action='store_true'); ap.add_argument('--run',action='store_true'); ap.add_argument('--summarize',action='store_true'); ap.add_argument('--root',type=Path); args=ap.parse_args()
    if args.root: ROOT=args.root.resolve()
    assert sum([args.prepare,args.run,args.summarize])==1
    if args.prepare: prepare()
    if args.run: run()
    if args.summarize: summarize()
