from pathlib import Path
import json,hashlib,shutil,datetime,sys
import numpy as np
import pandas as pd
R=P=PREP=None
T='alk braf egfr erbb2 kdr met ache dpp4 ptgs2 hdac1 pde4d'.split()
ROUTES=['source_exact','source_scaffold','temporal_exact','temporal_scaffold']
MODELS=[('rf',42),('rf',53),('rf',67),('ridge',0),('tanimoto1nn',0)]
K=['panel','view','scenario','target','route']
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v): Path(p).write_text(json.dumps(v,indent=2,ensure_ascii=False),encoding='utf8')
def sig(f): return hashlib.sha256(f[['record_uid','smiles','target']].to_csv(index=False,lineterminator='\r\n').encode()).hexdigest()
def csv(p): return pd.read_csv(p,float_precision='round_trip')
def contexts():
 for t in T:
  panel='kinase' if t in T[:6] else 'nonkinase'
  prep=PREP/('current_prepared' if panel=='kinase' else 'external_prepared')
  man=pd.read_csv(prep/'split_manifest.csv')
  for route in ROUTES:
   key=dict(zip(K,[panel,'main','revised',t,route])); paths={}
   for part in ['train','test']:
    row=man[(man.view=='main')&(man.scenario=='revised')&(man.target==t)&(man.route==route)&(man.partition==part)].iloc[0]
    path=prep/'data/main'/str(row.corrected_csv).replace('\\','/').split('/')[-1]
    paths[part]=(path,row)
   preds={(m,s):P/'results/predictions'/panel/f'main_revised_{t}_{route}_{m}_{s}.csv' for m,s in MODELS}
   yield key,prep,paths,preds

def init():
 R.mkdir(exist_ok=True); (R/'analysis').mkdir(exist_ok=True); (R/'checkpoints').mkdir(exist_ok=True)
 if (R/'ANALYSIS_PLAN.json').exists(): return
 plan=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),panels=['kinase','nonkinase'],view='main',scenario='revised',targets=T,routes=ROUTES,contexts=44,prediction_files=220,models=MODELS,contrasts=[['temporal_exact','source_exact'],['source_scaffold','source_exact'],['temporal_scaffold','temporal_exact']],bins=[0,.2,.4,.6,.8,1],bin_closure='left-closed, final includes 1; empty n=0 metrics blank',weighting='equal contexts; rows equally weighted within context; no molecule-equal sensitivity',rf_aggregation='mean row squared loss per seed, then metric-wise median over seeds',fingerprint=dict(radius=2,fpSize=2048,includeChirality=False),nearest='ordered training, first maximum; RDKit BulkTanimotoSimilarity; empty union retains original similarity 1',prohibitions='no fitting, tuning, bootstrap, selection, significance, causal attribution or external mutation',label_std_ddof=0)
 save(R/'ANALYSIS_PLAN.json',plan)
 hashes={}
 for key,prep,paths,preds in contexts():
  for p in [prep/x for x in ['split_manifest.csv','model_config.json','decision_record.json']]+[v[0] for v in paths.values()]+list(preds.values()): hashes[str(p)]=sha(p)
 for p in (P/'results/analysis').glob('*.csv'): hashes[str(p)]=sha(p)
 save(R/'INPUT_SHA256.json',hashes)
 save(R/'WORK_STATUS.json',dict(state='plan_and_hashes_frozen_before_computation',plan_sha256=sha(R/'ANALYSIS_PLAN.json'),inputs_sha256=sha(R/'INPUT_SHA256.json')))
 save(R/'commands.json',dict(historical_windows_command_provenance=True,compute='E:/paper/kinase_followup_execution/.venv-win/Scripts/python.exe analysis_pipeline.py',verify='E:/paper/kinase_followup_execution/.venv-win/Scripts/python.exe verify_diagnostics.py',resume='Same compute command reuses hash-bound completed similarity checkpoints'))

def run():
 init()
 from rdkit import Chem,DataStructs,rdBase
 from rdkit.Chem import rdFingerprintGenerator
 assert rdBase.rdkitVersion=='2025.09.2'
 gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
 binding=dict(plan=sha(R/'ANALYSIS_PLAN.json'),inputs=sha(R/'INPUT_SHA256.json'))
 frozen=csv(P/'results/analysis/target_metrics.csv'); rows=[]; bins=[]; allsim=[]; metricmax=0
 for key,prep,paths,preds in contexts():
  name='_'.join(key.values()); frames={}
  for part,(path,rec) in paths.items():
   assert sha(path)==rec.corrected_sha256
   # Original ordered signatures use default pandas parsing.
   assert sig(pd.read_csv(path))==rec.corrected_ordered_signature
   # Preserve the original runner's source-label parser; round-trip applies to saved outputs.
   frames[part]=pd.read_csv(path); assert len(frames[part])==rec.corrected_rows and frames[part].record_uid.is_unique
  tr,te=frames['train'],frames['test']; y=te.target.to_numpy(); mean=tr.target.mean()
  fs={}
  for ms,p in preds.items():
   f=csv(p); assert list(f.record_uid)==list(te.record_uid) and list(f.smiles)==list(te.smiles)
   assert np.array_equal(f.target.to_numpy(),y); assert np.max(abs(f.training_mean-mean))<1e-12
   fs[ms]=f
  cp=R/'checkpoints'/f'{name}.csv'; done=cp.with_suffix('.json')
  if done.exists():
   d=json.loads(done.read_text()); assert d['binding']==binding and d['sha256']==sha(cp); sim=csv(cp)
  else:
   tf=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in tr.smiles]; cache={}; out=[]
   for uid,s in zip(te.record_uid,te.smiles):
    if s not in cache:
     fp=gen.GetFingerprint(Chem.MolFromSmiles(s)); a=np.array(DataStructs.BulkTanimotoSimilarity(fp,tf))
     # Original runner defines empty union as 1; RDKit empty/empty returns 1 in pinned version.
     if fp.GetNumOnBits()==0:
      for j,v in enumerate(tf):
       if v.GetNumOnBits()==0: a[j]=1.0
     j=int(a.argmax()); cache[s]=(float(a[j]),j,fp.GetNumOnBits()==0)
    val,j,zero=cache[s]; out.append(dict(**key,record_uid=uid,max_tanimoto=val,nearest_train_record_uid=tr.record_uid.iloc[j],nearest_label=tr.target.iloc[j],test_zero_fingerprint=zero,nearest_train_zero_fingerprint=tf[j].GetNumOnBits()==0))
   sim=pd.DataFrame(out); sim.to_csv(cp,index=False)
   save(done,dict(binding=binding,sha256=sha(cp),zero_training_fingerprints=sum(f.GetNumOnBits()==0 for f in tf),zero_test_fingerprints=int(sim.test_zero_fingerprint.sum())))
  assert list(sim.record_uid)==list(te.record_uid)
  # Structural checkpoints do not depend on label parsing. Resolve labels from original ordered inputs.
  lookup=tr.set_index('record_uid').target
  sim=sim.copy(); sim['nearest_label']=sim.nearest_train_record_uid.map(lookup)
  assert sim.nearest_label.notna().all()
  assert np.max(abs(sim.nearest_label-fs[('tanimoto1nn',0)].prediction))<=1e-12
  allsim.append(sim); base=float(np.mean((y-mean)**2)); var=float(np.var(y)); shift=float((y.mean()-mean)**2); assert abs(base-var-shift)<1e-12
  row=dict(**key,train_n=len(tr),test_n=len(te),train_unique_molecules=tr.smiles.nunique(),test_unique_molecules=te.smiles.nunique(),train_documents=tr.document_chembl_id.nunique(),test_documents=te.document_chembl_id.nunique(),train_label_mean=mean,test_label_mean=y.mean(),train_label_std=tr.target.std(ddof=0),test_label_std=np.std(y),test_variance=var,squared_mean_shift=shift,baseline_mse=base,similarity_median=sim.max_tanimoto.median(),similarity_p25=sim.max_tanimoto.quantile(.25),similarity_p75=sim.max_tanimoto.quantile(.75),similarity_ge_08_fraction=(sim.max_tanimoto>=.8).mean())
  bi=np.minimum(np.searchsorted([0,.2,.4,.6,.8,1],sim.max_tanimoto,side='right')-1,4)
  for model in ['rf','ridge','tanimoto1nn']:
   losses=[(y-f.prediction.to_numpy())**2 for (m,s),f in fs.items() if m==model]
   mse=float(np.median([v.mean() for v in losses])); skill=float(np.median([1-v.mean()/base for v in losses])); row[model+'_mse']=mse; row[model+'_skill']=skill
   fr=frozen.copy()
   for k,v in dict(**key,model=model).items(): fr=fr[fr[k]==v]
   assert len(fr)==1
   for col,val in [('mse_model',mse),('mse_mean',base),('skill',skill)]:
    err=abs(fr.iloc[0][col]-val); metricmax=max(metricmax,err); assert err<=1e-12,(key,col,err)
   for b in range(5):
    mask=bi==b; n=int(mask.sum()); bm=float(np.mean((y[mask]-mean)**2)) if n else np.nan
    mm=float(np.median([v[mask].mean() for v in losses])) if n else np.nan
    bins.append(dict(**key,model=model,bin_index=b,bin_lower=[0,.2,.4,.6,.8,1][b],bin_upper=[0,.2,.4,.6,.8,1][b+1],n=n,model_mse=mm,baseline_mse=bm,skill=1-mm/bm if n and bm>0 else np.nan))
  rows.append(row)
  save(R/'WORK_STATUS.json',dict(state='computing',completed_contexts=len(rows),binding=binding,last_context=name)); print(name,'complete',flush=True)
 ctx=pd.DataFrame(rows); ctx.to_csv(R/'analysis/context_diagnostics.csv',index=False); pd.concat(allsim).to_csv(R/'analysis/row_similarity.csv',index=False); pd.DataFrame(bins).to_csv(R/'analysis/fixed_bin_diagnostics.csv',index=False)
 pairs=[]
 for t in T:
  for route,ref in json.loads((R/'ANALYSIS_PLAN.json').read_text())['contrasts']:
   a=ctx[(ctx.target==t)&(ctx.route==route)].iloc[0]; b=ctx[(ctx.target==t)&(ctx.route==ref)].iloc[0]
   z=dict(panel=a.panel,view='main',scenario='revised',target=t,route=route,reference_route=ref,train_n_ratio=a.train_n/b.train_n,test_n_ratio=a.test_n/b.test_n,similarity_change=a.similarity_median-b.similarity_median,baseline_mse_route=a.baseline_mse,baseline_mse_reference=b.baseline_mse,baseline_log_change=np.log(a.baseline_mse/b.baseline_mse))
   for m in ['rf','ridge','tanimoto1nn']:
    z.update({m+'_mse_route':a[m+'_mse'],m+'_mse_reference':b[m+'_mse'],m+'_skill_route':a[m+'_skill'],m+'_skill_reference':b[m+'_skill'],m+'_skill_change':a[m+'_skill']-b[m+'_skill'],m+'_mse_log_change':np.log(a[m+'_mse']/b[m+'_mse'])})
    z[m+'_relative_error_log_change']=z[m+'_mse_log_change']-z['baseline_log_change']; z[m+'_log_identity_rhs']=np.log((1-a[m+'_skill'])/(1-b[m+'_skill']))
    assert abs(z[m+'_relative_error_log_change']-z[m+'_log_identity_rhs'])<1e-12
   pairs.append(z)
 pd.DataFrame(pairs).to_csv(R/'analysis/target_route_diagnostics.csv',index=False)
 save(R/'analysis/computation_checks.json',dict(contexts=44,prediction_files=220,metric_max_abs_error=metricmax,binding=binding))
 save(R/'WORK_STATUS.json',dict(state='analysis_complete',completed_contexts=44,binding=binding))
if __name__=='__main__':
 import argparse
 ap=argparse.ArgumentParser(description='Frozen-prediction diagnostics: no fitting, tuning or bootstrap')
 ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--prepared',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--resume',action='store_true')
 args=ap.parse_args();P=args.bundle.resolve();PREP=args.prepared.resolve();R=args.output.resolve()
 assert R!=P and R!=PREP and P not in R.parents and PREP not in R.parents
 if R.exists() and not args.resume:raise FileExistsError('Use a new output or explicitly resume hash-bound similarity checkpoints')
 import rdkit,sklearn
 assert (np.__version__,pd.__version__,rdkit.__version__,sklearn.__version__)==('2.0.2','2.3.3','2025.09.2','1.6.1')
 run()

