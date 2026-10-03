"""Independent checks of frozen-prediction diagnostics; no fitting or resampling."""
from pathlib import Path
import json, hashlib, re
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

R=P=PREP=None
T='alk braf egfr erbb2 kdr met ache dpp4 ptgs2 hdac1 pde4d'.split()
K=['panel','view','scenario','target','route']
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def csv(p): return pd.read_csv(p,float_precision='round_trip')
def save(p,d): p.write_text(json.dumps(d,indent=2),encoding='utf8')

def main():
    frozen=json.loads((R/'INPUT_SHA256.json').read_text())
    for p,h in frozen.items(): assert sha(Path(p))==h,p
    c=csv(R/'analysis/context_diagnostics.csv'); s=csv(R/'analysis/row_similarity.csv')
    b=csv(R/'analysis/fixed_bin_diagnostics.csv'); p=csv(R/'analysis/target_route_diagnostics.csv')
    assert len(c)==44 and len(p)==33 and len(b)==660 and len(s)==51900
    assert not c.duplicated(K).any() and not s.duplicated(K+['record_uid']).any()
    assert set(c.target)==set(T) and c.groupby('target').size().eq(4).all()
    assert s.max_tanimoto.between(0,1).all()
    gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
    cache={}; sampled=0; pred_count=0; max_nn_error=0.; max_bin_error=0.
    for q in c.itertuples():
        key={k:getattr(q,k) for k in K}; mask=np.ones(len(s),bool); bm=np.ones(len(b),bool)
        for k,v in key.items(): mask&=s[k].eq(v); bm&=b[k].eq(v)
        ss=s.loc[mask].reset_index(drop=True); bb=b.loc[bm]
        prep=PREP/('current_prepared' if q.panel=='kinase' else 'external_prepared')
        man=pd.read_csv(prep/'split_manifest.csv'); frames={}
        for part in ['train','test']:
            z=man.query('view==@q.view and scenario==@q.scenario and target==@q.target and route==@q.route and partition==@part').iloc[0]
            path=prep/'data/main'/str(z.corrected_csv).replace('\\','/').rsplit('/',1)[-1]
            assert sha(path)==z.corrected_sha256;frames[part]=pd.read_csv(path)
        tr,te=frames['train'],frames['test']; assert len(ss)==len(te)==q.test_n and len(tr)==q.train_n
        assert list(ss.record_uid)==list(te.record_uid)
        yy=te.target.to_numpy(); mu=tr.target.mean(); base=(yy-mu)**2
        assert abs(base.mean()-q.baseline_mse)<1e-12
        assert abs(q.baseline_mse-q.test_variance-q.squared_mean_shift)<1e-12
        labels=tr.set_index('record_uid').target.reindex(ss.nearest_train_record_uid).to_numpy()
        assert np.isfinite(labels).all() and np.max(abs(labels-ss.nearest_label))<1e-12
        for model,seeds in [('rf',[42,53,67]),('ridge',[0]),('tanimoto1nn',[0])]:
            losses=[]
            for seed in seeds:
                f=csv(P/'results/predictions'/q.panel/f'main_revised_{q.target}_{q.route}_{model}_{seed}.csv')
                assert list(f.record_uid)==list(te.record_uid) and np.array_equal(f.target,yy)
                if model=='tanimoto1nn':
                    err=float(np.max(abs(f.prediction-labels)));max_nn_error=max(max_nn_error,err);assert err<1e-12
                losses.append((yy-f.prediction.to_numpy())**2);pred_count+=1
            mse=np.median([v.mean() for v in losses]); assert abs(mse-getattr(q,model+'_mse'))<1e-12
            assert abs(1-mse/base.mean()-getattr(q,model+'_skill'))<1e-12
            bz=bb[bb.model.eq(model)]; assert len(bz)==5 and bz.n.sum()==len(te)
            for z in bz.itertuples():
                mm=(ss.max_tanimoto>=z.bin_lower)&((ss.max_tanimoto<z.bin_upper) if z.bin_index<4 else (ss.max_tanimoto<=1))
                assert int(mm.sum())==z.n
                if z.n==0: assert pd.isna(z.model_mse) and pd.isna(z.baseline_mse) and pd.isna(z.skill)
                else:
                    actual=np.median([v[mm].mean() for v in losses]);e=abs(actual-z.model_mse);max_bin_error=max(max_bin_error,e);assert e<1e-12
                    assert abs(base[mm].mean()-z.baseline_mse)<1e-12
        # Independent integer bit-count implementation, including first-maximum tie rule.
        token=hashlib.sha256('\n'.join(tr.smiles).encode()).hexdigest()
        if token not in cache:
            mat=np.zeros((len(tr),2048),dtype=np.int32)
            for i,sm in enumerate(tr.smiles):mat[i,list(gen.GetFingerprint(Chem.MolFromSmiles(sm)).GetOnBits())]=1
            cache[token]=(mat,mat.sum(axis=1))
        mat,n=cache[token]
        for i in sorted(set([0,len(te)//2,len(te)-1])):
            bits=list(gen.GetFingerprint(Chem.MolFromSmiles(te.smiles.iloc[i])).GetOnBits())
            inter=mat[:,bits].sum(axis=1);union=n+len(bits)-inter
            sim=np.divide(inter,union,out=np.ones(len(tr),float),where=union!=0);j=int(np.argmax(sim))
            assert abs(sim[j]-ss.max_tanimoto.iloc[i])<1e-12
            assert tr.record_uid.iloc[j]==ss.nearest_train_record_uid.iloc[i];sampled+=1
    for model in ['rf','ridge','tanimoto1nn']:
        a=p[model+'_mse_log_change']-p.baseline_log_change
        rhs=np.log((1-p[model+'_skill_route'])/(1-p[model+'_skill_reference']))
        assert np.max(abs(a-rhs))<1e-12 and np.max(abs(a-p[model+'_relative_error_log_change']))<1e-12
    fm=csv(R/'analysis/Figure_6_data.csv');tm=csv(R/'analysis/Table_S21_data.csv')
    contrast=p[p.route.eq('temporal_exact')&p.reference_route.eq('source_exact')].set_index('target').loc[T]
    assert list(fm.target)==T and list(tm.target)==T
    for col in ['rf_mse_log_change','baseline_log_change']:
        assert np.array_equal(fm[col],contrast[col])
    assert np.array_equal(fm.relative_error_log_change,contrast.rf_relative_error_log_change)
    for col in tm.columns[1:]:assert np.array_equal(tm[col],contrast[col])
    result=dict(passed=True,protected_input_files=len(frozen),contexts=44,prediction_files=pred_count,row_contexts=len(s),bin_rows=660,empty_bin_rows=int(b.n.eq(0).sum()),paired_contrasts=33,independent_nearest_checks=sampled,nearest_label_max_abs_difference=max_nn_error,bin_max_abs_difference=max_bin_error,zero_test_fingerprints=int(s.test_zero_fingerprint.sum()),no_fitting=True,no_bootstrap=True)
    save(R/'VERIFICATION.json',result);print(json.dumps(result))
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--prepared',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    P=a.bundle.resolve();PREP=a.prepared.resolve();R=a.output.resolve();main()

