"""Portable input/prediction/results audit; no original paths, fitting or model files."""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main(root):
    manifest=json.loads((root/'SNAPSHOT_MANIFEST.json').read_text())
    for p,h in manifest.items(): assert sha(root/p)==h,p
    m=pd.read_csv(root/'manifest.csv',float_precision='round_trip'); r=pd.read_csv(root/'analysis/run_metrics.csv',float_precision='round_trip'); cache={}; checked={}
    def csv(p):
        if p not in cache: cache[p]=pd.read_csv(root/p,float_precision='round_trip')
        return cache[p]
    for x in m.drop_duplicates('prediction_file').itertuples():
        f=csv(x.prediction_file); te=csv(x.test_file); tr=csv(x.train_file)
        assert list(f.record_uid)==list(te.record_uid) and np.array_equal(f.target.to_numpy(),te.target.to_numpy())
        assert len(tr)==x.train_n and len(te)==x.test_n
        for k in ['record_uid','smiles','document_chembl_id']: assert not set(te[k])&set(tr[k])
        loss=np.mean((te.target-f.prediction)**2); base=np.mean((te.target-tr.target.mean())**2); var=np.var(te.target)
        checked[x.prediction_file]=(loss,base,var)
    maxerr=0.
    for x in r.itertuples():
        a,b,v=checked[x.prediction_file]
        for lhs,rhs in [(x.mse,a),(x.baseline_mse,b),(x.test_variance,v),(x.training_mean_skill,1-a/b),(x.normalized_mse,a/v)]:
            err=abs(lhs-rhs); maxerr=max(maxerr,err); assert err<=1e-12
    out=dict(status='PASS',verified_snapshot_files=len(manifest),unique_predictions=len(checked),logical_prediction_references=len(r),metrics_max_abs_error=maxerr,model_fits_executed=False,original_absolute_paths_accessed=False)
    print(json.dumps(out,indent=2))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parent); args=ap.parse_args(); main(args.root)
