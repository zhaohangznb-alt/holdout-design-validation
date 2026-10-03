"""Verify public bundle bytes and recompute metrics from published predictions."""
from pathlib import Path
import hashlib,json,math
import pandas as pd
import numpy as np
HERE=Path(__file__).resolve().parent
def main():
    manifest=json.loads((HERE/'FILE_MANIFEST.json').read_text(encoding='utf8'))
    for name,h in manifest.items():assert hashlib.sha256((HERE/name).read_bytes()).hexdigest()==h,name
    df=pd.read_csv(HERE/'results/analysis/run_metrics.csv');cases=set()
    for r in df.itertuples():
        key=(r.panel,r.view,r.scenario,r.target,r.route,r.model,r.seed);assert key not in cases;cases.add(key)
        p=HERE/'results/predictions'/r.panel/r.prediction_file;assert hashlib.sha256(p.read_bytes()).hexdigest()==r.prediction_sha256
        f=pd.read_csv(p);y=f.target.to_numpy(float);pred=f.prediction.to_numpy(float);b=f.training_mean.to_numpy(float)
        mse=math.fsum(map(float,(y-pred)**2))/len(f);base=math.fsum(map(float,(y-b)**2))/len(f)
        assert math.isclose(mse,r.mse_model,rel_tol=1e-11,abs_tol=1e-11)
        assert math.isclose(base,r.mse_mean,rel_tol=1e-11,abs_tol=1e-11)
        assert math.isclose(1-mse/base,r.skill,rel_tol=1e-11,abs_tol=1e-11)
        for name,val in [('mae_model',math.fsum(map(float,abs(y-pred)))/len(y)),('mae_mean',math.fsum(map(float,abs(y-b)))/len(y)),('rmse_model',math.sqrt(mse)),('rmse_mean',math.sqrt(base))]:
            assert math.isclose(getattr(r,name),val,rel_tol=1e-11,abs_tol=1e-11)
    print(json.dumps(dict(files_verified=len(manifest),prediction_cases_verified=len(cases),expected_prediction_cases=680,complete=len(cases)==680)))
    assert len(cases)==680
if __name__=='__main__':main()
