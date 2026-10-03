"""Independent gates for frozen inputs and completed fresh-model runs."""
from pathlib import Path
import argparse, itertools, json, math, sys
import numpy as np
import pandas as pd
import runner

def check_inputs(root):
    m,s,frames,hashes=runner.load_inputs(root)
    contexts=sorted({k[:-1] for k in frames})
    integrity=[]
    for k in contexts:
        p={part:frames[k+(part,)] for part in ['train','valid','test']}
        row=dict(zip(runner.KEY,k))
        for col,tag in [('record_uid','id'),('smiles','exact'),('murcko_scaffold','scaffold'),('document_chembl_id','document')]:
            row[tag+'_overlap']=sum(len(set(p[a][col])&set(p[b][col])) for a,b in itertools.combinations(p,2))
        assert row['id_overlap']==row['exact_overlap']==0, row
        if 'scaffold' in k[-1]: assert row['scaffold_overlap']==0, row
        if k[-1].startswith('source'): assert row['document_overlap']==0, row
        else:
            assert p['train'].document_year.max()<p['valid'].document_year.min()
            assert p['valid'].document_year.max()<p['test'].document_year.min()
        integrity.append(row)
    assert len(contexts)==96
    assert len({runner.sig(frames[k+('train',)]) for k in contexts})==26
    main=[k for k in contexts if k[0]=='main']
    assert len(main)==48 and len({runner.sig(frames[k+('train',)]) for k in main})==24
    from rdkit import DataStructs
    bits=np.array([[0,0,0,0],[1,1,0,0],[1,1,0,0],[0,0,1,1]],dtype=np.uint8)
    bv=[DataStructs.CreateFromBitString(''.join(map(str,x))) for x in bits]
    for q in [bits[0],bits[1],np.array([1,0,1,0],dtype=np.uint8)]:
        sims=DataStructs.BulkTanimotoSimilarity(DataStructs.CreateFromBitString(''.join(map(str,q))),bv)
        idx,val=runner.first_nn(q,bits)
        assert idx==int(np.argmax(sims)) and math.isclose(val,max(sims))
    assert runner.first_nn(bits[1],bits)[0]==1 # first tied original row
    v=runner.metrics(np.array([1.,3.]),np.array([1.,2.]),2.)
    assert v['mse_model']==.5 and v['mse_mean']==1 and v['skill']==.5 and v['D']==.5
    assert runner.metrics(np.ones(2),np.ones(2),1.)['skill'] is None
    return m,s,frames,hashes,integrity

def check_run(root,out,m,frames,hashes,protocol='PROTOCOL.md'):
    binding=json.loads((out/'binding.json').read_text())
    assert binding['input_hashes']==hashes
    assert binding['runner_sha256']==runner.sha(runner.__file__)
    assert binding['protocol_sha256']==runner.sha(runner.HERE/protocol)
    entries=[json.loads(p.read_text()) for p in out.glob('*.complete.json')]
    assert entries
    seen=set(); checks=[]
    for r in entries:
        k=tuple(r[x] for x in runner.KEY); ident=k+(r['model'],r['seed'])
        assert ident not in seen; seen.add(ident)
        f=pd.read_csv(out/r['prediction_file']); t=frames[k+('test',)]; tr=frames[k+('train',)]
        assert runner.sha(out/r['prediction_file'])==r['prediction_sha256']
        assert runner.sig(f)==runner.sig(t)==r['test_signature']
        assert f.record_uid.tolist()==t.record_uid.tolist()
        mean=math.fsum(map(float,tr.target))/len(tr)
        assert np.allclose(f.training_mean,mean,rtol=0,atol=1e-12)
        y=t.target.to_numpy(float); pred=f.prediction.to_numpy(float)
        assert np.isfinite(pred).all()
        model_mse=math.fsum(map(float,(y-pred)**2))/len(y)
        base_mse=math.fsum(map(float,(y-mean)**2))/len(y)
        for col,expected in [('squared_model_error',(y-pred)**2),('squared_baseline_error',(y-mean)**2),('absolute_model_error',abs(y-pred)),('absolute_baseline_error',abs(y-mean))]:
            assert np.allclose(f[col],expected,rtol=1e-12,atol=1e-12)
        for name,val in [('mse_model',model_mse),('mse_mean',base_mse),('D',base_mse-model_mse),('skill',1-model_mse/base_mse),('mae_model',math.fsum(map(float,abs(y-pred)))/len(y)),('mae_mean',math.fsum(map(float,abs(y-mean)))/len(y)),('rmse_model',math.sqrt(model_mse)),('rmse_mean',math.sqrt(base_mse))]:
            assert math.isclose(r[name],val,rel_tol=1e-11,abs_tol=1e-11),(ident,name,r[name],val)
        ledger=json.loads((out/(r['job_key']+'.fit.json')).read_text())
        assert ledger['state']=='complete' and ledger['identity']['training_signature']==runner.sig(tr)
        assert runner.sha(out/(r['job_key']+'.joblib'))==ledger['model_sha256']
        assert ledger['actual_fit_executed']==(r['model']!='tanimoto1nn')
        if r['model']=='rf':
            cfg=json.loads((root/'model_config.json').read_text())
            assert ledger['identity']['config']==dict(cfg['rf'],n_jobs=2,random_state=r['seed'])
        elif r['model']=='ridge':
            assert ledger['identity']['config']==dict(alpha=1.0,fit_intercept=True,solver='lsqr',tol=1e-4,max_iter=10000)
        else:
            assert ledger['identity']['config']==dict(tie='first_original_train_row',zero_union_similarity=1)
        checks.append(dict(zip(runner.KEY,k),model=r['model'],seed=r['seed'],test_rows=len(t),metrics_verified=True))
    for k in {tuple(r[x] for x in runner.KEY) for r in entries}:
        assert len({r['test_signature'] for r in entries if tuple(r[x] for x in runner.KEY)==k})==1
    expected={tuple(getattr(r,k) for k in runner.KEY)+(mod,seed) for r in m.drop_duplicates(runner.KEY).itertuples() for mod in ['rf','ridge','tanimoto1nn'] for seed in ([42,53,67] if mod=='rf' else [0])}
    return dict(verified_cases=len(checks),complete_expected_cases=seen==expected,expected_cases=len(expected),missing_cases=len(expected-seen),unique_model_artifacts=len(list(out.glob('*.fit.json'))))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input-root',type=Path,default=runner.DEFAULT);ap.add_argument('--run',type=Path);ap.add_argument('--save',type=Path);args=ap.parse_args()
    m,s,f,h,integrity=check_inputs(args.input_root)
    result=dict(input_entries_verified=len(s),input_files_hashed=len(h),contexts=96,unique_training_inputs=26,split_integrity_passed=True,nn_ties_and_zero_union_passed=True,metric_reference_checks_passed=True)
    if args.run: result.update(check_run(args.input_root,args.run,m,f,h))
    if args.save:
        args.save.resolve().relative_to(runner.HERE)
        args.save.write_text(json.dumps(result,indent=2),encoding='utf8')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
