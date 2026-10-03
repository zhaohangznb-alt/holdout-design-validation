"""Offline eight-case replay of the sealed bundle on Linux; no bootstrap."""
from pathlib import Path
import argparse, datetime, hashlib, json, os, platform, subprocess, sys, time
KEY=['view','scenario','target','route']
def read(p): return json.loads(p.read_text(encoding='utf8'))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x): p.write_text(json.dumps(x,indent=2),encoding='utf8')
def preflight(bundle):
    import numpy as np,pandas as pd,rdkit,sklearn,scipy,joblib
    cfg=read(bundle/'inputs/current_reference/model_config.json')
    versions=dict(numpy=np.__version__,pandas=pd.__version__,rdkit=rdkit.__version__,sklearn=sklearn.__version__)
    assert all(cfg['environment'][k]==v for k,v in versions.items()),versions
    assert sys.version_info[:3]==(3,9,10),sys.version
    assert scipy.__version__=='1.13.1' and joblib.__version__=='1.5.3'
    assert cfg['seeds']==[42,53,67] and cfg['fingerprint']==dict(radius=2,fpSize=2048,includeChirality=False)
    for panel,target in [('kinase','alk'),('nonkinase','ache')]:
        for model,seeds in [('rf',[42,53,67]),('tanimoto1nn',[0])]:
            for seed in seeds:
                assert (bundle/'results/predictions'/panel/f'main_revised_{target}_source_exact_{model}_{seed}.csv').is_file()
    return dict(versions=versions,python=sys.version,platform=platform.platform(),scipy=scipy.__version__,joblib=joblib.__version__)
def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--bundle',type=Path,required=True);ap.add_argument('--output',type=Path);ap.add_argument('--preflight-only',action='store_true');args=ap.parse_args()
    bundle=args.bundle.resolve();env=preflight(bundle)
    if args.preflight_only: print(json.dumps(dict(preflight='PASS',model_execution=False,**env)));return
    assert platform.system()=='Linux','This entry records Linux execution; use the Windows entry for Windows'
    if args.output is None:ap.error('--output required')
    out=args.output.resolve();out.relative_to(bundle);assert out!=bundle
    out.mkdir(parents=True,exist_ok=False);(out/'logs').mkdir();commands=[]
    # Preserve the archive. Explicit CRLF matches the Windows ordered-signature
    # serialization; it does not change source CSV bytes, row order or labels.
    code_changes=[]
    for source in sorted(bundle.glob('*.py')):
        data=source.read_text(encoding='utf8')
        if source.name in ['runner.py','external_runner.py']:
            old="f[['record_uid','smiles','target']].to_csv(index=False).encode()"
            new="f[['record_uid','smiles','target']].to_csv(index=False,lineterminator='\\r\\n').encode()"
            assert data.count(old)==1
            data=data.replace(old,new)
        if source.name=='prepare_external.py':
            old='.to_csv(dest,index=False)'
            assert data.count(old)==1
            data=data.replace(old,".to_csv(dest,index=False,lineterminator='\\r\\n')")
        dest=out/source.name;dest.write_text(data,encoding='utf8')
        code_changes.append(dict(file=source.name,archived_sha256=sha(source),runtime_sha256=sha(dest),ordered_signature_crlf=source.name in ['runner.py','external_runner.py']))
    for source in bundle.glob('*PROTOCOL.md'):
        (out/source.name).write_bytes(source.read_bytes())
    os.symlink(bundle/'inputs',out/'inputs',target_is_directory=True)
    save(out/'COMPATIBILITY_ADAPTER.json',dict(changes=code_changes,scientific_change=False,detail='Explicit CRLF only for ordered-input hash serialization; all 408 literal input hashes remain mandatory'))
    plan=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),cases=[dict(panel=p,target=t,view='main',scenario='revised',route='source_exact',model=m,seed=s) for p,t in [('kinase','alk'),('nonkinase','ache')] for m,seeds in [('rf',[42,53,67]),('tanimoto1nn',[0])] for s in seeds],prediction_tolerance=1e-12,environment=env,network_used_for_preparation_or_models=False)
    save(out/'CASE_PLAN.json',plan)
    def step(script,*params):
        cmd=[sys.executable,'-B',str(out/script),*map(str,params)];log=out/'logs'/('%02d_%s.log'%(len(commands)+1,Path(script).stem));entry=dict(command=cmd,state='running',started=time.time(),log=log.relative_to(out).as_posix());commands.append(entry);save(out/'commands.json',commands)
        with log.open('w',encoding='utf8') as stream:p=subprocess.run(cmd,cwd=bundle,stdout=stream,stderr=subprocess.STDOUT,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
        entry.update(state='complete' if p.returncode==0 else 'failed',finished=time.time(),returncode=p.returncode,log_sha256=sha(log));save(out/'commands.json',commands)
        if p.returncode:
            print(log.read_text(encoding='utf8'),flush=True)
            for nested in sorted((out/'fresh/logs').glob('*.log')):
                print(nested.name,nested.read_text(encoding='utf8'),flush=True)
            raise RuntimeError(str(log))
    step('reproduce.py','--output',out/'fresh','--prepare-only')
    import numpy as np,pandas as pd
    input_checks=[]
    for panel,sub,ref in [('kinase','current_prepared',bundle/'inputs/current_reference'),('nonkinase','external_prepared',Path(__file__).resolve().parent)]:
        manifest=pd.read_csv(ref/('split_manifest.csv' if panel=='kinase' else 'external_split_manifest.csv'))
        for r in manifest.itertuples():
            filename=str(r.corrected_csv).replace('\\','/').rsplit('/',1)[-1];f=out/'fresh'/sub/'data'/r.view/filename
            assert sha(f)==r.corrected_sha256,(str(f),sha(f),r.corrected_sha256)
            input_checks.append(dict(panel=panel,path=f.relative_to(out).as_posix(),sha256=sha(f),archived_sha256=r.corrected_sha256))
    assert len(input_checks)==408;save(out/'input_checks.json',input_checks)
    predictions=[];fits=indices=0
    for panel,target,script,sub,nrows in [('kinase','alk','runner.py','current_prepared',239),('nonkinase','ache','external_runner.py','external_prepared',921)]:
        manifest=pd.read_csv(out/'fresh'/sub/'prediction_manifest.csv');first=manifest[manifest.view.eq('main')].drop_duplicates(KEY).sort_values(KEY).iloc[0]
        assert tuple(first[k] for k in KEY)==('main','revised',target,'source_exact')
        for model,limit in [('rf',3),('tanimoto1nn',1)]:
            dest=out/(panel+'_'+model);extra=['--include-paired'] if panel=='kinase' else []
            step(script,'--input-root',out/'fresh'/sub,'--output',dest,'--model',model,'--limit-cases',limit,*extra)
            dones=list(dest.glob('*.complete.json'));assert len(dones)==limit and not (dest/'RUNNING.lock').exists()
            assert len(list(dest.glob('command_*.json')))==1
            for done in dones:
                d=read(done);assert tuple(d[k] for k in KEY)==('main','revised',target,'source_exact');f=dest/d['prediction_file'];ref=bundle/'results/predictions'/panel/f.name
                a=pd.read_csv(f,float_precision='round_trip');b=pd.read_csv(ref,float_precision='round_trip')
                assert list(a.columns)==list(b.columns) and len(a)==len(b)==nrows
                numeric=list(a.select_dtypes(include='number').columns);other=[k for k in a.columns if k not in numeric]
                assert a[other].equals(b[other])
                errors={k:float(np.max(np.abs(a[k]-b[k]))) for k in numeric};assert all(v<=1e-12 for v in errors.values()),errors
                if model=='tanimoto1nn':assert errors['prediction']==0
                ledger=read(dest/(d['job_key']+'.fit.json'));assert ledger['state']=='complete';assert ledger['actual_fit_executed']==(model=='rf') and ledger['index_built']==(model=='tanimoto1nn')
                assert ledger['started']>=datetime.datetime.fromisoformat(plan['created_utc']).timestamp()
                if model=='rf':fits+=1
                else:indices+=1
                predictions.append(dict(panel=panel,case=f.name,rows=nrows,prediction_max_abs_difference=errors['prediction'],all_numeric_errors=errors,byte_identical=sha(f)==sha(ref),fresh_sha256=sha(f),archived_sha256=sha(ref),ledger=done.parent.joinpath(d['job_key']+'.fit.json').relative_to(out).as_posix()))
    assert fits==6 and indices==2 and len(predictions)==8
    save(out/'LINUX_REPLAY_COMPLETE.json',dict(passed=True,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),input_tables_byte_identical=408,rf_fresh_fits=fits,nn_indices=indices,prediction_cases=predictions,prediction_tolerance=1e-12,environment=env,bootstrap_executed=False,full_680_retraining=False))
    print('LINUX_COMPLETION_JSON '+json.dumps(read(out/'LINUX_REPLAY_COMPLETE.json'),separators=(',',':')),flush=True)
    print('PASS: Linux offline 408 input tables, 6 RF fits, 2 NN indices',flush=True)
if __name__=='__main__':main()
