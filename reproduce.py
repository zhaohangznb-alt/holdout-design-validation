"""Offline package-relative preparation -> fresh training -> verification -> analysis."""
from pathlib import Path
import argparse,json,subprocess,sys,time,datetime,platform,hashlib
HERE=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--prepare-only',action='store_true');ap.add_argument('--quick-check',action='store_true');args=ap.parse_args()
    out=args.output.resolve();out.relative_to(HERE);out.mkdir(parents=True,exist_ok=False);(out/'logs').mkdir()
    records=[]
    def step(script,*params):
        command=[sys.executable,'-B',str(HERE/script),*map(str,params)]
        log=out/'logs'/('%02d_%s.log'%(len(records)+1,Path(script).stem));started=time.time()
        print('STEP',script,flush=True)
        with log.open('w',encoding='utf8') as stream:proc=subprocess.run(command,cwd=HERE,stdout=stream,stderr=subprocess.STDOUT)
        records.append(dict(command=command,returncode=proc.returncode,elapsed_seconds=time.time()-started,log=log.relative_to(out).as_posix(),log_sha256=hashlib.sha256(log.read_bytes()).hexdigest()))
        (out/'steps.json').write_text(json.dumps(records,indent=2),encoding='utf8')
        if proc.returncode:raise RuntimeError('Failed: '+script+'; see '+str(log))
    ref=HERE/'inputs/current_reference';origin=HERE/'inputs/current_origin';raw=HERE/'inputs/external_raw'
    current=out/'current_prepared';external=out/'external_prepared'
    step('prepare_current.py','--origin',origin,'--reference',ref,'--output',current)
    step('validate_stage1.py','--input-root',current,'--save',out/'current_input_check.json')
    step('prepare_external.py','--raw',raw,'--model-config',ref/'model_config.json','--output',external)
    step('validate_external.py','--input-root',external,'--save',out/'external_input_check.json')
    if not args.prepare_only:
        extra=['--model','ridge','--limit-cases','1'] if args.quick_check else []
        step('runner.py','--input-root',current,'--output',out/'current_models','--include-paired',*extra)
        step('validate_stage1.py','--input-root',current,'--run',out/'current_models','--save',out/'current_run_check.json')
        step('external_runner.py','--input-root',external,'--output',out/'external_models',*extra)
        step('validate_external.py','--input-root',external,'--run',out/'external_models','--save',out/'external_run_check.json')
        if not args.quick_check:step('analyze_results.py','--current',out/'current_models','--external',out/'external_models','--output',out/'analysis')
    (out/'PIPELINE_COMPLETE.json').write_text(json.dumps(dict(completed=True,scope='preparation_only' if args.prepare_only else 'two_case_smoke_check' if args.quick_check else 'full_fresh_training_and_analysis',network_used=False,platform=platform.platform(),python=sys.version,steps=len(records),completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()),indent=2),encoding='utf8')
    print('PIPELINE COMPLETE',flush=True)
if __name__=='__main__':main()
