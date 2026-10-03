"""Start fresh matching fits in a NEW directory, copying only required prepared references."""
from pathlib import Path
import argparse,shutil,json,sys,subprocess
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    src=Path(__file__).resolve().parent;dst=args.output.resolve()
    assert not dst.exists(),'Use a new, empty destination; no completed work will be overwritten'
    dst.mkdir(parents=True)
    for d in ['models','analysis','checkpoints','ledgers','predictions']: (dst/d).mkdir()
    shutil.copytree(src/'inputs',dst/'inputs')
    for name in ['ANALYSIS_PLAN.json','MODEL_CONFIG.json','SOURCE_HASHES.json','BINDING.json','INPUT_HASHES.json','PREPARED.json','POOL_COUNTS.csv','manifest.csv','experiment.py']:
        shutil.copy2(src/name,dst/name)
    copied=0
    for p in (src/'ledgers').glob('*.complete.json'):
        e=json.loads(p.read_text(encoding='utf8'))
        if e['reused_frozen_prediction']:
            shutil.copy2(p,dst/'ledgers'/p.name);shutil.copy2(src/e['prediction_file'],dst/e['prediction_file']);copied+=1
    assert copied==110
    subprocess.run([sys.executable,str(src/'experiment.py'),'--root',str(dst),'--run'],check=True)
    subprocess.run([sys.executable,str(src/'experiment.py'),'--root',str(dst),'--summarize'],check=True)
    print('Fresh execution complete:',dst)
if __name__=='__main__': main()
