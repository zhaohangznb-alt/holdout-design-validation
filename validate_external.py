"""Verify selected external identities, frozen splits, and independent run metrics."""
import argparse,itertools,json
from pathlib import Path
import external_runner as er
import validate_stage1 as v
from fetch_external import CANDIDATES

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input-root',type=Path,default=er.DEFAULT);ap.add_argument('--run',type=Path);ap.add_argument('--save',type=Path);args=ap.parse_args()
    plan=json.loads((args.input_root/'decision_record.json').read_text(encoding='utf8'))
    assert not plan['selection_used_model_scores'] and plan['protocol_sha256']==er.sha(er.HERE/'EXTERNAL_PROTOCOL.md')
    m,s,frames,hashes=er.load_inputs(args.input_root)
    keys={k[:-1] for k in frames}
    for k in keys:
        p={part:frames[k+(part,)] for part in ['train','valid','test']}
        for col in ['smiles','record_uid']+(['murcko_scaffold'] if 'scaffold' in k[-1] else [])+(['document_chembl_id'] if k[-1].startswith('source') else []):
            assert all(not(set(p[a][col])&set(p[b][col])) for a,b in itertools.combinations(p,2))
        if k[-1].startswith('temporal'):
            assert p['train'].document_year.max()<p['valid'].document_year.min()
            assert p['valid'].document_year.max()<p['test'].document_year.min()
        if k[1]=='single_protein':assert all(x.operational_eligible.all() for x in p.values())
    result=dict(selected_targets=plan['selected_targets'],contexts=len(keys),input_manifest_entries=len(s),unique_training_inputs=len({er.sig(frames[k+('train',)]) for k in keys}),split_integrity_passed=True,selection_score_blind=True,protocol_hash_verified=True)
    if args.run:
        v.runner=er;result.update(v.check_run(args.input_root,args.run,m,frames,hashes,protocol='EXTERNAL_PROTOCOL.md'))
    if args.save:
        args.save.resolve().relative_to(er.HERE);er.save(args.save,result)
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
