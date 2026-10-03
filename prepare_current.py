"""Rebuild all current retained cohorts from frozen original partitions + policy.

Uses no predictions or fitted models. No training labels or structures replaced.
"""
from pathlib import Path
import argparse,csv,json,itertools,shutil
import pandas as pd
from runner import HERE,KEY,sha,sig,save,DEFAULT
TARGETS={'alk':'CHEMBL4247','braf':'CHEMBL5145','egfr':'CHEMBL203','erbb2':'CHEMBL1824','kdr':'CHEMBL279','met':'CHEMBL3717'}
BAD={'egfr':{'CHEMBL5736995','CHEMBL5733991','CHEMBL5734361','CHEMBL5737313','CHEMBL5737314'},'erbb2':{'CHEMBL3883068'},'braf':{'CHEMBL5730741'}}
BIND={'CHEMBL5737276','CHEMBL5737277','CHEMBL5737278','CHEMBL5737279','CHEMBL5737282','CHEMBL5737283'}
OVERRIDE={'CHEMBL3888984','CHEMBL3888985'}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--origin',type=Path,default=DEFAULT.parent);ap.add_argument('--reference',type=Path,default=DEFAULT);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
    origin=args.origin;ref=args.reference;out=args.output.resolve();out.relative_to(HERE);out.mkdir(parents=True,exist_ok=False)
    pol=pd.read_csv(ref/'policy631.csv');common=set(pol[pol.scope.eq('common_source_eligibility')].context_uid);strict=set(pol[pol.scope.eq('main_single_protein_only')].context_uid)
    assert len(common)==600 and len(strict)==31
    q=set(json.loads((origin/'policy/policy_A.json').read_text(encoding='utf8'))['quarantine_ids']);assert len(q)==35
    api={};input_hashes={}
    for p in sorted((origin/'metadata').glob('assays_*.json')):
        api.update({a['assay_chembl_id']:a for a in json.loads(p.read_text(encoding='utf8'))['assays']});input_hashes['origin/'+p.relative_to(origin).as_posix()]=sha(p)
    allowed={}
    for t,tid in TARGETS.items():
        p=origin/'metadata'/('archived_assays_'+t+'.csv');input_hashes['origin/'+p.relative_to(origin).as_posix()]=sha(p)
        old=pd.read_csv(p).drop_duplicates().set_index('assay_chembl_id');assert old.index.is_unique
        allowed[t]={aid for aid in old.index if aid in api and api[aid].get('relationship_type')=='D' and api[aid].get('target_chembl_id')==tid and old.loc[aid,'assay_type']=='B' and (old.loc[aid,'bao_format']=='BAO_0000357' or aid in OVERRIDE) and aid not in BIND}
    m=pd.read_csv(ref/'split_manifest.csv');audit=[];new=[]
    for r in m.itertuples():
        src=origin/'frozen_inputs'/r.target/r.route/(r.partition+'.csv');input_hashes['origin/'+src.relative_to(origin).as_posix()]=sha(src)
        f=pd.read_csv(src);assert set(f.assay_chembl_id).issubset(api)
        keep=~f.record_uid.isin(common)&~f.assay_chembl_id.isin(BAD.get(r.target,set()))
        if r.scenario=='single_protein':
            keep &= f.assay_chembl_id.isin(allowed[r.target])
            if r.view=='main':keep &= ~f.record_uid.isin(strict)
        elif r.target=='braf' and r.partition=='test':keep &= ~f.record_uid.isin(q)
        ids=set(f.loc[keep,'record_uid']);dest=out/'data'/r.view/(r.scenario+'_'+r.target+'_'+r.route+'_'+r.partition+'.csv');dest.parent.mkdir(parents=True,exist_ok=True)
        with src.open(encoding='utf8',newline='') as stream:lines=list(csv.reader(stream))
        uid=lines[0].index('record_uid')
        with dest.open('w',encoding='utf8',newline='') as stream:
            w=csv.writer(stream,lineterminator='\n');w.writerow(lines[0]);w.writerows(line for line in lines[1:] if line[uid] in ids)
        g=pd.read_csv(dest);assert len(g)==r.corrected_rows and sig(g)==r.corrected_ordered_signature
        z=r._asdict();z.pop('Index');z['corrected_csv']=dest.name;z['corrected_sha256']=sha(dest);new.append(z)
        audit.append(dict(zip(KEY,[r.view,r.scenario,r.target,r.route]),partition=r.partition,original_rows=len(f),retained_rows=len(g),ordered_signature_match=True,literal_csv_hash_match=sha(dest)==r.corrected_sha256))
    pd.DataFrame(new).to_csv(out/'split_manifest.csv',index=False)
    for name in ['prediction_manifest.csv','model_config.json','decision_record.json','policy631.csv']:
        shutil.copyfile(ref/name,out/name);input_hashes['reference/'+name]=sha(ref/name)
    input_hashes['reference/split_manifest.csv']=sha(ref/'split_manifest.csv');input_hashes['origin/policy/policy_A.json']=sha(origin/'policy/policy_A.json')
    pd.DataFrame(audit).to_csv(out/'preparation_audit.csv',index=False)
    save(out/'PREPARATION_COMPLETE.json',dict(tables=len(audit),all_ordered_signatures_match=True,literal_csv_hash_matches=sum(x['literal_csv_hash_match'] for x in audit),data_values_replaced=0,structures_replaced=0,partitions_reassigned=0,predictions_used=False,code_sha256=sha(__file__),input_hashes=input_hashes))
    print(json.dumps(dict(tables=len(audit),ordered_signature_matches=len(audit),literal_csv_hash_matches=sum(x['literal_csv_hash_match'] for x in audit))))
if __name__=='__main__':main()
