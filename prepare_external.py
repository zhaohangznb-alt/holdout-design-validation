"""Offline cleaning and blind target screening from frozen API JSON."""
from pathlib import Path
import argparse,json,hashlib,itertools,math
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from fetch_external import CANDIDATES
from runner import HERE,KEY,sha,sig,save

def read_raw(raw):
    assert json.loads((raw/'RETRIEVAL_COMPLETE.json').read_text(encoding='utf8'))['version']=='ChEMBL_37'
    hashes={}
    for meta in sorted(raw.rglob('*.request.json')):
        p=meta.with_name(meta.name.replace('.request.json','.json'))
        assert sha(p)==json.loads(meta.read_text(encoding='utf8'))['sha256']
        hashes[p.relative_to(raw).as_posix()]=sha(p)
    return hashes
def clean(raw,symbol,tid):
    records=[];assays={}
    for p in sorted((raw/symbol).glob('activities_*.json')):
        if p.name.endswith(('.request.json','.error.json')):continue
        records.extend(json.loads(p.read_text(encoding='utf8'))['activities'])
    for p in sorted((raw/symbol).glob('assays_*.json')):
        if p.name.endswith(('.request.json','.error.json')):continue
        assays.update({a['assay_chembl_id']:a for a in json.loads(p.read_text(encoding='utf8'))['assays']})
    f=pd.DataFrame(records);assert f.activity_id.is_unique
    flow=[]
    def keep(name,mask):
        nonlocal f
        n=len(f);f=f.loc[mask.loc[f.index]].copy();flow.append(dict(target=symbol,step=name,before=n,after=len(f),removed=n-len(f)))
    keep('matching_human_target',f.target_chembl_id.eq(tid)&f.target_organism.eq('Homo sapiens'))
    keep('exact_standard_IC50',f.standard_type.eq('IC50')&f.standard_relation.eq('=')&pd.to_numeric(f.standard_flag).eq(1))
    keep('no_database_quality_flags',pd.to_numeric(f.potential_duplicate,errors='coerce').fillna(0).eq(0)&f.data_validity_comment.isna())
    val=pd.to_numeric(f.standard_value,errors='coerce');label=pd.to_numeric(f.pchembl_value,errors='coerce')
    keep('finite_positive_nM_and_pchembl',f.standard_units.eq('nM')&np.isfinite(val)&val.gt(0)&np.isfinite(label))
    f['target']=pd.to_numeric(f.pchembl_value);f['document_year']=pd.to_numeric(f.document_year,errors='coerce')
    required=['activity_id','canonical_smiles','assay_chembl_id','document_chembl_id','document_year']
    keep('complete_identifiers_year_structure',f[required].notna().all(axis=1)&np.isfinite(f.document_year)&f.document_year.eq(f.document_year.round()))
    cache={}
    def chem(sm):
        if sm not in cache:
            mol=Chem.MolFromSmiles(sm)
            if mol is None:cache[sm]=None
            else:
                c=Chem.MolToSmiles(mol,canonical=True);sc=MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
                cache[sm]=(c,sc or 'acyclic:'+c)
        return cache[sm]
    z=f.canonical_smiles.map(chem);keep('rdkit_parse',z.notna())
    f['smiles']=[chem(sm)[0] for sm in f.canonical_smiles];f['murcko_scaffold']=[chem(sm)[1] for sm in f.canonical_smiles]
    f['document_year']=f.document_year.astype(int);f['source_uid']='chembl37_activity_'+f.activity_id.astype(str)
    f['operational_eligible']=f.assay_chembl_id.map(lambda aid:assays[aid].get('relationship_type')=='D' and assays[aid].get('target_chembl_id')==tid and assays[aid].get('assay_type')=='B' and assays[aid].get('bao_format')=='BAO_0000357')
    cols=['smiles','assay_chembl_id','document_chembl_id','document_year','murcko_scaffold']
    g=f.groupby(cols,as_index=False,sort=True).agg(target=('target','median'),n_source_records=('source_uid','size'),source_record_uids=('source_uid',lambda v:'|'.join(sorted(set(v)))),operational_eligible=('operational_eligible','all'))
    g['record_uid']=['ext37_'+hashlib.sha256((tid+'|'+json.dumps(list(x),ensure_ascii=False)).encode()).hexdigest()[:24] for x in g[cols].itertuples(index=False,name=None)]
    assert g.record_uid.is_unique
    flow.append(dict(target=symbol,step='aggregate_context_median',before=len(f),after=len(g),removed=len(f)-len(g)))
    summary=dict(target=symbol,target_id=tid,retrieved_activities=len(records),cleaned_activities=len(f),contexts=len(g),molecules=g.smiles.nunique(),documents=g.document_chembl_id.nunique(),first_year=int(g.document_year.min()),last_year=int(g.document_year.max()))
    return g,flow,summary
def split(g):
    docs=np.array(sorted(g.document_chembl_id.unique()));docs=docs[np.random.default_rng(42).permutation(len(docs))]
    nt=min(max(1,round(.7*len(docs))),len(docs)-2);nv=min(max(1,round(.15*len(docs))),len(docs)-nt-1)
    source={p:g[g.document_chembl_id.isin(d)].copy() for p,d in zip(['train','valid','test'],[docs[:nt],docs[nt:nt+nv],docs[nt+nv:]])}
    counts=g.groupby('document_year').size().sort_index(ascending=False);testyears=[];total=len(g);acc=0
    for y,n in counts.items():
        testyears.append(y);acc+=n
        if acc>=.2*total:break
    validyears=[];acc=0
    for y,n in counts[counts.index<min(testyears)].items():
        validyears.append(y);acc+=n
        if acc>=.15*total:break
    assert validyears
    trainend=int(min(validyears)-1);validend=int(max(validyears))
    temporal={'train':g[g.document_year.le(trainend)].copy(),'valid':g[g.document_year.gt(trainend)&g.document_year.le(validend)].copy(),'test':g[g.document_year.gt(validend)].copy()}
    def filter_parts(parts,scaffold):
        out={};sm=set();sc=set()
        for p in ['train','valid','test']:
            x=parts[p];x=x[~x.smiles.isin(sm)]
            if scaffold:x=x[~x.murcko_scaffold.isin(sc)]
            out[p]=x.reset_index(drop=True);sm.update(x.smiles);sc.update(x.murcko_scaffold)
        for col in ['smiles']+(['murcko_scaffold'] if scaffold else []):
            assert all(not(set(out[a][col])&set(out[b][col])) for a,b in itertools.combinations(out,2))
        return out
    routes={prefix+'_'+suffix:filter_parts(parts,suffix=='scaffold') for prefix,parts in [('source',source),('temporal',temporal)] for suffix in ['exact','scaffold']}
    return routes,dict(training_through=trainend,validation_through=validend)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--raw',type=Path,default=HERE/'external/raw');ap.add_argument('--output',type=Path,default=HERE/'external/prepared_v2');ap.add_argument('--model-config',type=Path,default=HERE/'../claude_manuscript_polish_20261002_v1/final_delivery_20261003_v1/03_reproducibility/policy631_execution/model_config.json');args=ap.parse_args()
    out=args.output.resolve();out.relative_to(HERE);out.mkdir(parents=True,exist_ok=False)
    hashes=read_raw(args.raw); assert json.loads((args.raw/'RETRIEVAL_COMPLETE.json').read_text(encoding='utf8'))['protocol_sha256']==sha(HERE/'EXTERNAL_PROTOCOL.md')
    selected=[];summaries=[];allflow=[];splits=[];prediction=[];cutoffs={};unavailable=[]
    for symbol,(tid,accession) in CANDIDATES.items():
        g,flow,summary=clean(args.raw,symbol,tid);allflow.extend(flow);routes,cut=split(g);cutoffs[symbol]=cut
        reasons=[]
        for check,ok in [('cleaned_records>=1500',summary['cleaned_activities']>=1500),('molecules>=1000',summary['molecules']>=1000),('documents>=20',summary['documents']>=20),('year_span>=10',summary['last_year']-summary['first_year']>=10)]:
            if not ok:reasons.append(check)
        for r in ['source_exact','temporal_exact']:
            summary[r+'_test_contexts']=len(routes[r]['test'])
            if len(routes[r]['test'])<100:reasons.append(r+'_test>=100')
        summary.update(selected=not reasons,failed_criteria=';'.join(reasons));summaries.append(summary)
        (out/'screened_contexts').mkdir(exist_ok=True);g.to_csv(out/'screened_contexts'/(symbol+'.csv'),index=False)
        if reasons:continue
        selected.append(symbol)
        for scenario in ['revised','single_protein']:
            for route,parts in routes.items():
                p={n:(x if scenario=='revised' else x[x.operational_eligible]).copy() for n,x in parts.items()}
                if not all(len(x)>0 for x in p.values()):
                    unavailable.append(dict(target=symbol,scenario=scenario,route=route,reason='empty partition'));continue
                for part,x in p.items():
                    dest=out/'data/main'/(scenario+'_'+symbol+'_'+route+'_'+part+'.csv');dest.parent.mkdir(parents=True,exist_ok=True);x.to_csv(dest,index=False)
                    # Read back serialized floats: the frozen CSV is the model input.
                    p[part]=pd.read_csv(dest)
                    splits.append(dict(view='main',scenario=scenario,target=symbol,route=route,partition=part,corrected_csv=dest.name,corrected_sha256=sha(dest),corrected_rows=len(x),corrected_ordered_signature=sig(p[part]),below_100_test=part=='test' and len(x)<100))
                for seed in [42,53,67]:prediction.append(dict(view='main',scenario=scenario,target=symbol,route=route,seed=seed,training_signature=sig(p['train']),test_rows=len(p['test'])))
    assert selected and prediction
    pd.DataFrame(summaries).to_csv(out/'candidate_screen.csv',index=False);pd.DataFrame(allflow).to_csv(out/'cleaning_flow.csv',index=False)
    pd.DataFrame(splits).to_csv(out/'split_manifest.csv',index=False);pd.DataFrame(prediction).to_csv(out/'prediction_manifest.csv',index=False)
    cfg=json.loads(args.model_config.read_text(encoding='utf8'))
    cfg.pop('verified_script_sources',None);save(out/'model_config.json',cfg)
    save(out/'decision_record.json',dict(selected_targets=selected,candidate_screen='candidate_screen.csv',unavailable=unavailable,cutoffs=cutoffs,selection_used_model_scores=False,protocol_sha256=sha(HERE/'EXTERNAL_PROTOCOL.md'),raw_response_hashes=hashes,prepared_code_sha256=sha(__file__),contexts=len(prediction)//3))
    print(pd.DataFrame(summaries).to_string(index=False))
if __name__=='__main__':main()
