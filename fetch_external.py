"""Read-only, cached official ChEMBL retrieval. No score-dependent selection."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import hashlib,json,urllib.request,urllib.parse,datetime,time
HERE=Path(__file__).resolve().parent
RAW=HERE/'external'/'raw'; BASE='https://www.ebi.ac.uk/chembl/api/data/'
CANDIDATES={'ache':('CHEMBL220','P22303'),'dpp4':('CHEMBL284','P27487'),'ptgs2':('CHEMBL230','P35354'),'ca2':('CHEMBL205','P00918'),'hdac1':('CHEMBL325','Q13547'),'pde4d':('CHEMBL288','Q08499')}
def get(ep,dest):
    dest.parent.mkdir(parents=True,exist_ok=True); meta=dest.with_suffix('.request.json'); url=BASE+ep
    if dest.exists():
        m=json.loads(meta.read_text(encoding='utf8'));assert m['url']==url and m['sha256']==hashlib.sha256(dest.read_bytes()).hexdigest();return json.loads(dest.read_text(encoding='utf8'))
    started=datetime.datetime.now(datetime.timezone.utc).isoformat()
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Algorithm-provenance-audit/1.0'}),timeout=90) as r:body=r.read();status=r.status
    except Exception as e:
        (dest.with_suffix('.error.json')).write_text(json.dumps(dict(url=url,time=started,error=repr(e),automatic_retry=False),indent=2),encoding='utf8');raise
    obj=json.loads(body); temp=dest.with_suffix('.tmp');temp.write_bytes(body);temp.replace(dest)
    meta.write_text(json.dumps(dict(url=url,time=started,status=status,sha256=hashlib.sha256(body).hexdigest()),indent=2),encoding='utf8')
    return obj
def retrieve(item):
    symbol,(tid,accession)=item; d=RAW/symbol
    target=get('target/'+tid+'.json',d/'target.json')
    assert target['organism']=='Homo sapiens' and target['target_type']=='SINGLE PROTEIN'
    assert [c['accession'] for c in target['target_components']]==[accession]
    fields='activity_id,record_id,molecule_chembl_id,canonical_smiles,assay_chembl_id,document_chembl_id,document_year,pchembl_value,potential_duplicate,data_validity_comment,standard_flag,standard_relation,standard_type,standard_value,standard_units,target_chembl_id,target_organism,assay_type,bao_format,assay_description,activity_comment'
    offset=0; rows=[]; total=None
    while True:
        ep='activity.json?'+urllib.parse.urlencode(dict(target_chembl_id=tid,standard_type='IC50',standard_relation='=',standard_flag=1,limit=1000,offset=offset,order_by='activity_id',only=fields))
        obj=get(ep,d/('activities_%06d.json'%offset)); page=obj['activities']; n=obj['page_meta']['total_count']
        if total is None:total=n
        assert n==total
        rows.extend(page);print(symbol,'activities',len(rows),'/',total,flush=True)
        if not obj['page_meta']['next']:break
        assert page;offset+=len(page)
    assert len(rows)==total and len({r['activity_id'] for r in rows})==total
    aids=sorted({r['assay_chembl_id'] for r in rows if r['assay_chembl_id']})
    got={}
    for i in range(0,len(aids),50):
        batch=aids[i:i+50];ep='assay.json?'+urllib.parse.urlencode(dict(assay_chembl_id__in=','.join(batch),limit=1000,only='assay_chembl_id,target_chembl_id,relationship_type,confidence_score,assay_type,bao_format,description'))
        obj=get(ep,d/('assays_%06d.json'%i));assert obj['page_meta']['next'] is None
        got.update({a['assay_chembl_id']:a for a in obj['assays']})
        assert set(batch).issubset(got)
    completion=dict(target=tid,accession=accession,activity_records=len(rows),assays=len(got),complete=True)
    (d/'completion.json').write_text(json.dumps(completion,indent=2),encoding='utf8')
    print(symbol,'retrieval complete',flush=True)
    return completion
def main():
    s=get('status.json',RAW/'status_before.json');assert s['chembl_db_version']=='ChEMBL_37'
    # Two concurrent read streams; each stream paginates serially without retries.
    with ThreadPoolExecutor(max_workers=2) as ex: results=list(ex.map(retrieve,CANDIDATES.items()))
    e=get('status.json',RAW/'status_after.json');assert e['chembl_db_version']==s['chembl_db_version']
    (RAW/'RETRIEVAL_COMPLETE.json').write_text(json.dumps(dict(version=s['chembl_db_version'],targets=results,protocol_sha256=hashlib.sha256((HERE/'EXTERNAL_PROTOCOL.md').read_bytes()).hexdigest()),indent=2),encoding='utf8')
if __name__=='__main__':main()
