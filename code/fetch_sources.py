"""Retrieve and verify the immutable byte expectations in source_manifest.tsv.

Dynamic literature HTML may no longer reproduce the original snapshot hash.
Such a mismatch is a failure, never silently adopted as equivalent evidence.
"""
from pathlib import Path
import argparse,csv,hashlib,json,os,time,urllib.request,datetime
ROOT=Path(__file__).resolve().parents[1]
from data_paths import resolve_data_root,data_path,checked_data_destination,DEFAULT_DATA_ROOT

def sha256(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def safe_path(root,relative):
    return data_path(relative,root)

def records(manifest):
    with manifest.open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def download(row,path,data_root=None):
    path=checked_data_destination(path,data_root)
    path.parent.mkdir(parents=True,exist_ok=True)
    partial=checked_data_destination(path.with_suffix(path.suffix+'.part'),data_root)
    expected=int(row['bytes'])
    if expected>=16*1024*1024:
        start=partial.stat().st_size if partial.exists() else 0
        if start>expected:raise ValueError('Partial file larger than expected; inspect manually: '+str(partial))
        with partial.open('ab') as f:
            for low in range(start,expected,4*1024*1024):
                high=min(expected-1,low+4*1024*1024-1)
                for attempt in range(5):
                    try:
                        req=urllib.request.Request(row['url'],headers={'Range':f'bytes={low}-{high}'})
                        with urllib.request.urlopen(req,timeout=90) as response:
                            if response.status!=206:raise ValueError('Server did not honor HTTP byte range')
                            payload=response.read()
                        if len(payload)!=high-low+1:raise ValueError('Incomplete byte range')
                        f.write(payload);f.flush();break
                    except Exception:
                        if attempt==4:raise
                        time.sleep(min(8,2*(attempt+1)))
    else:
        with urllib.request.urlopen(row['url'],timeout=90) as response,partial.open('wb') as f:
            while block:=response.read(1024*1024):f.write(block)
    actual=sha256(partial)
    if partial.stat().st_size!=expected or actual!=row['sha256']:
        raise ValueError(f"Downloaded bytes differ from frozen source: {row['local_path']}; retained .part for inspection")
    partial.replace(path)

def verify_sources(data_root,manifest,fetch=False,all_snapshots=False):
    data_root=resolve_data_root(data_root,require_exists=not fetch,for_write=fetch)
    details=[];failures=[]
    for row in records(manifest):
        path=safe_path(data_root,row['local_path'])
        required=row['required_for_analysis'].lower()=='true' or all_snapshots
        if not path.exists() and fetch and required:
            print('Retrieving',row['local_path'],flush=True)
            download(row,path,data_root)
        if not path.exists():
            status='missing_required' if required else 'missing_optional_snapshot'
            if required:failures.append(row['local_path'])
            details.append({'path':row['local_path'],'status':status,'required':required});continue
        actual=sha256(path);valid=path.stat().st_size==int(row['bytes']) and actual==row['sha256']
        details.append({'path':row['local_path'],'status':'verified' if valid else 'checksum_mismatch','required':required,'sha256':actual,'bytes':path.stat().st_size})
        if not valid:failures.append(row['local_path'])
    result={'verified_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'data_root':str(data_root),'manifest_sha256':sha256(manifest),'status':'PASS' if not failures else 'FAIL','verified_files':sum(r['status']=='verified' for r in details),'missing_optional_snapshots':sum(r['status']=='missing_optional_snapshot' for r in details),'failed_paths':failures,'files':details}
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-root',type=Path,default=None,help='Absolute external data root; overrides STRIKE007_CACHE and external user policy')
    p.add_argument('--fetch',action='store_true',help='Retrieve missing required originals; refuse mismatched bytes')
    p.add_argument('--all-source-snapshots',action='store_true',help='Require/retrieve even optional documentary snapshots')
    p.add_argument('--receipt',type=Path,help='Optional new JSON file outside the repository; otherwise print the summary only')
    a=p.parse_args()
    if a.receipt is not None:
        receipt_root=resolve_data_root(a.receipt.parent,for_write=True)
        a.receipt=checked_data_destination(a.receipt,receipt_root)
        if a.receipt.exists():raise FileExistsError('Receipt already exists: '+str(a.receipt))
    r=verify_sources(resolve_data_root(a.data_root,require_exists=not a.fetch,for_write=a.fetch),ROOT/'source_manifest.tsv',a.fetch,a.all_source_snapshots)
    if a.receipt is not None:
        a.receipt.parent.mkdir(parents=True,exist_ok=True)
        with a.receipt.open('x',encoding='utf-8') as stream:json.dump(r,stream,indent=2)
    print(json.dumps({k:v for k,v in r.items() if k!='files'},indent=2))
    if r['status']!='PASS':raise SystemExit(2)
