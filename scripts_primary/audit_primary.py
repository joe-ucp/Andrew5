"""Audit submitted GEO sparse counts without estimating treatment-score effects.

Run after acquire_primary.py. One sample at a time bounds memory. Reusable load_sample
returns raw gene-by-cell CSR counts and preserves exact feature/barcode indexing.
"""
from pathlib import Path
import gzip,csv,json,sys,time,hashlib,os
import numpy as np
import pandas as pd
from scipy import io,sparse
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from data_paths import data_path,checked_data_destination,configure_data_root_cli
MARKERS=['INS','GCG','SST','PPY','CHGA','CHGB','ARX','MAFB','TTR','IRX2','PCSK1','PCSK2','IAPP','PDX1','MAFA','NKX6-1','SMOC1','MKI67','TOP2A','KRT19','PRSS1','PRSS2','CPA1','AMY2A','COL1A1','PTPRC']
def load_sample(row,cache=True):
    feat=pd.read_csv(data_path(row['features_path'],must_exist=True),sep='\t',header=None,names=['ensembl','symbol','feature_type'])
    bc=pd.read_csv(data_path(row['barcodes_path'],must_exist=True),sep='\t',header=None,names=['barcode'])
    npz=data_path('data_primary/sparse_cache/'+row['accession']+'.npz',for_write=cache)
    if npz.exists():mat=sparse.load_npz(npz)
    else:
        with gzip.open(data_path(row['matrix_path'],must_exist=True),'rb') as f:mat=io.mmread(f).tocsr()
        mat.sum_duplicates();mat.sort_indices()
        if cache:
            npz.parent.mkdir(exist_ok=True)
            partial=checked_data_destination(npz.with_suffix('.partial.npz'))
            sparse.save_npz(partial,mat,compressed=True)
            partial.replace(npz)
    if mat.shape!=(len(feat),len(bc)):raise ValueError('Feature/barcode dimension mismatch')
    return mat,feat,bc
def audit(row):
    QC=data_path('data_primary/audit',for_write=True)
    QC.mkdir(parents=True,exist_ok=True)
    t=time.time();x,feat,bc=load_sample(row)
    lib=np.asarray(x.sum(axis=0)).ravel()
    detected=np.asarray((x>0).sum(axis=0)).ravel()
    mt=feat.symbol.str.startswith('MT-').to_numpy()
    mt_count=np.asarray(x[mt].sum(axis=0)).ravel()
    mt_fraction=mt_count/np.maximum(lib,1)
    keep=(lib>=1000)&(detected>=400)&(mt_fraction<.2)
    metrics=pd.DataFrame({'barcode':bc.barcode,'UMI':lib,'detected_genes':detected,'mt_fraction':mt_fraction,'published_QC_pass':keep})
    marker_rows=[]
    for g in MARKERS:
        ix=np.flatnonzero(feat.symbol.to_numpy()==g)
        v=np.asarray(x[ix].sum(axis=0)).ravel() if len(ix) else np.zeros(len(bc))
        metrics[g]=v
        marker_rows.append({'accession':row['accession'],'gene':g,'feature_matches':len(ix),'total_UMI':int(v.sum()),'cells_detected':int((v>0).sum()),'QC_cells_detected':int(((v>0)&keep).sum())})
    metrics.to_csv(data_path('data_primary/audit/'+row['accession']+'_cell_metrics.tsv.gz',for_write=True),sep='\t',index=False)
    with gzip.open(data_path(row['matrix_path'],must_exist=True),'rt') as f:head=[next(f).rstrip() for _ in range(3)]
    record={'accession':row['accession'],'donor':row['donor'],'treatment_code':row['treatment_code'],'features':x.shape[0],'barcodes':x.shape[1],'nnz':x.nnz,'sum_counts':int(x.sum()),'minimum_count':int(x.data.min()),'maximum_count':int(x.data.max()),'negative_entries':int((x.data<0).sum()),'noninteger_entries':int((x.data!=np.floor(x.data)).sum()),'unique_barcodes':bc.barcode.nunique(),'unique_ensembl':feat.ensembl.nunique(),'unique_symbols':feat.symbol.nunique(),'zero_library_cells':int((lib==0).sum()),'published_QC_pass_cells':int(keep.sum()),'below_1000_UMI_cells':int((lib<1000).sum()),'below_400_gene_cells':int((detected<400).sum()),'mt_ge_20pct_cells':int((mt_fraction>=.2).sum()),'library_q10':float(np.quantile(lib,.1)),'library_median':float(np.median(lib)),'library_q90':float(np.quantile(lib,.9)),'QC_library_median':float(np.median(lib[keep])),'matrix_header':' | '.join(head),'audit_seconds':time.time()-t}
    return record,marker_rows
def main():
    QC=data_path('data_primary/audit',for_write=True)
    rows=pd.read_csv(ROOT/'primary_sample_manifest.tsv',sep='\t').to_dict('records')
    results=[];markers=[]
    for row in rows:
        data_path(row['matrix_path'],must_exist=True)
        res,mark=audit(row);results.append(res);markers.extend(mark)
        pd.DataFrame(results).to_csv(data_path('data_primary/audit/matrix_audit.tsv',for_write=True),sep='\t',index=False)
        pd.DataFrame(markers).to_csv(data_path('data_primary/audit/marker_audit.tsv',for_write=True),sep='\t',index=False)
        print(json.dumps(res),flush=True)
if __name__=='__main__':
    configure_data_root_cli(__doc__)
    main()
