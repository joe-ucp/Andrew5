"""Stream GEO's dense submitted expression layer to sparse CSR; do not assume counts."""
from pathlib import Path
import gzip,csv,json,re,time,os,sys
import numpy as np
import pandas as pd
from scipy import sparse
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from data_paths import data_path,configure_data_root_cli

def main():
    path=data_path('data_primary/GSE217837/GSE217837_scRNA_count_matrix.csv.gz',must_exist=True)
    values=[];indices=[];indptr=[0];genes=[]
    noninteger=negative=nonfinite=0;maximum=0.;total=0.;t=time.time()
    with gzip.open(path,'rt') as f:
        header=next(csv.reader([next(f)]));barcodes=header[1:]
        n=len(barcodes)
        for i,line in enumerate(f):
            name,nums=line.rstrip().split(',',1)
            name=name.strip('"')
            v=np.fromstring(nums,sep=',',dtype=np.float64)
            if v.size!=n:raise ValueError(f'Row length {name}: {v.size}!={n}')
            noninteger+=int((v!=np.floor(v)).sum());negative+=int((v<0).sum());nonfinite+=int((~np.isfinite(v)).sum())
            maximum=max(maximum,float(v.max()));total+=float(v.sum())
            ix=np.flatnonzero(v)
            indices.append(ix.astype(np.int32));values.append(v[ix]);indptr.append(indptr[-1]+len(ix));genes.append(name)
            if i%5000==0:print('processed genes',i,flush=True)
    x=sparse.csr_matrix((np.concatenate(values),np.concatenate(indices),np.array(indptr,dtype=np.int64)),shape=(len(genes),n))
    sparse.save_npz(data_path('data_primary/GSE217837/submitted_expression.npz',for_write=True),x)
    pd.Series(genes,name='symbol').to_csv(data_path('data_primary/GSE217837/features.tsv',for_write=True),sep='\t',index=False)
    pd.DataFrame({'barcode':barcodes,'barcode_suffix':[b.rsplit('_',1)[-1] for b in barcodes]}).to_csv(data_path('data_primary/GSE217837/barcodes.tsv',for_write=True),sep='\t',index=False)
    row={'matrix':'GSE217837_scRNA_count_matrix.csv.gz','features':x.shape[0],'cells':x.shape[1],'nnz':x.nnz,'noninteger_entries':noninteger,'negative_entries':negative,'nonfinite_entries':nonfinite,'max_value':maximum,'sum_values':total,'unique_symbols':len(set(genes)),'unique_barcodes':len(set(barcodes)),'barcode_suffix_counts':pd.Series([b.rsplit('_',1)[-1] for b in barcodes]).value_counts().to_dict(),'declared_layer':'filename says count_matrix; GEO sample description says integrated matrix','verified_layer':'noninteger submitted expression, normalization not established by matrix alone','donor_mapping':'barcode suffix groups do not establish donor identity or mapping to GSM without additional provenance','seconds':time.time()-t}
    (data_path('data_primary/GSE217837/structural_audit.json',for_write=True)).write_text(json.dumps(row,indent=2),encoding='utf-8')
    print(json.dumps(row),flush=True)

if __name__=='__main__':
    configure_data_root_cli(__doc__)
    main()
