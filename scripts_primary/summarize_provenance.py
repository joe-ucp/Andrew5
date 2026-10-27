from pathlib import Path
import re,json,sys
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from data_paths import resolve_data_root,data_path,configure_data_root_cli

def main():
    DATA_ROOT=resolve_data_root(require_exists=True)
    md=data_path('data_primary/metadata',must_exist=True)
    rows=[]
    for p in sorted(md.glob('GSM672*.soft.txt')):
        p=data_path('data_primary/metadata/'+p.name,must_exist=True)
        txt=p.read_text(encoding='utf-8')
        def v(k):
            m=re.search(r'^'+re.escape(k)+r' = ([^\r\n]*)',txt,re.M)
            return m.group(1) if m else ''
        rows.append({'series':'GSE217837','accession':p.name.split('.')[0],'title':v('!Sample_title'),'condition':'; '.join(re.findall(r'!Sample_characteristics_ch1 = ([^\r\n]*)',txt)),'donor_id':'not specified in sample metadata','assay':v('!Sample_title'),'SMOC1_perturbation':'not annotated','metadata_path':str(p.relative_to(DATA_ROOT))})
    pd.DataFrame(rows).to_csv(data_path('data_primary/GSE217837_sample_manifest.tsv',for_write=True),sep='\t',index=False)
    audit=pd.read_csv(data_path('data_primary/audit/matrix_audit.tsv',must_exist=True),sep='\t')
    print(audit[['accession','donor','treatment_code','barcodes','published_QC_pass_cells','library_median','QC_library_median']].to_string(index=False))
    print(json.dumps({'samples':len(audit),'total_barcodes':int(audit.barcodes.sum()),'QC_pass_barcodes':int(audit.published_QC_pass_cells.sum()),'counts_sum':int(audit.sum_counts.sum()),'nnz_sum':int(audit.nnz.sum())}))

if __name__=='__main__':
    configure_data_root_cli(__doc__)
    main()
