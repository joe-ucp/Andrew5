"""Report declared null ranges and operating points without changing frozen gate."""
from data_paths import evidence_path
import json
import numpy as np
import pandas as pd
from reference_analysis import EV,ROOT,references,auc

def main():
    s=pd.read_csv(evidence_path('calibration_simulations.tsv.gz'),sep='\t',float_precision='round_trip')
    out=[]
    for depth,df in s.groupby('depth'):
        pos=df[df.model=='PM'].score.to_numpy()
        for name,models in {'alpha_only':['alpha'],'low_mixtures':['alpha','mixture_0','mixture_0.01','mixture_0.05','mixture_0.1'],'stress_inclusive':[x for x in df.model.unique() if x!='PM']}.items():
            null=df[df.model.isin(models)].score.to_numpy();th=np.quantile(null,.95)
            sens=float((pos>th).mean());maxfpr=max(float((df[df.model==m].score>th).mean()) for m in models)
            out.append({'depth':depth,'null_range':name,'threshold_q95':th,'PM_sensitivity':sens,'pooled_null_FPR':float((null>th).mean()),'maximum_individual_null_FPR':maxfpr,'PM_exact_zero_fraction':float((pos==0).mean()),'null_exact_zero_fraction':float((null==0).mean()),'minimum_auc':min(auc(pos,df[df.model==m].score.to_numpy()) for m in models),'scope':'same_frozen_score_additional_reporting_not_new_gate'})
    pd.DataFrame(out).to_csv(EV/'calibration_operating_points.tsv',sep='\t',index=False)
    cfg=json.loads((EV/'frozen_score.json').read_text());bulk,sc=references();genes=cfg['up']+cfg['down'];rows=[]
    for prep,d in sc.items():
        for gene in genes:
            rows.append({'preparation':prep,'gene':gene,'direction':'up' if gene in cfg['up'] else 'down','fraction_cells_detected':float((d.loc[gene]>0).mean()),'total_counts':int(d.loc[gene].sum()),'bulk_four_pair_median_CPM':float((bulk.loc[gene,[c for _,a,b in cfg['pairs'] for c in (a,b)]]/bulk[[c for _,a,b in cfg['pairs'] for c in (a,b)]].sum()*1e6).median())})
    pd.DataFrame(rows).to_csv(EV/'reference_feature_capture.tsv',sep='\t',index=False)
    native=pd.read_csv(evidence_path('reference_cell_scores.tsv'),sep='\t',float_precision='round_trip')
    ss=[]
    for prep,d in sc.items():
        n=(d.loc[genes]>0).sum(axis=0)
        ss.append({'preparation':prep,'n_cells':d.shape[1],'median_total_counts':float(d.sum().median()),'median_detected_score_genes':float(n.median()),'mean_detected_score_genes':float(n.mean()),'zero_score_fraction':float((native[native.preparation==prep].score==0).mean())})
    pd.DataFrame(ss).to_csv(EV/'reference_capture_summary.tsv',sep='\t',index=False)
    print(pd.DataFrame(out).to_string(index=False))

if __name__ == '__main__':
    import argparse
    from pathlib import Path
    from public_runtime import execute_stage
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    execute_stage('calibration', args.data_root, args.output_dir)
