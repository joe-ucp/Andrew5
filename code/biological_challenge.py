"""Frozen-score observational stress test; not an adverse-lineage or function assay.

Uses released noninteger expression without count resampling or inferred donor IDs.
Feature-missing zero-fill and shared-universe reranking are transport diagnostics.
"""
from data_paths import evidence_path, prepare_original_results
from pathlib import Path
import json,hashlib,os
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata
from reference_analysis import auc
from data_paths import data_path,configure_data_root_cli
ROOT=Path(__file__).resolve().parents[1]
EV=ROOT/'evidence'
MODULE='SMOC1 PLCE1 PAPPA2 ZNF331 ALDH1A1 SLC30A8 BTG2 TM4SF4 NR4A1 PCSK2'.split()

def sparse_rank_scores(x,up,down):
    values=[]
    for j in range(0,x.shape[1],128):
        twice_ranks=np.rint(2*rankdata(x[:,j:j+128].toarray(),axis=0,method='average')).astype(np.int64)
        numerator=twice_ranks[up].sum(axis=0)*len(down)-twice_ranks[down].sum(axis=0)*len(up)
        values.extend((numerator/(2*x.shape[0]*len(up)*len(down))).tolist())
    return np.asarray(values)

def descriptive_correlation(x,y):
    if len(x)<3 or np.ptp(x)==0 or np.ptp(y)==0:return None
    return float(np.corrcoef(rankdata(x),rankdata(y))[0,1])

def compare(high,low,label):
    row={'comparison':label,'high_n':len(high),'low_n':len(low),'high_mean':None,'low_mean':None,'mean_difference':None,'auc_high_vs_low':None,'high_fraction_in_low_10_90_interval':None,'high_median':None,'low_median':None,'status':'not_evaluable'}
    if len(high) and len(low):
        lo,hi=np.quantile(low,[.1,.9])
        row.update({'high_mean':float(high.mean()),'low_mean':float(low.mean()),'mean_difference':float(high.mean()-low.mean()),'auc_high_vs_low':auc(high,low),'high_fraction_in_low_10_90_interval':float(((high>=lo)&(high<=hi)).mean()),'high_median':float(np.median(high)),'low_median':float(np.median(low)),'status':'descriptive_only'})
    return row

def main():
    prepare_original_results()
    freeze=EV/'frozen_score.json'
    if not freeze.exists():raise RuntimeError('Frozen external score required; no feature fitting here')
    conf=json.loads(freeze.read_text())
    x=sparse.load_npz(data_path('data_primary/GSE217837/submitted_expression.npz',must_exist=True)).tocsr()
    genes=pd.read_csv(data_path('data_primary/GSE217837/features.tsv',must_exist=True),sep='\t')['symbol'].astype(str).tolist()
    cells=pd.read_csv(data_path('data_primary/GSE217837/barcodes.tsv',must_exist=True),sep='\t')
    idx={g:i for i,g in enumerate(genes)}
    if len(idx)!=len(genes):raise ValueError('Challenge feature symbols not unique')
    universe=conf['universe']; ui={g:i for i,g in enumerate(universe)}
    present=[g for g in universe if g in idx];missing=[g for g in universe if g not in idx]
    # The main operational transport diagnostic preserves the exact frozen denominator.
    # Missing submitted features are NOT asserted to be biological zero.
    mapper=sparse.csr_matrix((np.ones(len(present)),([ui[g] for g in present],[idx[g] for g in present])),shape=(len(universe),len(genes)))
    transported=mapper@x
    conversion=sparse_rank_scores(transported,[ui[g] for g in conf['up']],[ui[g] for g in conf['down']])
    pi={g:i for i,g in enumerate(present)}
    shared=sparse_rank_scores(x[[idx[g] for g in present]],[pi[g] for g in conf['up'] if g in pi],[pi[g] for g in conf['down'] if g in pi])
    module_present=[g for g in MODULE if g in idx]
    module_values=[]
    # Separately defined module denominator is the full released assay universe.
    for j in range(0,x.shape[1],128):
        ranks=(rankdata(x[:,j:j+128].toarray(),axis=0,method='average')-.5)/x.shape[0]
        module_values.extend(ranks[[idx[g] for g in module_present]].mean(axis=0).tolist())
    module=np.asarray(module_values)
    hormones=['INS','GCG','SST','PPY']
    if any(g not in idx for g in hormones+['SMOC1']):raise ValueError('Required hormone/SMOC1 feature missing')
    h=x[[idx[g] for g in hormones]].toarray()
    maxi=h.max(axis=0);arg=h.argmax(axis=0);unique=(h==maxi).sum(axis=0)==1
    strata=np.full(x.shape[1],'other_or_ambiguous',dtype=object)
    for i,name in enumerate(['beta_like','alpha_like','delta_like','PP_like']):strata[(arg==i)&(maxi>0)&unique]=name
    smoc=x[idx['SMOC1']].toarray().ravel();beta=strata=='beta_like'
    q25,q75=np.quantile(smoc[beta],[.25,.75]) if beta.any() else (np.nan,np.nan)
    high=beta&(smoc>=q75)&(smoc>0);low=beta&(smoc<=q25)
    overlap=high&low;high[overlap]=False;low[overlap]=False
    groups=np.where(high,'SMOC1_high',np.where(low,'SMOC1_low','not_extreme'))
    rows=[]
    for variant,score in [('frozen_universe_zero_fill_transport',conversion),('shared_feature_rerank_sensitivity',shared),('ten_gene_module',module)]:
        rows.append(compare(score[high],score[low],variant))
    pd.DataFrame(rows).to_csv(EV/'biological_challenge_results.tsv',sep='\t',index=False)
    summaries=[]
    for stratum in sorted(set(strata)):
        mask=strata==stratum
        summaries.append({'stratum':stratum,'cells':int(mask.sum()),'conversion_mean':float(conversion[mask].mean()),'conversion_median':float(np.median(conversion[mask])),'SMOC1_positive_cells':int((smoc[mask]>0).sum()),'module_mean':float(module[mask].mean()),'spearman_conversion_module':descriptive_correlation(conversion[mask],module[mask]),'spearman_conversion_SMOC1':descriptive_correlation(conversion[mask],smoc[mask])})
    pd.DataFrame(summaries).to_csv(EV/'biological_challenge_strata.tsv',sep='\t',index=False)
    cells=cells.assign(stratum=strata,SMOC1_expression=smoc,SMOC1_beta_group=groups,conversion_frozen_universe=conversion,conversion_shared_universe=shared,ten_gene_module=module)
    cells.to_csv(evidence_path('biological_challenge_cell_scores.tsv.gz'),sep='\t',index=False)
    status={'classification':'observational_stress_test_only','source_matrix_sha256':'7431fbad604e919b955a3b9f328c992be889a3b96326eb621c3d4dd50ce2b363','frozen_score_sha256':hashlib.sha256(freeze.read_bytes()).hexdigest(),'n_cells':x.shape[1],'n_submitted_features':x.shape[0],'frozen_universe_n':len(universe),'available_universe_n':len(present),'missing_universe_genes':missing,'available_up_n':sum(g in idx for g in conf['up']),'available_down_n':sum(g in idx for g in conf['down']),'missing_up_genes':[g for g in conf['up'] if g not in idx],'missing_down_genes':[g for g in conf['down'] if g not in idx],'module_genes':MODULE,'available_module_genes':module_present,'module_score_feature_overlap':sorted(set(MODULE)&(set(conf['up'])|set(conf['down']))),'module_denominator':'all released challenge features, mean within-cell fractional rank','beta_like_n':int(beta.sum()),'SMOC1_q25':float(q25),'SMOC1_q75':float(q75),'SMOC1_high_n':int(high.sum()),'SMOC1_low_n':int(low.sum()),'removed_tied_group_overlap_n':int(overlap.sum()),'comparisons':rows,'beta_conversion_module_spearman':descriptive_correlation(conversion[beta],module[beta]),'beta_conversion_SMOC1_spearman':descriptive_correlation(conversion[beta],smoc[beta]),'adverse_lineage_specificity_requirement':'UNRESOLVED: no adverse lineage labels, intervention assignment, functional assay, or verified donor mapping in this expression layer','limitations':['Noninteger submitted expression is not raw counts; no depth/UMI inference or multinomial resampling.','Barcode suffixes retained as grouping metadata only, not asserted donor labels.','Dominant hormone stratum describes measured state, not ancestry.','SMOC1 association and published trajectory module are not validated measures of dedifferentiation or dysfunction.','Feature overlap mechanically couples module and conversion scores; reported explicitly.','Cell-level distributions are descriptive; no inferential p-values or population confidence intervals.','Missing-feature zero fill is an operational transport convention, not measured biological absence; shared-feature sensitivity changes the rank denominator.']}
    status['SMOC1_quartile_thresholds_tied']=bool(q25==q75)
    status['SMOC1_realized_group_interpretation']='Positive versus zero expression because both quartile cutoffs are zero; not balanced quartiles' if q25==q75==0 else 'Prespecified expression-threshold groups, retaining all ties except high-low overlap'
    status['signed_rank_implementation']='Exact half-rank integer numerator divided once; tied-null scores equal exactly zero'
    (EV/'biological_challenge_status.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in status.items() if k!='missing_universe_genes'},indent=2))

if __name__ == '__main__':
    import argparse
    from pathlib import Path
    from public_runtime import execute_stage
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    execute_stage('challenge', args.data_root, args.output_dir)
