"""Independent Baron author-labelled state/count-thinning validation, not conversion truth."""
from pathlib import Path
import argparse,json,os
import numpy as np
import pandas as pd
from scipy import sparse
from data_paths import data_path,resolve_data_root
from revision_benchmark import EV,METHODS,ROOT,load_protocol,sha,score_methods,corr
from reference_analysis import auc
from public_runtime import frozen_score, scientific_file, parameters_sha256
from revision_score import components_sparse,thin_sparse_counts
from rank_metric_comparison import signed_log_error_radii

def main():
    plan=load_protocol();conf=frozen_score()
    program_path=scientific_file('evidence/revision_01/benchmark_programs.json')
    locked=json.loads(program_path.read_text());programs=locked['programs']
    if (EV/'independent_status.json').exists():raise FileExistsError('Independent analysis already completed')
    bulk=data_path('revision_01/independent',for_write=True);bulk.mkdir(parents=True,exist_ok=True)
    registry=json.loads(scientific_file('evidence/revision_01/baron_input_manifest.json').read_text())
    reports=[];summaries=[];native_scores=[];thinrows=[];allcoverage=[];missing_universe=[];native_components=[];component_draws=[]
    for donor,record in enumerate(sorted(registry,key=lambda x:x['name']),1):
        path=data_path('revision_01/sources/'+record['name'],must_exist=True)
        if sha(path)!=record['sha256']:raise ValueError('Independent source hash mismatch')
        print('Independent donor',donor,path.name,flush=True)
        frame=pd.read_csv(path,index_col=0);labels=frame['assigned_cluster'].astype(str)
        if not frame.index.is_unique:raise ValueError('Original independent cell IDs not unique')
        count=frame.iloc[:,2:];arr=count.to_numpy()
        if not np.isfinite(arr).all() or (arr<0).any() or (arr!=np.floor(arr)).any():raise ValueError('Independent matrix is not nonnegative integer counts')
        total_all=arr.sum(axis=1);use=labels.isin(['alpha','beta']).to_numpy()
        totals=total_all[use];celltypes=labels.to_numpy()[use];ids=frame.index.to_numpy()[use].astype(str)
        original_genes=conf['universe'];x=sparse.csc_matrix(count.loc[use].reindex(columns=original_genes,fill_value=0).to_numpy().T)
        represented=[g for g in original_genes if g in count.columns]
        missing_universe.extend({'donor':donor,'gene':g} for g in original_genes if g not in count.columns)
        shared_idx=pd.Index(original_genes).get_indexer(represented)
        modes={'original_universe':(x,original_genes),'shared_features':(x[shared_idx],represented)}
        for program,feature in programs.items():
            for side in ['up','down']:
                missing=[g for g in feature[side] if g not in count.columns]
                allcoverage.append({'donor':donor,'program':program,'block':side,'n_selected':len(feature[side]),'n_represented':len(feature[side])-len(missing),'fraction_represented':1-len(missing)/len(feature[side]),'missing_genes':';'.join(missing),'original_universe_n':len(original_genes),'represented_universe_n':len(represented),'qualified_membership_95pct':1-len(missing)/len(feature[side])>=.95})
        for mode,(matrix,genes) in modes.items():
            lookup={g:i for i,g in enumerate(genes)};up=[lookup[g] for g in conf['up'] if g in lookup];down=[lookup[g] for g in conf['down'] if g in lookup]
            values=score_methods(matrix,totals,genes,programs);comp=components_sparse(matrix,up,down)
            log_roundoff_radius=signed_log_error_radii(matrix,totals,genes,programs)
            native_components.append(pd.DataFrame({'donor':donor,'mode':mode,'cell_id':ids,'cell_type':celltypes,'total_counts':totals,**comp}))
            np.savez_compressed(bulk/f'donor{donor}_{mode}_scores.npz',scores=values,cell_id=ids,cell_type=celltypes,total_counts=totals,programs=np.array(list(programs)),methods=np.array(METHODS),log_roundoff_radius=log_roundoff_radius,log_state_is_beta=(celltypes=='beta'))
            for pi,program in enumerate(programs):
                for mi,method in enumerate(METHODS):
                    v=values[pi,mi];a=v[celltypes=='alpha'];b=v[celltypes=='beta'];threshold=float(np.quantile(a,.95))
                    summaries.append({'donor':donor,'mode':mode,'program':program,'method':method,'alpha_n':len(a),'beta_n':len(b),'beta_vs_alpha_AUROC':auc(b,a),'alpha_native_q95':threshold,'achieved_alpha_FPR':float(np.mean(a>threshold)),'beta_detection_fraction':float(np.mean(b>threshold)),'alpha_score_depth_r':corr(a,totals[celltypes=='alpha']),'beta_score_depth_r':corr(b,totals[celltypes=='beta']),'scope':'author cell-state labels; not conversion truth'})
            native_scores.append((mode,values,comp))
        # Count draws use the original universe plus all outside mass. Shared-feature
        # sensitivities use exactly the same molecules with missing zero rows removed.
        native_by_mode={m:(v,c) for m,v,c in native_scores[-2:]}
        for depth in plan['independent_validation']['thinning_depths']:
            eligible=totals>=depth;types=celltypes[eligible]
            for replicate,seed in enumerate(plan['independent_validation']['seeds']):
                rng=np.random.default_rng(np.random.SeedSequence([seed,donor,depth]))
                thinned=thin_sparse_counts(x[:,eligible],totals[eligible],depth,rng)
                for mode,(_,genes) in modes.items():
                    matrix=thinned if mode=='original_universe' else thinned[shared_idx]
                    vals=score_methods(matrix,np.full(eligible.sum(),depth),genes,programs)
                    lookup={g:i for i,g in enumerate(genes)};up=[lookup[g] for g in conf['up'] if g in lookup];down=[lookup[g] for g in conf['down'] if g in lookup]
                    components=components_sparse(matrix,up,down);native_v,native_c=native_by_mode[mode]
                    for celltype in ['alpha','beta']:
                        inside=types==celltype;original=(celltypes==celltype)&eligible;count_n=int(inside.sum());valid=count_n>=20
                        base={'donor':donor,'mode':mode,'cell_type':celltype,'depth':depth,'replicate':replicate,'seed':seed,'population_n':int((celltypes==celltype).sum()),'eligible_n':count_n,'excluded_n':int((celltypes==celltype).sum())-count_n,'evaluable':valid}
                        for pi,program in enumerate(programs):
                            for mi,method in enumerate(METHODS):
                                before=float(native_v[pi,mi,original].mean()) if valid else np.nan
                                after=float(vals[pi,mi,inside].mean()) if valid else np.nan
                                thinrows.append({**base,'program':program,'method':method,'native_eligible_mean':before,'thinned_mean':after,'thinning_change':after-before})
                        for component in ['score','occupancy_component','occupancy_up','occupancy_down','positive_rank_residual','k_up','k_down']:
                            before=float(native_c[component][original].mean()) if valid else np.nan
                            after=float(components[component][inside].mean()) if valid else np.nan
                            component_draws.append({**base,'component':component,'native_eligible_mean':before,'thinned_mean':after,'thinning_change':after-before})
            print('  depth',depth,'completed',flush=True)
        reports.append({'donor':donor,'file':record['name'],'sha256':record['sha256'],'all_cells':len(frame),'alpha_beta_cells':int(use.sum()),'unique_original_cell_ids':bool(frame.index.is_unique),'unique_genes':bool(count.columns.is_unique),'integer_nonnegative':True})
    cells=pd.concat(native_components,ignore_index=True);cells.to_csv(bulk/'baron_native_components.tsv.gz',sep='\t',index=False)
    native=[]
    for (donor,mode,label),g in cells.groupby(['donor','mode','cell_type']):
        for component in ['score','occupancy_component','occupancy_up','occupancy_down','positive_rank_residual','k_up','k_down']:
            native.append({'donor':donor,'mode':mode,'cell_type':label,'component':component,'n':len(g),'mean':float(g[component].mean()),'median_depth':float(g.total_counts.median()),'spearman_depth':corr(g[component],g.total_counts)})
    pd.DataFrame(native).to_csv(EV/'independent_native_components.tsv',sep='\t',index=False)
    pd.DataFrame(summaries).to_csv(EV/'independent_state_benchmark.tsv',sep='\t',index=False)
    pd.DataFrame(allcoverage).to_csv(EV/'independent_feature_coverage.tsv',sep='\t',index=False)
    pd.DataFrame(missing_universe,columns=['donor','gene']).to_csv(EV/'independent_missing_universe.tsv',sep='\t',index=False)
    draws=pd.DataFrame(thinrows);componentdraws=pd.DataFrame(component_draws)
    draws.to_csv(bulk/'baron_method_thinning_draws.tsv.gz',sep='\t',index=False)
    componentdraws.to_csv(bulk/'baron_component_thinning_draws.tsv.gz',sep='\t',index=False)
    for frame,name,lastkeys in [(draws,'independent_thinning',['program','method']),(componentdraws,'independent_component_thinning',['component'])]:
        keys=['donor','mode','cell_type','depth']+lastkeys
        summary=frame.groupby(keys,dropna=False).agg(population_n=('population_n','first'),eligible_n=('eligible_n','first'),excluded_n=('excluded_n','first'),native_eligible_mean=('native_eligible_mean','first'),thinned_mean=('thinned_mean','mean'),thinning_change=('thinning_change','mean'),technical_repeat_SD=('thinning_change','std'),valid_repeats=('thinning_change','count')).reset_index()
        summary.to_csv(EV/(name+'.tsv'),sep='\t',index=False)
    status={'status':'COMPLETE','scientific_parameters_sha256':parameters_sha256(),'program_sha256':sha(program_path),'source_reports':reports,'all_source_cells':sum(r['all_cells'] for r in reports),'alpha_beta_cells':sum(r['alpha_beta_cells'] for r in reports),'true_thinning_replicates':5,'cell_identity':'original cell_id; duplicate bare barcodes never merged','limitations':['Published celltype labels are state labels, not conversion or lineage truth','Four source donors include donor1 age17 and donor4 with type2diabetes; not four healthy adults','Submitted human files8569 cells differ from reported filtered8629','Post-review specified independent measurement benchmark; no fitting to independent labels or treatment effects']}
    (EV/'independent_status.json').write_text(json.dumps(status,indent=2),encoding='utf-8')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    from public_runtime import execute_stage
    execute_stage('independent', args.data_root, args.output_dir)
