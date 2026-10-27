"""Fixed feature/scoring benchmark; results use configured output locations."""
from data_paths import evidence_path
from public_runtime import protocol, frozen_score, parameters_sha256, scientific_file
from pathlib import Path
import argparse,datetime,hashlib,json,sys
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata,spearmanr
from data_paths import data_path,resolve_data_root
from reference_analysis import references,auc,ROOT
from primary_analysis import map_universe,load_sample
from revision_score import components_sparse,thin_sparse_counts

EV=ROOT/'evidence/revision_01'
METHODS=['signed_rank','signed_log1p_CPT10K','ucell_form_1500']

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def load_protocol():
    return protocol('revision_01')

def freeze_programs():
    plan=load_protocol();conf=frozen_score()
    table=pd.read_csv(evidence_path('reference_feature_selection.tsv'),sep='\t',float_precision='round_trip').set_index('gene').loc[conf['universe']]
    programs={'original':{'up':conf['up'],'down':conf['down'],'omitted_donor':None}}
    cols=[p[0]+'_log2_CPM_change' for p in conf['pairs']]
    rows=[]
    for omitted in range(4):
        donor=conf['pairs'][omitted][0];values=table[[c for i,c in enumerate(cols) if i!=omitted]].to_numpy()
        med=np.median(values,axis=1);genes=table.index.to_numpy()
        up=sorted(np.flatnonzero((values>0).all(axis=1)),key=lambda i:(-med[i],genes[i]))[:50]
        down=sorted(np.flatnonzero((values<0).all(axis=1)),key=lambda i:(med[i],genes[i]))[:50]
        name='omit_'+donor;programs[name]={'up':genes[up].tolist(),'down':genes[down].tolist(),'omitted_donor':donor}
        for side in ['up','down']:
            a=set(conf[side]);b=set(programs[name][side]);rows.append({'program':name,'side':side,'selected_n':len(b),'overlap_original':len(a&b),'jaccard_original':len(a&b)/len(a|b)})
    result={'scientific_parameters_sha256':parameters_sha256(),'original_score_sha256':sha(scientific_file('evidence/frozen_score.json')),
            'reference_table_sha256':sha(evidence_path('reference_feature_selection.tsv')),'programs':programs,
            'universe_policy':'original frozen universe; LODO is conditional feature stability, not unbiased held-out validation',
            'scorers':METHODS,'ucell_exact':plan['benchmark']['ucell_exact']}
    path=EV/'benchmark_programs.json'
    if path.exists():
        if json.loads(path.read_text())!=result:raise RuntimeError('Locked benchmark programs changed')
    else:
        path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        (EV/'benchmark_program_freeze.json').write_text(json.dumps({'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'program_sha256':sha(path),'program_construction_precedes_treatment_scoring_in_this_execution':True,'historical_treatment_outcomes_already_known':True},indent=2),encoding='utf-8')
    pd.DataFrame(rows).to_csv(EV/'benchmark_feature_stability.tsv',sep='\t',index=False)
    membership=[]
    for side in ['up','down']:
        for gene in sorted(set().union(*(set(v[side]) for v in programs.values()))):
            membership.append({'gene':gene,'side':side,'original':gene in conf[side],'LODO_inclusion_count':sum(gene in v[side] for k,v in programs.items() if k!='original')})
    pd.DataFrame(membership).to_csv(EV/'benchmark_gene_membership.tsv',sep='\t',index=False)
    return conf,programs

def score_methods(x,totals,genes,programs):
    """All methods use the same universe; UCell-form uses the published 2026 normalization."""
    x=sparse.csc_matrix(x,copy=True);x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
    if not np.isfinite(x.data).all() or (x.data<0).any():raise ValueError('Nonnegative finite counts required')
    totals=np.asarray(totals);G=x.shape[0];names=list(programs);P=len(names);lookup={g:i for i,g in enumerate(genes)}
    up=np.zeros((P,G),dtype=np.int64);down=up.copy();nu=[];nd=[]
    for i,p in enumerate(programs.values()):
        u=[lookup[g] for g in p['up'] if g in lookup];v=[lookup[g] for g in p['down'] if g in lookup]
        if not u or not v:raise ValueError('Empty represented score block')
        up[i,u]=1;down[i,v]=1;nu.append(len(u));nd.append(len(v))
    nu=np.array(nu);nd=np.array(nd);iw=up*nd[:,None]-down*nu[:,None]
    weights=up/nu[:,None]-down/nd[:,None]
    result=np.empty((P,3,x.shape[1]));R=1500
    valid_u=(G>=R)&(nu<R)&(nd<R)
    for j in range(x.shape[1]):
        a,b=x.indptr[j:j+2];idx=x.indices[a:b];v=x.data[a:b];K=len(v);z=G-K
        twice=np.rint(2*rankdata(v,method='average')).astype(np.int64)
        result[:,0,j]=(iw[:,idx]@(twice+z-1))/(2*G*nu*nd)
        result[:,1,j]=weights[:,idx]@np.log1p(10000*v/max(float(totals[j]),1))
        baseline=min((G+K+1)/2,R);delta=np.minimum(K+1-twice/2,R)-baseline
        rank_up=nu*baseline+up[:,idx]@delta;rank_down=nd*baseline+down[:,idx]@delta
        u=(nu*R-rank_up)/(nu*R-nu*(nu+1)/2);d=(nd*R-rank_down)/(nd*R-nd*(nd+1)/2)
        result[:,2,j]=np.where(valid_u,np.maximum(0,u-d),np.nan)
    return result

def corr(a,b):
    return float(spearmanr(a,b).statistic) if len(a)>1 and np.ptp(a)>0 and np.ptp(b)>0 else np.nan

def run_primary():
    conf,programs=freeze_programs();outdir=data_path('revision_01/benchmark',for_write=True);outdir.mkdir(parents=True,exist_ok=True)
    if (EV/'benchmark_primary_status.json').exists():raise FileExistsError('Benchmark completed; no implicit overwrite')
    _,ref=references();cal=[];thresholds={};refrecords=[]
    for prep,frame in ref.items():
        x=sparse.csc_matrix(frame.reindex(conf['universe'],fill_value=0).to_numpy());totals=frame.sum(axis=0).to_numpy()
        values=score_methods(x,totals,conf['universe'],programs)
        for pi,program in enumerate(programs):
            for mi,method in enumerate(METHODS):
                v=values[pi,mi];refrecords.append(pd.DataFrame({'program':program,'method':method,'preparation':prep,'cell_id':frame.columns,'score':v,'total_counts':totals}))
    refcells=pd.concat(refrecords,ignore_index=True);refcells.to_csv(outdir/'reference_benchmark_cells.tsv.gz',sep='\t',index=False)
    for (program,method),g in refcells.groupby(['program','method']):
        neg=g.loc[g.preparation=='alpha','score'].to_numpy();pos=g.loc[g.preparation=='PM','score'].to_numpy();t=float(np.quantile(neg,.95));thresholds[program,method]=t
        cal.append({'program':program,'method':method,'native_PM_vs_alpha_AUROC':auc(pos,neg),'alpha_q95':t,'alpha_FPR':float(np.mean(neg>t)),'PM_preparation_detection_fraction':float(np.mean(pos>t)),'alpha_n':len(neg),'PM_n':len(pos)})
    pd.DataFrame(cal).to_csv(EV/'benchmark_reference.tsv',sep='\t',index=False)
    original=pd.read_csv(evidence_path('primary_cell_scores.tsv.gz'),sep='\t',float_precision='round_trip')
    manifest=pd.read_csv(ROOT/'primary_sample_manifest.tsv',sep='\t');manifest=manifest[manifest.treatment_code.isin(['D','H_GLP'])]
    results=[];checks=[];inputs=[]
    for row in manifest.to_dict('records'):
        print('Benchmark primary',row['accession'],flush=True)
        x,feat,bc=load_sample(row);old=original[(original.accession==row['accession'])&(original.annotation=='alpha_related')]
        if old.empty:
            old=original[(original.accession==row['accession'])&(original.annotation=='alpha_like')]
        if old.empty:raise RuntimeError('No original alpha population')
        ix=pd.Index(bc.barcode).get_indexer(old.barcode)
        if (ix<0).any():raise ValueError('Unmatched original cell barcode')
        totals=np.asarray(x[:,ix].sum(axis=0)).ravel();mapped=map_universe(x[:,ix],feat,conf['universe'])
        vals=score_methods(mapped,totals,conf['universe'],programs)
        maximum=float(np.max(np.abs(vals[0,0]-old.score.to_numpy())))
        if maximum>1e-15:raise AssertionError('Frozen primary score reproduction failed')
        checks.append({'accession':row['accession'],'n':len(ix),'maximum_score_error':maximum})
        np.savez_compressed(outdir/(row['accession']+'_scores.npz'),scores=vals,barcodes=old.barcode.to_numpy(dtype=str),total_counts=totals,programs=np.array(list(programs)),methods=np.array(METHODS))
        for pi,program in enumerate(programs):
            for mi,method in enumerate(METHODS):
                v=vals[pi,mi];results.append({'accession':row['accession'],'donor':row['donor'],'treatment':row['treatment_code'],'program':program,'method':method,'n':len(v),'mean':float(v.mean()),'median':float(np.median(v)),'spearman_score_depth':corr(v,totals),'fraction_above_reference_q95':float(np.mean(v>thresholds[program,method]))})
        for rel in [row['features_path'],row['barcodes_path'],'data_primary/sparse_cache/'+row['accession']+'.npz']:
            inputs.append({'path':rel,'sha256':sha(data_path(rel,must_exist=True))})
    table=pd.DataFrame(results);table.to_csv(EV/'benchmark_primary_samples.tsv',sep='\t',index=False)
    d=table[table.treatment=='D'].set_index(['donor','program','method']);t=table[table.treatment=='H_GLP'].set_index(['donor','program','method']);t=t.loc[d.index]
    contrasts=d.index.to_frame(index=False);contrasts['control_mean']=d['mean'].values;contrasts['treated_mean']=t['mean'].values;contrasts['delta']=t['mean'].values-d['mean'].values
    contrasts['control_n']=d.n.values;contrasts['treated_n']=t.n.values;contrasts.to_csv(EV/'benchmark_primary_contrasts.tsv',sep='\t',index=False)
    (EV/'benchmark_primary_status.json').write_text(json.dumps({'status':'COMPLETE','scientific_parameters_sha256':parameters_sha256(),'program_sha256':sha(EV/'benchmark_programs.json'),'all_program_method_combinations_reported':len(programs)*3,'original_score_checks':checks,'input_hashes':inputs,'scientific_scope':'Post-review alternative features/scorers; descriptive transfer benchmark, no new confirmatory biological test'},indent=2),encoding='utf-8')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    from public_runtime import execute_stage
    execute_stage('benchmark', args.data_root, args.output_dir)
