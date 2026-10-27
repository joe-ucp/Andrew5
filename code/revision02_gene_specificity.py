"""Revision02: fixed-gene interpretation and baseline-matched program controls.

No treatment sample is read until the baseline-only candidate/program files are
saved and hashed. Historical treatment effects were already known; this is an
internally frozen post-review descriptive extension, not outcome blinding.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os
import sys
import numpy as np
import pandas as pd
from data_paths import data_path,evidence_path,resolve_data_root,ENV_NAME
from primary_analysis import map_universe
from public_runtime import protocol, frozen_score, parameters_sha256, scientific_file, verify_registered_data
from revision02_gene_math import (summarize_genes,program_numerators,donor_contrasts,
                                 equal_donor_mean,match_candidates,draw_matched_programs)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'evidence/revision_02'
sys.path.insert(0,str(ROOT/'scripts_primary'))
from audit_primary import load_sample
DONORS=['D1','D2','D3','D4']
LOCKED_MATCHING={'seed':8202,'programs':500,'maximum_attempts':50000,
                'pbar_caliper':.025,'Lbar_caliper':.5,'nearest_candidates':200,
                'minimum_candidates':10,'donor_block_detection_caliper':.025,
                'donor_block_log_abundance_caliper':.25,'abundance_offset':.01,
                'up_genes':50,'down_genes':50}


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def write_json(path,value):
    with Path(path).open('x',encoding='utf-8') as f:json.dump(value,f,indent=2);f.write('\n')


def load_protocol():
    plan, config = protocol('revision_02'), frozen_score()
    if len(config['up'])!=50 or len(config['down'])!=50 or len(set(config['up']+config['down']))!=100:
        raise RuntimeError('The unchanged 50/50 original program is required')
    # This literal design is also retained in the program-freeze receipt. The
    # published parameters must contain this exact matching configuration.
    declared_plan=plan.get('matched_programs',{})
    declared={'seed':declared_plan.get('rng_seed'),'programs':declared_plan.get('n_programs'),
              'maximum_attempts':declared_plan.get('max_attempts'),
              'pbar_caliper':declared_plan.get('slot_calipers',{}).get('absolute_pbar'),
              'Lbar_caliper':declared_plan.get('slot_calipers',{}).get('absolute_mean_log2_abundance'),
              'nearest_candidates':declared_plan.get('candidate_limit'),
              'minimum_candidates':declared_plan.get('minimum_candidates_per_slot'),
              'donor_block_detection_caliper':declared_plan.get('acceptance_each_donor_and_block',{}).get('mean_detection_abs_mismatch'),
              'donor_block_log_abundance_caliper':declared_plan.get('acceptance_each_donor_and_block',{}).get('mean_log2_abundance_abs_mismatch'),
              'abundance_offset':.01,'up_genes':declared_plan.get('up_n'),'down_genes':declared_plan.get('down_n')}
    for key,value in LOCKED_MATCHING.items():
        if declared.get(key)!=value:raise RuntimeError('Revision02 matching plan mismatch: '+key)
    if not all(declared_plan.get(key) is True for key in ['exclude_original100','unique_within_program','cross_program_reuse']):
        raise RuntimeError('Revision02 candidate exclusion/reuse policy mismatch')
    return plan,config


def sample_summary(row,universe,inputs):
    """Hash-verified original sparse counts with fixed revision01 alpha membership."""
    sources=pd.read_csv(ROOT/'source_manifest.tsv',sep='\t').set_index('local_path')
    for key in ['matrix_path','features_path','barcodes_path']:
        relative=row[key].replace('\\','/');path=data_path(relative,must_exist=True);actual=sha(path)
        if actual!=sources.loc[relative,'sha256']:raise ValueError('Source checksum mismatch: '+relative)
        inputs[relative]=actual
    cache_relative='data_primary/sparse_cache/'+row['accession']+'.npz'
    cache=data_path(cache_relative)
    if cache.exists():
        cache=verify_registered_data(cache_relative)
        inputs[cache_relative]=sha(cache)
    cell_relative='revision_01/'+row['accession']+'_cell_components.tsv.gz'
    cell_path=verify_registered_data(cell_relative)
    actual=sha(cell_path)
    inputs[cell_relative]=actual
    cells=pd.read_csv(cell_path,sep='\t',usecols=['accession','donor','treatment','barcode','source_cell_index','library_counts','alpha_broad'])
    if not (cells.accession.eq(row['accession'])&cells.donor.eq(row['donor'])&cells.treatment.eq(row['treatment_code'])).all():raise ValueError('Population sample metadata mismatch')
    cells=cells[cells.alpha_broad].copy()
    counts,features,barcodes=load_sample(row,cache=False)
    ix=cells.source_cell_index.to_numpy(dtype=int)
    if len(np.unique(ix))!=len(ix) or not np.array_equal(barcodes.barcode.iloc[ix],cells.barcode):raise ValueError('Original source-index/barcode mapping differs')
    totals=np.asarray(counts.sum(axis=0)).ravel()[ix]
    if not np.array_equal(totals,cells.library_counts):raise ValueError('Original full-library totals differ')
    mapped=map_universe(counts[:,ix],features,universe)
    result=summarize_genes(mapped,totals)
    result['source_index_sha256']=np.asarray(hashlib.sha256(ix.astype('<i8').tobytes()).hexdigest())
    return result


def gene_rows(summary,row,genes,original):
    sign=np.r_[np.ones(50),-np.ones(50)]
    frame=pd.DataFrame({'gene':np.asarray(genes)[original],'block':['up']*50+['down']*50,
                        'donor':row['donor'],'treatment':row['treatment_code'],'accession':row['accession'],
                        'n_cells':int(summary['n'])})
    for key in ['detected_n','detection_fraction','mean_CP10K','positive_mean_raw','positive_mean_CP10K','positive_mean_log1p_CP10K']:
        frame[key]=summary[key][original]
    frame['signed_occupancy']=sign*summary['detected_n'][original]/(100*int(summary['n']))
    frame['signed_residual']=sign*summary['centered_rank_sum'][original]/(100*int(summary['G'])*int(summary['n']))
    frame['signed_zero_baseline_rank']=frame.signed_occupancy+frame.signed_residual
    return frame


def gene_tables(samples):
    index=['donor','gene','block'];control=samples[samples.treatment=='D'].set_index(index)
    treated=samples[samples.treatment=='H_GLP'].set_index(index).loc[control.index]
    contrast=control.index.to_frame(index=False)
    for key in ['n_cells','detected_n','detection_fraction','mean_CP10K','positive_mean_raw',
                'positive_mean_CP10K','positive_mean_log1p_CP10K','signed_occupancy','signed_residual','signed_zero_baseline_rank']:
        contrast['control_'+key]=control[key].to_numpy();contrast['treated_'+key]=treated[key].to_numpy()
        contrast['delta_'+key]=treated[key].to_numpy()-control[key].to_numpy()
    rows=[]
    for (gene,block),group in contrast.groupby(['gene','block'],sort=True):
        values=group.delta_signed_occupancy.to_numpy()
        rows.append({'gene':gene,'block':block,'mean_signed_occupancy_delta':values.mean(),
                     'mean_absolute_occupancy_delta':np.abs(values).mean(),
                     'positive_occupancy_donors':int((values>0).sum()),'negative_occupancy_donors':int((values<0).sum()),
                     'zero_occupancy_donors':int((values==0).sum()),
                     'mean_signed_residual_delta':group.delta_signed_residual.mean(),
                     'mean_detection_fraction_delta':group.delta_detection_fraction.mean(),
                     'mean_positive_raw_delta':group.delta_positive_mean_raw.mean(),
                     'mean_positive_CP10K_delta':group.delta_positive_mean_CP10K.mean(),
                     'mean_positive_log1p_CP10K_delta':group.delta_positive_mean_log1p_CP10K.mean(),
                     'positive_expression_evaluable_donors':int(group.delta_positive_mean_CP10K.notna().sum()),
                     'positive_CP10K_increase_donors':int((group.delta_positive_mean_CP10K>0).sum()),
                     'positive_CP10K_decrease_donors':int((group.delta_positive_mean_CP10K<0).sum())})
    summary=pd.DataFrame(rows).sort_values(['block','mean_absolute_occupancy_delta','gene'],ascending=[True,False,True])
    summary['within_block_display_rank']=summary.groupby('block',sort=False).cumcount()+1
    summary['display_top10']=summary.within_block_display_rank<=10
    ranked=summary.set_index('gene').within_block_display_rank
    contrast['within_block_display_rank']=contrast.gene.map(ranked)
    concentrations=[]
    for block,group in summary.groupby('block',sort=True):
        ordered=group.sort_values('within_block_display_rank').gene
        sets=[('equal_donor_mean_absolute',group.set_index('gene').loc[ordered].mean_absolute_occupancy_delta.to_numpy())]
        for donor in DONORS:
            values=contrast[(contrast.donor==donor)&(contrast.block==block)].set_index('gene').loc[ordered].delta_signed_occupancy.to_numpy()
            sets.append((donor,np.abs(values)))
        for donor,absolute in sets:
            total=absolute.sum();shares=absolute/total if total>0 else np.full(50,np.nan)
            concentrations.append({'donor':donor,'block':block,'absolute_occupancy_total':total,
                                   'fixed_top5_absolute_fraction':np.sum(shares[:5]),
                                   'fixed_top10_absolute_fraction':np.sum(shares[:10]),
                                   'herfindahl_absolute_concentration':np.sum(shares**2),
                                   'ranking':'equal-donor mean absolute occupancy contribution; gene-symbol tie break'})
    return contrast,summary,pd.DataFrame(concentrations)


def summarize_matched(contrast_sets):
    rows=[];draws=[]
    for component in ['score','occupancy','residual']:
        entries=[(donor,contrast[component]) for donor,contrast in zip(DONORS,contrast_sets)]
        entries.append(('equal_donor_mean',equal_donor_mean(contrast_sets,component)))
        for donor,item in entries:
            integers=item['numerators'];denom=int(item['denominator'])
            values=np.array([int(x)/denom for x in integers]);original=values[0];random=values[1:]
            below=sum(int(x)<int(integers[0]) for x in integers[1:]);equal=sum(int(x)==int(integers[0]) for x in integers[1:])
            rows.append({'donor':donor,'component':component,'original_delta':original,'matched_programs':len(random),
                         'matched_mean':random.mean(),'matched_SD':random.std(ddof=1),'matched_minimum':random.min(),
                         'matched_q025':np.quantile(random,.025,method='linear'),'matched_q05':np.quantile(random,.05,method='linear'),
                         'matched_median':np.median(random),'matched_q95':np.quantile(random,.95,method='linear'),
                         'matched_q975':np.quantile(random,.975,method='linear'),'matched_maximum':random.max(),
                         'descriptive_midrank_percentile':(below+.5*equal)/len(random),
                         'n_less':below,'n_equal':equal,'n_greater':len(random)-below-equal,
                         'descriptive_rank_interval_lower':below/len(random),
                         'descriptive_rank_interval_upper':(below+equal)/len(random),
                         'quantile_method':'linear','not_a_randomization_p_value':True})
            for program,(numerator,value) in enumerate(zip(integers,values)):
                draws.append({'donor':donor,'component':component,'program':program,'original':program==0,
                              'delta':value,'exact_numerator':str(int(numerator)),'exact_denominator':str(denom)})
    return pd.DataFrame(rows),pd.DataFrame(draws)


def main():
    plan,config=load_protocol()
    bulk=data_path('revision_02/gene_specificity',for_write=True)
    owned=['gene_sample_summaries.tsv','gene_donor_contrasts.tsv','gene_contribution_summary.tsv',
           'gene_concentration.tsv','matched_candidate_support.tsv','matched_baseline_quality.tsv',
           'matched_program_comparison.tsv','gene_program_freeze.json','gene_specificity_status.json']
    for name in owned:
        if (OUT/name).exists():raise FileExistsError('No implicit overwrite: '+str(OUT/name))
    if bulk.exists():raise FileExistsError('Gene-specificity bulk directory already exists; explicit resume/review required')
    OUT.mkdir(parents=True,exist_ok=True);bulk.mkdir(parents=True)
    universe=config['universe'];index={gene:i for i,gene in enumerate(universe)}
    original=np.array([index[g] for g in config['up']+config['down']],dtype=int)
    manifest=pd.read_csv(ROOT/'primary_sample_manifest.tsv',sep='\t')
    manifest=manifest[manifest.treatment_code.isin(['D','H_GLP'])].sort_values(['donor','treatment_code'])
    if len(manifest)!=8 or sorted(manifest.donor.unique())!=DONORS or not all(
        len(group)==2 and set(group.treatment_code)=={'D','H_GLP'} for _,group in manifest.groupby('donor')):
        raise ValueError('Four paired donors required')
    summaries={};inputs={};sample_rows=[]
    # Strict first phase: vehicle samples only, before any new treated summaries.
    vehicle=manifest[manifest.treatment_code=='D'].sort_values('donor')
    for row in vehicle.to_dict('records'):
        summary=sample_summary(row,universe,inputs);summaries[(row['donor'],'D')]=summary
        np.savez_compressed(bulk/(row['accession']+'_gene_summaries.npz'),**summary)
        sample_rows.append(gene_rows(summary,row,universe,original))
        print('Vehicle baseline complete:',row['accession'],flush=True)
    p=np.stack([summaries[(d,'D')]['detection_fraction'] for d in DONORS])
    L=np.stack([np.log2(summaries[(d,'D')]['mean_CP10K']+.01) for d in DONORS])
    np.savez_compressed(bulk/'vehicle_matching_covariates.npz',genes=np.asarray(universe),p=p,L=L)
    candidates,support=match_candidates(universe,original,p,L)
    pd.DataFrame(support).to_csv(OUT/'matched_candidate_support.tsv',sep='\t',index=False)
    candidate_rows=[{'slot':slot,'original_gene':universe[original[slot]],'block':'up' if slot<50 else 'down',
                     'candidate_rank':rank+1,'candidate_gene':universe[int(gene)]}
                    for slot,values in enumerate(candidates) for rank,gene in enumerate(values)]
    pd.DataFrame(candidate_rows).to_csv(bulk/'candidate_genes.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
    programs,quality,matching=draw_matched_programs(original,candidates,p,L,seed=8202,target=500,max_attempts=50000)
    np.savez_compressed(bulk/'matched_programs.npz',programs=programs,genes=np.asarray(universe),original=original)
    membership=[{'program':j+1,'slot':slot,'block':'up' if slot<50 else 'down','matched_gene':universe[int(g)],
                 'original_gene':universe[original[slot]]} for j,program in enumerate(programs) for slot,g in enumerate(program)]
    pd.DataFrame(membership,columns=['program','slot','block','matched_gene','original_gene']).to_csv(bulk/'matched_program_membership.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
    quality_frame=pd.DataFrame(quality,columns=['program','donor','block','original_detection','matched_detection',
                    'delta_detection','original_log_abundance','matched_log_abundance','delta_log_abundance'])
    quality_frame.to_csv(bulk/'accepted_matching_quality.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
    quality_rows=[]
    if len(quality_frame):
        for (donor,block),group in quality_frame.groupby(['donor','block'],sort=True):
            quality_rows.append({'donor':donor,'block':block,'accepted_programs':len(group),
                                 'original_detection':group.original_detection.iloc[0],
                                 'original_log_abundance':group.original_log_abundance.iloc[0],
                                 'mean_detection_mismatch':group.delta_detection.mean(),
                                 'maximum_absolute_detection_mismatch':group.delta_detection.abs().max(),
                                 'mean_log_abundance_mismatch':group.delta_log_abundance.mean(),
                                 'maximum_absolute_log_abundance_mismatch':group.delta_log_abundance.abs().max()})
    pd.DataFrame(quality_rows,columns=['donor','block','accepted_programs','original_detection','original_log_abundance',
                 'mean_detection_mismatch','maximum_absolute_detection_mismatch','mean_log_abundance_mismatch',
                 'maximum_absolute_log_abundance_mismatch']).to_csv(OUT/'matched_baseline_quality.tsv',sep='\t',index=False)
    freeze={'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scientific_parameters_sha256':parameters_sha256(),
            'original_score_sha256':sha(evidence_path('frozen_score.json')),'matching_design':LOCKED_MATCHING,
            'matching_status':matching,'treatment_samples_read_by_this_analysis':False,'historical_treatment_effects_already_known':True,
            'vehicle_input_sha256':dict(inputs),
            'candidate_and_program_sha256':{path.relative_to(data_path('revision_02')).as_posix():sha(path)
                    for path in sorted(bulk.glob('*')) if path.is_file()},
            'code_sha256':{name:sha(ROOT/'code'/name) for name in ['revision02_gene_specificity.py','revision02_gene_math.py']}}
    write_json(OUT/'gene_program_freeze.json',freeze)
    print('Baseline-only program freeze:',matching,flush=True)
    # A failed match is retained. Gene interpretation can still finish, but no
    # partial-set treatment percentile or relaxed matching is produced.
    for row in manifest[manifest.treatment_code=='H_GLP'].sort_values('donor').to_dict('records'):
        summary=sample_summary(row,universe,inputs);summaries[(row['donor'],'H_GLP')]=summary
        np.savez_compressed(bulk/(row['accession']+'_gene_summaries.npz'),**summary)
        sample_rows.append(gene_rows(summary,row,universe,original))
        print('Treated gene summary complete:',row['accession'],flush=True)
    samples=pd.concat(sample_rows,ignore_index=True)
    contrasts,genes,concentration=gene_tables(samples)
    samples.to_csv(OUT/'gene_sample_summaries.tsv',sep='\t',index=False)
    contrasts.to_csv(OUT/'gene_donor_contrasts.tsv',sep='\t',index=False)
    genes.to_csv(OUT/'gene_contribution_summary.tsv',sep='\t',index=False)
    concentration.to_csv(OUT/'gene_concentration.tsv',sep='\t',index=False)
    selected=np.vstack([original,programs]) if matching['status']=='COMPLETE' else original[None,:]
    pairs=[donor_contrasts(summaries[(d,'D')],summaries[(d,'H_GLP')],selected) for d in DONORS]
    if matching['status']=='COMPLETE':
        comparison,draws=summarize_matched(pairs)
        comparison.to_csv(OUT/'matched_program_comparison.tsv',sep='\t',index=False)
        draws.to_csv(bulk/'matched_program_contrasts.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
    else:
        pd.DataFrame([{'status':matching['status'],'specificity_percentile_reported':False,
                       'reason':'Fixed baseline matching requirements were not fulfilled; no relaxation or partial-program inference.'}]).to_csv(OUT/'matched_program_comparison.tsv',sep='\t',index=False)
    previous=pd.read_csv(scientific_file('evidence/revision_01/observed_components.tsv'),sep='\t',float_precision='round_trip')
    checks=[]
    for (donor,treatment),summary in summaries.items():
        native=int(program_numerators(summary,original[None,:])['score'][0])/(100*int(summary['G'])*int(summary['n']))
        expected=previous[(previous.donor==donor)&(previous.treatment==treatment)&(previous.population=='alpha_broad')&(previous.component=='score')]['mean'].iloc[0]
        error=abs(native-expected)
        if error>1e-13:raise AssertionError('Original frozen sample mean did not reproduce')
        checks.append({'donor':donor,'treatment':treatment,'original_mean_error':error})
    for donor,pair in zip(DONORS,pairs):
        group=contrasts[contrasts.donor==donor]
        exact=int(pair['score']['numerators'][0])/int(pair['score']['denominator'])
        if abs((group.delta_signed_occupancy+group.delta_signed_residual).sum()-exact)>1e-13:
            raise AssertionError('Signed gene components failed to sum to original donor score contrast')
    status={'status':'COMPLETE' if matching['status']=='COMPLETE' else 'GENE_ANALYSIS_COMPLETE_MATCHED_CONTROL_FAILED',
            'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scientific_parameters_sha256':parameters_sha256(),
            'program_freeze_sha256':sha(OUT/'gene_program_freeze.json'),
            'matching':matching,'all_original_genes_reported':100,'paired_donors':DONORS,'population':'fixed broad-alpha expression proxy',
            'sample_mean_reproduction':checks,'source_and_membership_sha256':inputs,
            'compact_output_sha256':{name:sha(OUT/name) for name in owned if (OUT/name).exists()},
            'bulk_output_sha256':{p.relative_to(data_path('revision_02')).as_posix():sha(p) for p in sorted(bulk.glob('*')) if p.is_file()},
            'limits':['Conditional-positive cells differ across arms; no same-cell retention or lineage claim.',
                      'Occupancy coefficient 1/100 applies because both fixed blocks contain exactly50 genes.',
                      'Gene residuals use centered positive ranks; no arbitrary allocation of changing zero baseline to individual genes.',
                      'Matched programs are descriptive controls, not treatment-randomization nulls or calibrated p-values.',
                      'Matching marginal baseline detection/abundance does not reproduce gene correlations or original source-selection history.',
                      'Historical treatment effects were known; candidate/program freeze precedes only the new program treatment outcomes.',
                      'Published curated cell labels and corrected-layer identity remain unresolved.']}
    write_json(OUT/'gene_specificity_status.json',status)
    print(json.dumps({k:v for k,v in status.items() if k not in ['source_and_membership_sha256','compact_output_sha256','bulk_output_sha256']},indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    from public_runtime import execute_stage
    execute_stage('genes', args.data_root, args.output_dir)
