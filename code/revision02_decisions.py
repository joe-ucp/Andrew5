"""Frozen source-threshold transport and paired measurement diagnostics (Revision 02).

Full per-cell/per-draw data stay at the external data root. Source thresholds
are saved before target classifications. Reuses exactly Revision 01 draw seeds.
"""
import argparse
import datetime
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from data_paths import data_path, resolve_data_root
from public_runtime import scientific_file, verify_registered_data, parameters_sha256
from primary_analysis import load_sample, map_universe
from reference_analysis import references
from revision_score import components_sparse, thin_sparse_counts
from revision02_math import METHODS, score_bundle, zero_centered_components, q95, decision_metrics, state_auc
from revision02_common import ROOT, EV, OLD, load_protocol, sha, external, write_json


def verify_old(relative):
    return verify_registered_data(relative)


def freeze_source(conf):
    source = data_path('data_reference/GSE123844_RAW.tar', must_exist=True)
    source_registry = pd.read_csv(ROOT/'source_manifest.tsv', sep='\t')
    match = source_registry[source_registry.local_path.str.replace('\\','/',regex=False) == 'data_reference/GSE123844_RAW.tar']
    if len(match) != 1 or match.sha256.iloc[0] != sha(source):
        raise ValueError('Reference source identity failed')
    _, preparations = references()
    programs = {'original': conf}
    values = {}
    for label, frame in preparations.items():
        x = sparse.csc_matrix(frame.reindex(conf['universe'], fill_value=0).to_numpy())
        values[label] = score_bundle(x, frame.sum(axis=0).to_numpy(), conf['universe'], programs)[0]
    thresholds = {method: q95(values['alpha'][i]) for i, method in enumerate(METHODS)}
    freeze = {'scientific_parameters_sha256':parameters_sha256(), 'source_sha256':sha(source),
              'program_sha256':sha(scientific_file('evidence/frozen_score.json')), 'thresholds':thresholds,
              'quantile':'numpy linear q95', 'comparison':'strict score > threshold',
              'alpha_n':values['alpha'].shape[1], 'source_threshold_saved_before_target_classification_in_this_execution':True,'historical_target_outcomes_already_known':True}
    write_json(EV/'source_threshold_freeze.json', freeze)
    rows = []
    for label, vals in values.items():
        for mi, method in enumerate(METHODS):
            rows.append({'preparation':label,'method':method,
                         **decision_metrics(vals[mi], thresholds[method])})
    pd.DataFrame(rows).to_csv(EV/'source_operating_points.tsv', sep='\t', index=False)
    return thresholds


def target_rows(donor, cohort, depth, replicate, seed, values, native, comps, native_comps,
                types, population_types, source_thresholds, native_thresholds, depth_thresholds):
    rows = []
    for mi, method in enumerate(METHODS):
        auroc, native_auc = state_auc(values[mi], types), state_auc(native[mi], types)
        rules = {'source_fixed': (source_thresholds[method], source_thresholds[method]),
                 'target_native_fixed': (native_thresholds[method], native_thresholds[method])}
        if depth:
            rules['target_depth_recalibrated'] = (depth_thresholds[method], native_thresholds[method])
        for rule, (threshold, prior_threshold) in rules.items():
            for label in ['alpha','beta']:
                use = types == label
                n = int(use.sum()); valid = n >= 20
                metrics = decision_metrics(values[mi,use],threshold,native[mi,use],prior_threshold)
                native_metrics = decision_metrics(native[mi,use],prior_threshold)
                base = {'donor':donor,'cohort':cohort,'depth':depth,'replicate':replicate,'seed':seed,
                        'method':method,'rule':rule,'cell_type':label,'population_n':int((population_types==label).sum()),
                        'retained_fraction':n/int((population_types==label).sum()),
                        'native_threshold':prior_threshold,'native_positive_fraction':native_metrics['positive_fraction'],
                        'native_tie_n':native_metrics['tie_n'],
                        'AUROC':auroc,'native_same_eligible_AUROC':native_auc,
                        'score_mean':float(values[mi,use].mean()) if valid else np.nan,
                        'native_score_mean':float(native[mi,use].mean()) if valid else np.nan,
                        'ucell_zero_floor_fraction':float((values[2,use]==0).mean()) if valid else np.nan,
                        'native_ucell_zero_floor_fraction':float((native[2,use]==0).mean()) if valid else np.nan,
                        'decision_scope':'author-labelled mature state, not cell conversion'}
                for key in ['k_up','k_down','occupancy_component','positive_rank_residual','zero_rank_baseline']:
                    base[key+'_mean'] = float(comps[key][use].mean()) if valid else np.nan
                    base['native_'+key+'_mean'] = float(native_comps[key][use].mean()) if valid else np.nan
                rows.append({**base, **metrics})
    return rows


def transport():
    plan, conf = load_protocol(); spec = plan['threshold_transport']
    if (EV/'threshold_status.json').exists() or (EV/'source_threshold_freeze.json').exists():
        raise FileExistsError('No implicit rerun or overwrite of threshold outcomes')
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    source_thresholds = freeze_source(conf)
    rows, checks, target_threshold_rows, native_states = [], [], [], []
    genes = conf['universe']; lookup = {g:i for i,g in enumerate(genes)}
    up, down = [lookup[g] for g in conf['up']], [lookup[g] for g in conf['down']]
    registry = json.loads((scientific_file('evidence/revision_01/baron_input_manifest.json')).read_text())
    for donor, record in enumerate(sorted(registry,key=lambda r:r['name']),1):
        path = data_path('revision_01/sources/'+record['name'],must_exist=True)
        if sha(path) != record['sha256']: raise ValueError('Independent source hash mismatch')
        frame = pd.read_csv(path,index_col=0)
        if not frame.index.is_unique: raise ValueError('Nonunique original cell IDs')
        count = frame.iloc[:,2:]; raw = count.to_numpy()
        if not np.isfinite(raw).all() or (raw<0).any() or (raw!=np.floor(raw)).any():
            raise ValueError('Nonnegative integer counts required')
        use = frame.assigned_cluster.isin(['alpha','beta']).to_numpy()
        types = frame.assigned_cluster.to_numpy(dtype=str)[use]
        ids = frame.index.to_numpy(dtype=str)[use]
        totals = raw.sum(1)[use].astype(np.int64)
        x = sparse.csc_matrix(count.loc[use].reindex(columns=genes,fill_value=0).to_numpy().T)
        native = score_bundle(x,totals,genes,{'original':conf})[0]
        comp = zero_centered_components(components_sparse(x,up,down),len(genes))
        old = np.load(verify_old(f'revision_01/independent/donor{donor}_original_universe_scores.npz'))
        # Read numerical arrays only: historical state labels were object dtype.
        error = float(np.max(np.abs(native[:3]-old['scores'][0])))
        if error > 1e-15: raise AssertionError('Native independent score reproduction failed')
        native_thresholds = {method:q95(native[i,types=='alpha']) for i,method in enumerate(METHODS)}
        target_threshold_rows.extend({'donor':donor,'method':method,'native_alpha_q95':t,'native_alpha_n':int((types=='alpha').sum())} for method,t in native_thresholds.items())
        write_json(external(f'threshold_transport/donor{donor}_native_thresholds.json'),
                   {'donor':donor,'thresholds':native_thresholds,'rule':'all native author-labelled alpha; reused unchanged across all depths'})
        np.savez_compressed(external(f'threshold_transport/donor{donor}_native.npz'),
                            cell_id=ids,cell_type=types,total_counts=totals,scores=native,
                            methods=np.asarray(METHODS),**comp)
        native_states.append(pd.DataFrame({'donor':donor,'cell_id':ids,'cell_type':types,'total_counts':totals,**comp}))
        rows.extend(target_rows(donor,'all_native',0,-1,-1,native,native,comp,comp,types,types,
                                source_thresholds,native_thresholds,{}))
        for depth in spec['depths']:
            eligible = totals >= depth
            subtypes = types[eligible]; subnative = native[:,eligible]
            subcomp = {k:v[eligible] for k,v in comp.items()}
            for replicate,seed in enumerate(spec['seeds']):
                rng = np.random.default_rng(np.random.SeedSequence([seed,donor,depth]))
                thin = thin_sparse_counts(x[:,eligible],totals[eligible],depth,rng)
                vals = score_bundle(thin,np.full(eligible.sum(),depth),genes,{'original':conf})[0]
                parts = zero_centered_components(components_sparse(thin,up,down),len(genes))
                recal = {method:q95(vals[i,subtypes=='alpha']) for i,method in enumerate(METHODS)}
                np.savez_compressed(external(f'threshold_transport/donor{donor}_depth{depth}_seed{seed}.npz'),
                                    cell_id=ids[eligible],cell_type=subtypes,original_total_counts=totals[eligible],
                                    scores=vals,methods=np.asarray(METHODS),**parts)
                for cohort,mask in [('depth_eligible',np.ones(eligible.sum(),dtype=bool)),
                                    ('fixed_10000',totals[eligible]>=10000)]:
                    rows.extend(target_rows(donor,cohort,depth,replicate,seed,vals[:,mask],subnative[:,mask],
                                {k:v[mask] for k,v in parts.items()},{k:v[mask] for k,v in subcomp.items()},
                                subtypes[mask],types,source_thresholds,native_thresholds,recal))
            print('Threshold transport donor',donor,'depth',depth,'complete',flush=True)
        checks.append({'donor':donor,'source_sha256':sha(path),'native_score_maximum_error':error})
    draws = pd.DataFrame(rows)
    draws.to_csv(external('threshold_transport/decision_draws.tsv.gz'),sep='\t',index=False)
    pd.concat(native_states,ignore_index=True).to_csv(external('threshold_transport/native_components.tsv.gz'),sep='\t',index=False)
    pd.DataFrame(target_threshold_rows).to_csv(EV/'target_native_thresholds.tsv',sep='\t',index=False)
    keys = ['donor','cohort','depth','method','rule','cell_type']
    numeric = [c for c in draws.select_dtypes(include=np.number).columns if c not in keys+['replicate','seed']]
    summary = draws.groupby(keys,sort=True)[numeric].mean().reset_index()
    for metric in ['positive_fraction','negative_to_positive_fraction','positive_to_negative_fraction','score_mean','AUROC']:
        summary[metric+'_technical_SD'] = draws.groupby(keys,sort=True)[metric].std().to_numpy()
    summary['draws'] = draws.groupby(keys,sort=True).size().to_numpy()
    summary['evaluable'] = draws.groupby(keys,sort=True).evaluable.all().to_numpy()
    summary['suppression_reason'] = draws.groupby(keys,sort=True).reason.agg(lambda s:';'.join(sorted(set(s)-{''}))).to_numpy()
    summary['AUROC_reason'] = np.where(summary.AUROC.notna(),'','fewer_than_20_cells_in_at_least_one_state')
    summary.to_csv(EV/'threshold_summary.tsv',sep='\t',index=False)
    # Same seeds and molecules: verify the original three method means against Rev01.
    previous = pd.read_csv(scientific_file('evidence/revision_01/independent_thinning.tsv'),sep='\t',float_precision='round_trip')
    previous = previous[(previous['mode']=='original_universe')&(previous.program=='original')]
    check = summary[(summary.cohort=='depth_eligible')&(summary.rule=='source_fixed')&(summary.method.isin(METHODS[:3]))]
    joined = previous.merge(check,on=['donor','depth','method','cell_type'],validate='one_to_one')
    error = float(np.nanmax(np.abs(joined.thinned_mean-joined.score_mean)))
    if len(joined)!=len(previous) or error>1e-14:
        raise AssertionError('Same-seed thinning did not reproduce prior method means')
    write_json(EV/'threshold_status.json',{'status':'COMPLETE','started_utc':started,
               'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'scientific_parameters_sha256':parameters_sha256(),'source_threshold_freeze_sha256':sha(EV/'source_threshold_freeze.json'),
               'source_checks':checks,'prior_thinning_comparisons':len(joined),'prior_thinning_maximum_error':error,
               'draw_summary_rows':len(draws),'code_sha256':{n:sha(ROOT/'code'/n) for n in ['revision02_decisions.py','revision02_math.py','revision_score.py']},
               'limits':['State labels do not label conversion','Depth recalibration is resubstitution on the evaluated alpha cells',
                         'Five draws quantify conditional sampling variation, not new biological replicates',
                         'Fixed high-depth cohort is selected; source threshold remains unchanged']})


def primary():
    plan, conf = load_protocol()
    if (EV/'primary_diagnostics_status.json').exists(): raise FileExistsError('Primary diagnostics already complete')
    manifest = pd.read_csv(ROOT/'primary_sample_manifest.tsv',sep='\t')
    manifest = manifest[manifest.treatment_code.isin(['D','H_GLP'])].sort_values(['donor','treatment_code'])
    baseline, fixed, floors, checks = [], [], [], []
    programs = json.loads((scientific_file('evidence/revision_01/benchmark_programs.json')).read_text())['programs']
    for row in manifest.to_dict('records'):
        accession = row['accession']; base = {'accession':accession,'donor':row['donor'],'treatment':row['treatment_code']}
        cells = pd.read_csv(verify_old(f'revision_01/{accession}_cell_components.tsv.gz'),sep='\t',float_precision='round_trip')
        cells = cells[cells.alpha_broad].copy()
        comp = zero_centered_components({k:cells[k].to_numpy() for k in ['score','k_positive','k_up','k_down',
                        'occupancy_up','occupancy_down','positive_rank_up','positive_rank_down',
                        'occupancy_component','positive_rank_residual','up_rank_contribution','down_rank_contribution']},len(conf['universe']))
        for k,v in comp.items(): baseline.append({**base,'component':k,'n':len(cells),'mean':float(v.mean())})
        same_ids = cells.loc[cells.library_counts>=10000,'source_cell_index'].to_numpy()
        for depth in plan['primary_fixed_cohort']['depths']:
            cache = np.load(verify_old(f'revision_01/{accession}_depth{depth}_thinning.npz'))
            use = cache['library_counts']>=10000
            if not np.array_equal(cache['source_cell_index'][use],same_ids):
                raise AssertionError('Fixed-cohort identity changed across depths')
            native = cells.loc[cells.library_counts>=10000]
            for replicate,seed in enumerate(cache['seeds']):
                fixed.append({**base,'depth':depth,'replicate':replicate,'seed':int(seed),'population_n':len(cells),
                              'eligible_n':len(native),'retained_fraction':len(native)/len(cells),
                              **{f'native_{k}':float(native[k].mean()) for k in ['score','k_up','k_down','occupancy_component','positive_rank_residual']},
                              **{f'thinned_{k}':float(cache[k][replicate,use].mean()) for k in ['score','k_up','k_down','occupancy_component','positive_rank_residual']}})
        # Reuse the original count cache with its frozen identities, and the
        # original broad population. This computes only the supplementary raw clip.
        for relative in [row['features_path'], row['barcodes_path'], 'data_primary/sparse_cache/'+accession+'.npz']:
            verify_registered_data(relative)
        x,features,barcodes = load_sample(row)
        positions = pd.Index(barcodes.barcode).get_indexer(cells.barcode)
        if (positions<0).any(): raise ValueError('Original cells missing')
        totals = np.asarray(x[:,positions].sum(0)).ravel()
        if not np.array_equal(totals,cells.library_counts): raise ValueError('Original cell depths changed')
        mapped = map_universe(x[:,positions],features,conf['universe'])
        values = score_bundle(mapped,totals,conf['universe'],programs)
        prior = np.load(verify_old(f'revision_01/benchmark/{accession}_scores.npz'))
        order = pd.Index(prior['barcodes']).get_indexer(cells.barcode)
        error = float(np.max(np.abs(values[:,:3]-prior['scores'][:,:,order])))
        if error>1e-15: raise AssertionError('Locked comparator reproduction failed')
        np.savez_compressed(external(f'primary_diagnostics/{accession}_scores.npz'),scores=values,
                            barcode=cells.barcode.to_numpy(dtype=str),total_counts=totals,
                            programs=np.asarray(list(programs)),methods=np.asarray(METHODS))
        for pi,program in enumerate(programs):
            floors.append({**base,'program':program,'n':len(cells),'floor_n':int((values[pi,2]==0).sum()),
                           'floor_fraction':float((values[pi,2]==0).mean()),'clipped_mean':float(values[pi,2].mean()),
                           'unclipped_mean':float(values[pi,3].mean())})
        checks.append({**base,'locked_score_maximum_error':error,'fixed_cohort_n':len(same_ids)})
        print('Primary baseline/floor/fixed-cohort',accession,'complete',flush=True)
    baseline = pd.DataFrame(baseline); fixed = pd.DataFrame(fixed); floors = pd.DataFrame(floors)
    baseline.to_csv(EV/'zero_baseline_samples.tsv',sep='\t',index=False)
    floors.to_csv(EV/'ucell_floor_samples.tsv',sep='\t',index=False)
    fixed.to_csv(external('primary_diagnostics/fixed_cohort_draws.tsv.gz'),sep='\t',index=False)
    for frame,name,keys,metrics in [
        (baseline,'zero_baseline_contrasts',['donor','component'],['mean']),
        (floors,'ucell_floor_contrasts',['donor','program'],['floor_fraction','clipped_mean','unclipped_mean']),
        (fixed,'fixed_cohort_contrasts',['donor','depth','replicate','seed'],
         [c for c in fixed if c.startswith(('native_','thinned_'))]+['eligible_n','retained_fraction'])]:
        control = frame[frame.treatment=='D'].set_index(keys)
        treated = frame[frame.treatment=='H_GLP'].set_index(keys).loc[control.index]
        result = control.index.to_frame(index=False)
        for metric in metrics:
            result['control_'+metric] = control[metric].to_numpy()
            result['treated_'+metric] = treated[metric].to_numpy()
            result['delta_'+metric] = treated[metric].to_numpy()-control[metric].to_numpy()
        if name=='fixed_cohort_contrasts':
            result.to_csv(external('primary_diagnostics/fixed_cohort_contrasts.tsv.gz'),sep='\t',index=False)
            summary = result.groupby(['donor','depth']).mean(numeric_only=True).drop(columns=['replicate','seed']).reset_index()
            summary['delta_thinned_score_technical_SD'] = result.groupby(['donor','depth']).delta_thinned_score.std().to_numpy()
            summary.to_csv(EV/'fixed_cohort_summary.tsv',sep='\t',index=False)
        else: result.to_csv(EV/(name+'.tsv'),sep='\t',index=False)
    write_json(EV/'primary_diagnostics_status.json',{'status':'COMPLETE','scientific_parameters_sha256':parameters_sha256(),
               'checks':checks,'code_sha256':{n:sha(ROOT/'code'/n) for n in ['revision02_decisions.py','revision02_math.py']},
               'fixed_cohort':'N>=10000 per arm; same source IDs at all four depths; five existing draws',
               'zero_baseline':'b=(G-K)/(2G) cancels in the complete signed score',
               'unclipped_diagnostic':'Supplementary explanation; locked clipped comparator unchanged'})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['transport', 'primary'])
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    from public_runtime import execute_stage
    execute_stage('transport' if args.stage == 'transport' else 'primary-diagnostics', args.data_root, args.output_dir)
