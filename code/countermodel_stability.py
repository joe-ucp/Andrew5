"""Numerical stability of declared nuisance scenarios with regenerated profile pools.

This is a post-calibration Monte Carlo verification extension. It does not refit
the frozen score, estimate actual contamination, or validate biological ancestry.
Ten independent 400-cell pools and twenty group draws per pool retain variation
that conditioning all draws on one synthetic pool would omit. The zero-RNA null
shares its exact pool with the baseline so its conditional expected effect is zero.
"""
from data_paths import evidence_path, prepare_original_results, logical_evidence_key
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
from primary_analysis import (ROOT,EV,map_universe,prob_profiles,simulate_scores)
from audit_primary import load_sample

BATCHES=10
POOL_N=400
GROUP_DRAWS=20
FRACTIONS=[0,.01,.05,.10,.20,.50]
DOUBLET_FRACTIONS=[.01,.05,.10]
CYCLE_WEIGHTS=[0,.25,.5,.75,1]

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    prepare_original_results()
    calibration=json.loads((EV/'calibration_status.json').read_text())
    if calibration['technical_calibration_pass']:
        raise RuntimeError('This extension is scoped to gate-failed diagnostics; re-review required.')
    config=json.loads((EV/'frozen_score.json').read_text())
    universe=config['universe']; index={g:i for i,g in enumerate(universe)}
    weights=np.zeros(len(universe))
    weights[[index[g] for g in config['up']]]=1/len(config['up'])
    weights[[index[g] for g in config['down']]]=-1/len(config['down'])
    original=EV/'countermodel_results.tsv'
    original_hash=sha(original)
    manifest=pd.read_csv(ROOT/'primary_sample_manifest.tsv',sep='\t')
    cell_scores=pd.read_csv(evidence_path('primary_cell_scores.tsv.gz'),sep='\t',float_precision='round_trip')
    primary=pd.read_csv(EV/'primary_results.tsv',sep='\t',float_precision='round_trip')
    targets=primary[(primary.contrast=='H_GLP-D')&(primary.annotation=='alpha_related')&(primary.sensitivity=='primary')].set_index('donor')
    batches=[]; group_draws=[]; pools=[]; donor_checks=[]
    for _,sample in manifest[manifest.treatment_code=='D'].sort_values('donor').iterrows():
        donor=sample.donor; donor_number=int(donor.removeprefix('D'))
        x,features,barcodes=load_sample(sample.to_dict())
        saved=cell_scores[cell_scores.accession==sample.accession].set_index('barcode')
        ix=np.flatnonzero(barcodes.barcode.isin(saved.index).to_numpy())
        aligned=saved.loc[barcodes.barcode.iloc[ix]]
        if len(aligned)!=len(saved):raise ValueError('Retained barcode mismatch')
        xu=map_universe(x[:,ix],features,universe)
        profiles=prob_profiles(xu,aligned.library_counts.to_numpy())
        am=aligned.annotation.eq('alpha_related').to_numpy()
        bm=aligned.annotation.eq('beta_like').to_numpy()
        alpha=profiles[:,am]; beta=profiles[:,bm]
        depths=aligned.library_counts.to_numpy()[am]
        cycle=aligned.cycle_fraction.to_numpy()[am]
        quartiles=np.quantile(cycle,[.25,.75])
        low=cycle<=quartiles[0]
        high=(cycle>=quartiles[1])&(cycle>0)&~low
        cycle_valid=bool((cycle>0).sum()>=20 and low.sum()>=20 and high.sum()>=20)
        target=targets.loc[donor];n0=int(target.n_control);n1=int(target.n_treated)
        if len(depths)!=n0:raise ValueError('Vehicle group-count mismatch')
        donor_checks.append(dict(donor=donor,vehicle_alpha_n=n0,treated_alpha_n=n1,vehicle_beta_n=int(bm.sum()),cycle_positive_n=int((cycle>0).sum()),cycle_low_n=int(low.sum()),cycle_high_n=int(high.sum()),cycle_valid=cycle_valid,depth_distribution='vehicle alpha',seed_rule='7400 + numeric_donor*100 + batch_index; batch_index0..9'))
        del x,xu,profiles
        for batch in range(BATCHES):
            seed=7400+donor_number*100+batch
            rng=np.random.default_rng(seed)
            baseline=simulate_scores(alpha,depths,weights,rng,n=POOL_N)
            scenarios={'beta_RNA_0':baseline.copy()}
            for fraction in FRACTIONS[1:]:
                scenarios[f'beta_RNA_{fraction:g}']=simulate_scores(alpha,depths,weights,rng,n=POOL_N,beta=beta,frac=fraction)
            for fraction in DOUBLET_FRACTIONS:
                values=baseline.copy()
                selected=rng.choice(POOL_N,int(round(fraction*POOL_N)),replace=False)
                values[selected]=rng.choice(scenarios['beta_RNA_0.5'],len(selected))
                scenarios[f'doublet_fraction_{fraction:g}']=values
            if cycle_valid:
                lowpool=simulate_scores(alpha[:,low],depths,weights,rng,n=POOL_N)
                highpool=simulate_scores(alpha[:,high],depths,weights,rng,n=POOL_N)
                for fraction in CYCLE_WEIGHTS:
                    n_high=int(POOL_N*fraction)
                    scenarios[f'cycle_high_weight_{fraction:g}']=np.r_[rng.choice(lowpool,POOL_N-n_high),rng.choice(highpool,n_high)]
            for j,score in enumerate(baseline):
                pools.append(dict(donor=donor,batch=batch,seed=seed,model='baseline',synthetic_cell=j,score=float(score)))
            for model,values in scenarios.items():
                conditional_effect=float(values.mean()-baseline.mean())
                if model=='beta_RNA_0' and conditional_effect!=0:
                    raise AssertionError('Zero-contamination model must share baseline pool')
                draws=np.array([rng.choice(values,n1).mean()-rng.choice(baseline,n0).mean() for _ in range(GROUP_DRAWS)])
                batches.append(dict(donor=donor,contrast='H_GLP-D',batch=batch,seed=seed,model=model,pool_n=POOL_N,group_draws=GROUP_DRAWS,actual_n_control=n0,actual_n_treated=n1,pool_baseline_mean=float(baseline.mean()),pool_scenario_mean=float(values.mean()),conditional_pool_mean_effect=conditional_effect,group_mean_effect=float(draws.mean()),group_sd=float(draws.std(ddof=1)),conditional_group_q95=float(np.quantile(draws,.95)),observed_delta=float(target.delta),cycle_calibrated=cycle_valid))
                for j,value in enumerate(draws):group_draws.append(dict(donor=donor,batch=batch,seed=seed,model=model,draw=j,group_mean_delta=float(value)))
                for j,value in enumerate(values):pools.append(dict(donor=donor,batch=batch,seed=seed,model=model,synthetic_cell=j,score=float(value)))
            print(f'{donor} batch {batch+1}/{BATCHES}: {len(scenarios)} fixed scenarios',flush=True)
        del alpha,beta
    batch_table=pd.DataFrame(batches);draw_table=pd.DataFrame(group_draws)
    records=[]
    for (donor,model),part in batch_table.groupby(['donor','model'],sort=True):
        draw=draw_table[(draw_table.donor==donor)&(draw_table.model==model)].group_mean_delta.to_numpy()
        effects=part.conditional_pool_mean_effect.to_numpy()
        observed=float(targets.loc[donor].delta)
        pooled_q95=float(np.quantile(draw,.95));max_batch_q95=float(part.conditional_group_q95.max())
        records.append(dict(donor=donor,contrast='H_GLP-D',model=model,independent_pool_batches=len(part),synthetic_cells_per_pool=POOL_N,total_group_draws=len(draw),conditional_effect_mean=float(effects.mean()),conditional_effect_min=float(effects.min()),conditional_effect_max=float(effects.max()),conditional_effect_MC_SD=float(effects.std(ddof=1)),group_delta_mean=float(draw.mean()),group_delta_MC_SD=float(draw.std(ddof=1)),group_delta_min=float(draw.min()),group_delta_max=float(draw.max()),pooled_MC_q05=float(np.quantile(draw,.05)),pooled_MC_q95=pooled_q95,maximum_batch_conditional_q95=max_batch_q95,observed_delta=observed,exceeds_pooled_MC_q95=bool(observed>pooled_q95),exceeds_maximum_batch_q95=bool(observed>max_batch_q95),interpretation='Monte Carlo stability envelope; conditional model sensitivity, not empirical contamination or confidence bound'))
    summary=pd.DataFrame(records)
    batch_table.to_csv(evidence_path('countermodel_stability_batches.tsv'),sep='\t',index=False)
    draw_table.to_csv(evidence_path('countermodel_stability_group_draws.tsv'),sep='\t',index=False)
    pd.DataFrame(pools).to_csv(evidence_path('countermodel_stability_pools.tsv.gz'),sep='\t',index=False,compression={'method':'gzip','mtime':0})
    summary.to_csv(EV/'countermodel_stability_summary.tsv',sep='\t',index=False)
    pd.DataFrame(donor_checks).to_csv(EV/'countermodel_stability_donor_checks.tsv',sep='\t',index=False)
    if sha(original)!=original_hash:raise AssertionError('Provisional original was modified')
    status=dict(status='EXECUTED_DIAGNOSTIC_MONTE_CARLO_STABILITY',original_countermodel_status='PROVISIONAL_SINGLE_POOL; retained unchanged',original_countermodel_sha256=original_hash,scope='H_GLP-D; fixed observations and score; ten independently regenerated pools; no measured contamination/biological CI',batches=BATCHES,pool_n=POOL_N,group_draws_per_batch=GROUP_DRAWS,summary_rows=len(summary),batch_rows=len(batch_table),group_draw_rows=len(draw_table),pool_score_rows=len(pools),zero_contamination_conditional_effect_max_abs=float(summary[summary.model=='beta_RNA_0'].conditional_effect_max.abs().max()),source_hashes={(p.relative_to(ROOT).as_posix() if p.is_relative_to(ROOT) else logical_evidence_key(p)):sha(p) for p in [EV/'frozen_score.json',EV/'primary_results.tsv',evidence_path('primary_cell_scores.tsv.gz'),ROOT/'code/countermodel_stability.py',ROOT/'code/primary_analysis.py']},policy='All declared models and batches retained; q95 pooled across all200 draws and realized extrema are numerical sensitivity diagnostics, not guaranteed conservative statistical coverage.')
    (EV/'countermodel_stability_status.json').write_text(json.dumps(status,indent=2)+'\n')
    print(json.dumps(status,indent=2))

if __name__ == '__main__':
    import argparse
    from pathlib import Path
    from public_runtime import execute_stage
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    execute_stage('countermodels', args.data_root, args.output_dir)
