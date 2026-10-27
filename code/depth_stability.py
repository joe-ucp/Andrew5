"""After-outcome numerical stability check of fixed 2000-count resampling.

Exactly 20 independent observation-model realizations; fixed original alpha cell
sets and frozen score. This does not fit a score, change the hypothesis, or add
biological replication. Uses the original primary-analysis implementation.
"""
from data_paths import evidence_path, prepare_original_results
from pathlib import Path
import json
import hashlib
import datetime
import numpy as np
import pandas as pd
from scipy import sparse
from primary_analysis import ROOT, EV, load_sample, map_universe, score_sparse, prob_profiles


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resample_fixed_cells(profiles, weights, rng, depth=2000):
    """Multinomial on nonzero support, exactly equivalent to the full vector law."""
    n = profiles.shape[1]
    outside = profiles.shape[0]-1
    values = []
    indices = []
    pointers = [0]
    for cell in range(n):
        start, stop = profiles.indptr[cell:cell+2]
        idx = profiles.indices[start:stop]
        p = np.maximum(profiles.data[start:stop], 0)
        p /= p.sum()
        counts = rng.multinomial(depth, p)
        take = (counts > 0) & (idx < outside)
        values.append(counts[take])
        indices.append(idx[take])
        pointers.append(pointers[-1]+int(take.sum()))
    sampled = sparse.csc_matrix((np.concatenate(values), np.concatenate(indices),
                                 np.asarray(pointers)), shape=(outside, n))
    return score_sparse(sampled, weights)


def main():
    prepare_original_results()
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    config_path = EV/'frozen_score.json'
    source_path = evidence_path('primary_cell_scores.tsv.gz')
    cfg = json.loads(config_path.read_text())
    universe = cfg['universe']
    feature_ix = {g:i for i,g in enumerate(universe)}
    weights = np.zeros(len(universe))
    weights[[feature_ix[g] for g in cfg['up']]] = 1/len(cfg['up'])
    weights[[feature_ix[g] for g in cfg['down']]] = -1/len(cfg['down'])
    cells = pd.read_csv(source_path, sep='\t', float_precision='round_trip')
    fixed = cells[(cells.annotation=='alpha_related') & cells.treatment.isin(['D','H_GLP'])].copy()
    if fixed.duplicated(['accession','barcode']).any():
        raise ValueError('Duplicate primary cell keys.')
    manifest = pd.read_csv(ROOT/'primary_sample_manifest.tsv',sep='\t')
    manifest = manifest[manifest.treatment_code.isin(['D','H_GLP'])].sort_values(['donor','treatment_code'])
    if len(manifest)!=8 or manifest.donor.nunique()!=4:
        raise ValueError('Expected exactly four paired donors.')
    seeds = list(range(7200,7220))
    generators = [np.random.default_rng(seed) for seed in seeds]
    records = []
    source_records = []
    for row in manifest.to_dict('records'):
        sub = fixed[fixed.accession==row['accession']]
        x, feat, barcodes = load_sample(row)
        barcode_ix = pd.Index(barcodes.barcode)
        if not barcode_ix.is_unique:
            raise ValueError('Nonunique submitted barcode.')
        indices = barcode_ix.get_indexer(sub.barcode)
        if (indices<0).any():
            raise ValueError('Primary cell key missing in source matrix.')
        selected = x[:,indices].tocsc()
        totals = np.asarray(selected.sum(axis=0)).ravel()
        if not np.array_equal(totals,sub.library_counts.to_numpy()):
            raise ValueError('Original cell depths differ from source matrix.')
        mapped = map_universe(selected,feat,universe)
        native = score_sparse(mapped,weights)
        if not np.allclose(native,sub.score.to_numpy(),rtol=0,atol=1e-14):
            raise ValueError('Original frozen scores do not reproduce.')
        profiles = prob_profiles(mapped,totals)
        for rep,(seed,rng) in enumerate(zip(seeds,generators)):
            score = resample_fixed_cells(profiles,weights,rng)
            records.append({'donor':row['donor'],'treatment':row['treatment_code'],
                            'accession':row['accession'],'replicate':rep,'seed':seed,
                            'depth':2000,'fixed_cells':len(sub),'mean_score':float(score.mean())})
        source_records.append({'accession':row['accession'],'fixed_alpha_cells':len(sub),
                               'native_score_max_abs_reproduction_error':float(np.max(np.abs(native-sub.score.to_numpy())))})
        print(row['accession'],row['donor'],row['treatment_code'],len(sub),'cells; 20 repeats complete',flush=True)
        del x,selected,mapped,profiles
    samples = pd.DataFrame(records)
    wide = samples.pivot(index=['donor','replicate','seed','depth'],columns='treatment',values='mean_score').reset_index()
    wide['delta'] = wide['H_GLP']-wide['D']
    wide['positive'] = wide.delta>0
    wide['interpretation'] = 'after_outcome_MC_stability_fixed_cells_not_biological_replication'
    summaries = []
    for donor,subset in wide.groupby('donor'):
        summaries.append({'donor':donor,'depth':2000,'MC_replicates':len(subset),
                          'mean_delta':float(subset.delta.mean()),'SD_delta':float(subset.delta.std(ddof=1)),
                          'minimum_delta':float(subset.delta.min()),'maximum_delta':float(subset.delta.max()),
                          'positive_fraction':float(subset.positive.mean()),
                          'fixed_control_cells':int(samples[(samples.donor==donor)&(samples.treatment=='D')].fixed_cells.iloc[0]),
                          'fixed_treated_cells':int(samples[(samples.donor==donor)&(samples.treatment=='H_GLP')].fixed_cells.iloc[0]),
                          'interpretation':'MC_variability_conditional_on_original_cells'})
    samples.to_csv(evidence_path('depth_stability_sample_means.tsv'),sep='\t',index=False)
    wide.to_csv(evidence_path('depth_stability_replicates.tsv'),sep='\t',index=False)
    summary = pd.DataFrame(summaries)
    summary.to_csv(EV/'depth_stability_summary.tsv',sep='\t',index=False)
    status = {'status':'COMPLETE','started_utc':started,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'design':'After-outcome numerical stability verification; no new hypothesis or endpoint selection',
              'seeds':seeds,'depth':2000,'replicates':20,'cells_reselected':False,
              'sampling_law':'Independent multinomial at fixed depth using empirical whole-library probabilities; zero-probability categories omitted computationally; outside-universe mass retained',
              'primary_cell_scores_sha256':sha(source_path),'frozen_score_sha256':sha(config_path),
              'sample_manifest_sha256':sha(ROOT/'primary_sample_manifest.tsv'),
              'limitations':['Monte Carlo repeats are not independent donors or new biological observations',
                             'Empirical read resampling does not model unseen molecules, capture differences, or systematic library effects',
                             'Samples share fixed observed profiles; interpretation is conditional numerical stability only'],
              'sources':source_records,'summary':summaries}
    (EV/'depth_stability_status.json').write_text(json.dumps(status,indent=2)+'\n',encoding='utf-8')
    print(summary.to_string(index=False),flush=True)


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    from public_runtime import execute_stage
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    execute_stage('depth', args.data_root, args.output_dir)
