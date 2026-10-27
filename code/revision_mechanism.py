"""Revision 01: raw-count rank components and genuine finite-library thinning.

Execution uses the published scientific parameters and immutable score identity.
Scientific tables are new files in evidence/revision_01; per-cell material is
kept under the external data root's revision_01 directory. No score is fitted.
"""
from data_paths import evidence_path
from pathlib import Path
import argparse
import datetime
import gzip
import hashlib
import json
import os
import time
import numpy as np
import pandas as pd
from scipy import io
from scipy.stats import rankdata
from data_paths import data_path, resolve_data_root, ENV_NAME
from primary_analysis import map_universe, full_gene
from reference_analysis import CELL_CYCLE
from revision_score import components_sparse, thin_sparse_counts

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'evidence' / 'revision_01'
from public_runtime import protocol, frozen_score, parameters_sha256, scientific_file
HORMONES = ['INS', 'GCG', 'SST', 'PPY', 'GHRL']
POPULATIONS = ['alpha_broad', 'alpha_strict', 'alpha_cycle_high', 'beta',
               'mixed_endocrine', 'other_unassigned_endocrine']
COMPONENTS = ['score', 'up_rank_contribution', 'down_rank_contribution',
              'occupancy_component', 'occupancy_up', 'occupancy_down',
              'positive_rank_residual', 'positive_rank_up', 'positive_rank_down',
              'k_positive', 'k_up', 'k_down', 'residual_bound']


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_raw_sample(row):
    """Read deposited processed counts directly, preserving feature/barcode order."""
    feature_path = data_path(row['features_path'], must_exist=True)
    barcode_path = data_path(row['barcodes_path'], must_exist=True)
    matrix_path = data_path(row['matrix_path'], must_exist=True)
    features = pd.read_csv(feature_path, sep='\t', header=None,
                           names=['ensembl', 'symbol', 'feature_type'])
    barcodes = pd.read_csv(barcode_path, sep='\t', header=None, names=['barcode'])
    with gzip.open(matrix_path, 'rb') as handle:
        counts = io.mmread(handle).tocsr()
    counts.sum_duplicates()
    counts.sort_indices()
    if counts.shape != (len(features), len(barcodes)) or not barcodes.barcode.is_unique:
        raise ValueError('Source feature/barcode mapping is not one-to-one in ' + row['accession'])
    if not np.isfinite(counts.data).all() or (counts.data < 0).any() or (counts.data != np.floor(counts.data)).any():
        raise ValueError('Submitted matrix does not contain nonnegative integer counts.')
    return counts, features, barcodes


def select_populations(counts, features, barcodes, row, cycle_threshold=None):
    """Reproduce original QC and declare expression proxies without cell annotations."""
    total = np.asarray(counts.sum(axis=0)).ravel().astype(np.int64)
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    mitochondrial = np.asarray(counts[features.symbol.str.startswith('MT-').to_numpy()].sum(axis=0)).ravel()
    mt_fraction = mitochondrial / np.maximum(total, 1)
    qc = (total >= 1000) & (detected >= 400) & (mt_fraction < .2)
    markers = {gene: full_gene(counts, features, gene) for gene in HORMONES + ['ARX', 'IRX2', 'SMOC1']}
    alpha = qc & (markers['GCG'] >= 3)
    beta = qc & (markers['INS'] >= 3)
    for gene in HORMONES:
        if gene != 'GCG':
            alpha &= markers['GCG'] > markers[gene]
        if gene != 'INS':
            beta &= markers['INS'] > markers[gene]
    cycle_counts = np.asarray(counts[np.flatnonzero(features.symbol.isin(CELL_CYCLE))].sum(axis=0)).ravel()
    cycle_fraction = cycle_counts / np.maximum(total, 1)
    endocrine = qc & (sum(markers[gene] for gene in HORMONES) >= 3)
    strict = alpha & (markers['GCG'] >= 10) & ((markers['ARX'] + markers['IRX2']) > 0)
    if row['treatment_code'] == 'D':
        positive = cycle_fraction[alpha & (cycle_fraction > 0)]
        cycle_threshold = float(np.quantile(positive, .75, method='linear')) if len(positive) >= 20 else None
    high = alpha & (cycle_fraction > 0)
    high &= False if cycle_threshold is None else cycle_fraction >= cycle_threshold
    hormones = np.vstack([markers[gene] for gene in HORMONES])
    mixed = qc & ((hormones >= 3).sum(axis=0) >= 2) & (hormones.max(axis=0) <= .8 * hormones.sum(axis=0))
    masks = dict(alpha_broad=alpha, alpha_cycle_high=high, alpha_strict=strict,
                 mixed_endocrine=mixed, other_unassigned_endocrine=endocrine & ~alpha & ~beta, beta=beta)
    selected = np.flatnonzero(endocrine)
    frame = pd.DataFrame({'accession': row['accession'], 'donor': row['donor'],
                          'treatment': row['treatment_code'], 'barcode': barcodes.barcode.iloc[selected].to_numpy(),
                          'source_cell_index': selected, 'library_counts': total[selected],
                          'detected_genes': detected[selected], 'mt_fraction': mt_fraction[selected],
                          'cycle_fraction': cycle_fraction[selected]})
    for gene, values in markers.items():
        frame[gene] = values[selected]
    for population, mask in masks.items():
        frame[population] = mask[selected]
    audit = {'accession': row['accession'], 'donor': row['donor'], 'treatment': row['treatment_code'],
             'submitted_cells': len(total), 'original_QC_cells': int(qc.sum()),
             'endocrine_proxy_cells': int(endocrine.sum()),
             'vehicle_positive_cycle_q75': cycle_threshold,
             'cycle_threshold_status': 'evaluable' if cycle_threshold is not None else 'fewer_than_20_positive_vehicle_alpha_cycle_values',
             'vehicle_positive_alpha_cycle_cells': int((alpha & (cycle_fraction > 0)).sum()) if row['treatment_code'] == 'D' else None}
    return frame, selected, cycle_threshold, audit


def reproduce_original_cells(frame, historical):
    selected = frame[frame.alpha_broad | frame.beta]
    prior = historical[historical.accession == frame.accession.iloc[0]]
    if set(selected.barcode) != set(prior.barcode):
        raise ValueError('Original alpha/beta cell keys did not reproduce.')
    prior = prior.set_index('barcode').loc[selected.barcode]
    if not np.array_equal(selected.library_counts, prior.library_counts):
        raise ValueError('Original library totals did not reproduce.')
    expected = np.where(selected.alpha_broad, 'alpha_related', 'beta_like')
    if not np.array_equal(expected, prior.annotation):
        raise ValueError('Original alpha/beta expression-proxy membership changed.')
    error = float(np.max(np.abs(selected.score.to_numpy() - prior.score.to_numpy()))) if len(selected) else 0.
    if error > 1e-14:
        raise ValueError('Original frozen scores did not reproduce.')
    return {'historical_selected_cells': len(prior), 'maximum_original_score_error': error}


def correlation(x, y):
    x, y = np.asarray(x), np.asarray(y)
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def observed_tables(frame, audit):
    summaries, relationships, bins = [], [], []
    for population in POPULATIONS:
        sub = frame[frame[population]]
        base = {key: audit[key] for key in ['accession', 'donor', 'treatment']}
        base.update(population=population, n=len(sub), QC_denominator=audit['original_QC_cells'],
                    endocrine_denominator=audit['endocrine_proxy_cells'],
                    fraction_of_QC=len(sub) / audit['original_QC_cells'] if audit['original_QC_cells'] else np.nan,
                    status='descriptive_proxy' if len(sub) else 'empty_population')
        for component in COMPONENTS:
            summaries.append({**base, 'component': component,
                              'mean': float(sub[component].mean()) if len(sub) else np.nan,
                              'median': float(sub[component].median()) if len(sub) else np.nan})
        for predictor in ['library_counts', 'detected_genes', 'k_positive', 'score_genes_detected']:
            for outcome in ['score', 'occupancy_component', 'positive_rank_residual']:
                relationships.append({**base, 'predictor': predictor, 'outcome': outcome,
                                      'spearman_r': correlation(sub[predictor], sub[outcome]),
                                      'inference': 'within_sample_descriptive_no_cell_level_p_value'})
        if len(sub):
            # Within-sample/population quantile bins; tied edges may yield fewer than ten.
            assignments = pd.qcut(sub.library_counts, 10, labels=False, duplicates='drop')
            for decile in sorted(assignments.dropna().unique()):
                in_bin = sub[assignments == decile]
                bins.append({**base, 'depth_quantile_bin': int(decile) + 1,
                             'number_of_nonempty_bins': int(assignments.nunique()),
                             'minimum_library_counts': in_bin.library_counts.min(),
                             'maximum_library_counts': in_bin.library_counts.max(),
                             'bin_n': len(in_bin), 'median_library_counts': in_bin.library_counts.median(),
                             **{'mean_' + name: in_bin[name].mean() for name in ['score', 'occupancy_component', 'positive_rank_residual', 'k_up', 'k_down']}})
    return summaries, relationships, bins


def thinning_tables(mapped, frame, row, up, down, depths, seeds, populations, minimum_cells):
    records = []
    for depth in depths:
        eligible = (frame.library_counts.to_numpy() >= depth) & frame.alpha_broad.to_numpy()
        positions = np.flatnonzero(eligible)
        source = mapped[:, positions]
        totals = frame.library_counts.to_numpy()[positions]
        native = frame.iloc[positions]
        per_cell = {name: [] for name in COMPONENTS}
        for replicate, seed in enumerate(seeds):
            rng = np.random.default_rng(np.random.SeedSequence([int(seed), int(row['accession'][3:]), int(depth)]))
            sampled = thin_sparse_counts(source, totals, depth, rng)
            components = components_sparse(sampled, up, down)
            for name in COMPONENTS:
                per_cell[name].append(components[name])
            for population in populations:
                original_mask = frame[population].to_numpy()
                use = native[population].to_numpy()
                record = {'accession': row['accession'], 'donor': row['donor'],
                          'treatment': row['treatment_code'], 'population': population,
                          'depth': depth, 'replicate': replicate, 'seed': seed,
                          'population_n': int(original_mask.sum()), 'eligible_n': int(use.sum()),
                          'excluded_below_depth': int(original_mask.sum() - use.sum()),
                          'retained_fraction': float(use.sum() / original_mask.sum()) if original_mask.any() else np.nan,
                          'minimum_cells_report': minimum_cells,
                          'status': 'conditional_on_eligible_original_cells' if use.sum() >= minimum_cells else 'fewer_than_minimum_eligible_cells'}
                for component in COMPONENTS:
                    record['original_' + component] = float(native.loc[use, component].mean()) if use.sum() >= minimum_cells else np.nan
                    record['thinned_' + component] = float(components[component][use].mean()) if use.sum() >= minimum_cells else np.nan
                records.append(record)
        cache = data_path(f"revision_01/{row['accession']}_depth{depth}_thinning.npz", for_write=True)
        if cache.exists():
            raise FileExistsError('Refusing to overwrite revision cache: ' + str(cache))
        np.savez_compressed(cache, source_cell_index=native.source_cell_index.to_numpy(),
                            barcode=native.barcode.to_numpy(dtype=str), seeds=np.asarray(seeds),
                            library_counts=totals, **{name: np.asarray(value) for name, value in per_cell.items()})
        print(row['accession'], 'depth', depth, 'eligible alpha proxies', len(positions), len(seeds), 'draws complete', flush=True)
    return records


def paired_contrasts(table, keys, value_columns):
    control = table[table.treatment == 'D'].set_index(keys)
    treated = table[table.treatment == 'H_GLP'].set_index(keys)
    if not control.index.equals(treated.index):
        control, treated = control.align(treated, join='outer', axis=0)
    output = control.index.to_frame(index=False)
    for name in value_columns:
        output['control_' + name] = control[name].to_numpy()
        output['treated_' + name] = treated[name].to_numpy()
        output['delta_' + name] = treated[name].to_numpy() - control[name].to_numpy()
    output['interpretation'] = 'paired_donor_descriptive_H_GLP_minus_vehicle'
    return output


def load_protocol():
    """Verify the registered scientific parameters and score, without private records."""
    return protocol('revision_01'), frozen_score()


def main():
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    plan, config = load_protocol()
    thinning_plan = plan['thinning']
    depths, seeds = thinning_plan['depths'], thinning_plan['seeds']
    thinning_populations = thinning_plan['populations']
    minimum_cells = thinning_plan['minimum_cells_report']
    if depths != [1000, 2000, 5000, 10000] or seeds != list(range(7600, 7605)) or thinning_plan['replicates'] != len(seeds):
        raise ValueError('The frozen four-depth, five-seed thinning design is required.')
    if thinning_populations != ['alpha_broad', 'alpha_strict', 'alpha_cycle_high'] or minimum_cells != 20:
        raise ValueError('Thinning population or minimum cell requirements differ from the published scientific parameters.')
    if (OUT / 'mechanism_status.json').exists():
        raise FileExistsError('Revision mechanism outputs already exist; no automatic overwrite.')
    universe = config['universe']
    index = {gene: position for position, gene in enumerate(universe)}
    up = [index[gene] for gene in config['up']]
    down = [index[gene] for gene in config['down']]
    manifest = pd.read_csv(ROOT / 'primary_sample_manifest.tsv', sep='\t')
    manifest = manifest[manifest.treatment_code.isin(['D', 'H_GLP'])].sort_values(['donor', 'treatment_code'])
    if len(manifest) != 8 or manifest.donor.nunique() != 4 or not (manifest.groupby('donor').treatment_code.nunique() == 2).all():
        raise ValueError('Exactly four donors with D and H_GLP are required.')
    historical = pd.read_csv(evidence_path('primary_cell_scores.tsv.gz'), sep='\t', float_precision='round_trip')
    if historical.duplicated(['accession', 'barcode']).any():
        raise ValueError('Historical primary cell keys are not unique.')
    source_manifest = pd.read_csv(ROOT / 'source_manifest.tsv', sep='\t')
    source_hashes = dict(zip(source_manifest.local_path.str.replace('\\', '/', regex=False), source_manifest.sha256))
    # Check all owned destinations before calculation so a partial earlier run
    # cannot silently replace results or trigger avoidable recomputation.
    table_names = ['population_audit', 'observed_components', 'observed_component_contrasts',
                   'observed_depth_relationships', 'observed_depth_bins', 'thinning_sample_summaries',
                   'thinning_donor_contrasts', 'thinning_stability_summary']
    for name in table_names:
        if (OUT / (name + '.tsv')).exists():
            raise FileExistsError('Revision table already exists: ' + name)
    for accession in manifest.accession:
        cache_names = [f'{accession}_cell_components.tsv.gz'] + [f'{accession}_depth{depth}_thinning.npz' for depth in depths]
        for name in cache_names:
            if data_path('revision_01/' + name, for_write=True).exists():
                raise FileExistsError('Revision cache already exists: ' + name)
    external = data_path('revision_01', for_write=True)
    external.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    thresholds, audits, summaries, relationships, bins, draws = {}, [], [], [], [], []
    for row in manifest.to_dict('records'):
        sample_start = time.monotonic()
        verified_hashes = {}
        for key in ['matrix_path', 'features_path', 'barcodes_path']:
            relative = row[key].replace('\\', '/')
            digest = sha(data_path(relative, must_exist=True))
            if source_hashes.get(relative) != digest:
                raise ValueError('Deposited source checksum differs from committed manifest: ' + relative)
            verified_hashes[key] = digest
        counts, features, barcodes = read_raw_sample(row)
        cells, selected, threshold, audit = select_populations(counts, features, barcodes, row, thresholds.get(row['donor']))
        if row['treatment_code'] == 'D':
            thresholds[row['donor']] = threshold
        mapped = map_universe(counts[:, selected], features, universe)
        components = components_sparse(mapped, up, down)
        for name, values in components.items():
            cells[name] = values
        cells['score_genes_detected'] = cells.k_up + cells.k_down
        audit.update(reproduce_original_cells(cells, historical))
        audit['matrix_sha256'] = verified_hashes['matrix_path']
        audit['features_sha256'] = verified_hashes['features_path']
        audit['barcodes_sha256'] = verified_hashes['barcodes_path']
        cache = data_path(f"revision_01/{row['accession']}_cell_components.tsv.gz", for_write=True)
        if cache.exists():
            raise FileExistsError('Refusing to overwrite revision cache: ' + str(cache))
        cells.to_csv(cache, sep='\t', index=False)
        sample_summaries, sample_relationships, sample_bins = observed_tables(cells, audit)
        summaries.extend(sample_summaries)
        relationships.extend(sample_relationships)
        bins.extend(sample_bins)
        draws.extend(thinning_tables(mapped, cells, row, up, down, depths, seeds, thinning_populations, minimum_cells))
        audit['elapsed_seconds'] = time.monotonic() - sample_start
        audits.append(audit)
        print(row['accession'], row['donor'], row['treatment_code'], 'complete', flush=True)
        del counts, mapped, cells
    observed = pd.DataFrame(summaries)
    thinning = pd.DataFrame(draws)
    observed_contrasts = paired_contrasts(observed, ['donor', 'population', 'component'], ['mean', 'n'])
    thinning_contrasts = paired_contrasts(thinning, ['donor', 'population', 'depth', 'replicate', 'seed'],
                                         ['eligible_n', 'retained_fraction'] +
                                         ['original_' + name for name in COMPONENTS] +
                                         ['thinned_' + name for name in COMPONENTS])
    thinning_summary = []
    for keys, group in thinning_contrasts.groupby(['donor', 'population', 'depth'], sort=True):
        values = group.delta_thinned_score.dropna()
        thinning_summary.append(dict(zip(['donor', 'population', 'depth'], keys)) |
                                {'evaluable_replicates': len(values),
                                 'mean_delta': values.mean(), 'SD_delta': values.std(ddof=1),
                                 'minimum_delta': values.min(), 'maximum_delta': values.max(),
                                 'positive_fraction': float((values > 0).mean()) if len(values) else np.nan,
                                 'original_same_eligible_delta': group.delta_original_score.iloc[0],
                                 'control_eligible_n': group.control_eligible_n.iloc[0],
                                 'treated_eligible_n': group.treated_eligible_n.iloc[0],
                                 'interpretation': 'MC_variation_conditional_on_observed_eligible_cells_not_biological_replication'})
    tables = {'population_audit': pd.DataFrame(audits), 'observed_components': observed,
              'observed_component_contrasts': observed_contrasts,
              'observed_depth_relationships': pd.DataFrame(relationships), 'observed_depth_bins': pd.DataFrame(bins),
              'thinning_sample_summaries': thinning, 'thinning_donor_contrasts': thinning_contrasts,
              'thinning_stability_summary': pd.DataFrame(thinning_summary)}
    for name, table in tables.items():
        path = OUT / (name + '.tsv')
        if path.exists():
            raise FileExistsError('Refusing to overwrite revision table: ' + str(path))
        table.to_csv(path, sep='\t', index=False)
    status = {'status': 'COMPLETE', 'started_utc': started,
              'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'scientific_parameters_sha256': parameters_sha256(),
              'frozen_score_sha256': sha(scientific_file('evidence/frozen_score.json')),
              'code_sha256': {name: sha(ROOT / 'code' / name) for name in ['revision_mechanism.py', 'revision_score.py']},
              'sampling_law': 'multivariate hypergeometric without replacement; exact outside-universe category',
              'depths': depths, 'seeds': seeds, 'rng': 'PCG64 SeedSequence([seed,numeric_GSM,depth]); marginals',
              'raw_count_input': 'direct released MatrixMarket processed-count files; not FASTQ',
              'population_labels': 'expression proxies; not published cell annotations',
              'cycle_thresholds_by_donor': thresholds,
              'per_cell_cache_root': str(external), 'tables': {name: len(table) for name, table in tables.items()},
              'limitations': ['After-outcome revision analyses, not a new preregistered discovery',
                              'Thinning conditions on observed molecules and cannot restore unseen expression',
                              'Depth eligibility changes cells; original scores on exactly the same cells are reported',
                              'Monte Carlo draws are not independent biological replicates',
                              'Expression proxies do not establish lineage, cell function, or published annotation']}
    (OUT / 'mechanism_status.json').write_text(json.dumps(status, indent=2) + '\n', encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    from public_runtime import execute_stage
    execute_stage('mechanism', args.data_root, args.output_dir)
