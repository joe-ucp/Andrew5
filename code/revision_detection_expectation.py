"""Exact conditional detection expectations for the existing finite-library draws.

This numerical addendum changes no program, population, scorer, or random draw.
Only occupancy has an analytic expectation here; rank residuals are not predicted.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from data_paths import data_path, resolve_data_root, ENV_NAME

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'evidence' / 'revision_01'
from public_runtime import parameters_sha256, scientific_file
COMPONENTS = ['occupancy_component', 'occupancy_up', 'occupancy_down', 'k_up', 'k_down']


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def integer_array(values, label):
    values = np.asarray(values)
    if values.dtype.kind not in 'iuf' or values.dtype.kind == 'b':
        raise ValueError(label + ' must contain numeric nonnegative integers, not booleans.')
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
        raise ValueError(label + ' must contain finite nonnegative integers.')
    if (values >= 2**63).any():
        raise ValueError(label + ' exceeds the supported signed 64-bit integer range.')
    return values.astype(np.int64)


def detection_probability(total, count, depth):
    """P(gene detected) under uniform sampling without replacement.

    Inputs broadcast. The empty library at depth zero is supported. Zero count,
    zero depth, full depth, and forced detection are handled without log-PMF
    boundary ambiguities. Singleton successes/draws have exact simple ratios.
    Remaining cases use -expm1(log P(no detection)), avoiding 1-exp cancellation.
    """
    total, count, depth = np.broadcast_arrays(
        integer_array(total, 'Total'), integer_array(count, 'Gene count'),
        integer_array(depth, 'Depth'))
    if (count > total).any() or (depth > total).any():
        raise ValueError('Gene count and depth must each lie between zero and total.')
    probability = np.zeros(total.shape, dtype=float)
    positive = (count > 0) & (depth > 0)
    forced = positive & (depth > total - count)
    probability[forced] = 1.
    interior = positive & ~forced
    singleton = interior & ((count == 1) | (depth == 1))
    probability[singleton] = (count[singleton].astype(float) * depth[singleton]
                              / total[singleton])
    general = interior & ~singleton
    if general.any():
        log_zero = hypergeom.logpmf(0, total[general], count[general], depth[general])
        probability[general] = -np.expm1(log_zero)
    if (not np.isfinite(probability).all() or (probability < -1e-10).any()
            or (probability > 1 + 1e-10).any()):
        raise ArithmeticError('Hypergeometric detection probability is invalid.')
    return np.clip(probability, 0., 1.)


def expected_components(counts_up, counts_down, totals, depth):
    """Expected detection components, genes by cells; no independence assumed."""
    upper = integer_array(counts_up, 'Up-block counts')
    lower = integer_array(counts_down, 'Down-block counts')
    totals = integer_array(totals, 'Totals')
    if upper.ndim != 2 or lower.ndim != 2 or not len(upper) or not len(lower):
        raise ValueError('Gene blocks must be nonempty two-dimensional arrays.')
    if totals.ndim != 1 or upper.shape[1] != len(totals) or lower.shape[1] != len(totals):
        raise ValueError('Gene-block cell dimensions and totals do not agree.')
    # Blocks are disjoint categories. Outside-block molecules remain in totals.
    if (upper.sum(axis=0) + lower.sum(axis=0) > totals).any():
        raise ValueError('Disjoint gene-block counts exceed the whole library.')
    kup = detection_probability(totals[None, :], upper, depth).sum(axis=0)
    kdown = detection_probability(totals[None, :], lower, depth).sum(axis=0)
    up = .5 * kup / len(upper)
    down = -.5 * kdown / len(lower)
    return dict(occupancy_component=up + down, occupancy_up=up,
                occupancy_down=down, k_up=kup, k_down=kdown)


def native_components(upper, lower):
    kup, kdown = (upper > 0).sum(axis=0), (lower > 0).sum(axis=0)
    up, down = .5 * kup / len(upper), -.5 * kdown / len(lower)
    return dict(occupancy_component=up + down, occupancy_up=up,
                occupancy_down=down, k_up=kup, k_down=kdown)


def checked_mc_rows(group, expected_seeds, eligible_n, native, components=COMPONENTS):
    if len(group) != len(expected_seeds) or sorted(group.seed.astype(int)) != list(expected_seeds):
        raise ValueError('Existing Monte Carlo draws do not match the locked seeds.')
    if not (group.eligible_n == eligible_n).all():
        raise ValueError('Existing draw population differs from analytic eligibility.')
    for component in components:
        actual = group['original_' + component].to_numpy()
        target = native[component].mean() if eligible_n >= 20 else np.nan
        if not np.allclose(actual, target, atol=2e-14, rtol=0, equal_nan=True):
            raise ValueError('Existing draw native component does not reproduce: ' + component)


def summary_rows(base, expected, native, draws):
    result = []
    valid = base['eligible_n'] >= 20
    for component in COMPONENTS:
        prediction = float(expected[component].mean()) if valid else np.nan
        before = float(native[component].mean()) if valid else np.nan
        values = np.asarray(draws[component], dtype=float)
        values = values[np.isfinite(values)]
        if valid and len(values) != 5:
            raise ValueError('All five existing technical draws are required.')
        if not valid and len(values):
            raise ValueError('Below-minimum group unexpectedly reports a draw mean.')
        mean = float(values.mean()) if len(values) else np.nan
        sd = float(values.std(ddof=1)) if len(values) > 1 else np.nan
        error = mean - prediction
        result.append({**base, 'component': component,
                       'status': 'conditional_on_eligible_observed_cells' if valid else 'fewer_than_20_eligible_cells',
                       'native_eligible_mean': before, 'expected_thinned_mean': prediction,
                       'expected_thinning_change': prediction - before,
                       'mc_thinned_mean': mean, 'mc_repeats': len(values),
                       'mc_repeat_sd': sd, 'mc_mean_mcse': sd / np.sqrt(len(values)) if len(values) else np.nan,
                       'mc_minimum': float(values.min()) if len(values) else np.nan,
                       'mc_maximum': float(values.max()) if len(values) else np.nan,
                       'mc_mean_minus_expected': error,
                       'mc_mean_error_over_repeat_sd': error / sd if sd > 0 else np.nan,
                       'inference': 'numerical_MC_check_not_biological_replication_or_causal_depth_effect'})
    return result


def save_cells(bulk, filename, ids, totals, predicted):
    target = bulk / filename
    if target.exists():
        raise FileExistsError('Refusing to replace analytic per-cell output: ' + str(target))
    pd.DataFrame({'cell_id': ids, 'total_counts': totals, **predicted}).to_csv(target, sep='\t', index=False)
    return {'path': str(target), 'sha256': sha(target), 'cells': len(ids)}


def primary_predictions(config, plan, bulk):
    from primary_analysis import load_sample, map_universe
    manifest = pd.read_csv(ROOT / 'primary_sample_manifest.tsv', sep='\t')
    manifest = manifest[manifest.treatment_code.isin(['D', 'H_GLP'])].sort_values(['donor', 'treatment_code'])
    if len(manifest) != 8 or manifest.donor.nunique() != 4:
        raise ValueError('Expected eight fixed primary libraries from four donors.')
    mc_table = pd.read_csv(scientific_file('evidence/revision_01/thinning_sample_summaries.tsv'), sep='\t', float_precision='round_trip')
    idx = {gene: i for i, gene in enumerate(config['universe'])}
    up, down = [idx[g] for g in config['up']], [idx[g] for g in config['down']]
    rows, outputs, inputs = [], [], []
    for row in manifest.to_dict('records'):
        accession = row['accession']
        path = data_path(f'revision_01/{accession}_cell_components.tsv.gz', must_exist=True)
        frame = pd.read_csv(path, sep='\t', float_precision='round_trip')
        frame = frame.loc[frame.alpha_broad].reset_index(drop=True)
        if frame.barcode.duplicated().any():
            raise ValueError('Primary barcode keys are not unique.')
        x, features, barcodes = load_sample(row)
        selected = frame.source_cell_index.to_numpy(dtype=np.int64)
        if not np.array_equal(barcodes.barcode.iloc[selected].to_numpy(), frame.barcode):
            raise ValueError('Fixed primary cell indices and barcodes differ.')
        total = np.asarray(x[:, selected].sum(axis=0)).ravel().astype(np.int64)
        if not np.array_equal(total, frame.library_counts):
            raise ValueError('Primary whole-library totals differ from fixed cell table.')
        mapped = map_universe(x[:, selected], features, config['universe'])
        upper, lower = mapped[up].toarray(), mapped[down].toarray()
        native = native_components(upper, lower)
        for component in COMPONENTS:
            if not np.allclose(native[component], frame[component], atol=2e-15, rtol=0):
                raise ValueError('Native primary detection component mismatch.')
        inputs.append({'accession': accession, 'components_sha256': sha(path),
                       'sparse_cache_sha256': sha(data_path(f'data_primary/sparse_cache/{accession}.npz', must_exist=True))})
        for depth in plan['thinning']['depths']:
            eligible = total >= depth
            expected = expected_components(upper[:, eligible], lower[:, eligible], total[eligible], depth)
            before = {k: v[eligible] for k, v in native.items()}
            cache = data_path(f'revision_01/{accession}_depth{depth}_thinning.npz', must_exist=True)
            group = mc_table[(mc_table.accession == accession) & (mc_table.depth == depth)
                             & (mc_table.population == 'alpha_broad')].sort_values('replicate')
            checked_mc_rows(group, plan['thinning']['seeds'], int(eligible.sum()), before)
            with np.load(cache) as z:
                if (not np.array_equal(z['barcode'], frame.barcode.to_numpy()[eligible])
                        or not np.array_equal(z['source_cell_index'], selected[eligible])
                        or not np.array_equal(z['library_counts'], total[eligible])
                        or not np.array_equal(z['seeds'], plan['thinning']['seeds'])):
                    raise ValueError('Cached draws do not use the exact eligible primary cells.')
                for component in COMPONENTS:
                    from_cells = z[component].mean(axis=1) if eligible.sum() >= 20 else np.full(5, np.nan)
                    if not np.allclose(from_cells, group['thinned_' + component], atol=2e-14, rtol=0, equal_nan=True):
                        raise ValueError('Cached draw components disagree with existing summary.')
            base = dict(dataset='GSE266427', accession=accession, donor=row['donor'],
                        treatment=row['treatment_code'], population='alpha_broad', mode='original_universe',
                        depth=depth, population_n=len(total), eligible_n=int(eligible.sum()),
                        excluded_n=int((~eligible).sum()))
            rows += summary_rows(base, expected, before, {k: group['thinned_' + k] for k in COMPONENTS})
            outputs.append(save_cells(bulk, f'{accession}_alpha_broad_depth{depth}.tsv.gz',
                                      frame.barcode.to_numpy()[eligible], total[eligible], expected))
        print('Analytic detection complete', accession, flush=True)
        del x, mapped, upper, lower
    return rows, outputs, inputs


def independent_predictions(config, plan, bulk):
    registry = json.loads((scientific_file('evidence/revision_01/baron_input_manifest.json')).read_text())
    draws_path = data_path('revision_01/independent/baron_component_thinning_draws.tsv.gz', must_exist=True)
    draws = pd.read_csv(draws_path, sep='\t', float_precision='round_trip')
    draws = draws[draws['mode'] == 'original_universe']
    report = pd.read_csv(scientific_file('evidence/revision_01/independent_component_thinning.tsv'), sep='\t', float_precision='round_trip')
    report = report[report['mode'] == 'original_universe']
    cell_path = data_path('revision_01/independent/baron_native_components.tsv.gz', must_exist=True)
    fixed = pd.read_csv(cell_path, sep='\t', float_precision='round_trip')
    fixed = fixed[fixed['mode'] == 'original_universe']
    rows, outputs, inputs = [], [], []
    for donor, record in enumerate(sorted(registry, key=lambda r: r['name']), 1):
        path = data_path('revision_01/sources/' + record['name'], must_exist=True)
        if sha(path) != record['sha256']:
            raise ValueError('Baron original count-file hash mismatch.')
        frame = pd.read_csv(path, index_col=0)
        if not frame.index.is_unique:
            raise ValueError('Original Baron cell IDs must be unique.')
        count = frame.iloc[:, 2:]
        total = count.sum(axis=1).to_numpy(dtype=np.int64)
        # Zero-fill absent frozen genes: this is original-universe mode only.
        upper = count.reindex(columns=config['up'], fill_value=0).to_numpy().T
        lower = count.reindex(columns=config['down'], fill_value=0).to_numpy().T
        native = native_components(upper, lower)
        original_ids = frame.index.to_numpy(dtype=str)
        for label in ['alpha', 'beta']:
            membership = frame.assigned_cluster.to_numpy() == label
            matching = fixed[(fixed.donor == donor) & (fixed.cell_type == label)]
            if not np.array_equal(original_ids[membership], matching.cell_id):
                raise ValueError('Baron fixed author-labelled cell IDs differ.')
            if not np.array_equal(total[membership], matching.total_counts):
                raise ValueError('Baron whole-library totals differ.')
            for component in COMPONENTS:
                if not np.allclose(native[component][membership], matching[component], atol=2e-15, rtol=0):
                    raise ValueError('Baron native detection component mismatch.')
            for depth in plan['independent_validation']['thinning_depths']:
                eligible = membership & (total >= depth)
                expected = expected_components(upper[:, eligible], lower[:, eligible], total[eligible], depth)
                before = {k: v[eligible] for k, v in native.items()}
                group = draws[(draws.donor == donor) & (draws.cell_type == label) & (draws.depth == depth)]
                values = {}
                for component in COMPONENTS:
                    subset = group[group.component == component].sort_values('replicate')
                    check = subset.rename(columns={'native_eligible_mean': 'original_' + component})
                    checked_mc_rows(check, plan['independent_validation']['seeds'], int(eligible.sum()), before, [component])
                    existing = report[(report.donor == donor) & (report.cell_type == label)
                                      & (report.depth == depth) & (report.component == component)]
                    if len(existing) != 1 or not np.allclose(existing.thinned_mean.iloc[0], subset.thinned_mean.mean(),
                                                           atol=2e-14, rtol=0, equal_nan=True):
                        raise ValueError('Baron draws differ from the published revision summary.')
                    values[component] = subset.thinned_mean.to_numpy()
                base = dict(dataset='GSE84133', accession=record['name'].split('_')[0], donor=f'human{donor}',
                            treatment='unperturbed_reference', population=label, mode='original_universe',
                            depth=depth, population_n=int(membership.sum()), eligible_n=int(eligible.sum()),
                            excluded_n=int((membership & ~eligible).sum()))
                rows += summary_rows(base, expected, before, values)
                outputs.append(save_cells(bulk, f'baron{donor}_{label}_depth{depth}.tsv.gz',
                                          original_ids[eligible], total[eligible], expected))
        inputs.append({'donor': donor, 'path': str(path), 'sha256': record['sha256']})
        print('Analytic detection complete Baron donor', donor, flush=True)
    return rows, outputs, inputs


def main():
    from revision_mechanism import load_protocol
    plan, config = load_protocol()
    if plan['thinning']['depths'] != [1000, 2000, 5000, 10000] or plan['independent_validation']['thinning_depths'] != [1000, 2000, 5000, 10000]:
        raise ValueError('Only the published four depths are supported.')
    target = OUT / 'detection_expectation_samples.tsv'
    status_path = OUT / 'detection_expectation_status.json'
    if target.exists() or status_path.exists():
        raise FileExistsError('Detection-expectation outputs already exist.')
    bulk = data_path('revision_01/detection_expectation', for_write=True)
    if bulk.exists():
        raise FileExistsError('Detection-expectation bulk output already exists.')
    bulk.mkdir(parents=True)
    rows, outputs, sources = primary_predictions(config, plan, bulk)
    second, more_outputs, more_sources = independent_predictions(config, plan, bulk)
    table = pd.DataFrame(rows + second)
    table.to_csv(target, sep='\t', index=False)
    occ = table[(table.component == 'occupancy_component') & (table.mc_repeats == 5)]
    diagnostics = []
    for dataset, subset in occ.groupby('dataset'):
        diagnostics.append({'dataset': dataset, 'evaluable_sample_population_depth_groups': len(subset),
                            'maximum_absolute_MC_mean_error': float(subset.mc_mean_minus_expected.abs().max()),
                            'median_absolute_MC_mean_error': float(subset.mc_mean_minus_expected.abs().median()),
                            'maximum_absolute_error_over_repeat_SD': float(subset.mc_mean_error_over_repeat_sd.abs().max()),
                            'scope': 'descriptive_numerical_check_no_threshold_selected_or_test_of_biology'})
    status = {'status': 'COMPLETE', 'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'scientific_parameters_sha256': parameters_sha256(), 'score_sha256': sha(scientific_file('evidence/frozen_score.json')),
              'code_sha256': sha(Path(__file__)), 'table_sha256': sha(target),
              'table_rows': len(table), 'components': COMPONENTS, 'new_draws': 0,
              'new_scorers_or_refits': 0, 'conditional_expectation': '1-C(N-c,d)/C(N,d)',
              'evaluation': 'float64 hypergeom.logpmf with -expm1; exact boundary/singleton simplifications',
              'occupancy_diagnostics': diagnostics,
              'source_inputs': sources + more_sources,
              'existing_primary_summary_sha256': sha(scientific_file('evidence/revision_01/thinning_sample_summaries.tsv')),
              'existing_independent_summary_sha256': sha(scientific_file('evidence/revision_01/independent_component_thinning.tsv')),
              'existing_independent_draws_sha256': sha(data_path('revision_01/independent/baron_component_thinning_draws.tsv.gz', must_exist=True)),
              'per_cell_outputs': outputs + more_outputs,
              'limitations': ['Exact conditional expectation of detection occupancy only, not of the rank residual or total score.',
                              'Conditions on observed finite libraries, fixed cell membership, and locked depth eligibility.',
                              'Five existing draw means estimate Monte Carlo variability, not donor uncertainty.',
                              'No causal attribution to sequencing depth, conversion, ancestry, or function.']}
    status_path.write_text(json.dumps(status, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'status': status['status'], 'table_rows': len(table), 'diagnostics': diagnostics}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    from public_runtime import execute_stage
    execute_stage('detection', args.data_root, args.output_dir)
