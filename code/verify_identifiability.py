"""Construct and verify a conditional observation-equivalence witness.

This is a nonparametric identifiability demonstration, not an estimated biological
lineage model. Uses observed score strata and preserves all input row counts.
Only the Python standard library is required.
"""
from __future__ import annotations
from data_paths import evidence_path, resolve_data_root, ENV_NAME
import os
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def choose_field(names, supplied, candidates, required=True):
    if supplied:
        if supplied not in names:
            raise ValueError(f"Column {supplied!r} not present: {names}")
        return supplied
    for item in candidates:
        if item in names:
            return item
    if required:
        raise ValueError(f"Need an explicit column; candidate names {candidates}; actual {names}")
    return None


def open_text(path):
    return gzip.open(path, 'rt', encoding='utf-8', newline='') if path.suffix == '.gz' else path.open(encoding='utf-8', newline='')


def construct(args):
    source = args.input.resolve()
    with open_text(source) as handle:
        reader = csv.DictReader(handle, delimiter=args.delimiter)
        names = reader.fieldnames or []
        donor = choose_field(names, args.donor_column, ['donor', 'donor_id'])
        treatment = choose_field(names, args.treatment_column, ['treatment', 'treatment_code'])
        score = choose_field(names, args.score_column, ['score', 'conversion_score', 'primary_score'])
        annotation = choose_field(names, args.annotation_column, ['annotation', 'cell_type', 'celltype', 'state'], required=False)
        if args.input_scope != 'alpha_related_cells' and not annotation:
            raise ValueError('Mixed/all-QC scope requires an explicit annotation column.')
        counts = collections.Counter()
        total = 0
        for row in reader:
            value = float(row[score])
            if not math.isfinite(value) or not -1.0000000001 <= value <= 1.0000000001:
                raise ValueError(f'Expected finite rank-score in [-1,1], found {value}.')
            label = row[annotation] if annotation else 'alpha_related'
            # Quarter-width bins, with +1 retained in the last bin. Split at 0.
            bin_id = min(7, max(0, math.floor((value + 1) / .25)))
            counts[(row[donor], row[treatment], label, bin_id)] += 1
            total += 1
    sample_totals = collections.Counter()
    for (d, t, a, b), n in counts.items():
        sample_totals[d, t] += n
    fields = ['donor', 'treatment', 'annotation', 'score_lower_inclusive', 'score_upper',
              'observed_n', 'sample_captured_n', 'observed_probability',
              'recruitment_stable_alpha_n', 'recruitment_alpha_transition_n',
              'recruitment_beta_loss_n', 'loss_stable_alpha_n', 'loss_alpha_transition_n',
              'loss_beta_loss_n', 'other_observed_n', 'likelihood_ratio']
    rows = []
    for (d, t, label, b), n in sorted(counts.items()):
        lower = -1 + b*.25
        upper = lower + .25
        alpha = args.input_scope == 'alpha_related_cells' or label.lower() in args.alpha_labels
        # Deliberate witness choice, not an inferred rate. Only nonnegative-score bins.
        switched = (n + 1)//2 if alpha and lower >= 0 else 0
        stable = n - switched if alpha else 0
        other = 0 if alpha else n
        rows.append(dict(zip(fields, [d, t, label, lower, upper, n, sample_totals[d, t],
                                     n/sample_totals[d, t], stable, switched, 0,
                                     stable, 0, switched, other, 1.0])))
    if not any(row['recruitment_alpha_transition_n'] for row in rows):
        raise ValueError('No nonnegative-score alpha-related stratum: no nontrivial specified witness can be generated.')
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out = args.output_dir/'identifiability_witness.csv'
    with out.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    # Reopen the delivered witness; verify integer conservation independently
    # of the construction variables and compute the likelihood both ways.
    with out.open(encoding='utf-8', newline='') as handle:
        retained = list(csv.DictReader(handle))
    reconstructed = collections.Counter()
    rec_ll = 0.0
    loss_ll = 0.0
    transition_n = 0
    for row in retained:
        n = int(row['observed_n'])
        denominator = int(row['sample_captured_n'])
        rec_n = sum(int(row[k]) for k in ['recruitment_stable_alpha_n', 'recruitment_alpha_transition_n', 'recruitment_beta_loss_n', 'other_observed_n'])
        loss_n = sum(int(row[k]) for k in ['loss_stable_alpha_n', 'loss_alpha_transition_n', 'loss_beta_loss_n', 'other_observed_n'])
        if rec_n != n or loss_n != n:
            raise AssertionError('Witness lost an observed cell.')
        if int(row['recruitment_alpha_transition_n']) != int(row['loss_beta_loss_n']):
            raise AssertionError('Witness changed selected-cell support.')
        if not math.isclose(float(row['observed_probability']), n/denominator, abs_tol=1e-14):
            raise AssertionError('Witness changed the observed probability.')
        rec_ll += n * math.log(rec_n/denominator)
        loss_ll += n * math.log(loss_n/denominator)
        transition_n += int(row['recruitment_alpha_transition_n'])
        reconstructed[row['donor'], row['treatment']] += n
    if reconstructed != sample_totals or sum(reconstructed.values()) != total:
        raise AssertionError('Sample totals changed.')
    if rec_ll != loss_ll:
        raise AssertionError('Observed log-likelihoods differ.')
    result = {
        'status': 'VERIFIED_CONDITIONAL_OBSERVATIONAL_EQUIVALENCE',
        'input': str(source.relative_to(ROOT)) if source.is_relative_to(ROOT) else source.name,
        'input_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'witness_sha256': hashlib.sha256(out.read_bytes()).hexdigest(),
        'input_scope': args.input_scope,
        'columns': {'donor': donor, 'treatment': treatment, 'score': score, 'annotation': annotation},
        'observed_rows': total, 'samples': len(sample_totals), 'nonempty_strata': len(retained),
        'arbitrarily_reassigned_latent_labels': transition_n,
        'score_edges': [-1, -.75, -.5, -.25, 0, .25, .5, .75, 1],
        'bin_convention': 'left-closed right-open except last bin includes 1',
        'recruitment_conditional_log_likelihood_without_common_multinomial_constant': rec_ll,
        'loss_conditional_log_likelihood_without_common_multinomial_constant': loss_ll,
        'log_likelihood_difference': rec_ll-loss_ll,
        'all_observed_probabilities_positive': True,
        'likelihood_ratio': 1.0,
        'absolute_tissue_counts_measured': False,
        'assumptions': ['ancestry and historical transition unobserved',
                        'no externally fixed lineage-specific emission laws or transition bounds',
                        'likelihood conditional on capture and input selection',
                        'arbitrary hidden labels are witnesses, not biological estimates'],
        'claim': 'These observations alone cannot select between these latent histories in the declared unrestricted model.',
        'nonclaim': 'No equal biological plausibility, measured transition rate, function, or cure is established.'
    }
    (args.output_dir/'identifiability_verification.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--input-scope', choices=['all_qc_cells', 'alpha_related_cells', 'selected_scored_cells'], required=True)
    parser.add_argument('--output-dir', type=Path, required=True, help='New output directory outside the repository')
    parser.add_argument('--delimiter', default='\t')
    for field in ['donor', 'treatment', 'score', 'annotation']:
        parser.add_argument('--'+field+'-column')
    parser.add_argument('--alpha-labels', nargs='+', default=['alpha', 'alpha_related', 'alpha-related', 'alpha_like', 'alpha-like'])
    arguments = parser.parse_args()
    arguments.output_dir = resolve_data_root(arguments.output_dir, for_write=True)
    for name in ['identifiability_witness.csv', 'identifiability_verification.json']:
        if (arguments.output_dir/name).exists():
            raise FileExistsError('Output already exists: '+str(arguments.output_dir/name))
    if arguments.data_root is not None or arguments.input is None:
        os.environ[ENV_NAME] = str(resolve_data_root(arguments.data_root, require_exists=True))
    if arguments.input is None:
        arguments.input = evidence_path('primary_cell_scores.tsv.gz', must_exist=True)
    arguments.alpha_labels = [x.lower() for x in arguments.alpha_labels]
    construct(arguments)
