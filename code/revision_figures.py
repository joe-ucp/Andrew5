"""Publication figures from retained revision tables; no scientific reanalysis."""
from data_paths import evidence_path, logical_evidence_key
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os
import textwrap
import numpy as np
import pandas as pd
from data_paths import data_path, resolve_data_root, ENV_NAME

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / 'evidence'
REV = EV / 'revision_01'
OUT = ROOT / 'manuscript' / 'figures'
DONORS = ['D1', 'D2', 'D3', 'D4']
COLORS = {'D': '#355C7D', 'H_GLP': '#C96F32', 'up': '#35826A', 'down': '#8D4E8A',
          'residual': '#B5A044', 'neutral': '#63727E'}
DONOR_COLORS = ['#355C7D', '#C96F32', '#35826A', '#8D4E8A']
METHODS = ['signed_rank', 'signed_log1p_CPT10K', 'ucell_form_1500']
METHOD_LABELS = ['Signed rank', 'Signed log1p', 'UCell-form']
PROGRAMS = ['original', 'omit_IIDP1697', 'omit_IIDP1738', 'omit_IIDP1770', 'omit_IIDP1776']
PROGRAM_LABELS = ['Original', 'Omit 1697', 'Omit 1738', 'Omit 1770', 'Omit 1776']
PROGRAM_COLORS = ['#172D3D', '#C96F32', '#35826A', '#8D4E8A', '#729BA5']
records = {}


def read(name, revision=True):
    return pd.read_csv(REV / name if revision else evidence_path(name, must_exist=True), sep='\t')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def panel(ax, label, title):
    ax.set_title(label + '  ' + textwrap.fill(title, 42), loc='left', fontweight='bold', pad=10)
    ax.grid(axis='y', color='#E6E9EC', linewidth=.6, zorder=0)


def save(fig, name, sources, note):
    from matplotlib.text import Text
    for item in fig.findobj(Text):
        item.set_fontsize(max(7, item.get_fontsize() * .8))
    fig.text(.02, .012, textwrap.fill(note, 132), fontsize=6.8, color='#435460', va='bottom')
    OUT.mkdir(parents=True, exist_ok=True)
    png, pdf = OUT / (name + '.png'), OUT / (name + '.pdf')
    fig.savefig(png, dpi=220, bbox_inches='tight', facecolor='white', metadata={'Software': 'STRIKE007 revision figures'})
    fig.savefig(pdf, bbox_inches='tight', facecolor='white', metadata={'Creator': 'STRIKE007 revision figures', 'CreationDate': None, 'ModDate': None})
    plt.close(fig)
    records[name] = {'sources': [{'path': logical_evidence_key(path), 'sha256': sha(path)} for path in sources],
                     'png_sha256': sha(png), 'pdf_sha256': sha(pdf), 'code_sha256': sha(Path(__file__)),
                     'note': note, 'visual_review': 'PENDING'}
    print('Rendered', name, flush=True)


def figure1():
    points = read('operating_points.tsv')
    cells = read('reference_cell_scores.tsv', revision=False)
    capture = read('reference_capture_summary.tsv', revision=False)
    native = points[(points.depth == 'native') & (points.comparator == 'alpha_only')].iloc[0]
    fig, axes = plt.subplots(2, 2, figsize=(6.65, 6.2))
    a, b, c, d = axes.ravel()
    positive = cells[cells.preparation == 'PM'].score.to_numpy()
    control = cells[cells.preparation == 'alpha'].score.to_numpy()
    thresholds = np.r_[np.inf, np.sort(np.unique(np.r_[positive, control]))[::-1], -np.inf]
    fpr = np.asarray([(control > threshold).mean() for threshold in thresholds])
    tpr = np.asarray([(positive > threshold).mean() for threshold in thresholds])
    a.step(fpr, tpr, where='post', color=COLORS['D'], lw=2, label=f"AUROC = {native.empirical_AUROC:.3f}")
    a.plot([0, 1], [0, 1], color='#A3ADB5', ls=':', lw=1)
    a.scatter([native.achieved_comparator_FPR], [native.positive_preparation_detection_fraction], color=COLORS['H_GLP'], s=60, zorder=4)
    a.set(xlim=(0, 1), ylim=(0, 1.02), xlabel='Empirical alpha false-positive fraction', ylabel='PM-preparation detection fraction')
    a.legend(frameon=False, loc='lower right')
    panel(a, 'A', 'Native preparation ROC')
    values = [native.achieved_comparator_FPR, native.positive_preparation_detection_fraction]
    b.barh([1, 0], values, color=[COLORS['D'], COLORS['H_GLP']], height=.5, zorder=3)
    for y, value in zip([1, 0], values):
        b.text(value + .025, y, f'{100 * value:.1f}%', va='center', fontsize=12, fontweight='bold')
    b.set(yticks=[1, 0], yticklabels=['Alpha FPR', 'PM detection'], xlim=(0, .72), ylim=(-.6, 1.6), xlabel='Observed fraction at alpha q95')
    b.text(.96, .95, 'Threshold = 0; score > threshold\n158 alpha cells; 860 PM cells', transform=b.transAxes, ha='right', va='top', fontsize=9)
    panel(b, 'B', 'Achieved native operating point')
    for comparator, label, color, marker in [('alpha_only', 'Alpha control', COLORS['D'], 'o'),
                                            ('low_mixtures', '0-10% beta-RNA mixtures', COLORS['H_GLP'], 's'),
                                            ('stress_inclusive', 'Including 20-50% stress', COLORS['down'], '^')]:
        subset = points[(points.depth != 'native') & (points.comparator == comparator)].copy()
        subset['depth_numeric'] = pd.to_numeric(subset.depth)
        subset = subset.sort_values('depth_numeric')
        c.plot(subset.depth_numeric, subset.positive_preparation_detection_fraction, marker=marker, color=color, label=label, lw=1.7,
               markerfacecolor='none', markersize={'o': 6, 's': 7, '^': 8}[marker], markeredgewidth=1.4)
    c.axhline(native.positive_preparation_detection_fraction, color='#65747F', ls=':', lw=1.3, label='Native PM detection: 48.5%')
    c.set(xlabel='Simulated count depth (multinomial)', ylabel='PM detection at comparator q95', ylim=(0, .70))
    c.legend(frameon=False, fontsize=8.5, loc='upper left')
    panel(c, 'C', 'Depth simulations, separate from native cells')
    order = ['alpha', 'PM', 'beta', 'graft']
    cap = capture.set_index('preparation').loc[order]
    bars = d.bar(np.arange(4), cap.median_detected_score_genes, color=[COLORS['D'], COLORS['H_GLP'], COLORS['up'], '#9B919F'], width=.6, zorder=3)
    for bar, value in zip(bars, cap.median_detected_score_genes):
        d.text(bar.get_x() + bar.get_width() / 2, value + .45, f'{value:g}', ha='center', fontsize=11)
    d.set(xticks=np.arange(4), xticklabels=['Alpha', 'PM', 'Beta', 'Graft*'], ylim=(0, 19), ylabel='Median detected score genes (of 100)')
    panel(d, 'D', 'Native feature capture')
    fig.subplots_adjust(left=.11, right=.98, bottom=.15, top=.93, wspace=.40, hspace=.64)
    save(fig, 'rev1_figure1_reference', [REV / 'operating_points.tsv', evidence_path('reference_cell_scores.tsv'), EV / 'reference_capture_summary.tsv'],
         'Preparation detection is not sensitivity for true conversion. *Graft preparation was author-flagged for quality; shown descriptively.')


def figure2():
    depth = read('observed_depth_summary.tsv')
    bins = read('observed_depth_bins.tsv')
    correlations = read('observed_depth_relationships.tsv')
    contrasts = read('observed_component_contrasts.tsv')
    fig = plt.figure(figsize=(6.65, 8.4))
    grid = fig.add_gridspec(3, 2, height_ratios=[1.05, 1, 1])
    a, b = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
    positions = np.arange(4)
    for shift, treatment, label in [(-.18, 'D', 'Vehicle'), (.18, 'H_GLP', 'H + GLP-1')]:
        subset = depth.set_index('group').loc[[f'{donor}:{treatment}:alpha_related' for donor in DONORS]]
        values = subset.depth_median.to_numpy()
        a.bar(positions + shift, values, width=.34, color=COLORS[treatment], label=label, zorder=3)
        for x, value in zip(positions + shift, values):
            a.text(x, value + 300, f'{value / 1000:.1f}k', ha='center', fontsize=8.5)
    a.set(xticks=positions, xticklabels=DONORS, ylim=(0, 20000), ylabel='Median whole-library counts')
    a.legend(frameon=False, ncol=2, fontsize=9)
    panel(a, 'A', 'Observed broad-alpha depth')
    populations = ['alpha_broad', 'alpha_strict', 'alpha_cycle_high', 'beta', 'mixed_endocrine', 'other_unassigned_endocrine']
    labels = ['Broad\nalpha', 'Strict\nalpha', 'Cycle-\nhigh', 'Beta', 'Mixed', 'Other']
    native = contrasts[contrasts.component == 'score']
    for donor, offset, color, marker in zip(DONORS, [-.21, -.07, .07, .21], DONOR_COLORS, ['o', 's', '^', 'D']):
        values = native[native.donor == donor].set_index('population').loc[populations].delta_mean
        b.scatter(np.arange(6) + offset, values, color=color, marker=marker, s=32, label=donor, zorder=3)
    b.axhline(0, color='#7C8790', lw=.8)
    b.set(xticks=np.arange(6), xticklabels=labels, ylim=(-.02, .058), ylabel='Native mean score contrast')
    b.legend(frameon=False, ncol=4, fontsize=8.5, loc='upper left')
    panel(b, 'B', 'Population-specific native contrasts')
    for number, donor in enumerate(DONORS):
        ax = fig.add_subplot(grid[1 + number // 2, number % 2])
        for treatment, label in [('D', 'Vehicle'), ('H_GLP', 'H + GLP-1')]:
            sub = bins[(bins.donor == donor) & (bins.treatment == treatment) & (bins.population == 'alpha_broad')].sort_values('median_library_counts')
            rho = correlations[(correlations.donor == donor) & (correlations.treatment == treatment) & (correlations.population == 'alpha_broad') & (correlations.predictor == 'library_counts') & (correlations.outcome == 'score')].spearman_r.iloc[0]
            ax.plot(sub.median_library_counts, sub.mean_score, marker='o', ms=3.5, lw=1.5, color=COLORS[treatment], label=f'{label}; rho={rho:.2f}')
        ax.set(xscale='log', xlim=(1000, 40000), xticks=[1000, 5000, 10000, 30000], xticklabels=['1k', '5k', '10k', '30k'],
               ylim=(-.20, -.005), xlabel='Median counts in depth decile')
        ax.minorticks_off()
        if number == 0:
            ax.set_ylabel('Mean frozen score')
        ax.legend(frameon=False, fontsize=7.5, loc='lower left')
        panel(ax, chr(ord('C') + number), donor + ': observed deciles')
    fig.subplots_adjust(left=.11, right=.98, bottom=.12, top=.94, hspace=.65, wspace=.38)
    save(fig, 'rev1_figure2_depth_population', [REV / name for name in ['observed_depth_summary.tsv', 'observed_depth_bins.tsv', 'observed_depth_relationships.tsv', 'observed_component_contrasts.tsv']],
         'Contrasts are H + GLP-1 minus vehicle. Expression proxies are not published annotations. Mixed may overlap alpha/beta; Other is unassigned endocrine. Deciles and rho are descriptive.')


def figure3():
    contrasts = read('observed_component_contrasts.tsv')
    native = contrasts[contrasts.population == 'alpha_broad'].pivot(index='donor', columns='component', values='delta_mean').loc[DONORS]
    thinning = read('thinning_stability_summary.tsv')
    thinning = thinning[thinning.population == 'alpha_broad']
    fig = plt.figure(figsize=(6.65, 8.2))
    grid = fig.add_gridspec(3, 2, height_ratios=[1, 1.08, 1.08])
    a, b = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
    x = np.arange(4)
    for offset, component, label, color in [(-.18, 'up_rank_contribution', 'Signed up-rank component', COLORS['up']),
                                             (.18, 'down_rank_contribution', 'Signed down-rank component', COLORS['down'])]:
        a.bar(x + offset, native[component], width=.34, color=color, label=label, zorder=3)
    a.scatter(x, native.score, color='#152D3C', marker='D', s=38, label='Total score contrast', zorder=4)
    a.axhline(0, color='#7C8790', lw=.8)
    a.set(xticks=x, xticklabels=DONORS, ylim=(-.012, .058), ylabel='H + GLP-1 minus vehicle')
    a.legend(frameon=False, fontsize=8.3)
    panel(a, 'A', 'Signed fractional-rank components')
    for offset, component, label, color in [(-.24, 'occupancy_up', 'Up detection / 100', COLORS['up']),
                                             (0, 'occupancy_down', 'Minus down detection / 100', COLORS['down']),
                                             (.24, 'positive_rank_residual', 'Positive-rank remainder', COLORS['residual'])]:
        b.bar(x + offset, native[component], width=.22, color=color, label=label, zorder=3)
    b.scatter(x, native.score, color='#152D3C', marker='D', s=38, zorder=4)
    b.axhline(0, color='#7C8790', lw=.8)
    b.set(xticks=x, xticklabels=DONORS, ylim=(-.009, .045), ylabel='H + GLP-1 minus vehicle')
    b.legend(frameon=False, fontsize=8.3)
    panel(b, 'B', 'Exact detection imbalance plus rank remainder')
    for number, donor in enumerate(DONORS):
        ax = fig.add_subplot(grid[1 + number // 2, number % 2])
        sub = thinning[thinning.donor == donor].sort_values('depth')
        ax.plot(sub.depth, sub.original_same_eligible_delta, marker='s', ms=4, color=COLORS['neutral'], ls='--', label='Native, same eligible cells')
        ax.errorbar(sub.depth, sub.mean_delta, yerr=sub.SD_delta, marker='o', ms=4, color=COLORS['D'], capsize=3, label='True thinning: mean +/- SD')
        ax.axhline(0, color='#7C8790', lw=.8)
        labels = [f'{int(row.depth / 1000)}k\n{int(row.control_eligible_n)}:{int(row.treated_eligible_n)}' for row in sub.itertuples()]
        ax.set(xscale='log', xticks=sub.depth, xticklabels=labels, ylim=(-.010, .032), xlabel='Counts; vehicle:treated n')
        ax.minorticks_off()
        ax.tick_params(axis='x', labelsize=7.3)
        if number == 0:
            ax.set_ylabel('Mean score contrast')
            ax.legend(frameon=False, fontsize=7.4, loc='upper left')
        panel(ax, chr(ord('C') + number), donor + ': without replacement')
    fig.subplots_adjust(left=.11, right=.98, bottom=.13, top=.94, hspace=.65, wspace=.38)
    save(fig, 'rev1_figure3_mechanism_thinning', [REV / 'observed_component_contrasts.tsv', REV / 'thinning_stability_summary.tsv'],
         'Frozen score = (detected up - detected down) / 100 + rank remainder. Five technical draws; fixed population masks; no upward draws or new biological replication.')


def figure4():
    stability = read('benchmark_feature_stability.tsv')
    reference = read('benchmark_reference.tsv')
    primary = read('benchmark_primary_contrasts.tsv')
    fig, axes = plt.subplots(3, 2, figsize=(6.65, 8.4))
    a, b, c = axes[0, 0], axes[0, 1], axes[1, 0]
    for offset, side, color in [(-.17, 'up', COLORS['up']), (.17, 'down', COLORS['down'])]:
        sub = stability[stability.side == side].set_index('program').loc[PROGRAMS[1:]]
        a.bar(np.arange(4) + offset, sub.jaccard_original, width=.32, color=color, label=side.capitalize(), zorder=3)
    a.set(xticks=np.arange(4), xticklabels=['1697', '1738', '1770', '1776'], ylim=(0, 1), xlabel='Training donor omitted', ylabel='Side-specific Jaccard with original')
    a.legend(frameon=False, fontsize=9, ncol=2)
    panel(a, 'A', 'Feature selection stability')
    for ax, field, letter, title, ylabel in [(b, 'native_PM_vs_alpha_AUROC', 'B', 'Native preparation ranking', 'PM vs alpha AUROC'),
                                            (c, 'PM_preparation_detection_fraction', 'C', 'Own alpha-q95 operating point', 'PM-preparation detection')]:
        for j, (program, color) in enumerate(zip(PROGRAMS, PROGRAM_COLORS)):
            sub = reference[reference.program == program].set_index('method').loc[METHODS]
            ax.scatter(np.arange(3) + (j - 2) * .10, sub[field], color=color, s=45 if j == 0 else 27,
                       marker='D' if j == 0 else 'o', zorder=4 if j == 0 else 3)
        ax.set(xticks=np.arange(3), xticklabels=METHOD_LABELS, ylim=(.5, 1) if letter == 'B' else (0, .8), ylabel=ylabel)
        ax.tick_params(axis='x', labelsize=8.3)
        panel(ax, letter, title)
    primary_axes = [axes[1, 1], axes[2, 0], axes[2, 1]]
    for number, method in enumerate(METHODS):
        ax = primary_axes[number]
        for j, (program, label, color) in enumerate(zip(PROGRAMS, PROGRAM_LABELS, PROGRAM_COLORS)):
            sub = primary[(primary.program == program) & (primary.method == method)].set_index('donor').loc[DONORS]
            ax.scatter(np.arange(4) + (j - 2) * .11, sub.delta, color=color, s=45 if j == 0 else 27,
                       marker='D' if j == 0 else 'o', label=label, zorder=4 if j == 0 else 3)
        ax.axhline(0, color='#7C8790', lw=.8)
        ax.set(xticks=np.arange(4), xticklabels=DONORS, ylabel='H + GLP-1 minus vehicle')
        if number == 2:
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        panel(ax, chr(ord('D') + number), METHOD_LABELS[number] + ': donor estimates')
    handles, labels = primary_axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=5, frameon=False, loc='lower center', bbox_to_anchor=(.5, .058), fontsize=9)
    fig.subplots_adjust(left=.11, right=.98, bottom=.16, top=.94, hspace=.65, wspace=.40)
    save(fig, 'rev1_figure4_benchmark', [REV / name for name in ['benchmark_feature_stability.tsv', 'benchmark_reference.tsv', 'benchmark_primary_contrasts.tsv']],
         'LODO is training-feature stability, not independent validation. Each scorer has its own units and q95; all 5 programs and 3 methods are shown without winner selection.')


def figure5():
    native = read('independent_state_benchmark.tsv')
    native = native[native['mode'] == 'original_universe']
    thinning = read('independent_thinning.tsv')
    thinning = thinning[(thinning['mode'] == 'original_universe') & (thinning.program == 'original') & (thinning.method == 'signed_rank')]
    components = read('independent_component_thinning.tsv')
    components = components[(components['mode'] == 'original_universe') & (components.cell_type == 'alpha') & (components.depth == 1000)]
    fig, axes = plt.subplots(3, 2, figsize=(6.65, 8.4))
    state_axes = [axes[0, 0], axes[0, 1], axes[1, 0]]
    lower = np.floor(native.beta_vs_alpha_AUROC.min() * 10) / 10
    for number, (ax, method) in enumerate(zip(state_axes, METHODS)):
        for j, (program, color, label) in enumerate(zip(PROGRAMS, PROGRAM_COLORS, PROGRAM_LABELS)):
            sub = native[(native.program == program) & (native.method == method)].set_index('donor').loc[[1, 2, 3, 4]]
            ax.scatter(np.arange(4) + (j - 2) * .11, sub.beta_vs_alpha_AUROC, color=color,
                       marker='D' if j == 0 else 'o', s=45 if j == 0 else 27, label=label, zorder=3)
        ax.set(xticks=np.arange(4), xticklabels=['1', '2', '3', '4'], ylim=(lower, 1.005),
               xlabel='Baron donor', ylabel='Author beta vs alpha state AUROC')
        panel(ax, chr(ord('A') + number), METHOD_LABELS[number] + ': independent state')
    for ax, cell_type, letter in [(axes[1, 1], 'alpha', 'D'), (axes[2, 0], 'beta', 'E')]:
        for donor, color, marker in zip([1, 2, 3, 4], DONOR_COLORS, ['o', 's', '^', 'D']):
            sub = thinning[(thinning.donor == donor) & (thinning.cell_type == cell_type)].sort_values('depth')
            ax.errorbar(sub.depth, sub.thinning_change, yerr=sub.technical_repeat_SD, marker=marker, color=color,
                        ms=3.3, lw=1.3, capsize=2, label='Donor ' + str(donor))
        ax.axhline(0, color='#7C8790', lw=.8)
        ax.set(xscale='log', xticks=[1000, 2000, 5000, 10000], xticklabels=['1k', '2k', '5k', '10k'],
               xlabel='Target counts without replacement', ylabel='Thinned minus native score')
        ax.minorticks_off()
        if cell_type == 'alpha':
            ax.legend(frameon=False, ncol=2, fontsize=8, loc='upper right')
        panel(ax, letter, cell_type.capitalize() + ': original rank, same cells')
    ax = axes[2, 1]
    values = components.pivot(index='donor', columns='component', values='thinning_change').loc[[1, 2, 3, 4]]
    for offset, name, label, color in [(-.24, 'occupancy_up', 'Up detection', COLORS['up']),
                                      (0, 'occupancy_down', 'Minus down detection', COLORS['down']),
                                      (.24, 'positive_rank_residual', 'Rank remainder', COLORS['residual'])]:
        ax.bar(np.arange(4) + offset, values[name], width=.22, color=color, label=label, zorder=3)
    ax.scatter(np.arange(4), values.score, marker='D', color='#172D3D', s=32, zorder=4)
    ax.axhline(0, color='#7C8790', lw=.8)
    ax.set(xticks=np.arange(4), xticklabels=['1', '2', '3', '4'], ylim=(-.01, .085), xlabel='Baron donor', ylabel='Thinned minus native component')
    ax.legend(frameon=False, fontsize=7.5, loc='upper right')
    panel(ax, 'F', 'Alpha: components at 1000 counts')
    handles, labels = state_axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=5, frameon=False, loc='lower center', bbox_to_anchor=(.5, .058), fontsize=9)
    fig.subplots_adjust(left=.11, right=.98, bottom=.17, top=.94, hspace=.68, wspace=.40)
    save(fig, 'rev1_figure5_independent', [REV / name for name in ['independent_state_benchmark.tsv', 'independent_thinning.tsv', 'independent_component_thinning.tsv', 'independent_feature_coverage.tsv']],
         'State discrimination is not conversion validation. Original universe: 50/50 up, 49/50 down represented; absent features zero-filled. Five draws on the same eligible cells; points omitted if n < 20.')


def supplement_nuisance():
    original = read('countermodel_stability_summary.tsv', revision=False)
    original = original[original.contrast == 'H_GLP-D']
    precision = read('nuisance_precision_summary.tsv')
    drift = read('nuisance_pool_drift.tsv')
    coupling = read('nuisance_coupling_differences.tsv')
    final = precision[(precision.pool_n == 6400) & (precision.group_draws_per_batch == 1000)].copy()
    fig, axes = plt.subplots(3, 2, figsize=(6.65, 8.5))
    a, b, c, d, e, f = axes.ravel()
    x = np.arange(4)
    families = [(['beta_RNA_0', 'beta_RNA_0.01', 'beta_RNA_0.05', 'beta_RNA_0.1'], '0-10% beta RNA', COLORS['H_GLP']),
                (['doublet_fraction_0.01', 'doublet_fraction_0.05', 'doublet_fraction_0.1'], '1-10% doublets', COLORS['up']),
                ([name for name in original.model.unique() if name.startswith('cycle_')], 'Cycling weights 0-1', COLORS['down'])]
    for offset, (models, label, color) in zip([-.18, 0, .18], families):
        values = original[original.model.isin(models)].groupby('donor').pooled_MC_q95.max().loc[DONORS]
        a.scatter(x + offset, values, color=color, marker='o', s=27, label=label, zorder=3)
    observed = original.groupby('donor').observed_delta.first().loc[DONORS]
    a.scatter(x, observed, marker='D', color='#172D3D', s=42, label='Observed native contrast', zorder=4)
    a.set(xticks=x, xticklabels=DONORS, ylim=(0, .046), ylabel='Maximum pooled q95 / native contrast')
    a.legend(frameon=False, fontsize=7.5, loc='upper left')
    panel(a, 'A', 'Original, fixed family ranges')
    coupling_names = [('independent_profile_depth', 'Independent profile-depth', COLORS['neutral']),
                      ('paired_profile_own_depth', 'Paired profile-own depth', COLORS['up'])]
    for offset, (key, label, color) in zip([-.10, .10], coupling_names):
        sub = final[final.coupling == key].groupby('donor').pooled_q95.max().loc[DONORS]
        b.scatter(x + offset, sub, color=color, marker='s', s=30, label=label, zorder=3)
    b.scatter(x, observed, marker='D', color='#172D3D', s=42, label='Observed native contrast', zorder=4)
    b.set(xticks=x, xticklabels=DONORS, ylim=(0, .046), ylabel='Maximum of 5 scenario q95 values')
    b.legend(frameon=False, fontsize=7.2, loc='upper left')
    panel(b, 'B', 'Revised representative scenarios')
    maximum_drift = drift.absolute_q95_difference.max() * 1.12
    for ax, (key, title, _) in zip([c, d], coupling_names):
        subset = drift[drift.coupling == key]
        for donor, color in zip(DONORS, DONOR_COLORS):
            for draws, style, marker in [(200, ':', 's'), (1000, '-', 'o')]:
                sub = subset[(subset.donor == donor) & (subset.group_draws_per_batch == draws)].groupby('pool_n').absolute_q95_difference.max().sort_index()
                ax.plot(sub.index, sub.to_numpy(), color=color, ls=style, marker=marker, ms=3, lw=1.2,
                        label=donor if draws == 1000 else None)
        ax.set(xscale='log', xticks=[400, 1600, 6400], xticklabels=['400', '1600', '6400'],
               ylim=(0, maximum_drift), xlabel='Synthetic cells per pool', ylabel='Largest absolute q95 drift')
        ax.minorticks_off()
        ax.text(.97, .96, '200 draws: dotted\n1000 draws: solid', transform=ax.transAxes, ha='right', va='top', fontsize=7)
        if ax is c:
            ax.legend(frameon=False, ncol=2, fontsize=7.2, loc='upper left')
        panel(ax, 'C' if ax is c else 'D', title)
    models = ['beta_RNA_0', 'beta_RNA_0.1', 'beta_RNA_0.2', 'doublet_fraction_0.1', 'cycle_high_vs_low']
    labels = ['Beta\n0%', 'Beta\n10%', 'Beta\n20%', 'Doublet\n10%', 'Cycle\nhigh-low']
    comparisons = coupling[(coupling.pool_n == 6400) & (coupling.group_draws_per_batch == 1000)]
    for donor, offset, color, marker in zip(DONORS, [-.21, -.07, .07, .21], DONOR_COLORS, ['o', 's', '^', 'D']):
        values = comparisons[comparisons.donor == donor].set_index('model').loc[models].paired_minus_independent_q95
        e.scatter(np.arange(5) + offset, values, color=color, marker=marker, s=25, zorder=3)
    e.axhline(0, color='#7C8790', lw=.8)
    e.set(xticks=np.arange(5), xticklabels=labels, ylabel='Paired minus independent q95')
    panel(e, 'E', 'Dependence assumptions change q95')
    final['batch_q95_range'] = final.batch_q95_max - final.batch_q95_min
    for offset, (key, label, color) in zip([-.18, .18], coupling_names):
        spread = final[final.coupling == key].groupby('donor').batch_q95_range.max().loc[DONORS]
        f.bar(x + offset, spread, width=.33, color=color, label=label, zorder=3)
    f.set(xticks=x, xticklabels=DONORS, ylabel='Largest three-batch q95 range')
    panel(f, 'F', 'Residual Monte Carlo variability')
    fig.subplots_adjust(left=.11, right=.98, bottom=.13, top=.94, hspace=.65, wspace=.40)
    save(fig, 'rev1_supplement_nuisance', [EV / 'countermodel_stability_summary.tsv'] + [REV / name for name in ['nuisance_precision_summary.tsv', 'nuisance_pool_drift.tsv', 'nuisance_coupling_differences.tsv']],
         'B,E,F: pools 6400, 1000 group draws per batch, 3 batches. C,D: maximum over the 5 scenarios relative to that setting. Conditional MC diagnostics are not calibrated contamination bounds or convergence guarantees.')


def main():
    global plt
    runtime = data_path('revision_01/figure_runtime', for_write=True)
    os.environ['MPLCONFIGDIR'] = str(runtime)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plotting
    plt = plotting
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.titlesize': 10.5,
                         'axes.labelsize': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.axisbelow': True, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    for function in [figure1, figure2, figure3, figure4, figure5, supplement_nuisance]:
        function()
    status = {'updated_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'figures': records,
              'pending': [],
              'scope': 'Plots of retained tables; no new fitted scores or statistical tests'}
    (REV / 'figures_status.json').write_text(json.dumps(status, indent=2) + '\n', encoding='utf-8')


if __name__=='__main__':
    import argparse
    from public_runtime import execute_stage
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();execute_stage('figures',args.data_root,args.output_dir)
