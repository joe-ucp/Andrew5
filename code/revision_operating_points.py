"""Revision 01: independent score algebra and retained-result operating points.

This does not alter the frozen program or run a new treatment benchmark.
Empirical cell-level fractions describe preparations, not verified conversion.
Full ROC coordinates go to the external D: data root; compact summaries stay here.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from data_paths import data_path, resolve_data_root, evidence_path, logical_evidence_key

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "evidence"


def _vectors(pos, neg):
    pos, neg = np.asarray(pos, dtype=float), np.asarray(neg, dtype=float)
    if pos.ndim != 1 or neg.ndim != 1 or not len(pos) or not len(neg):
        raise ValueError("Require nonempty one-dimensional positive and comparator arrays")
    if not np.isfinite(pos).all() or not np.isfinite(neg).all():
        raise ValueError("Scores must be finite")
    return pos, neg


def rank_score_components(x, up, down):
    """Exact occupancy plus centered-positive-rank decomposition for zero ties.

    For G genes and k positive counts, average positive-only rank r+:
      s = .5*(kU/m-kD/n)
          + [sum_U+(r+-(k+1)/2)/m - sum_D+(r+-(k+1)/2)/n]/G.
    The second term is retained, never assumed negligible.
    """
    x = np.asarray(x)
    up, down = np.asarray(up, dtype=int), np.asarray(down, dtype=int)
    if x.ndim != 2 or x.shape[0] == 0 or not np.isfinite(x).all() or (x < 0).any():
        raise ValueError("Require finite nonnegative genes-by-cells array")
    if not len(up) or not len(down) or len(np.unique(up)) != len(up) or len(np.unique(down)) != len(down):
        raise ValueError("Require nonempty unique up and down blocks")
    if set(up) & set(down) or min(up.min(), down.min()) < 0 or max(up.max(), down.max()) >= x.shape[0]:
        raise ValueError("Blocks must be disjoint valid indices")
    G, cells = x.shape
    m, n = len(up), len(down)
    detected = x > 0
    k = detected.sum(axis=0)
    ku, kd = detected[up].sum(axis=0), detected[down].sum(axis=0)
    occupancy = .5 * (ku / m - kd / n)
    residual = np.zeros(cells)
    for j in range(cells):
        positive = np.flatnonzero(detected[:, j])
        if not len(positive):
            continue
        centered = np.zeros(G)
        centered[positive] = rankdata(x[positive, j], method="average") - (k[j] + 1) / 2
        residual[j] = (centered[up].sum() / m - centered[down].sum() / n) / G
    ranks = rankdata(x, method="average", axis=0)
    # General definition, independently evaluated without production integer shortcut.
    score = (ranks[up].mean(axis=0) - ranks[down].mean(axis=0)) / G
    bound = np.maximum(k - 1, 0) / (2 * G) * (ku / m + kd / n)
    return dict(score=score, occupancy_component=occupancy,
                positive_rank_residual=residual, k_positive=k, k_up=ku,
                k_down=kd, residual_bound=bound)


def empirical_auc(pos, neg):
    pos, neg = _vectors(pos, neg)
    ranks = rankdata(np.r_[pos, neg], method="average")
    return float((ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) /
                 (len(pos) * len(neg)))


def threshold_metrics(pos, neg, threshold):
    """Strict greater-than rule with ties disclosed; no random tie breaking."""
    pos, neg = _vectors(pos, neg)
    threshold = float(threshold)
    if np.isnan(threshold):
        raise ValueError("Threshold cannot be NaN")
    return dict(threshold=threshold, positive_preparation_n=len(pos),
                comparator_n=len(neg),
                positive_preparation_detection_fraction=float(np.mean(pos > threshold)),
                achieved_comparator_FPR=float(np.mean(neg > threshold)),
                positive_tie_fraction=float(np.mean(pos == threshold)),
                comparator_tie_fraction=float(np.mean(neg == threshold)),
                positive_detection_if_greater_equal=float(np.mean(pos >= threshold)),
                comparator_FPR_if_greater_equal=float(np.mean(neg >= threshold)))


def strict_roc(pos, neg):
    """All attainable deterministic score>threshold states, including all-positive."""
    pos, neg = _vectors(pos, neg)
    thresholds = np.r_[-np.inf, np.unique(np.r_[pos, neg])]
    ps, ns = np.sort(pos), np.sort(neg)
    positive = len(pos) - np.searchsorted(ps, thresholds, side="right")
    negative = len(neg) - np.searchsorted(ns, thresholds, side="right")
    return pd.DataFrame(dict(threshold=thresholds,
                             positive_preparation_detection_fraction=positive / len(pos),
                             achieved_comparator_FPR=negative / len(neg),
                             positive_n_above=positive, comparator_n_above=negative))


def spearman_descriptive(x, y):
    """Observed association only; intentionally omit cell-level p-values."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if len(x) < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan
    return float(spearmanr(x, y).statistic)


def main():
    external = resolve_data_root(require_exists=True)
    output = EVIDENCE / "revision_01"
    bulk = data_path("revision_01/operating_points", external, for_write=True)
    output.mkdir(exist_ok=True)
    bulk.mkdir(parents=True, exist_ok=True)
    paths = [evidence_path(name, data_root=external, must_exist=True) for name in ["reference_cell_scores.tsv",
             "calibration_simulations.tsv.gz", "primary_cell_scores.tsv.gz",
             "frozen_score.json"]]
    before = {logical_evidence_key(p, data_root=external): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    native = pd.read_csv(paths[0], sep="\t", float_precision="round_trip")
    sim = pd.read_csv(paths[1], sep="\t", float_precision="round_trip")
    primary = pd.read_csv(paths[2], sep="\t", float_precision="round_trip")
    cfg = json.loads(paths[3].read_text())
    cohorts = []
    pos = native[native.preparation == "PM"].score.to_numpy()
    for name, labels in [("alpha_only", ["alpha"]), ("beta_only", ["beta"]),
                         ("alpha_and_beta", ["alpha", "beta"])]:
        neg = native[native.preparation.isin(labels)].score.to_numpy()
        cohorts.append(("native", name, pos, neg))
    for depth, part in sim.groupby("depth"):
        pos = part[part.model == "PM"].score.to_numpy()
        ranges = {"alpha_only": ["alpha"],
                  "low_mixtures": ["alpha", "mixture_0", "mixture_0.01", "mixture_0.05", "mixture_0.1"],
                  "stress_inclusive": [x for x in part.model.unique() if x != "PM"]}
        for name, models in ranges.items():
            cohorts.append((str(depth), name, pos, part[part.model.isin(models)].score.to_numpy()))
    points, roc_summaries, curves = [], [], []
    for depth, comparator, pos, neg in cohorts:
        row = threshold_metrics(pos, neg, np.quantile(neg, .95, method="linear"))
        row.update(depth=depth, comparator=comparator, threshold_rule="empirical q95 linear; score strictly greater",
                   empirical_AUROC=empirical_auc(pos, neg),
                   positive_exact_zero_fraction=float(np.mean(pos == 0)),
                   comparator_exact_zero_fraction=float(np.mean(neg == 0)),
                   scope="captured preparation cells/simulated profiles; not cell-level true conversion")
        points.append(row)
        roc = strict_roc(pos, neg)
        roc.insert(0, "comparator", comparator)
        roc.insert(0, "depth", depth)
        curves.append(roc)
        for cap in [0., .01, .05, .10, .20]:
            feasible = roc[roc.achieved_comparator_FPR <= cap]
            best = feasible.sort_values(["positive_preparation_detection_fraction", "achieved_comparator_FPR", "threshold"],
                                        ascending=[False, True, True]).iloc[0]
            roc_summaries.append(dict(depth=depth, comparator=comparator, empirical_FPR_cap=cap,
                                     threshold=float(best.threshold),
                                     achieved_comparator_FPR=float(best.achieved_comparator_FPR),
                                     maximum_attainable_positive_preparation_detection_fraction=float(best.positive_preparation_detection_fraction),
                                     scope="descriptive empirical ROC frontier; not held-out threshold performance"))
    pd.DataFrame(points).to_csv(output / "operating_points.tsv", sep="\t", index=False)
    pd.DataFrame(roc_summaries).to_csv(output / "roc_frontiers.tsv", sep="\t", index=False)
    pd.concat(curves, ignore_index=True).to_csv(bulk / "all_roc_coordinates.tsv.gz", sep="\t", index=False,
                                               compression={"method": "gzip", "mtime": 0})
    correlations, depths = [], []
    for prep, group in native.groupby("preparation"):
        correlations.append(dict(dataset="native_reference", grouping="preparation", group=prep,
                                 annotation="", donor="", treatment="", n=len(group), predictor="depth",
                                 outcome="score", spearman=spearman_descriptive(group.depth, group.score)))
        depths.append(dict(dataset="native_reference", group=prep, n=len(group),
                           depth_min=float(group.depth.min()), depth_q25=float(group.depth.quantile(.25)),
                           depth_median=float(group.depth.median()), depth_q75=float(group.depth.quantile(.75)),
                           depth_max=float(group.depth.max()), score_median=float(group.score.median()),
                           score_q25=float(group.score.quantile(.25)), score_q75=float(group.score.quantile(.75))))
    for (annotation, donor, treatment), group in primary.groupby(["annotation", "donor", "treatment"]):
        for predictor in ["library_counts", "detected_genes", "score_genes_detected"]:
            for outcome in ["score", "score_depth2000"]:
                correlations.append(dict(dataset="primary", grouping="within_donor_treatment_annotation",
                                         group=f"{donor}:{treatment}:{annotation}", annotation=annotation, donor=donor,
                                         treatment=treatment, n=len(group), predictor=predictor, outcome=outcome,
                                         spearman=spearman_descriptive(group[predictor], group[outcome])))
        depths.append(dict(dataset="primary", group=f"{donor}:{treatment}:{annotation}", n=len(group),
                           depth_min=float(group.library_counts.min()), depth_q25=float(group.library_counts.quantile(.25)),
                           depth_median=float(group.library_counts.median()), depth_q75=float(group.library_counts.quantile(.75)),
                           depth_max=float(group.library_counts.max()), score_median=float(group.score.median()),
                           score_q25=float(group.score.quantile(.25)), score_q75=float(group.score.quantile(.75))))
    for annotation, group in primary.groupby("annotation"):
        for predictor in ["library_counts", "detected_genes", "score_genes_detected"]:
            correlations.append(dict(dataset="primary", grouping="pooled_descriptive_mixes_donors_and_treatments",
                                     group=annotation, annotation=annotation, donor="pooled", treatment="pooled", n=len(group),
                                     predictor=predictor, outcome="score",
                                     spearman=spearman_descriptive(group[predictor], group.score)))
    pd.DataFrame(correlations).to_csv(output / "observed_spearman.tsv", sep="\t", index=False)
    pd.DataFrame(depths).to_csv(output / "observed_depth_summary.tsv", sep="\t", index=False)
    for p in paths:
        if hashlib.sha256(p.read_bytes()).hexdigest() != before[logical_evidence_key(p, data_root=external)]:
            raise AssertionError(f"Original input changed: {p}")
    status = dict(status="COMPLETED_RETAINED_RESULT_DIAGNOSTICS", original_evidence_unchanged=True,
                  inputs_sha256=before, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  universe_n=len(cfg["universe"]), up_n=len(cfg["up"]), down_n=len(cfg["down"]),
                  score_identity="s=.5*(kU/m-kD/n)+(sum_U+(r+-(k+1)/2)/m-sum_D+(r+-(k+1)/2)/n)/G",
                  approximation="Equal50 blocks: occupancy=(kU-kD)/100; residual is retained, not presumed negligible.",
                  threshold_convention="strict >; empirical q95 uses NumPy linear interpolation; ties separately reported",
                  operating_point_rows=len(points), roc_frontier_rows=len(roc_summaries), correlation_rows=len(correlations),
                  full_roc_path=str(bulk / "all_roc_coordinates.tsv.gz"),
                  inference_limits=["PM is a positive preparation without verified cell-level conversion truth",
                                    "80% criterion is the preserved analyst-defined rule, not a biological law",
                                    "Within-library correlations are observed associations without cell-level p-values",
                                    "ROC thresholds and their empirical performance use the same observations",
                                    "No score tuning, new treatment benchmark, donor-level uncertainty, or causal depth adjustment"])
    (output / "operating_points_status.json").write_text(json.dumps(status, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(status, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    from public_runtime import execute_stage
    execute_stage('operating-points', args.data_root, args.output_dir)
