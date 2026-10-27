"""Narrow numerical comparison for native signed-log rank statistics.

This module never changes scores or published results. Positive radii certify
only nonempty, integer-count cancellation: nd*up_count(c)=nu*down_count(c)
for every positive count c. Equal-count terms have a common evaluated log value,
so their mathematical weighted sum is zero even if the log implementation varies.
Other cells, including absent signatures, receive radius zero.

The budget gamma_(m+4)*L1/(1-gamma_(m+4)) covers float64 weight/product and dot
rounding, with slack for the nonnegative L1 reduction. It assumes ordinary IEEE
float64 arithmetic without overflow/underflow. It is not a biological uncertainty
interval or a general tolerance for arbitrary log-score changes.
"""
from collections import Counter
import math

import numpy as np
from scipy import sparse
from scipy.stats import rankdata, spearmanr


def log_score_error_radii(term_l1, nonzero_terms, exact_cancellation):
    """Return zero-only roundoff budgets; caller supplies proven cancellation."""
    l1, count, cancel = np.broadcast_arrays(np.asarray(term_l1, dtype=float),
                                           np.asarray(nonzero_terms),
                                           np.asarray(exact_cancellation, dtype=bool))
    if (not np.isfinite(l1).all() or (l1 < 0).any() or
            not np.isfinite(count).all() or (count < 0).any() or
            (count != np.floor(count)).any() or ((count == 0) & (l1 != 0)).any()):
        raise ValueError('Require finite nonnegative L1 and integer term counts')
    unit = np.finfo(np.float64).eps / 2
    steps = count.astype(float) + 4
    if (steps * unit >= .5).any():
        raise ValueError('Unsupported accumulation length')
    gamma = steps * unit / (1 - steps * unit)
    radius = gamma * l1 / (1 - gamma)
    active = cancel & (count > 0) & (l1 > 0)
    return np.where(active, np.nextafter(radius, np.inf), 0.)


def signed_log_error_radii(x, totals, genes, programs):
    """Compute (program, cell) certified-zero radii from the native count matrix.

    Input order and represented block denominators match score_methods. Raw
    integer counts and unique gene labels are required for the count proof.
    """
    x = sparse.csc_matrix(x, copy=True)
    x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
    genes = list(genes); totals = np.asarray(totals, dtype=float)
    if len(genes) != x.shape[0] or len(set(genes)) != len(genes):
        raise ValueError('Unique gene labels must match matrix rows')
    if (totals.shape != (x.shape[1],) or not np.isfinite(totals).all() or
            (totals < 0).any() or not np.isfinite(x.data).all() or
            (x.data < 0).any() or (x.data != np.floor(x.data)).any()):
        raise ValueError('Finite nonnegative integer counts and valid totals required')
    lookup = {g:i for i,g in enumerate(genes)}
    up = np.zeros((len(programs), len(genes)), dtype=bool)
    down = np.zeros_like(up)
    for pi, program in enumerate(programs.values()):
        u, d = program['up'], program['down']
        if len(set(u)) != len(u) or len(set(d)) != len(d) or set(u) & set(d):
            raise ValueError('Signature blocks must be unique and disjoint')
        up[pi, [lookup[g] for g in u if g in lookup]] = True
        down[pi, [lookup[g] for g in d if g in lookup]] = True
    nu, nd = up.sum(1), down.sum(1)
    if (nu == 0).any() or (nd == 0).any():
        raise ValueError('Both represented signature blocks must be nonempty')
    weights = up / nu[:, None] - down / nd[:, None]
    result = np.zeros((len(programs), x.shape[1]))
    for j in range(x.shape[1]):
        a, b = x.indptr[j:j+2]; idx, values = x.indices[a:b], x.data[a:b]
        if len(values) and totals[j] <= 0:
            raise ValueError('Positive counts require a positive whole-library total')
        if len(values) and np.max(values) > totals[j]:
            raise ValueError('A gene count cannot exceed the whole-library total')
        logs = np.log1p(10000 * values / max(float(totals[j]), 1))
        if not np.isfinite(logs).all():
            raise ValueError('Nonfinite normalized log term')
        for pi in range(len(programs)):
            um, dm = up[pi, idx], down[pi, idx]
            active = um | dm; m = int(active.sum())
            if not m:
                continue
            cu, cd = Counter(values[um].tolist()), Counter(values[dm].tolist())
            if any(int(nd[pi])*cu[c] != int(nu[pi])*cd[c] for c in cu.keys() | cd.keys()):
                continue
            # A platform with lane-dependent unequal evaluations for the same
            # argument is outside this common-term proof and fails explicitly.
            for c in cu.keys() | cd.keys():
                common = logs[active & (values == c)]
                if np.any(common != common[0]):
                    raise ValueError('Equal-count log values differ within a cell')
            l1 = math.fsum(abs(float(w)*float(v)) for w,v in zip(weights[pi,idx][active], logs[active]))
            result[pi,j] = log_score_error_radii(l1, m, True)
    return result


def _clusters(scores, radii):
    scores, radii = np.asarray(scores, dtype=float), np.asarray(radii, dtype=float)
    if (scores.ndim != 1 or radii.shape != scores.shape or not len(scores) or
            not np.isfinite(scores).all() or not np.isfinite(radii).all() or (radii < 0).any()):
        raise ValueError('Require equally sized finite score/radius vectors')
    uncertain = radii > 0
    if (np.abs(scores[uncertain]) > radii[uncertain]).any():
        raise ValueError('Fresh score exceeds its certified-zero roundoff budget')
    low = np.where(uncertain, -radii, scores)
    high = np.where(uncertain, radii, scores)
    order = np.argsort(low, kind='stable')
    groups, current, right = [], [], -np.inf
    for index in order:
        if current and low[index] > right:
            groups.append(np.asarray(current, dtype=int)); current = []
        current.append(int(index)); right = max(right, high[index]) if len(current) > 1 else high[index]
    groups.append(np.asarray(current, dtype=int))
    canonical = np.empty(len(scores)); start = 1
    variable = []
    for group in groups:
        canonical[group] = start + (len(group)-1)/2
        start += len(group)
        variable.append(bool(uncertain[group].any() and len(group) > 1))
    return scores, canonical, groups, variable


def spearman_envelope(scores, radii, depth):
    """Conservative rank-correlation interval over certified-zero tie splits.

    Canonical cluster midranks c and perturbations e satisfy c·e=0 after
    centering. Rearrangement bounds e·rank(depth); full within-cluster splitting
    bounds ||e||². Fixed-only clusters remain tied. Within an uncertain cluster,
    the full-splitting relaxation may overbound its fixed-point tie constraints.
    """
    scores, canonical, groups, variable = _clusters(scores, radii)
    depth = np.asarray(depth, dtype=float)
    if depth.shape != scores.shape or not np.isfinite(depth).all():
        raise ValueError('Finite depth vector must match scores')
    y = rankdata(depth, method='average'); y -= y.mean()
    c = canonical - canonical.mean()
    v0, vy = float(c @ c), float(y @ y)
    if v0 == 0 or vy == 0:
        raise ValueError('Degenerate rank variance cannot certify a metric mismatch')
    adjustment = 0.; extra_variance = 0.
    for group, changes in zip(groups, variable):
        if not changes:
            continue
        m = len(group)
        deviation = np.arange(1, m+1) - (m+1)/2
        adjustment += float(deviation @ np.sort(y[group]))
        extra_variance += m*(m*m-1)/12
    numerator = float(c @ y)
    minimum_denominator = math.sqrt(v0*vy)
    maximum_denominator = math.sqrt((v0+extra_variance)*vy)
    lo_num, hi_num = numerator-adjustment, numerator+adjustment
    lower = lo_num/(maximum_denominator if lo_num >= 0 else minimum_denominator)
    upper = hi_num/(minimum_denominator if hi_num >= 0 else maximum_denominator)
    return {'lower':max(-1., lower), 'upper':min(1., upper),
            'canonical':numerator/minimum_denominator,
            'uncertain_clusters':sum(variable), 'uncertain_cluster_cells':sum(len(g) for g,v in zip(groups,variable) if v),
            'covariance_adjustment_bound':adjustment, 'rank_variance_increment_bound':extra_variance}


def empirical_auc(scores, positive):
    scores, positive = np.asarray(scores, dtype=float), np.asarray(positive, dtype=bool)
    if scores.ndim != 1 or positive.shape != scores.shape or not np.isfinite(scores).all():
        raise ValueError('Finite scores and matching state labels required')
    p, n = int(positive.sum()), int((~positive).sum())
    if not p or not n:
        raise ValueError('Both states required')
    ranks = rankdata(scores, method='average')
    return float((ranks[positive].sum()-p*(p+1)/2)/(p*n))


def auc_envelope(scores, radii, positive):
    scores, canonical, groups, variable = _clusters(scores, radii)
    positive = np.asarray(positive, dtype=bool)
    if positive.shape != scores.shape:
        raise ValueError('State labels must match scores')
    p, n = int(positive.sum()), int((~positive).sum())
    if not p or not n:
        raise ValueError('Both states required')
    center = empirical_auc(canonical, positive)
    pairs = sum(int(positive[g].sum())*int((~positive[g]).sum()) for g,v in zip(groups,variable) if v)
    delta = .5*pairs/(p*n)
    return {'lower':max(0., center-delta), 'upper':min(1., center+delta),
            'canonical':center, 'uncertain_cross_state_pairs':pairs,
            'uncertain_clusters':sum(variable), 'uncertain_cluster_cells':sum(len(g) for g,v in zip(groups,variable) if v)}


def compare_rank_metric(retained, fresh, scores, radii, *, metric, covariate=None,
                        positive=None, rtol=1e-12, atol=2e-14):
    """Validate a narrowly eligible metric mismatch; never change either value.

    Integration must restrict this to native independent signed-log AUROC or
    score/depth Spearman columns. All other fields retain their usual comparison.
    """
    if metric == 'spearman':
        recomputed = float(spearmanr(scores, covariate).statistic)
        envelope = lambda: spearman_envelope(scores, radii, covariate)
    elif metric == 'auc':
        recomputed = empirical_auc(scores, positive)
        envelope = lambda: auc_envelope(scores, radii, positive)
    else:
        raise ValueError('Only spearman and auc metrics are eligible')
    result = {'accepted':False, 'metric':metric, 'retained':float(retained),
              'fresh':float(fresh), 'recomputed':recomputed}
    if not (np.isfinite(retained) and np.isfinite(fresh) and np.isfinite(recomputed)):
        return {**result, 'reason':'nonfinite_metric_requires_ordinary_comparison'}
    if not np.isclose(fresh, recomputed, rtol=rtol, atol=atol):
        return {**result, 'reason':'fresh_metric_does_not_recompute'}
    try:
        bounds = envelope()
    except ValueError as error:
        return {**result, 'reason':str(error)}
    tolerance = atol + rtol*abs(float(retained))
    fresh_tolerance = atol + rtol*abs(float(fresh))
    if not bounds['lower']-fresh_tolerance <= fresh <= bounds['upper']+fresh_tolerance:
        return {**result, **bounds, 'reason':'fresh_metric_outside_zero_roundoff_envelope'}
    accepted = bounds['lower']-tolerance <= retained <= bounds['upper']+tolerance
    return {**result, **bounds, 'accepted':bool(accepted),
            'reason':'certified_zero_roundoff_rank_envelope' if accepted else 'retained_metric_outside_zero_roundoff_envelope'}
