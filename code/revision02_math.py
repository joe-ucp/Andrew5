"""Pure decision and baseline diagnostics; the original three scorers are unchanged."""
import numpy as np
from scipy import sparse
from scipy.stats import rankdata
from revision_benchmark import score_methods
from reference_analysis import auc

METHODS = ['signed_rank', 'signed_log1p_CPT10K', 'ucell_form_1500',
           'unclipped_ucell_difference']


def zero_centered_components(comp, universe_n):
    """Remove the common rank of zeros from separately plotted signed blocks."""
    if universe_n <= 0:
        raise ValueError('Positive universe size required')
    result = dict(comp)
    result['zero_rank_baseline'] = (universe_n - np.asarray(comp['k_positive'])) / (2 * universe_n)
    result['centered_up'] = np.asarray(comp['occupancy_up']) + np.asarray(comp['positive_rank_up'])
    result['centered_down'] = np.asarray(comp['occupancy_down']) + np.asarray(comp['positive_rank_down'])
    for observed, expected in [
        (comp['up_rank_contribution'], result['zero_rank_baseline'] + result['centered_up']),
        (comp['down_rank_contribution'], -result['zero_rank_baseline'] + result['centered_down']),
        (comp['score'], result['centered_up'] + result['centered_down']),
    ]:
        if not np.allclose(observed, expected, rtol=0, atol=2e-15):
            raise AssertionError('Zero-rank baseline identity failed')
    return result


def unclipped_ucell(x, genes, programs, rank_cap=1500):
    """Exact block difference used by the locked comparator before its final clip.

    This is an explanatory diagnostic, not a claim of package equivalence.
    Missing program features follow the old scorer's represented-block policy.
    """
    x = sparse.csc_matrix(x, copy=True)
    x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
    if not np.isfinite(x.data).all() or (x.data < 0).any():
        raise ValueError('Finite nonnegative counts required')
    G, N = x.shape
    if len(genes) != G or len(set(genes)) != G:
        raise ValueError('Unique ordered universe genes required')
    lookup = {g: i for i, g in enumerate(genes)}
    up = np.zeros((len(programs), G), dtype=np.int64)
    down = up.copy()
    for i, program in enumerate(programs.values()):
        for target, key in [(up, 'up'), (down, 'down')]:
            target[i, [lookup[g] for g in program[key] if g in lookup]] = 1
    nu, nd = up.sum(1), down.sum(1)
    if (nu == 0).any() or (nd == 0).any():
        raise ValueError('Nonempty represented blocks required')
    valid = (G >= rank_cap) & (nu < rank_cap) & (nd < rank_cap)
    result = np.empty((len(programs), N))
    for j in range(N):
        a, b = x.indptr[j:j+2]
        idx, v = x.indices[a:b], x.data[a:b]
        K = len(v)
        twice = np.rint(2 * rankdata(v, method='average')).astype(np.int64)
        baseline = min((G + K + 1) / 2, rank_cap)
        delta = np.minimum(K + 1 - twice / 2, rank_cap) - baseline
        rank_up = nu * baseline + up[:, idx] @ delta
        rank_down = nd * baseline + down[:, idx] @ delta
        u = (nu * rank_cap - rank_up) / (nu * rank_cap - nu * (nu + 1) / 2)
        d = (nd * rank_cap - rank_down) / (nd * rank_cap - nd * (nd + 1) / 2)
        result[:, j] = np.where(valid, u-d, np.nan)
    return result


def score_bundle(x, totals, genes, programs):
    locked = score_methods(x, totals, genes, programs)
    raw = unclipped_ucell(x, genes, programs)
    if not np.allclose(np.maximum(raw, 0), locked[:, 2], rtol=0, atol=0, equal_nan=True):
        raise AssertionError('Supplementary diagnostic changed the locked clip')
    return np.concatenate([locked, raw[:, None]], axis=1)


def q95(values, minimum=20):
    values = np.asarray(values, dtype=float)
    if len(values) < minimum or not np.isfinite(values).all():
        return np.nan
    return float(np.quantile(values, .95, method='linear'))


def decision_metrics(values, threshold, native=None, native_threshold=None, minimum=20):
    """Strict threshold, ties retained, paired transitions on identical cells.

    Counts always survive small-group suppression. A recalibrated analysis may
    compare different native and thinned thresholds; callers must label that.
    negative_n is the nonpositive decision (<= threshold); tie_n is a subset.
    Transition fractions use all paired cells, not only cells at risk, as denominator.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or not np.isfinite(values).all():
        raise ValueError('One-dimensional finite score vector required')
    n = len(values)
    result = {'n': n, 'threshold': float(threshold), 'evaluable': n >= minimum and np.isfinite(threshold)}
    result['reason'] = '' if result['evaluable'] else ('threshold_unavailable' if not np.isfinite(threshold) else 'fewer_than_20_cells')
    for key in ['positive_n', 'negative_n', 'tie_n', 'positive_fraction', 'tie_fraction',
                'native_positive_n', 'negative_to_positive_n', 'positive_to_negative_n',
                'negative_to_positive_fraction', 'positive_to_negative_fraction']:
        result[key] = np.nan
    if not np.isfinite(threshold):
        return result
    positive = values > threshold
    result.update(positive_n=int(positive.sum()), negative_n=int((~positive).sum()),
                  tie_n=int((values == threshold).sum()))
    if result['evaluable']:
        result.update(positive_fraction=float(positive.mean()), tie_fraction=float((values == threshold).mean()))
    if native is not None:
        native = np.asarray(native, dtype=float)
        if native.shape != values.shape or not np.isfinite(native).all():
            raise ValueError('Finite paired native scores on the same cells required')
        native_threshold = threshold if native_threshold is None else native_threshold
        if np.isfinite(native_threshold):
            prior = native > native_threshold
            new = int((~prior & positive).sum()); lost = int((prior & ~positive).sum())
            result.update(native_positive_n=int(prior.sum()), negative_to_positive_n=new, positive_to_negative_n=lost)
            if result['evaluable']:
                result.update(negative_to_positive_fraction=new/n, positive_to_negative_fraction=lost/n)
    return result


def state_auc(scores, states, minimum=20):
    scores, states = np.asarray(scores), np.asarray(states)
    a, b = scores[states == 'alpha'], scores[states == 'beta']
    return auc(b, a) if len(a) >= minimum and len(b) >= minimum else np.nan
