"""Pure frozen-rank decomposition and finite-library thinning for revision 01.

No file access, feature selection, or fitting occurs in this module.
"""
import numpy as np
from scipy import sparse
from scipy.stats import rankdata


def _blocks(n_genes, up, down):
    up = np.asarray(up, dtype=np.int64)
    down = np.asarray(down, dtype=np.int64)
    if up.ndim != 1 or down.ndim != 1 or not len(up) or not len(down):
        raise ValueError('Score blocks must be nonempty one-dimensional indices.')
    combined = np.r_[up, down]
    if len(np.unique(combined)) != len(combined):
        raise ValueError('Score blocks must be unique and disjoint.')
    if combined.min() < 0 or combined.max() >= n_genes:
        raise ValueError('Score block index is outside the frozen universe.')
    labels = np.zeros(n_genes, dtype=np.int8)
    labels[up] = 1
    labels[down] = -1
    return up, down, labels


def components_sparse(x, up, down):
    """Return exact frozen score and its signed rank/detection components.

    G is the number of universe genes and K the number positive in a cell.
    kU/kD count detected genes in the two frozen blocks of sizes m/n.
    S = .5*(kU/m-kD/n) +
        [sum_U+(r+-(K+1)/2)/m-sum_D+(r+-(K+1)/2)/n]/G.
    Positive ranks r+ use average ties. For equal 50-gene blocks, occupancy
    is exactly (kU-kD)/100. Integer half-rank numerators preserve exact zeros.
    Signed up/down rank contributions use (full-universe rank minus .5)/G;
    their sum equals S. These terms are arithmetic, not biological pathways.
    """
    x = sparse.csc_matrix(x, copy=True)
    x.sum_duplicates()
    x.sort_indices()
    if x.shape[0] == 0 or not np.isfinite(x.data).all() or (x.data < 0).any():
        raise ValueError('The frozen universe must contain finite nonnegative expression.')
    G = x.shape[0]
    up, down, labels = _blocks(G, up, down)
    m, n = len(up), len(down)
    denominator = 2 * G * m * n
    names = ['score', 'up_rank_contribution', 'down_rank_contribution',
             'occupancy_component', 'occupancy_up', 'occupancy_down',
             'positive_rank_residual', 'positive_rank_up', 'positive_rank_down',
             'residual_bound', 'k_positive', 'k_up', 'k_down']
    result = {name: np.empty(x.shape[1], dtype=np.int64 if name.startswith('k_') else float)
              for name in names}
    for cell in range(x.shape[1]):
        start, stop = x.indptr[cell:cell+2]
        values, indices = x.data[start:stop], x.indices[start:stop]
        positive = values > 0
        values, indices = values[positive], indices[positive]
        K = len(values)
        z = G - K
        twice_positive_rank = np.rint(2 * rankdata(values, method='average')).astype(np.int64)
        upper, lower = labels[indices] == 1, labels[indices] == -1
        kU, kD = int(upper.sum()), int(lower.sum())
        increment = twice_positive_rank + z - 1
        centered = twice_positive_rank - K - 1
        up_increment = int(increment[upper].sum())
        down_increment = int(increment[lower].sum())
        score_numerator = up_increment * n - down_increment * m
        occupancy_numerator = G * (kU * n - kD * m)
        remainder_numerator = int(centered[upper].sum()) * n - int(centered[lower].sum()) * m
        if score_numerator != occupancy_numerator + remainder_numerator:
            raise AssertionError('Exact rank-decomposition identity failed.')
        result['score'][cell] = score_numerator / denominator
        result['up_rank_contribution'][cell] = (m * z + up_increment) / (2 * G * m)
        result['down_rank_contribution'][cell] = -(n * z + down_increment) / (2 * G * n)
        result['occupancy_component'][cell] = occupancy_numerator / denominator
        result['occupancy_up'][cell] = kU / (2 * m)
        result['occupancy_down'][cell] = -kD / (2 * n)
        result['positive_rank_residual'][cell] = remainder_numerator / denominator
        result['positive_rank_up'][cell] = int(centered[upper].sum()) / (2 * G * m)
        result['positive_rank_down'][cell] = -int(centered[lower].sum()) / (2 * G * n)
        result['residual_bound'][cell] = max(0, K - 1) / (2 * G) * (kU / m + kD / n)
        result['k_positive'][cell], result['k_up'][cell], result['k_down'][cell] = K, kU, kD
    return result


def _integer_counts(values, label):
    values = np.asarray(values)
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
        raise ValueError(label + ' must contain finite nonnegative integer counts.')
    if (values >= 10**9).any():
        raise ValueError(label + ' exceeds the supported hypergeometric count range.')
    return values.astype(np.int64)


def thin_sparse_counts(x, total_counts, depth, rng):
    """Draw without replacement from each eligible observed whole library.

    x contains frozen-universe counts. The difference between each observed
    full-library total and its universe sum is an outside-universe category.
    Aggregating that category preserves the exact marginal sampling law for
    all universe genes. Callers must pass only cells with total_counts>=depth.
    """
    if isinstance(depth, (bool, np.bool_)) or int(depth) != depth or depth < 0:
        raise ValueError('Thinning depth must be a nonnegative integer.')
    depth = int(depth)
    x = sparse.csc_matrix(x, copy=True)
    x.sum_duplicates()
    x.sort_indices()
    x.data = _integer_counts(x.data, 'Observed expression')
    totals = _integer_counts(total_counts, 'Whole-library totals')
    if totals.shape != (x.shape[1],):
        raise ValueError('One whole-library total is required for each cell.')
    if (totals < depth).any():
        raise ValueError('Upward draws are forbidden: every cell must have N>=depth.')
    values, indices, pointers = [], [], [0]
    for cell in range(x.shape[1]):
        start, stop = x.indptr[cell:cell+2]
        counts, gene_indices = x.data[start:stop], x.indices[start:stop]
        take = counts > 0
        counts, gene_indices = counts[take], gene_indices[take]
        outside = int(totals[cell]) - int(counts.sum())
        if outside < 0:
            raise ValueError('Universe counts exceed the observed whole-library total.')
        categories = np.r_[counts, outside]
        drawn = (categories.copy() if depth == totals[cell] else
                 rng.multivariate_hypergeometric(categories, depth, method='marginals'))
        if int(drawn.sum()) != depth or (drawn < 0).any() or (drawn > categories).any():
            raise AssertionError('Without-replacement thinning invariant failed.')
        selected = drawn[:-1] > 0
        values.append(drawn[:-1][selected])
        indices.append(gene_indices[selected])
        pointers.append(pointers[-1] + int(selected.sum()))
    return sparse.csc_matrix((np.concatenate(values) if values else np.array([], dtype=np.int64),
                              np.concatenate(indices) if indices else np.array([], dtype=np.int64),
                              np.asarray(pointers)), shape=x.shape)
