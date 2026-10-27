"""Scientific invariants using only synthetic fixtures; no external datasets."""
from pathlib import Path
import sys
import unittest
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from reference_analysis import scores_dense
from primary_analysis import score_sparse, map_universe, prob_profiles
from revision_operating_points import rank_score_components, threshold_metrics, strict_roc, empirical_auc
from revision_score import components_sparse, thin_sparse_counts
from revision_nuisance import sample_profile_pool, nested_doublet_mask, group_draws
from revision_benchmark import score_methods
from scipy.stats import rankdata


def weights_for(G, up, down):
    w = np.zeros(G)
    w[up] = 1 / len(up)
    w[down] = -1 / len(down)
    return w


class ScoreInvariants(unittest.TestCase):
    def check_score(self, x, up, down):
        parts = rank_score_components(x, up, down)
        dense = scores_dense(x, up, down)
        sp = score_sparse(sparse.csc_matrix(x), weights_for(len(x), up, down))
        np.testing.assert_array_equal(dense, sp)
        np.testing.assert_allclose(dense, parts["score"], atol=2e-16, rtol=0)
        np.testing.assert_allclose(dense, parts["occupancy_component"] + parts["positive_rank_residual"],
                                   atol=2e-16, rtol=0)
        self.assertTrue(np.all(np.abs(parts["positive_rank_residual"]) <= parts["residual_bound"] + 2e-16))
        independent = components_sparse(sparse.csc_matrix(x), up, down)
        np.testing.assert_array_equal(dense, independent["score"])
        for key in ["occupancy_component", "positive_rank_residual", "residual_bound"]:
            np.testing.assert_allclose(parts[key], independent[key], atol=2e-16, rtol=0)
        np.testing.assert_allclose(independent["up_rank_contribution"] + independent["down_rank_contribution"],
                                   dense, atol=2e-16, rtol=0)
        return parts

    def test_all_zero_and_positive_ties_unequal_blocks(self):
        x = np.array([[0, 1, 4], [0, 1, 4], [0, 3, 1], [0, 3, 1],
                      [0, 0, 4], [0, 0, 1], [0, 0, 0]], dtype=np.int64)
        parts = self.check_score(x, [0, 1, 4], [2, 3])
        self.assertEqual(parts["score"][0], 0)
        self.assertEqual(parts["occupancy_component"][0], 0)
        self.assertEqual(parts["positive_rank_residual"][0], 0)

    def test_equal50_occupancy_is_not_entire_score(self):
        x = np.zeros((150, 2), dtype=np.int64)
        x[:50, 0] = 1
        x[50:100, 0] = 2
        x[100:, 1] = 1
        parts = self.check_score(x, list(range(50)), list(range(50, 100)))
        np.testing.assert_array_equal(parts["occupancy_component"], (parts["k_up"] - parts["k_down"]) / 100)
        self.assertEqual(parts["occupancy_component"][0], 0)
        self.assertLess(parts["score"][0], 0)  # detected counts alone cannot explain this contrast.

    def test_many_sparse_count_fixtures(self):
        rng = np.random.default_rng(61701)
        for G, m, n in [(17, 3, 5), (103, 50, 50), (211, 17, 31)]:
            x = rng.integers(0, 8, size=(G, 23))
            x[rng.random(x.shape) < .84] = 0
            x[:, 0] = 0
            self.check_score(x, list(range(m)), list(range(m, m+n)))

    def test_explicit_sparse_zeros_do_not_change_score(self):
        x = sparse.csc_matrix((np.array([0, 4, 4, 0, 2]), np.array([0, 1, 3, 0, 2]),
                               np.array([0, 3, 5])), shape=(6, 2))
        w = weights_for(6, [0, 1], [2, 3, 4])
        np.testing.assert_array_equal(score_sparse(x, w), scores_dense(x.toarray(), [0, 1], [2, 3, 4]))

    def test_duplicate_gene_symbols_sum_before_ranking(self):
        x = sparse.csc_matrix([[1, 0], [2, 4], [0, 2], [3, 0], [7, 8]])
        features = pd.DataFrame({"symbol": ["A", "A", "B", "C", "outside"]})
        mapped = map_universe(x, features, ["A", "B", "C", "missing"])
        expected = np.array([[3, 4], [0, 2], [3, 0], [0, 0]])
        np.testing.assert_array_equal(mapped.toarray(), expected)
        self.check_score(expected, [0], [1, 2])

    def test_outside_transcripts_keep_full_library_mass(self):
        x = sparse.csc_matrix([[2, 0], [3, 4]])
        profiles = prob_profiles(x, np.array([10, 20]))
        np.testing.assert_allclose(profiles.toarray(), [[.2, 0], [.3, .2], [.5, .8]], atol=1e-15)
        np.testing.assert_allclose(profiles.sum(axis=0), np.ones((1, 2)), atol=1e-15)

    def test_invalid_blocks_and_negative_counts_rejected_by_decomposition(self):
        with self.assertRaises(ValueError):
            rank_score_components(np.ones((4, 2)), [0, 0], [1])
        with self.assertRaises(ValueError):
            rank_score_components(np.ones((4, 2)), [0], [0, 1])
        with self.assertRaises(ValueError):
            rank_score_components(-np.ones((4, 2)), [0], [1])

    def test_new_components_canonicalize_duplicate_sparse_entries(self):
        duplicate = sparse.csc_matrix((np.array([1, 2, 3, 0]), np.array([0, 0, 2, 3]),
                                       np.array([0, 4])), shape=(5, 1))
        actual = components_sparse(duplicate, [0, 1], [2])
        expected = rank_score_components(duplicate.toarray(), [0, 1], [2])
        np.testing.assert_allclose(actual["score"], expected["score"], atol=2e-16, rtol=0)
        np.testing.assert_array_equal(actual["k_positive"], [2])
        # Frozen fractional ranks are (r-.5)/G; all-zero means are exactly .5.
        zeros = components_sparse(sparse.csc_matrix((5, 1)), [0, 1], [2])
        self.assertEqual(zeros["up_rank_contribution"][0], .5)
        self.assertEqual(zeros["down_rank_contribution"][0], -.5)


class ExactThinningInvariants(unittest.TestCase):
    def test_full_library_and_zero_draw_endpoints(self):
        x = sparse.csc_matrix([[2, 0, 4], [1, 3, 0]])
        totals = np.array([7, 7, 7])
        full = thin_sparse_counts(x, totals, 7, np.random.default_rng(57))
        empty = thin_sparse_counts(x, totals, 0, np.random.default_rng(57))
        np.testing.assert_array_equal(full.toarray(), x.toarray())
        self.assertEqual(empty.nnz, 0)

    def test_observed_capacity_and_outside_mass_are_preserved(self):
        rng = np.random.default_rng(841)
        x = sparse.csc_matrix(rng.integers(0, 9, size=(13, 21)))
        totals = np.asarray(x.sum(axis=0)).ravel() + rng.integers(0, 20, size=21)
        depth = int(totals.min())
        drawn = thin_sparse_counts(x, totals, depth, rng)
        self.assertTrue((drawn.toarray() <= x.toarray()).all())
        inside = np.asarray(drawn.sum(axis=0)).ravel()
        outside_draw = depth - inside
        outside_observed = totals - np.asarray(x.sum(axis=0)).ravel()
        self.assertTrue((outside_draw >= 0).all())
        self.assertTrue((outside_draw <= outside_observed).all())
        self.assertTrue((drawn.data == np.floor(drawn.data)).all())

    def test_finite_population_variance_matches_hypergeometric_fixture(self):
        x = sparse.csc_matrix(np.full((1, 4000), 3, dtype=np.int64))
        drawn = thin_sparse_counts(x, np.full(4000, 10), 4, np.random.default_rng(57091)).toarray()[0]
        self.assertLess(abs(drawn.mean() - 1.2), .05)
        # Var=draw*p*(1-p)*(N-draw)/(N-1)=.56, rather than multinomial .84.
        self.assertLess(abs(drawn.var() - .56), .05)
        self.assertLessEqual(drawn.max(), 3)

    def test_upward_noninteger_and_inconsistent_totals_rejected(self):
        for x, totals, depth in [
            (sparse.csc_matrix([[2]]), [2], 3),
            (sparse.csc_matrix([[1.5]]), [3], 2),
            (sparse.csc_matrix([[4]]), [3], 2),
        ]:
            with self.assertRaises(ValueError):
                thin_sparse_counts(x, totals, depth, np.random.default_rng(0))


class OperatingPointInvariants(unittest.TestCase):
    def test_zero_threshold_ties_make_nominal_and_achieved_FPR_different(self):
        neg = np.array([-1] + [0] * 19)
        pos = np.array([-1, 0, 0, 1])
        threshold = np.quantile(neg, .95)
        got = threshold_metrics(pos, neg, threshold)
        self.assertEqual(threshold, 0)
        self.assertEqual(got["achieved_comparator_FPR"], 0)
        self.assertEqual(got["positive_preparation_detection_fraction"], .25)
        self.assertEqual(got["comparator_FPR_if_greater_equal"], .95)
        self.assertEqual(got["positive_detection_if_greater_equal"], .75)

    def test_roc_all_tied_and_perfect_and_reversed(self):
        for pos, neg, expected in [([0, 0], [0, 0, 0], .5), ([2, 3], [0, 1], 1), ([0, 1], [2, 3], 0)]:
            curve = strict_roc(pos, neg)
            self.assertEqual(empirical_auc(pos, neg), expected)
            area = np.trapezoid(curve.positive_preparation_detection_fraction.to_numpy()[::-1],
                                curve.achieved_comparator_FPR.to_numpy()[::-1])
            self.assertEqual(area, expected)
            self.assertEqual(curve.iloc[0].achieved_comparator_FPR, 1)
            self.assertEqual(curve.iloc[-1].achieved_comparator_FPR, 0)
            self.assertEqual(curve.iloc[-1].positive_preparation_detection_fraction, 0)

class NuisanceSamplingInvariants(unittest.TestCase):
    def test_nested_doublet_fraction_is_exact_at_every_declared_prefix(self):
        mask = nested_doublet_mask(6400, .1, np.random.default_rng(0))
        for size in [400, 1600, 6400]:
            self.assertEqual(mask[:size].sum(), size // 10)

    def test_paired_profile_draws_use_the_selected_profiles_own_depth(self):
        profiles = sparse.csc_matrix([[1., 0.], [0., 1.], [0., 0.]])
        own = np.array([2, 11])
        paired = sample_profile_pool(profiles, own, own, np.array([1., -1.]),
                                     np.random.default_rng(31), 200, True)
        np.testing.assert_array_equal(paired["observed_depth"], own[paired["profile_index"]])
        independent = sample_profile_pool(profiles, own, own, np.array([1., -1.]),
                                          np.random.default_rng(31), 200, False)
        self.assertTrue(np.any(independent["observed_depth"] != own[independent["profile_index"]]))

    def test_group_draw_prefixes_reproduce_without_changing_source_pools(self):
        scenario = np.array([-2., -1., 0., 1., 2.])
        baseline = np.array([-1., 0., 1.])
        small = group_draws(scenario, baseline, 7, 11, 8123, 200)
        large = group_draws(scenario, baseline, 7, 11, 8123, 1000)
        np.testing.assert_array_equal(small, large[:200])
        np.testing.assert_array_equal(scenario, [-2., -1., 0., 1., 2.])
        np.testing.assert_array_equal(baseline, [-1., 0., 1.])


class LockedBenchmarkInvariants(unittest.TestCase):
    def test_three_methods_match_independent_dense_formulas(self):
        rng = np.random.default_rng(984)
        G = 1603
        x = rng.integers(0, 7, size=(G, 9))
        x[rng.random(x.shape) < .94] = 0
        x[:, 0] = 0
        genes = [f"g{i}" for i in range(G)]
        programs = {"unequal": {"up": genes[:17], "down": genes[17:48]},
                    "equal": {"up": genes[50:100], "down": genes[100:150]},
                    "missing": {"up": genes[:13] + ["absent"], "down": genes[23:42]}}
        totals = x.sum(axis=0) + 500
        actual = score_methods(sparse.csc_matrix(x), totals, genes, programs)
        ascending = rankdata(x, axis=0, method="average")
        descending = np.minimum(rankdata(-x, axis=0, method="average"), 1500)
        logs = np.log1p(10000 * x / totals)
        for j, program in enumerate(programs.values()):
            u = [genes.index(g) for g in program["up"] if g in genes]
            d = [genes.index(g) for g in program["down"] if g in genes]
            expected_rank = (ascending[u].mean(axis=0)-ascending[d].mean(axis=0))/G
            expected_log = logs[u].mean(axis=0)-logs[d].mean(axis=0)
            def block(indices):
                n = len(indices)
                return (n*1500-descending[indices].sum(axis=0))/(n*1500-n*(n+1)/2)
            expected_ucell = np.maximum(0, block(u)-block(d))
            for index, expected in enumerate([expected_rank, expected_log, expected_ucell]):
                np.testing.assert_allclose(actual[j,index], expected, atol=3e-15, rtol=0)

    def test_capped_rank_non_evaluable_below_committed_universe_size(self):
        x = sparse.csc_matrix([[1, 0], [0, 2], [3, 1]])
        actual = score_methods(x, [4, 3], ["A","B","C"],
                               {"program":{"up":["A"],"down":["B"]}})
        self.assertTrue(np.isnan(actual[0,2]).all())


if __name__ == "__main__":
    unittest.main()
