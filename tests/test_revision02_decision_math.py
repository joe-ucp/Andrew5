"""Independent dense and boundary fixtures; no source data or new outcomes."""
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy import sparse
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from revision02_math import (zero_centered_components, unclipped_ucell,
    score_bundle, q95, decision_metrics, state_auc)
from revision_score import components_sparse
from revision_benchmark import score_methods


class BaselineAndScorers(unittest.TestCase):
    def test_unequal_blocks_zero_ties_and_full_detection_cancel_baseline(self):
        x = np.array([[0, 0, 1], [0, 4, 1], [0, 0, 1], [0, 1, 1],
                      [0, 4, 1], [0, 2, 1], [0, 0, 1]])
        comp = components_sparse(sparse.csc_matrix(x), [0, 1], [2, 3, 4])
        actual = zero_centered_components(comp, 7)
        dense = (rankdata(x, axis=0, method='average') - .5) / 7
        baseline = (7 - (x > 0).sum(axis=0)) / 14
        np.testing.assert_array_equal(actual['zero_rank_baseline'], baseline)
        np.testing.assert_allclose(actual['centered_up'], dense[[0, 1]].mean(axis=0) - baseline, rtol=0, atol=2e-16)
        np.testing.assert_allclose(actual['centered_down'], -dense[[2, 3, 4]].mean(axis=0) + baseline, rtol=0, atol=2e-16)
        self.assertEqual(actual['zero_rank_baseline'][0], .5)
        self.assertEqual(actual['centered_up'][0], 0)
        self.assertEqual(actual['centered_down'][0], 0)
        self.assertEqual(actual['zero_rank_baseline'][2], 0)
        np.testing.assert_allclose(actual['centered_up'] + actual['centered_down'], comp['score'], rtol=0, atol=2e-16)
        with self.assertRaises(ValueError):
            zero_centered_components(comp, 0)

    def test_unclipped_dense_definition_and_exact_locked_clip(self):
        G = 2000
        genes = [f'g{i}' for i in range(G)]
        programs = {
            'unequal': {'up': ['g0', 'g1', 'g3', 'missing'], 'down': ['g2', 'g4']},
            'balanced': {'up': genes[:50], 'down': genes[50:100]},
        }
        x = np.zeros((G, 4), dtype=int)
        x[[0, 1, 2, 3, 50, 51, 1999], 1] = [2, 2, 8, 1, 9, 9, 1]
        x[:, 2] = 1
        x[[4, 51, 99], 3] = [20, 20, 1]
        totals = x.sum(axis=0) + 30
        expected = []
        lookup = {gene: i for i, gene in enumerate(genes)}
        for program in programs.values():
            blocks = [[lookup[g] for g in program[side] if g in lookup] for side in ['up', 'down']]
            row = []
            for column in x.T:
                ranks = np.minimum(rankdata(-column, method='average'), 1500)
                scores = []
                for block in blocks:
                    n = len(block)
                    scores.append((n * 1500 - ranks[block].sum()) / (n * 1500 - n * (n + 1) / 2))
                row.append(scores[0] - scores[1])
            expected.append(row)
        actual = unclipped_ucell(sparse.csc_matrix(x), genes, programs)
        np.testing.assert_array_equal(actual, np.asarray(expected))
        bundled = score_bundle(x, totals, genes, programs)
        self.assertEqual(bundled.shape, (2, 4, 4))
        np.testing.assert_array_equal(bundled[:, :3], score_methods(x, totals, genes, programs))
        np.testing.assert_array_equal(np.maximum(bundled[:, 3], 0), bundled[:, 2])
        np.testing.assert_array_equal(bundled[:, 3], actual)


class ThresholdAndTransitions(unittest.TestCase):
    def test_linear_threshold_and_strict_zero_ties(self):
        values = np.tile([-1., 0., 1.], 7)
        result = decision_metrics(values, 0)
        self.assertEqual((result['positive_n'], result['negative_n'], result['tie_n']), (7, 14, 7))
        self.assertEqual(result['positive_fraction'], 1 / 3)
        self.assertEqual(result['tie_fraction'], 1 / 3)
        self.assertEqual(q95(np.arange(20)), 18.05)
        self.assertEqual(q95(np.zeros(20)), 0)

    def test_upward_movement_without_crossing_does_not_create_calls(self):
        native = np.full(20, -2.)
        values = np.full(20, -1.)
        result = decision_metrics(values, 0, native=native)
        self.assertTrue(np.all(values > native))
        self.assertEqual(result['positive_n'], 0)
        self.assertEqual(result['negative_to_positive_n'], 0)
        self.assertEqual(result['positive_to_negative_n'], 0)

    def test_paired_transition_counts_and_denominator_are_all_cells(self):
        values = np.tile([-1., 1., 0., -.5], 5)
        native = np.tile([-2., -1., 1., -.5], 5)
        result = decision_metrics(values, 0, native=native)
        self.assertEqual(result['native_positive_n'], 5)
        self.assertEqual(result['negative_to_positive_n'], 5)
        self.assertEqual(result['positive_to_negative_n'], 5)
        self.assertEqual(result['negative_to_positive_fraction'], .25)
        self.assertEqual(result['positive_to_negative_fraction'], .25)
        recalibrated = decision_metrics(np.zeros(20), -.5, native=np.zeros(20), native_threshold=0)
        self.assertEqual(recalibrated['negative_to_positive_n'], 20)

    def test_small_state_keeps_counts_but_not_fractions(self):
        result = decision_metrics(np.arange(19), 5, native=np.zeros(19))
        self.assertFalse(result['evaluable'])
        self.assertEqual(result['n'], 19)
        self.assertEqual(result['positive_n'], 13)
        self.assertEqual(result['tie_n'], 1)
        self.assertEqual(result['negative_to_positive_n'], 13)
        for key in ['positive_fraction', 'tie_fraction', 'negative_to_positive_fraction', 'positive_to_negative_fraction']:
            self.assertTrue(np.isnan(result[key]))
        self.assertEqual(result['reason'], 'fewer_than_20_cells')
        self.assertTrue(np.isnan(q95(np.arange(19))))

    def test_missing_threshold_and_misaligned_pairs_are_explicit(self):
        result = decision_metrics(np.zeros(20), np.nan)
        self.assertFalse(result['evaluable'])
        self.assertEqual(result['n'], 20)
        self.assertEqual(result['reason'], 'threshold_unavailable')
        self.assertTrue(np.isnan(result['positive_n']))
        self.assertTrue(np.isnan(q95(np.r_[np.arange(19), np.nan])))
        with self.assertRaises(ValueError):
            decision_metrics(np.zeros(20), 0, native=np.zeros(19))

    def test_auc_requires_both_states_and_preserves_ties(self):
        states = np.array(['alpha'] * 20 + ['beta'] * 20)
        self.assertEqual(state_auc(np.zeros(40), states), .5)
        scores = np.r_[np.zeros(20), np.ones(20)]
        self.assertEqual(state_auc(scores, states), 1)
        self.assertEqual(state_auc(-scores, states), 0)
        self.assertTrue(np.isnan(state_auc(scores[:-1], states[:-1])))


if __name__ == '__main__':
    unittest.main()
