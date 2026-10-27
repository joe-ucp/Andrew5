"""Independent synthetic checks of revision02 algebra and matched references."""
from fractions import Fraction
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy import sparse
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from revision02_gene_math import (summarize_genes, program_numerators,
    donor_contrasts, equal_donor_mean, match_candidates, draw_matched_programs)
from revision_score import components_sparse


class IndependentGeneAlgebra(unittest.TestCase):
    def setUp(self):
        self.counts = np.zeros((160, 4), dtype=np.int64)
        self.counts[[0, 1, 51, 80, 130], 1] = [1, 1, 3, 7, 3]
        self.counts[[3, 8, 55, 120], 2] = [6, 1, 6, 2]
        self.counts[:, 3] = 1
        self.totals = self.counts.sum(axis=0) + 20
        self.programs = np.array([np.arange(100), np.arange(20, 120)])

    def test_program_integer_numerators_equal_independent_dense_ranks(self):
        summary = summarize_genes(self.counts, self.totals)
        numerators = program_numerators(summary, self.programs)
        G, n = self.counts.shape
        for i, program in enumerate(self.programs):
            exact = Fraction(0)
            for column in self.counts.T:
                twice = np.rint(2 * rankdata(column, method='average')).astype(int)
                exact += Fraction(int(twice[program[:50]].sum() - twice[program[50:]].sum()), 100 * G * n)
            self.assertEqual(Fraction(int(numerators['score'][i]), 100 * G * n), exact)
            comp = components_sparse(sparse.csc_matrix(self.counts), program[:50], program[50:])
            for name, key in [('occupancy', 'occupancy_component'), ('residual', 'positive_rank_residual'), ('score', 'score')]:
                denominator = 100 * n * (1 if name == 'occupancy' else G)
                self.assertAlmostEqual(float(numerators[name][i] / denominator), float(comp[key].mean()), places=15)

    def test_signed_gene_contributions_add_and_conditional_summaries_are_explicit(self):
        counts = np.zeros((160, 3), dtype=int)
        counts[0] = [1, 0, 4]
        totals = np.array([10, 10, 20])
        summary = summarize_genes(counts, totals)
        self.assertEqual(summary['detected_n'][0], 2)
        self.assertAlmostEqual(summary['detection_fraction'][0], 2 / 3)
        self.assertEqual(summary['positive_mean_raw'][0], 2.5)
        self.assertEqual(summary['positive_mean_CP10K'][0], 1500)
        self.assertEqual(summary['mean_CP10K'][0], 1000)
        self.assertAlmostEqual(summary['positive_mean_log1p_CP10K'][0], (np.log1p(1000) + np.log1p(2000)) / 2)
        for key in ['positive_mean_raw', 'positive_mean_CP10K', 'positive_mean_log1p_CP10K']:
            self.assertTrue(np.isnan(summary[key][1]))
        program = np.arange(100)
        signs = np.r_[np.ones(50), -np.ones(50)]
        occupancy = signs * summary['detection_fraction'][program] / 100
        residual = signs * summary['centered_rank_sum'][program] / (100 * 160 * 3)
        dense = (rankdata(counts, axis=0, method='average') - .5) / 160
        self.assertAlmostEqual(float((occupancy + residual).sum()), float((dense[:50].mean(axis=0) - dense[50:100].mean(axis=0)).mean()), places=15)

    def test_duplicate_sparse_entries_count_one_detection(self):
        matrix = sparse.coo_matrix(([1, 3], ([0, 0], [0, 0])), shape=(160, 1))
        summary = summarize_genes(matrix, np.array([10]))
        self.assertEqual(summary['detected_n'][0], 1)
        self.assertEqual(summary['positive_mean_raw'][0], 4)

    def test_duplicated_empirical_population_gives_exact_zero_contrast(self):
        first = summarize_genes(self.counts, self.totals)
        second = summarize_genes(np.tile(self.counts, 2), np.tile(self.totals, 2))
        contrast = donor_contrasts(first, second, self.programs)
        for component in ['score', 'occupancy', 'residual']:
            self.assertTrue(np.array_equal(contrast[component]['numerators'], [0, 0]))
        aggregate = equal_donor_mean([
            {'score': {'numerators': np.array([1, 2]), 'denominator': 3}},
            {'score': {'numerators': np.array([-2, -4]), 'denominator': 6}},
        ], 'score')
        self.assertEqual(aggregate['numerators'].tolist(), [0, 0])


class IndependentMatchingDesign(unittest.TestCase):
    def fixture(self):
        genes = np.array([f'g{i:04}' for i in range(400)])
        original = np.arange(100)
        p = np.full((4, 400), .2)
        L = np.ones((4, 400))
        p[:, 300:350] = .3
        L[:, 350:] = 2
        return genes, original, p, L

    def test_calipers_exclude_original_and_programs_are_disjoint_reproducible(self):
        genes, original, p, L = self.fixture()
        candidates, support = match_candidates(genes, original, p, L)
        self.assertTrue(all(s['eligible_candidates'] == 200 for s in support))
        self.assertTrue(all(np.array_equal(c, np.arange(100, 300)) for c in candidates))
        programs, quality, status = draw_matched_programs(original, candidates, p, L, seed=42, target=3, max_attempts=3)
        repeated, _, _ = draw_matched_programs(original, candidates, p, L, seed=42, target=3, max_attempts=3)
        self.assertEqual(status['status'], 'COMPLETE')
        self.assertEqual(len(quality), 3 * 4 * 2)
        self.assertTrue(np.array_equal(programs, repeated))
        self.assertTrue(all(len(set(row)) == 100 and not set(row).intersection(original) for row in programs))
        self.assertTrue(all(abs(row['delta_detection']) <= .025 and abs(row['delta_log_abundance']) <= .25 for row in quality))

    def test_insufficient_support_and_donor_quality_fail_without_relaxation(self):
        genes, original, p, L = self.fixture()
        short = [np.arange(100, 109) for _ in range(100)]
        programs, quality, status = draw_matched_programs(original, short, p, L, target=3, max_attempts=4)
        self.assertEqual(status['status'], 'FAILED_CANDIDATE_SUPPORT')
        self.assertEqual(status['attempts'], 0)
        self.assertEqual(programs.shape, (0, 100))
        self.assertEqual(quality, [])
        p[:, 100:300] += np.array([.1, -.1, .1, -.1])[:, None]
        candidates, support = match_candidates(genes, original, p, L)
        self.assertTrue(all(s['support_pass'] for s in support))
        programs, _, status = draw_matched_programs(original, candidates, p, L, target=3, max_attempts=4)
        self.assertEqual(status['status'], 'FAILED_TARGET_NOT_REACHED')
        self.assertEqual(status['attempts'], 4)
        self.assertEqual(status['quality_rejections'], 4)
        self.assertEqual(programs.shape, (0, 100))


if __name__ == '__main__':
    unittest.main()
