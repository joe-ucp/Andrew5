"""Independent, synthetic checks of fixed-contribution subtraction semantics."""
from fractions import Fraction
from pathlib import Path
import unittest
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'code'))
from revision03_influence import calculate_influence


def fixture(rows, donor='synthetic'):
    return pd.DataFrame([
        {'donor': donor, 'gene': gene, 'block': block,
         'delta_signed_occupancy': occupancy,
         'delta_signed_residual': remainder,
         'delta_signed_zero_baseline_rank': float(Fraction(occupancy) + Fraction(remainder))}
        for gene, block, occupancy, remainder in rows
    ])


class InfluenceTests(unittest.TestCase):
    def example(self):
        # Binary-exact terms; A is positive while B opposes it, C/D tie.
        return fixture([('D', 'down', .125, 0.), ('B', 'down', -.25, -.125),
                        ('A', 'up', .25, .25), ('C', 'up', .0625, .0625)])

    def test_every_single_and_prefix_match_independent_exact_subtraction(self):
        frame = self.example()
        before = frame.copy(deep=True)
        single, top, selected, summary = calculate_influence(frame)
        terms = {row.gene: Fraction(row.delta_signed_occupancy) +
                 Fraction(row.delta_signed_residual) for row in frame.itertuples()}
        total = sum(terms.values(), Fraction())
        self.assertEqual(total, Fraction(3, 8))
        self.assertEqual(set(single.gene), set(terms))
        for row in single.itertuples():
            self.assertEqual(row.single_residual, float(total - terms[row.gene]))
            self.assertEqual(row.component_sum, float(total))
            self.assertEqual(row.single_residual,
                             row.remaining_occupancy + row.remaining_rank_remainder)
        self.assertEqual(selected.gene.tolist(), ['A', 'C', 'D'])
        for k, expected in [(1, Fraction(-1, 8)), (2, Fraction(-1, 4)),
                            (3, Fraction(-3, 8))]:
            row = top[top.k == k].iloc[0]
            self.assertEqual(row.topk_residual, float(expected))
            self.assertEqual(row.remaining_sign, 'negative')
        self.assertEqual(summary.iloc[0].min_single_removed_gene, 'A')
        self.assertEqual(summary.iloc[0].min_single_residual, -.125)
        self.assertEqual(int(summary.iloc[0].single_strict_sign_reversals), 1)
        # Deleting an opposing contribution must increase the remaining contrast.
        self.assertEqual(single.set_index('gene').loc['B', 'single_residual'], .75)
        pd.testing.assert_frame_equal(frame, before)

    def test_all_requested_k_rows_preserved_when_positive_support_is_insufficient(self):
        _, top, selected, summary = calculate_influence(self.example())
        self.assertEqual(top.k.tolist(), [1, 2, 3, 4, 5])
        self.assertEqual(len(selected), 3)
        for k in [4, 5]:
            row = top[top.k == k].iloc[0]
            self.assertEqual(row.status, 'insufficient_positive_terms')
            self.assertTrue(pd.isna(row.topk_residual))
            self.assertEqual(int(row.available_positive_terms), 3)
        self.assertFalse(bool(summary.iloc[0].top5_available))
        self.assertTrue(pd.isna(summary.iloc[0].top5_residual))

    def test_exact_ties_are_alphabetical_and_donor_selection_is_separate(self):
        a = self.example()
        b = fixture([('C', 'up', -.25, 0.), ('B', 'down', .125, 0.),
                     ('D', 'down', .25, .25), ('A', 'up', .125, 0.)], donor='other')
        source = pd.concat([a, b], ignore_index=True).iloc[::-1]
        _, top, selected, _ = calculate_influence(source)
        self.assertEqual(selected[selected.donor == 'synthetic'].gene.tolist(), ['A', 'C', 'D'])
        self.assertEqual(selected[selected.donor == 'other'].gene.tolist(), ['D', 'A', 'B'])
        self.assertEqual(top[(top.donor == 'other') & (top.k == 2)].iloc[0].selected_genes, 'D;A')

    def test_strict_sign_is_not_replaced_by_near_zero_flag(self):
        frame = fixture([('zero', 'down', 0., 0.), ('positive', 'up', 1e-12, 0.)])
        single, _, _, _ = calculate_influence(frame)
        at_boundary = single.set_index('gene').loc['zero']
        self.assertEqual(at_boundary.single_residual, 1e-12)
        self.assertEqual(at_boundary.remaining_sign, 'positive')
        self.assertFalse(bool(at_boundary.near_zero_1e_12))
        exact_zero = single.set_index('gene').loc['positive']
        self.assertEqual(exact_zero.remaining_sign, 'zero')
        self.assertTrue(bool(exact_zero.near_zero_1e_12))
        tiny = fixture([('zero', 'down', 0., 0.), ('positive', 'up', 5e-13, 0.)])
        small = calculate_influence(tiny)[0].set_index('gene').loc['zero']
        self.assertEqual(small.remaining_sign, 'positive')
        self.assertTrue(bool(small.near_zero_1e_12))

    def test_zero_and_negative_only_inputs_do_not_invent_positive_selections(self):
        for terms in [[('a', 'up', 0., 0.), ('b', 'down', 0., 0.)],
                      [('a', 'up', -.125, 0.), ('b', 'down', -.25, 0.)]]:
            _, top, selected, summary = calculate_influence(fixture(terms))
            self.assertEqual(top.k.tolist(), [1, 2, 3, 4, 5])
            self.assertTrue((top.status == 'insufficient_positive_terms').all())
            self.assertTrue(top.topk_residual.isna().all())
            self.assertEqual(len(selected), 0)
            self.assertEqual(int(summary.iloc[0].available_positive_terms), 0)

    def test_invalid_contributions_and_inconsistent_native_totals_fail(self):
        frame = self.example()
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            calculate_influence(pd.concat([frame, frame.iloc[[0]]], ignore_index=True))
        with self.assertRaisesRegex(ValueError, 'Missing'):
            calculate_influence(frame.drop(columns='delta_signed_residual'))
        bad = frame.copy(); bad.loc[0, 'block'] = 'unassigned'
        with self.assertRaisesRegex(ValueError, 'block'):
            calculate_influence(bad)
        bad = frame.copy(); bad.loc[0, 'delta_signed_occupancy'] = np.inf
        with self.assertRaisesRegex(ValueError, 'Finite'):
            calculate_influence(bad)
        bad = frame.copy(); bad.loc[0, 'delta_signed_zero_baseline_rank'] += 1e-10
        with self.assertRaisesRegex(AssertionError, 'Canonical'):
            calculate_influence(bad)
        with self.assertRaisesRegex(AssertionError, 'original donor'):
            calculate_influence(frame, native={'synthetic': .5})


if __name__ == '__main__':
    unittest.main()
