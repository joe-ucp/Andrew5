"""Numerical zero certificates and narrowly bounded rank-metric comparison."""
from pathlib import Path
import itertools
import sys
import unittest

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'code'))
from rank_metric_comparison import (log_score_error_radii, signed_log_error_radii,
                                    spearman_envelope, auc_envelope, empirical_auc,
                                    compare_rank_metric)
from revision_benchmark import score_methods


class ZeroCertificateTests(unittest.TestCase):
    def test_only_nonempty_certified_cancellations_receive_radius(self):
        radius = log_score_error_radii([0., 2., 2.], [0, 2, 2], [True, False, True])
        np.testing.assert_array_equal(radius[:2], [0., 0.])
        self.assertGreater(radius[2], 0.)
        self.assertLess(radius[2], 1e-13)

    def test_count_classes_not_just_equal_total_expression(self):
        # Column0 equal multisets; column1 equal sums with unequal multisets;
        # column2 absent; column3 unmatched detected counts.
        x = np.array([[2, 1, 0, 1], [3, 4, 0, 1], [3, 2, 0, 1], [2, 3, 0, 0]])
        genes = list('abcd'); programs = {'original':{'up':['a','b'],'down':['c','d']}}
        radii = signed_log_error_radii(x, [20]*4, genes, programs)
        self.assertGreater(radii[0,0], 0.)
        np.testing.assert_array_equal(radii[0,1:], [0., 0., 0.])
        values = score_methods(x, np.array([20]*4), genes, programs)[0,1]
        self.assertLessEqual(abs(values[0]), radii[0,0])

    def test_unequal_represented_blocks_use_integer_weighted_balance(self):
        x = np.array([[2,2],[2,2],[2,0]])
        p = {'p':{'up':['a'],'down':['b','c']}}
        r = signed_log_error_radii(x, [10,10], list('abc'), p)
        self.assertGreater(r[0,0],0.)
        self.assertEqual(r[0,1],0.)

    def test_many_term_residual_is_inside_certificate(self):
        rng = np.random.default_rng(9401)
        block = rng.integers(1,10,size=(50,20))
        other = np.stack([rng.permutation(block[:,j]) for j in range(20)],axis=1)
        x = np.vstack([block,other]); totals = x.sum(0)+100
        genes = [f'g{i}' for i in range(100)]
        programs = {'original':{'up':genes[:50],'down':genes[50:]}}
        r = signed_log_error_radii(x,totals,genes,programs)[0]
        v = score_methods(x,totals,genes,programs)[0,1]
        self.assertTrue((r>0).all())
        self.assertTrue((np.abs(v)<=r).all())

    def test_invalid_count_and_signature_inputs_fail(self):
        p = {'p':{'up':['a'],'down':['b']}}
        for x, totals, genes in [(np.array([[1.5],[1.5]]),[4],['a','b']),
                                  (np.array([[2],[2]]),[1],['a','b']),
                                  (np.array([[1],[1]]),[2],['a','a'])]:
            with self.subTest(x=x.tolist()), self.assertRaises(ValueError):
                signed_log_error_radii(x,totals,genes,p)
        with self.assertRaises(ValueError):log_score_error_radii([1.],[0],[True])


class MetricEnvelopeTests(unittest.TestCase):
    def setUp(self):
        self.radius = np.array([0.,1e-14,0.,1e-14,0.])
        self.fresh = np.array([-1.,-1e-16,0.,1e-16,1.])
        self.old = np.array([-1.,1e-16,0.,-1e-16,1.])
        self.depth = np.array([1.,5.,3.,2.,4.])
        self.positive = np.array([False,True,False,False,True])

    def test_admissible_residue_rank_changes_are_bounded(self):
        a = float(spearmanr(self.old,self.depth).statistic)
        b = float(spearmanr(self.fresh,self.depth).statistic)
        self.assertNotEqual(a,b)
        r = compare_rank_metric(a,b,self.fresh,self.radius,metric='spearman',covariate=self.depth)
        self.assertTrue(r['accepted'])
        r = compare_rank_metric(empirical_auc(self.old,self.positive),
                                empirical_auc(self.fresh,self.positive),
                                self.fresh,self.radius,metric='auc',positive=self.positive)
        self.assertTrue(r['accepted'])

    def test_envelopes_cover_all_small_cluster_weak_orderings(self):
        radii = np.array([0.,1e-14,1e-14,1e-14,0.])
        scores = np.array([-1.,0.,0.,0.,1.])
        s = spearman_envelope(scores,radii,self.depth)
        a = auc_envelope(scores,radii,self.positive)
        for middle in itertools.product([-1e-16,0.,1e-16],repeat=3):
            values = np.array([-1.,*middle,1.])
            rho = float(spearmanr(values,self.depth).statistic)
            auc = empirical_auc(values,self.positive)
            self.assertGreaterEqual(rho,s['lower']-1e-15)
            self.assertLessEqual(rho,s['upper']+1e-15)
            self.assertGreaterEqual(auc,a['lower']-1e-15)
            self.assertLessEqual(auc,a['upper']+1e-15)

    def test_fixed_absent_ties_have_zero_width(self):
        scores = np.array([-1.,0.,0.,0.,1.]); zero = np.zeros(5)
        for bounds in [spearman_envelope(scores,zero,self.depth),auc_envelope(scores,zero,self.positive)]:
            self.assertEqual(bounds['lower'],bounds['upper'])
            self.assertEqual(bounds['uncertain_clusters'],0)

    def test_real_non_tie_rank_change_is_rejected(self):
        scores = np.array([-2.,-1.,0.,1.,2.]); depth = np.arange(5.)
        changed = scores.copy(); changed[[0,4]] = changed[[4,0]]
        old = float(spearmanr(scores,depth).statistic)
        fresh = float(spearmanr(changed,depth).statistic)
        result = compare_rank_metric(old,fresh,changed,np.zeros(5),metric='spearman',covariate=depth)
        self.assertFalse(result['accepted'])

    def test_fresh_metric_must_recompute_even_when_bound_is_wide(self):
        actual = float(spearmanr(self.fresh,self.depth).statistic)
        result = compare_rank_metric(actual,actual+.01,self.fresh,self.radius,metric='spearman',covariate=self.depth)
        self.assertFalse(result['accepted'])
        self.assertEqual(result['reason'],'fresh_metric_does_not_recompute')

    def test_score_outside_zero_certificate_rejected(self):
        broken = self.fresh.copy(); broken[1]=1e-6
        with self.assertRaises(ValueError):spearman_envelope(broken,self.radius,self.depth)

    def test_degenerate_uncertain_rank_variance_is_not_blanket_acceptance(self):
        scores = np.array([-1e-16,0.,1e-16]); depth = np.arange(3.)
        actual = float(spearmanr(scores,depth).statistic)
        result = compare_rank_metric(0.,actual,scores,np.ones(3)*1e-14,metric='spearman',covariate=depth)
        self.assertFalse(result['accepted'])
        self.assertIn('Degenerate',result['reason'])


if __name__ == '__main__':unittest.main()
