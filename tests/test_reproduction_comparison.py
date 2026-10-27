"""The cancellation exception must not hide unrelated scientific changes."""
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from reproduce import compare_native_log_rank_metrics
from reference_analysis import auc
from revision_benchmark import corr


class NativeLogComparisonTests(unittest.TestCase):
    def fixture(self):
        positive=np.array([False,True,False,True,False,True])
        depths=np.array([10.,12.,20.,16.,30.,40.])
        canonical=np.array([-.8,0.,0.,0.,.2,.8])
        observed=canonical.copy();observed[1]=-1e-18;observed[3]=1e-18
        radii=np.array([0.,1e-16,0.,1e-16,0.,0.])
        bundle={'programs':np.array(['original']),
                'methods':np.array(['signed_rank','signed_log1p_CPT10K','ucell_form_1500']),
                'scores':np.array([[canonical,observed,canonical]]),
                'total_counts':depths,'log_roundoff_radius':radii[None,:],
                'log_state_is_beta':positive}
        def table(scores,method='signed_log1p_CPT10K'):
            return pd.DataFrame([{'donor':1,'mode':'original_universe','program':'original',
                'method':method,'alpha_n':3,'beta_n':3,
                'beta_vs_alpha_AUROC':auc(scores[positive],scores[~positive]),
                'alpha_score_depth_r':corr(scores[~positive],depths[~positive]),
                'beta_score_depth_r':corr(scores[positive],depths[positive]),
                'alpha_native_q95':.18,'achieved_alpha_FPR':1/3,
                'beta_detection_fraction':1/3,'scope':'synthetic state labels'}])
        return table(canonical),table(observed),bundle,table

    def test_certified_zero_rank_variation_is_explicitly_reported(self):
        old,new,bundle,_=self.fixture()
        old_before=old.copy(deep=True);new_before=new.copy(deep=True)
        audits=compare_native_log_rank_metrics(old,new,lambda *unused:bundle)
        self.assertGreater(len(audits),0)
        self.assertTrue(all(row['accepted'] for row in audits))
        pd.testing.assert_frame_equal(old,old_before)
        pd.testing.assert_frame_equal(new,new_before)

    def test_nonzero_rank_change_is_rejected(self):
        old,_,bundle,table=self.fixture()
        scores=bundle['scores'][0,1].copy();scores[0]=.9
        bundle['scores'][0,1]=scores
        with self.assertRaises(AssertionError):
            compare_native_log_rank_metrics(old,table(scores),lambda *unused:bundle)

    def test_rank_and_ucell_rows_have_no_exception(self):
        for method in ['signed_rank','ucell_form_1500']:
            old,new,bundle,_=self.fixture()
            old['method']=method;new['method']=method
            with self.subTest(method=method),self.assertRaises(AssertionError):
                compare_native_log_rank_metrics(old,new,lambda *unused:bundle)

    def test_thresholds_decision_fractions_and_counts_remain_strict(self):
        for column,change in [('alpha_native_q95',1e-8),('achieved_alpha_FPR',1e-8),
                              ('beta_detection_fraction',1e-8),('alpha_n',1)]:
            old,new,bundle,_=self.fixture();new[column]+=change
            with self.subTest(column=column),self.assertRaises((AssertionError,ValueError)):
                compare_native_log_rank_metrics(old,new,lambda *unused:bundle)

    def test_generated_metric_must_match_its_actual_cell_scores(self):
        old,new,bundle,_=self.fixture();new['beta_vs_alpha_AUROC']=.12345
        with self.assertRaises(AssertionError):
            compare_native_log_rank_metrics(old,new,lambda *unused:bundle)

    def test_forged_wide_radius_is_rejected(self):
        old,new,bundle,_=self.fixture();bundle['log_roundoff_radius'][0,1]=.1
        with self.assertRaises(ValueError):
            compare_native_log_rank_metrics(old,new,lambda *unused:bundle)


if __name__=='__main__':unittest.main()
