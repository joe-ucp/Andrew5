"""Optional retained-data regressions. No data resolution or reads on import."""
from pathlib import Path
import hashlib
import json
import os
import sys
import unittest

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))

@unittest.skipUnless(os.environ.get('STRIKE007_EXTERNAL_TESTS')=='1',
                     'Optional external regression: run code/run_tests.py --external-data --data-root ABSOLUTE_PATH')
class RetainedDataRegression(unittest.TestCase):
    def test_refined_zero_contamination_uses_identical_pool_means(self):
        from data_paths import evidence_path
        path=evidence_path('countermodel_stability_batches.tsv',must_exist=True)
        manifest=json.loads((ROOT/'evidence/original_bulk_manifest.json').read_text(encoding='utf-8-sig'))
        record=next(r for r in manifest['files'] if r['logical_path']=='evidence/'+path.name)
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),record['sha256'])
        batches=pd.read_csv(path,sep='\t',float_precision='round_trip')
        zero=batches[batches.model=='beta_RNA_0']
        self.assertEqual(len(zero),40)
        np.testing.assert_array_equal(zero.pool_baseline_mean,zero.pool_scenario_mean)
        np.testing.assert_array_equal(zero.conditional_pool_mean_effect,np.zeros(40))
        self.assertTrue((zero.group_sd>0).all())

if __name__=='__main__':unittest.main()
