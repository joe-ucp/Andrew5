"""Arbitrary-input diagnostic contracts and one entirely synthetic CLI run."""
from pathlib import Path
import json
import subprocess
import sys
import unittest

import numpy as np
import pandas as pd
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
from signature_diagnostic import diagnose, report_html, threshold_counts
from test_support import temporary_fixture


class DiagnosticTests(unittest.TestCase):
    def fixture(self):
        genes=['up','down','outside'];cells=['a','b','c','d']
        x=sparse.csc_matrix([[1,0,5,5],[5,0,1,0],[10,10,10,10]])
        signature={'up':['up'],'down':['down'],'universe':['up','down']+[f'other{i}' for i in range(1601)]}
        meta=pd.DataFrame({'cell_id':cells,'donor':'D1','state':['low','low','high','high'],
                           'treatment':['vehicle','exposure','vehicle','exposure'],'total_counts':100})
        return x,genes,cells,signature,meta

    def test_detection_directions_baseline_and_explicit_state_metrics(self):
        x,genes,cells,signature,meta=self.fixture()
        result,table=diagnose(x,genes,cells,signature,meta,thresholds={'signed_rank':0},
                              baseline_treatment='vehicle',comparison_treatment='exposure',
                              state_positive='high',state_negative='low',minimum_group_size=1)
        self.assertEqual(result['eligibility']['positive_total_cells'],4)
        low=next(r for r in result['treatment_contrasts'] if r['stratum']['state']=='low')
        self.assertTrue(low['directional_flags']['score_rises_while_down_detection_falls'])
        self.assertTrue(low['directional_flags']['score_rises_while_up_detection_falls'])
        self.assertEqual(table.loc[1,'signed_rank'],0)
        self.assertEqual(table.loc[1,'zero_rank_baseline'],.5)
        np.testing.assert_allclose(table.signed_rank,table.centered_up+table.centered_down,atol=2e-15,rtol=0)
        ranks=[r for r in result['state_discrimination'] if r['method']=='signed_rank']
        self.assertTrue(all(r['AUROC']==1 for r in ranks))
        self.assertEqual(table.loc[0,'total_counts'],100)

    def test_metadata_reordered_by_explicit_ids_and_duplicate_ids_rejected(self):
        x,genes,cells,signature,meta=self.fixture()
        result,table=diagnose(x,genes,cells,signature,meta.iloc[::-1])
        self.assertTrue(result['eligibility']['metadata_reordered_by_cell_id'])
        self.assertEqual(table.cell_id.tolist(),cells)
        bad=meta.copy();bad.loc[1,'cell_id']='a'
        with self.assertRaisesRegex(ValueError,'duplicate'):diagnose(x,genes,cells,signature,bad)
        bad=meta.copy();bad.loc[1,'cell_id']='unknown'
        with self.assertRaisesRegex(ValueError,'match'):diagnose(x,genes,cells,signature,bad)
        with self.assertRaisesRegex(ValueError,'duplicate'):diagnose(x,['up','up','outside'],cells,signature,meta)

    def test_missing_program_feature_zero_fills_original_denominator(self):
        x,genes,cells,signature,meta=self.fixture()
        signature['up']=['up','other0']
        result,table=diagnose(x,genes,cells,signature,meta)
        self.assertEqual(result['signature']['up_genes_absent_from_input'],['other0'])
        self.assertEqual(table.loc[0,'up_detection_fraction'],.5)
        self.assertEqual(table.loc[0,'occupancy_component'],-.25)

    def test_zero_totals_small_groups_and_unavailable_method_are_explicit(self):
        x,genes,cells,signature,meta=self.fixture()
        x=x.tolil();x[:,1]=0;x=x.tocsc();meta.loc[1,'total_counts']=0
        signature['universe']=['up','down']
        result,table=diagnose(x,genes,cells,signature,meta,thresholds={'signed_rank':0})
        self.assertEqual(result['eligibility']['excluded_zero_total_cells'],1)
        self.assertTrue(table.ucell_form_1500.isna().all())
        for group in result['groups']:
            self.assertFalse(group['rate_eligible'])
            self.assertIsNone(group['supplied_thresholds']['signed_rank']['positive_fraction'])
            self.assertEqual(group['method_evaluable_cells']['ucell_form_1500'],0)
        counts=threshold_counts([0,0,1,np.nan],0,20)
        self.assertEqual((counts['positive_n'],counts['not_positive_n'],counts['tie_n']),(1,2,2))
        self.assertEqual(counts['non_evaluable_cells'],1)
        self.assertIsNone(counts['positive_fraction'])

    def test_bad_totals_and_fractional_counts_are_rejected(self):
        x,genes,cells,signature,meta=self.fixture()
        bad=meta.copy();bad['total_counts']=1
        with self.assertRaisesRegex(ValueError,'totals'):diagnose(x,genes,cells,signature,bad)
        with self.assertRaisesRegex(ValueError,'integer'):diagnose(x*.5,genes,cells,signature,meta)
        with self.assertRaisesRegex(ValueError,'dimensions'):diagnose(x[:,:2],genes,cells,signature,meta)
        with self.assertRaisesRegex(ValueError,'method-to-number'):diagnose(x,genes,cells,signature,meta,thresholds=[])
        with self.assertRaisesRegex(ValueError,'finite numbers'):diagnose(x,genes,cells,signature,meta,thresholds={'signed_rank':True})
        with self.assertRaisesRegex(ValueError,'positive integer'):diagnose(x,genes,cells,signature,meta,minimum_group_size=True)

    def test_html_escapes_metadata_and_provenance(self):
        x,genes,cells,signature,meta=self.fixture();meta['donor']='<script>alert(1)</script>'
        result,_=diagnose(x,genes,cells,signature,meta)
        result['provenance']={'label':'<img src=x onerror=alert(1)>'}
        html=report_html(result)
        self.assertNotIn('<script>',html);self.assertNotIn('<img src=x',html)
        self.assertIn('&lt;script&gt;',html);self.assertIn('&lt;img src=x',html)

    def test_explicit_pair_ids_produce_threshold_transitions_without_silent_subset(self):
        x,genes,cells,signature,meta=self.fixture();meta['pair_id']=['low_pair','low_pair','high_pair','high_pair']
        result,_=diagnose(x,genes,cells,signature,meta,thresholds={'signed_rank':-.0001},
                          baseline_treatment='vehicle',comparison_treatment='exposure',minimum_group_size=1)
        low=next(r for r in result['treatment_contrasts'] if r['stratum']['state']=='low')
        self.assertTrue(low['paired_observations']['complete_matching'])
        self.assertEqual(low['paired_observations']['transitions']['signed_rank']['negative_to_positive_n'],1)
        meta.loc[1,'pair_id']='not_matched'
        result,_=diagnose(x,genes,cells,signature,meta,thresholds={'signed_rank':0},
                          baseline_treatment='vehicle',comparison_treatment='exposure',minimum_group_size=1)
        low=next(r for r in result['treatment_contrasts'] if r['stratum']['state']=='low')
        self.assertFalse(low['paired_observations']['complete_matching'])
        self.assertEqual(low['paired_observations']['transitions'],{})

    def test_synthetic_cli_writes_external_reports_and_refuses_overwrite(self):
        with temporary_fixture(prefix='strike007_diagnostic_') as folder:
            command=[sys.executable,'-B',str(ROOT/'code/signature_diagnostic.py'),'--example',
                     '--data-root',folder,'--output-prefix','reports/tiny','--threshold','0']
            run=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            directory=Path(folder)/'reports'
            result=json.loads((directory/'tiny.json').read_text())
            self.assertTrue(result['provenance']['synthetic_example'])
            self.assertEqual(result['eligibility']['input_cells'],96)
            self.assertTrue((directory/'tiny.html').is_file())
            self.assertTrue((directory/'tiny.cells.tsv').is_file())
            repeat=subprocess.run(command,cwd=ROOT,capture_output=True,text=True)
            self.assertNotEqual(repeat.returncode,0)
            self.assertIn('already exists',repeat.stderr)

if __name__=='__main__':unittest.main()
