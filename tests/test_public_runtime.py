"""Scientific execution must use verified inputs and protect published results."""
from pathlib import Path
import hashlib
import json
import os
import sys
import unittest
from unittest import mock

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'code'))
import public_runtime as runtime
import reproduce
from test_support import temporary_fixture


class PublicExecutionTests(unittest.TestCase):
    def test_required_external_output_rejects_repository(self):
        with self.assertRaises(ValueError):runtime.Execution(None,ROOT,'fixture')
        with self.assertRaises(ValueError):runtime.Execution(None,None,'fixture')

    def test_registered_data_and_generated_precedence(self):
        with temporary_fixture(prefix='public_runtime_') as folder:
            base=Path(folder);source=base/'source';source.mkdir()
            raw=source/'input.txt';raw.write_text('retained')
            expected=hashlib.sha256(raw.read_bytes()).hexdigest()
            records={'compact':{},'data':{'input.txt':{'sha256':expected,'kind':'derived_result'}}}
            with mock.patch.object(runtime,'registry',return_value=records):
                execution=runtime.Execution(source,base/'outputs','fixture')
                with execution:
                    self.assertEqual(runtime.verify_registered_data('input.txt'),raw)
                    generated=execution.output/'data/input.txt';generated.parent.mkdir();generated.write_text('new computation')
                    execution.inventory()
                    self.assertEqual(runtime.verify_registered_data('input.txt'),generated)
                    generated.write_text('tampered')
                    with self.assertRaises(ValueError):runtime.verify_registered_data('input.txt')
                self.assertEqual(raw.read_text(),'retained')

    def test_checksum_failure_and_read_only_inputs(self):
        with temporary_fixture(prefix='public_readonly_') as folder:
            base=Path(folder);source=base/'source';source.mkdir()
            raw=source/'input.txt';raw.write_text('retained')
            with mock.patch.object(runtime,'registry',return_value={'compact':{},'data':{'input.txt':{'sha256':'0'*64,'kind':'source'}}}):
                execution=runtime.Execution(source,base/'outputs','fixture')
                with execution:
                    with self.assertRaises(ValueError):runtime.verify_registered_data('input.txt')
                    with self.assertRaises(PermissionError):raw.write_text('overwrite')
                    with self.assertRaises(PermissionError):(ROOT/'should_never_exist.txt').write_text('overwrite')
                self.assertEqual(raw.read_text(),'retained')

    def test_complete_chain_requires_completed_dependencies(self):
        with temporary_fixture(prefix='public_chain_') as folder:
            execution=runtime.Execution(None,Path(folder)/'outputs','calibration')
            with execution,self.assertRaises(FileNotFoundError):reproduce.seed_inputs(execution,'calibration',True)

    def test_complete_chain_never_substitutes_retained_for_missing_generated_table(self):
        with temporary_fixture(prefix='public_missing_') as folder:
            output=Path(folder)/'output'
            first=runtime.Execution(None,output,'reference')
            with first:
                (output/'receipts').mkdir()
                (output/'receipts/reference.json').write_text(json.dumps({'status':'PASS'}))
                first.inventory()
            second=runtime.Execution(None,output,'calibration')
            with second:
                reproduce.seed_inputs(second,'calibration',True)
                with self.assertRaises(FileNotFoundError):runtime.scientific_file('evidence/calibration_results.tsv')

    def test_unregistered_preexisting_files_are_not_adopted(self):
        with temporary_fixture(prefix='public_stale_') as folder:
            output=Path(folder)/'output';output.mkdir()
            (output/'unexpected.tsv').write_text('stale')
            with self.assertRaises(ValueError):runtime.Execution(None,output,'fixture')

    def test_previous_inventory_hashes_are_never_refreshed_on_tampering(self):
        with temporary_fixture(prefix='public_hash_') as folder:
            output=Path(folder)/'output'
            first=runtime.Execution(None,output,'fixture')
            with first:
                (output/'first.tsv').write_text('original');first.inventory()
            second=runtime.Execution(None,output,'next')
            (output/'first.tsv').write_text('changed')
            with self.assertRaises(ValueError):second.inventory()
            with self.assertRaises(ValueError):runtime.Execution(None,output,'later')
            state=json.loads((output/'execution_inputs.json').read_text())
            self.assertEqual(state['first.tsv']['sha256'],hashlib.sha256(b'original').hexdigest())

    def test_unrelated_writes_and_symlink_escapes_are_rejected(self):
        with temporary_fixture(prefix='public_escape_') as folder:
            base=Path(folder);outside=base/'outside';outside.mkdir()
            execution=runtime.Execution(None,base/'output','fixture')
            with execution:
                with self.assertRaises(PermissionError):(outside/'file').write_text('bad')
            try:os.symlink(outside,execution.output/'escaped',target_is_directory=True)
            except OSError as error:self.skipTest('Symlink unavailable: '+str(error))
            with execution:
                with self.assertRaises(ValueError):execution.output_path('escaped/new.tsv')
                with self.assertRaises(PermissionError):(execution.output/'escaped/new.tsv').write_text('bad')

    def test_imported_aliases_restore_between_stage_contexts(self):
        import reference_analysis
        with temporary_fixture(prefix='public_alias_') as folder:
            base=Path(folder)
            for name in ['first','second','third']:
                execution=runtime.Execution(None,base/name,'fixture')
                with execution:
                    self.assertIs(reference_analysis.data_path.__self__,execution)
                    self.assertEqual(reference_analysis.EV,execution.output/'evidence')
                self.assertIs(reference_analysis.data_path,runtime._ORIGINAL['data_path'])
                self.assertEqual(reference_analysis.EV,ROOT/'evidence')

    def test_protocol_uses_scientific_configuration_only(self):
        actual=runtime.protocol('revision_01')
        self.assertEqual(actual['thinning']['seeds'],list(range(7600,7605)))
        self.assertEqual(len(runtime.frozen_score()['up']),50)
        self.assertNotIn('governance',actual)
        for name in ['revision_01','revision_02']:
            self.assertNotIn('revision_gate',json.dumps(runtime.protocol(name)))


if __name__=='__main__':unittest.main()
