"""Storage routing tests only: no research data, downloads, or scientific reruns."""
from pathlib import Path
import importlib
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

ROOT=Path(__file__).resolve().parents[1]
sys.dont_write_bytecode=True
sys.path[:0]=[str(ROOT/'code'),str(ROOT/'scripts_primary')]
import data_paths as paths
from test_support import temporary_fixture

class DataRootTests(unittest.TestCase):
    def setUp(self):
        self.temp=temporary_fixture(prefix='strike007_path_tests_')
        self.base=Path(self.temp.name).resolve()
        self.repo=self.base/'repository';self.repo.mkdir()
        self.external=self.base/'external';self.external.mkdir()
        self.override=self.base/'override';self.override.mkdir()
        self.patch=mock.patch.object(paths,'REPO_ROOT',self.repo);self.patch.start()
        self.env=mock.patch.dict(os.environ,{},clear=False);self.env.start()
        os.environ.pop(paths.ENV_NAME,None)
        self.policy=mock.patch.object(paths,'load_storage_policy',return_value={});self.policy.start()
    def tearDown(self):
        self.policy.stop();self.env.stop();self.patch.stop();self.temp.cleanup()
    def test_unconfigured_public_use_requires_explicit_external_root(self):
        self.assertIsNone(paths.DEFAULT_DATA_ROOT)
        with self.assertRaisesRegex(ValueError,'No external data root configured'):paths.resolve_data_root()
    def test_policy_default_is_overridden_by_environment_and_argument(self):
        with mock.patch.object(paths,'load_storage_policy',return_value={'data_root':str(self.external)}):
            self.assertEqual(paths.resolve_data_root(),self.external)
            os.environ[paths.ENV_NAME]=str(self.override)
            self.assertEqual(paths.resolve_data_root(),self.override)
            self.assertEqual(paths.resolve_data_root(self.external),self.external)
    def test_explicit_overrides_environment(self):
        os.environ[paths.ENV_NAME]=str(self.external)
        self.assertEqual(paths.resolve_data_root(require_exists=True),self.external)
        self.assertEqual(paths.resolve_data_root(self.override,require_exists=True),self.override)
    @unittest.skipUnless(os.name=='nt','Windows portable configuration')
    def test_portable_configuration_allows_other_external_drives(self):
        with mock.patch.object(paths,'load_storage_policy',return_value={}):
            self.assertEqual(paths.resolve_data_root(r'E:\ExternalResearch'),Path(r'E:\ExternalResearch'))
            with self.assertRaisesRegex(ValueError,'must not overlap'):
                paths.resolve_data_root(self.repo)
    def test_invalid_environment_never_falls_back(self):
        os.environ[paths.ENV_NAME]=str(self.repo/'cache')
        with self.assertRaisesRegex(ValueError,'must not overlap'):paths.resolve_data_root()
        self.assertFalse((self.repo/'cache').exists())
    def test_repository_and_nested_normalized_paths_rejected(self):
        for path in [self.repo,self.repo/'cache',self.repo/'code'/'..'/'cache',self.repo.parent]:
            with self.subTest(path=path),self.assertRaisesRegex(ValueError,'must not overlap'):
                paths.resolve_data_root(path)
    def test_empty_or_relative_override_rejected(self):
        for path in ['', 'relative_cache']:
            with self.subTest(path=path),self.assertRaises(ValueError):paths.resolve_data_root(path)
    @unittest.skipUnless(os.name=='nt','Windows D: policy')
    def test_non_D_windows_overrides_rejected(self):
        with mock.patch.object(paths,'load_storage_policy',return_value={'windows_allowed_drives':['D:']}):
            for path in [r'C:\SomeCache',r'E:\SomeCache',r'Z:\SomeCache',r'\\server\share\SomeCache']:
                with self.subTest(path=path),self.assertRaisesRegex(ValueError,'must live on D:'):
                    paths.resolve_data_root(path,for_write=True)
                os.environ[paths.ENV_NAME]=path
                with self.subTest(environment=path),self.assertRaisesRegex(ValueError,'must live on D:'):
                    paths.resolve_data_root(for_write=True)
    @unittest.skipUnless(os.name=='nt','Windows resolved-drive policy')
    def test_D_junction_resolving_to_C_is_rejected(self):
        original=Path.resolve
        selected=Path(r'D:\ExternalResearch\link')
        def simulated_resolve(path,*args,**kwargs):
            return Path(r'C:\OutsideCache') if path==selected else original(path,*args,**kwargs)
        with mock.patch.object(paths,'load_storage_policy',return_value={'windows_allowed_drives':['D:']}),mock.patch.object(Path,'resolve',simulated_resolve):
            with self.assertRaisesRegex(ValueError,'must live on D:'):paths.resolve_data_root(selected)
    def test_missing_read_root_fails_without_creation(self):
        missing=self.base/'missing'
        with self.assertRaisesRegex(FileNotFoundError,'No data will be created'):paths.resolve_data_root(missing,require_exists=True)
        self.assertFalse(missing.exists())
    @unittest.skipUnless(os.name=='nt','Windows drive availability check')
    def test_missing_drive_fails_without_fallback(self):
        original=Path.is_dir
        def unavailable_drive(path):return False if path==Path(self.base.anchor) else original(path)
        with mock.patch.object(Path,'is_dir',unavailable_drive):
            with self.assertRaisesRegex(FileNotFoundError,'drive/share is unavailable'):
                paths.resolve_data_root(self.base/'missing',for_write=True)
    def test_manifest_escape_and_absolute_paths_rejected(self):
        for rel in ['../escape','..\\escape',str(self.repo/'file'),'C:relative','C:/absolute']:
            with self.subTest(relative=rel),self.assertRaises(ValueError):paths.data_path(rel,self.external)
        self.assertEqual(paths.data_path('data_primary/example.mtx.gz',self.external),self.external/'data_primary/example.mtx.gz')
    def test_missing_file_is_a_clear_error(self):
        with self.assertRaisesRegex(FileNotFoundError,'Required external data is missing'):
            paths.data_path('data_primary/missing.gz',self.external,must_exist=True)
        self.assertFalse((self.external/'data_primary').exists())
    def test_direct_writer_destination_is_guarded(self):
        with self.assertRaisesRegex(ValueError,'escapes'):
            paths.checked_data_destination(self.repo/'data_primary/file',self.external)
    def test_resolved_directory_link_into_repo_rejected(self):
        link=self.external/'linked_repo'
        try:link.symlink_to(self.repo,target_is_directory=True)
        except OSError as error:self.skipTest(f'Symlink unavailable: {error}')
        with self.assertRaisesRegex(ValueError,'must not overlap'):paths.resolve_data_root(link)
        with self.assertRaisesRegex(ValueError,'escapes'):paths.data_path('linked_repo/data',self.external)
    def test_partial_file_link_cannot_target_repo(self):
        target=self.repo/'sentinel';target.write_text('unchanged')
        link=self.external/'example.gz.part'
        try:link.symlink_to(target)
        except OSError as error:self.skipTest(f'Symlink unavailable: {error}')
        with self.assertRaisesRegex(ValueError,'escapes'):paths.checked_data_destination(link,self.external)
        self.assertEqual(target.read_text(),'unchanged')
    def test_original_bulk_routes_externally_without_stale_fallback(self):
        old=self.repo/'evidence/primary_cell_scores.tsv.gz'
        old.parent.mkdir();old.write_bytes(b'legacy sentinel')
        target=paths.evidence_path(old.name,data_root=self.external)
        self.assertEqual(target,self.external/'original_results'/old.name)
        self.assertFalse(target.parent.exists())
        with self.assertRaisesRegex(FileNotFoundError,'Required external data is missing'):
            paths.evidence_path(old.name,data_root=self.external,must_exist=True)
        self.assertEqual(old.read_bytes(),b'legacy sentinel')
    def test_compact_evidence_stays_in_selected_repository(self):
        self.assertEqual(paths.evidence_path('revision_01/operating_points.tsv',repo_root=self.repo),
                         self.repo/'evidence/revision_01/operating_points.tsv')
        self.assertFalse((self.repo/'evidence').exists())
    def test_evidence_names_reject_traversal_and_absolute_paths(self):
        for name in ['../escape','..\\escape',str(self.repo/'absolute'),'C:relative','', '.']:
            with self.subTest(name=name),self.assertRaises(ValueError):
                paths.evidence_path(name,data_root=self.external)
    def test_logical_evidence_keys_do_not_claim_physical_repo_location(self):
        bulk=paths.evidence_path('reference_cell_scores.tsv',data_root=self.external)
        self.assertEqual(paths.logical_evidence_key(bulk,data_root=self.external),'evidence/reference_cell_scores.tsv')
        self.assertEqual(paths.logical_evidence_key(self.repo/'evidence/revision_01/result.tsv'),
                         'evidence/revision_01/result.tsv')
        with self.assertRaises(ValueError):paths.logical_evidence_key(self.external/'unregistered.tsv',data_root=self.external)
    def test_bulk_writer_prepares_only_valid_external_directory(self):
        with self.assertRaises(ValueError):paths.prepare_original_results(self.repo)
        self.assertFalse((self.repo/'original_results').exists())
        actual=paths.prepare_original_results(self.external)
        self.assertEqual(actual,self.external/'original_results')
        self.assertTrue(actual.is_dir())
    def test_original_manifest_hash_checks_and_explicit_gzip_wrapper_equivalence(self):
        import gzip,hashlib,json
        original=paths.prepare_original_results(self.external);rows=[]
        for name in sorted(paths.ORIGINAL_BULK_NAMES):
            decoded=b'column\nretained value\n';content=gzip.compress(decoded,mtime=1) if name.endswith('.gz') else decoded
            (original/name).write_bytes(content)
            rows.append({'logical_path':'evidence/'+name,'data_relative_path':'original_results/'+name,
                         'bytes':len(content),'sha256':hashlib.sha256(content).hexdigest(),
                         'decoded_sha256':hashlib.sha256(decoded).hexdigest()})
        (self.repo/'evidence').mkdir()
        (self.repo/'evidence/original_bulk_manifest.json').write_text(json.dumps({'files':rows}))
        self.assertTrue(all(r['original_bytes_identical'] for r in paths.verify_original_bulk_manifest(data_root=self.external)))
        path=original/'primary_cell_scores.tsv.gz';path.write_bytes(gzip.compress(decoded,mtime=2))
        with self.assertRaisesRegex(RuntimeError,'checksum mismatch'):paths.verify_original_bulk_manifest(data_root=self.external)
        checked=paths.verify_original_bulk_manifest(data_root=self.external,allow_reconstructed_gzip=True)
        self.assertEqual(sum(r['status']=='verified_decoded_gzip_wrapper_differs' for r in checked),1)
        path.write_bytes(gzip.compress(b'different science',mtime=2))
        with self.assertRaisesRegex(RuntimeError,'checksum mismatch'):
            paths.verify_original_bulk_manifest(data_root=self.external,allow_reconstructed_gzip=True)

class ExternalPolicyTests(unittest.TestCase):
    def setUp(self):
        self.temp=temporary_fixture(prefix='strike007_policy_tests_')
        self.base=Path(self.temp.name).resolve()
        self.repo=self.base/'repository';self.repo.mkdir()
        self.config=self.base/'config'/'storage.json';self.config.parent.mkdir()
        self.env=mock.patch.dict(os.environ,{},clear=False);self.env.start()
        os.environ.pop(paths.POLICY_ENV_NAME,None)
        self.patch=mock.patch.object(paths,'REPO_ROOT',self.repo);self.patch.start()
    def tearDown(self):
        self.patch.stop();self.env.stop();self.temp.cleanup()
    def test_explicit_policy_is_absolute_external_and_not_silently_missing(self):
        for invalid in ['', 'relative.json',str(self.repo/'policy.json')]:
            os.environ[paths.POLICY_ENV_NAME]=invalid
            with self.subTest(path=invalid),self.assertRaises(ValueError):paths.load_storage_policy()
        os.environ[paths.POLICY_ENV_NAME]=str(self.config)
        with self.assertRaisesRegex(FileNotFoundError,'Explicit storage policy is missing'):paths.load_storage_policy()
        self.config.write_text('{"windows_allowed_drives":["D:"]}',encoding='utf-8')
        self.assertEqual(paths.load_storage_policy()['windows_allowed_drives'],['D:'])
    def test_user_location_and_repository_policy_is_ignored(self):
        (self.repo/'.storage.local.json').write_text('{"data_root":"wrong"}',encoding='utf-8')
        if os.name=='nt':
            os.environ['LOCALAPPDATA']=str(self.base/'user_config')
            expected=self.base/'user_config/UCP/STRIKE007/storage.json'
        else:
            os.environ['XDG_CONFIG_HOME']=str(self.base/'user_config')
            expected=self.base/'user_config/ucp/strike007/storage.json'
        self.assertEqual(paths.storage_policy_path(),expected)
        self.assertEqual(paths.load_storage_policy(),{})
        self.assertFalse(expected.parent.exists())
    def test_malformed_policy_is_rejected_without_fallback(self):
        os.environ[paths.POLICY_ENV_NAME]=str(self.config)
        for invalid in ['{', '[]','{"data_root":""}', '{"windows_allowed_drives":"D:"}', '{"windows_allowed_drives":["D"]}']:
            self.config.write_text(invalid,encoding='utf-8')
            with self.subTest(policy=invalid),self.assertRaises(ValueError):paths.load_storage_policy()

class EntryPointTests(unittest.TestCase):
    def test_imports_do_not_download_or_create_directories(self):
        # Preload third-party libraries so the guard tests this repository's imports.
        import numpy,pandas,scipy,matplotlib.pyplot
        names=['fetch_sources','audit_primary','audit_challenge','reference_analysis','primary_analysis','biological_challenge','calibration_diagnostics','depth_stability','countermodel_stability','make_figures','public_runtime','reproduce','revision03_influence']
        with mock.patch.dict(os.environ,{paths.ENV_NAME:str(ROOT)}),mock.patch('urllib.request.urlopen',side_effect=AssertionError('Network access during import')),mock.patch.object(Path,'mkdir',side_effect=AssertionError('Directory creation during import')):
            for name in names:
                with self.subTest(module=name):importlib.import_module(name)
    def test_download_helpers_reject_repo_before_network(self):
        import fetch_sources
        with mock.patch('urllib.request.urlopen',side_effect=AssertionError('Network must not run')),mock.patch.object(Path,'mkdir',side_effect=AssertionError('mkdir must not run')):
            with self.assertRaises(ValueError):fetch_sources.download({},ROOT/'data_primary/forbidden.gz',ROOT)
    def test_invalid_cli_root_fails_before_outputs(self):
        sentinels=[ROOT/'evidence/reconstruction_receipt.json',ROOT/'evidence/source_verification.json']
        before={p:p.read_bytes() for p in sentinels if p.exists()}
        commands=[['code/fetch_sources.py','--fetch']]
        for command in commands:
            result=subprocess.run([sys.executable,'-B',str(ROOT/command[0]),*command[1:],'--data-root',str(ROOT)],cwd=ROOT,capture_output=True,text=True)
            with self.subTest(command=command):
                self.assertNotEqual(result.returncode,0)
                self.assertIn('must not overlap',result.stderr)
        for path,content in before.items():self.assertEqual(path.read_bytes(),content)
    def test_unconfigured_cli_fails_helpfully_without_loading_fixtures(self):
        with temporary_fixture(prefix='strike007_unconfigured_') as folder:
            env=dict(os.environ)
            env.pop(paths.ENV_NAME,None);env.pop(paths.POLICY_ENV_NAME,None)
            env['LOCALAPPDATA']=folder;env['XDG_CONFIG_HOME']=folder
            result=subprocess.run([sys.executable,'-B','-c',"import sys;sys.path.insert(0,'code');from data_paths import resolve_data_root;resolve_data_root()"],
                                  cwd=ROOT,env=env,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('No external data root configured',result.stderr)
            self.assertIn('--data-root',result.stderr)
            self.assertNotIn('FileNotFoundError',result.stderr)
    def test_all_owned_python_sources_compile(self):
        for folder in ['code','scripts_primary']:
            for path in (ROOT/folder).glob('*.py'):
                with self.subTest(path=path.name):compile(path.read_text(encoding='utf-8'),str(path),'exec')

if __name__=='__main__':unittest.main(verbosity=2)
