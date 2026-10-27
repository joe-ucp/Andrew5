"""Temporary synthetic fixtures honor local storage without hardcoded drives."""
from pathlib import Path
import os
import sys
import tempfile

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import data_paths

def temporary_fixture(prefix):
    explicit=os.environ.get('STRIKE007_TEST_TMPDIR')
    if explicit:
        parent=Path(explicit).expanduser()
        if not parent.is_absolute():raise ValueError('STRIKE007_TEST_TMPDIR must be absolute')
        parent.mkdir(parents=True,exist_ok=True)
    else:
        policy=data_paths.load_storage_policy()
        if policy.get('data_root'):
            parent=data_paths.data_path('test_tmp',for_write=True)
            parent.mkdir(parents=True,exist_ok=True)
        else:
            parent=None
    return tempfile.TemporaryDirectory(prefix=prefix,dir=parent)
