"""Portable synthetic tests; --external-data explicitly enables retained regression."""
from pathlib import Path
import argparse
import datetime
import importlib.metadata
import json
import os
import platform
import sys
import time
import unittest

ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--external-data',action='store_true',help='Also verify the retained, hash-registered zero-null fixture')
    parser.add_argument('--data-root',type=Path,help='Absolute external data root; required when no user policy/environment exists')
    parser.add_argument('--receipt',type=Path,help='Optional compact JSON result; no file is written by default')
    args=parser.parse_args()
    sys.dont_write_bytecode=True
    os.environ['PYTHONDONTWRITEBYTECODE']='1'
    if args.external_data or args.data_root is not None:
        from data_paths import resolve_data_root,ENV_NAME
        try:root=resolve_data_root(args.data_root,require_exists=args.external_data)
        except (ValueError,OSError) as error:parser.error(str(error))
        os.environ[ENV_NAME]=str(root)
    os.environ['STRIKE007_EXTERNAL_TESTS']='1' if args.external_data else '0'
    started=time.monotonic()
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.discover(str(ROOT/'tests')))
    packages={}
    for name in ['numpy','pandas','scipy','matplotlib']:
        try:packages[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:packages[name]='not installed'
    receipt={'status':'PASS' if result.wasSuccessful() else 'FAIL',
             'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
             'python':platform.python_version(),'platform':platform.platform(),
             'kernel':platform.release(),'packages':packages,'external_data_enabled':args.external_data,
             'tests_run':result.testsRun,'failure_count':len(result.failures),'error_count':len(result.errors),
             'skipped':[{'test':str(test),'reason':reason} for test,reason in result.skipped],
             'seconds':time.monotonic()-started}
    if args.receipt:
        args.receipt.parent.mkdir(parents=True,exist_ok=True)
        args.receipt.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,indent=2))
    return 0 if result.wasSuccessful() else 1

if __name__=='__main__':raise SystemExit(main())
