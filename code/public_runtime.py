"""Public scientific configuration and external-only execution context.

No approval records, historical run manifests or machine-specific defaults are
required. Retained inputs and newly generated inputs are explicitly distinguished.
"""
from pathlib import Path
import contextlib
import gzip
import hashlib
import json
import os
import sys

import data_paths as storage

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / 'config'
_ACTIVE = None
_AUDIT_INSTALLED = False
_ORIGINAL = {name:getattr(storage,name) for name in
             ['data_path','evidence_path','prepare_original_results','logical_evidence_key','checked_data_destination']}


def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()


def relative_path(value):
    value=str(value).replace('\\','/');path=Path(value)
    if path.is_absolute() or path.drive or '..' in path.parts or value in ('','.'):raise ValueError('A nonempty relative scientific path is required')
    return path


def registry():return json.loads((CONFIG/'scientific_inputs.json').read_text(encoding='utf-8'))


def parameters_sha256(name=None):return sha(CONFIG/'scientific_parameters.json')


def protocol(name):
    result=json.loads((CONFIG/'scientific_parameters.json').read_text(encoding='utf-8'))
    conf=ROOT/'evidence/frozen_score.json'
    if sha(conf)!=result['frozen_score_sha256']:raise ValueError('Registered scientific score has changed')
    if name not in result:raise ValueError('Unknown scientific protocol: '+name)
    return result[name]


def frozen_score():
    path=scientific_file('evidence/frozen_score.json')
    result=json.loads(path.read_text(encoding='utf-8'))
    if sha(path)!=json.loads((CONFIG/'scientific_parameters.json').read_text())['frozen_score_sha256']:
        raise ValueError('Scientific score identity mismatch')
    return result


def scientific_file(relative):
    rel=relative_path(relative).as_posix()
    if _ACTIVE and (_ACTIVE.output/rel).is_file():return _ACTIVE.generated_file(rel)
    if _ACTIVE and rel in _ACTIVE.required_generated:
        raise FileNotFoundError('Generated prerequisite is missing; retained fallback is disabled in the complete chain: '+rel)
    record=registry()['compact'].get(rel)
    if record is None:raise ValueError('Unregistered scientific compact input: '+rel)
    path=(ROOT/rel).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file():raise FileNotFoundError('Required scientific input is missing: '+rel)
    if sha(path)!=record['sha256']:raise ValueError('Scientific input checksum mismatch: '+rel)
    if _ACTIVE:_ACTIVE.reads[rel]='retained_compact'
    return path


def verify_registered_data(relative):
    rel=relative_path(relative).as_posix()
    if _ACTIVE:return _ACTIVE.verified_data(rel)
    path=_ORIGINAL['data_path'](rel,must_exist=True)
    record=registry()['data'].get(rel)
    if record is None:raise ValueError('Unregistered external scientific input: '+rel)
    _verify(path,record,rel)
    return path


def _verify(path,record,label):
    if not path.is_file():raise FileNotFoundError('Missing required scientific input: '+label)
    valid=sha(path)==record['sha256']
    if not valid and record.get('decoded_sha256') and path.suffix=='.gz':
        valid=hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest()==record['decoded_sha256']
    if not valid:raise ValueError('Scientific input checksum mismatch: '+label)


class Execution:
    """One same-code-tree execution; published inputs are read-only.

    The external directory contains generated data/evidence and small scientific
    metadata copies, never Python source or a second checkout.
    """
    def __init__(self,data_root,output_dir,stage,bulk_prefixes=(),bulk_names=()):
        if output_dir is None:raise ValueError('An explicit external --output-dir is required')
        self.output=storage.resolve_data_root(output_dir,for_write=True)
        self.source=storage.resolve_data_root(data_root,require_exists=True) if data_root is not None else None
        if self.source and (self.source==self.output or self.source.is_relative_to(self.output)):
            raise ValueError('Output must not replace or contain the source data root')
        self.stage=stage;self.prefixes=tuple(bulk_prefixes);self.names=set(bulk_names)
        self.reads={};self.verified={};self.records=registry();self.old=[]
        self.required_generated=set();self.required_data_prefixes=()
        self.output.mkdir(parents=True,exist_ok=True)
        self.state_path=self.output/'execution_inputs.json'
        self.generated=json.loads(self.state_path.read_text()) if self.state_path.exists() else {}
        self.initial=dict(self.generated)
        self.initial_paths=set()
        for path in self.output.rglob('*'):
            if not path.is_file() or path==self.state_path:continue
            relative=path.relative_to(self.output).as_posix()
            self.output_path(relative)
            if relative.startswith('.cache/'):continue
            self.initial_paths.add(relative)
            if relative not in self.initial:raise ValueError('Unregistered preexisting execution file: '+relative)
        self.verify_prior_outputs()

    def output_path(self,relative):
        path=(self.output/relative_path(relative)).resolve()
        if not path.is_relative_to(self.output) or path==self.output:
            raise ValueError('Generated path escapes the configured output directory: '+str(path))
        return path

    def verify_prior_outputs(self):
        for relative,record in self.initial.items():
            path=self.output_path(relative)
            if not path.is_file() or sha(path)!=record['sha256']:
                raise ValueError('Earlier generated output is missing or changed: '+relative)

    def owns(self,relative):
        return relative in self.names or any(relative==p.rstrip('/') or relative.startswith(p.rstrip('/')+'/') for p in self.prefixes)

    def generated_file(self,relative):
        path=self.output_path(relative)
        record=self.generated.get(relative)
        if record is None and relative in self.initial_paths:
            raise ValueError('Unregistered file in generated workspace: '+relative)
        # Current-stage outputs may be read before the end-of-stage inventory.
        if record and sha(path)!=record['sha256']:raise ValueError('Generated scientific input changed: '+relative)
        self.reads[relative]=record.get('origin','generated_in_output_workspace') if record else 'generated_in_output_workspace'
        return path

    def verified_data(self,relative):
        target=self.output_path('data/'+relative)
        if target.is_file():return self.generated_file('data/'+relative)
        if 'data/'+relative in self.required_generated or any(relative==p or relative.startswith(p+'/') for p in self.required_data_prefixes):
            raise FileNotFoundError('Generated prerequisite is missing; retained fallback is disabled in the complete chain: '+relative)
        if self.source is None:raise FileNotFoundError('This stage requires --data-root: '+relative)
        target=_ORIGINAL['data_path'](relative,self.source,must_exist=True)
        record=self.records['data'].get(relative)
        if record is None:raise ValueError('Unregistered external scientific input: '+relative)
        if relative not in self.verified:
            _verify(target,record,relative);self.verified[relative]=record['sha256']
        self.reads['data/'+relative]='registered_retained_external'
        return target

    def data_path(self,relative,data_root=None,*,must_exist=False,for_write=False):
        rel=relative_path(relative).as_posix();target=self.output_path('data/'+rel)
        # The original sparse loader marks a cache path writable before checking
        # whether it exists. Reuse a verified existing cache without modifying it.
        if for_write and rel.startswith('data_primary/sparse_cache/') and self.source:
            existing=self.source/rel
            if existing.is_file():return self.verified_data(rel)
        if must_exist:return self.verified_data(rel)
        if for_write or self.owns(rel):return target
        if target.exists():return target
        if self.source is None:return target
        return _ORIGINAL['data_path'](rel,self.source)

    def evidence_path(self,name,repo_root=None,data_root=None,*,must_exist=False,for_write=False):
        rel=relative_path(name).as_posix()
        if rel in storage.ORIGINAL_BULK_NAMES:
            datarel='original_results/'+rel
            if datarel in self.names or for_write:return self.output_path('data/'+datarel)
            return self.verified_data(datarel)
        if for_write:return self.output_path('evidence/'+rel)
        target=self.output_path('evidence/'+rel)
        if target.is_file():return self.generated_file('evidence/'+rel)
        return scientific_file('evidence/'+rel)

    def prepare_original_results(self,data_root=None):
        path=self.output/'data/original_results';path.mkdir(parents=True,exist_ok=True);return path

    def logical_evidence_key(self,path,repo_root=None,data_root=None):
        path=Path(path).resolve()
        for base in [self.output/'evidence',ROOT/'evidence']:
            if path.is_relative_to(base):return 'evidence/'+path.relative_to(base).as_posix()
        if path.name in storage.ORIGINAL_BULK_NAMES:return 'evidence/'+path.name
        raise ValueError('Not a registered evidence path: '+str(path))

    def checked_data_destination(self,path,data_root=None):
        path=Path(path).resolve()
        if path.is_relative_to(self.output) and path!=self.output:return path
        return _ORIGINAL['checked_data_destination'](path,data_root or self.source)

    def seed(self,relative):
        """Copy one small registered input, never result tables owned by a stage."""
        import shutil
        rel=relative_path(relative).as_posix();target=self.output_path(rel)
        if target.exists():return self.generated_file(rel)
        source=scientific_file(rel)
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,target)
        self.generated[rel]={'sha256':sha(target),'origin':'copied_registered_compact_input'}
        return target

    def inventory(self):
        self.verify_prior_outputs()
        for p in self.output.rglob('*'):
            if p.is_file() and p!=self.state_path and not p.name.endswith('.log') and '/.cache/' not in p.as_posix():
                rel=p.relative_to(self.output).as_posix()
                self.output_path(rel)
                if rel not in self.initial:self.generated[rel]={'sha256':sha(p),'origin':self.generated.get(rel,{}).get('origin','generated:'+self.stage)}
        self.state_path.write_text(json.dumps(self.generated,indent=2)+'\n',encoding='utf-8')

    def bind_modules(self):
        """Redirect imported path aliases and output globals; never numerical code."""
        replacements={name:getattr(self,name) for name in _ORIGINAL}
        for module in list(sys.modules.values()):
            filename=getattr(module,'__file__',None)
            if not filename:continue
            path=Path(filename).resolve()
            if not (path.is_relative_to(ROOT/'code') or path.is_relative_to(ROOT/'scripts_primary')):continue
            if module.__name__ in {'public_runtime','data_paths','reproduce'}:continue
            for name,new in replacements.items():
                current=getattr(module,name,None)
                if current in _ORIGINAL.values() or isinstance(getattr(current,'__self__',None),Execution):
                    self.old.append((module,name,_ORIGINAL[name]));setattr(module,name,new)
            for name in ['EV','OUT','OLD','EVIDENCE','FIG','REV']:
                value=getattr(module,name,None)
                if isinstance(value,Path) and value.is_relative_to(ROOT) and value!=ROOT:
                    relative=value.relative_to(ROOT)
                    if relative.parts[0] in {'evidence','manuscript'}:
                        self.old.append((module,name,value));setattr(module,name,self.output/relative)

    def __enter__(self):
        global _ACTIVE,_AUDIT_INSTALLED
        if _ACTIVE is not None:raise RuntimeError('Nested public execution is unsupported')
        _ACTIVE=self
        if not _AUDIT_INSTALLED:sys.addaudithook(_protect_inputs);_AUDIT_INSTALLED=True
        for name in _ORIGINAL:setattr(storage,name,getattr(self,name))
        self.bind_modules()
        return self

    def __exit__(self,*unused):
        global _ACTIVE
        for module,name,old in reversed(self.old):setattr(module,name,old)
        for module in list(sys.modules.values()):
            filename=getattr(module,'__file__',None)
            if not filename:continue
            path=Path(filename).resolve()
            if not (path.is_relative_to(ROOT/'code') or path.is_relative_to(ROOT/'scripts_primary')):continue
            if module.__name__ in {'public_runtime','data_paths','reproduce'}:continue
            for name,original in _ORIGINAL.items():
                if getattr(getattr(module,name,None),'__self__',None) is self:setattr(module,name,original)
            for name in ['EV','OUT','OLD','EVIDENCE','FIG','REV']:
                value=getattr(module,name,None)
                if isinstance(value,Path) and value.is_relative_to(self.output):
                    relative=value.relative_to(self.output)
                    if relative.parts and relative.parts[0] in {'evidence','manuscript'}:setattr(module,name,ROOT/relative)
        for name,old in _ORIGINAL.items():setattr(storage,name,old)
        _ACTIVE=None


def _protect_inputs(event,args):
    """Fail closed on accidental repository/source writes from scientific runners."""
    context=_ACTIVE
    if context is None:return
    path=None
    if event=='open' and isinstance(args[0],(str,bytes,os.PathLike)):
        mode=args[1];flags=args[2]
        writing=(isinstance(mode,str) and any(x in mode for x in 'wax+')) or (isinstance(flags,int) and flags & (os.O_WRONLY|os.O_RDWR|os.O_CREAT|os.O_TRUNC))
        if writing:path=args[0]
    elif event in {'os.remove','os.rmdir','os.mkdir','os.rename','os.chmod','os.truncate'} and args:path=args[0]
    elif event in {'os.symlink','os.link'} and len(args)>1:path=args[1]
    if path is None:return
    if not isinstance(path,(str,bytes,os.PathLike)):return
    target=Path(os.fsdecode(path)).resolve()
    if not target.is_relative_to(context.output):
        raise PermissionError('Scientific execution writes must stay inside the resolved external output directory: '+str(target))
    if target.relative_to(context.output).as_posix() in context.initial:
        raise PermissionError('Earlier generated outputs are immutable inputs: '+str(target))
    if event=='os.rename' and len(args)>1:
        target=Path(os.fsdecode(args[1])).resolve()
        if not target.is_relative_to(context.output):raise PermissionError('Rename destination escapes external output directory: '+str(target))
        if target.relative_to(context.output).as_posix() in context.initial:raise PermissionError('Cannot replace an earlier generated output: '+str(target))


def execute_stage(stage_name,data_root=None,output_dir=None):
    from reproduce import run_stage
    return run_stage(stage_name,data_root,output_dir)
