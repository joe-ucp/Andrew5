"""External research-data paths. Resolution never creates folders or downloads data."""
from pathlib import Path, PureWindowsPath
import argparse
import json
import os

REPO_ROOT=Path(__file__).resolve().parents[1]
POLICY_ENV_NAME='STRIKE007_STORAGE_POLICY'
# Kept for older CLI imports; there is no public machine-specific default.
DEFAULT_DATA_ROOT=None
ENV_NAME='STRIKE007_CACHE'
HINT='Set STRIKE007_CACHE or pass --data-root to an absolute external data folder, or configure an external user storage policy.'

def storage_policy_path():
    """Locate configuration outside the checkout without creating anything."""
    override=os.environ.get(POLICY_ENV_NAME)
    if override is not None:
        path=Path(override).expanduser()
        if not override.strip() or not path.is_absolute():
            raise ValueError(f'{POLICY_ENV_NAME} must name an absolute external JSON file')
    elif os.name=='nt':
        base=Path(os.environ.get('LOCALAPPDATA',Path.home()/'AppData'/'Local'))
        path=base/'UCP'/'STRIKE007'/'storage.json'
    else:
        base=Path(os.environ.get('XDG_CONFIG_HOME',Path.home()/'.config')).expanduser()
        path=base/'ucp'/'strike007'/'storage.json'
    if not path.is_absolute():
        raise ValueError(f'User configuration directory must be absolute: {path}')
    path=path.resolve()
    if path==REPO_ROOT or path.is_relative_to(REPO_ROOT):
        raise ValueError(f'Storage policy must be outside the code repository: {path}')
    return path

def load_storage_policy():
    """Read only at resolution time, so importing pure helpers needs no fixtures."""
    path=storage_policy_path()
    if not path.is_file():
        if POLICY_ENV_NAME in os.environ:
            raise FileNotFoundError(f'Explicit storage policy is missing: {path}')
        return {}
    try:
        policy=json.loads(path.read_text(encoding='utf-8-sig'))
    except (ValueError,OSError) as error:
        raise ValueError(f'Cannot read storage policy {path}: {error}') from error
    if not isinstance(policy,dict):raise ValueError('Storage policy must be a JSON object')
    drives=policy.get('windows_allowed_drives',[])
    if not isinstance(drives,list) or any(not isinstance(d,str) or len(d)!=2 or not d[0].isalpha() or d[1]!=':' for d in drives):
        raise ValueError('windows_allowed_drives must be a list of drive letters, for example ["D:"]')
    if 'data_root' in policy and (not isinstance(policy['data_root'],str) or not policy['data_root'].strip()):
        raise ValueError('Storage policy data_root must be a nonempty absolute path string')
    return policy

# Historical evidence/ keys stay stable; their bulk payloads live externally.
ORIGINAL_BULK_NAMES=frozenset({
    'primary_cell_scores.tsv.gz','biological_challenge_cell_scores.tsv.gz',
    'reference_cell_scores.tsv','calibration_simulations.tsv.gz',
    'countermodel_stability_pools.tsv.gz','countermodel_stability_group_draws.tsv',
    'countermodel_stability_batches.tsv','depth_stability_replicates.tsv',
    'depth_stability_sample_means.tsv','reference_feature_selection.tsv',
})

def resolve_data_root(value=None, *, require_exists=False, for_write=False, repo_root=None):
    """Explicit argument > environment > external policy; reject repository overlap.

    Existing junctions/symlinks are resolved before containment checks. Reads can
    require a populated folder; writers require an available drive/share but
    create directories only at their explicit operation, never during import.
    """
    policy=load_storage_policy()
    selected=value if value is not None else os.environ.get(ENV_NAME,policy.get('data_root'))
    if selected is None:raise ValueError(f'No external data root configured. {HINT}')
    if not str(selected).strip():raise ValueError(f'Empty external data root. {HINT}')
    root=Path(selected).expanduser()
    if not root.is_absolute():raise ValueError(f'External data root must be absolute: {root}. {HINT}')
    requested_drive=root.drive.upper()
    root=root.resolve();repo=Path(repo_root or REPO_ROOT).resolve()
    if root==repo or root.is_relative_to(repo) or repo.is_relative_to(root):
        raise ValueError(f'External data root must not overlap the code repository: {root}. {HINT}')
    allowed=tuple(d.upper() for d in policy.get('windows_allowed_drives',[]))
    if os.name=='nt' and allowed and (requested_drive not in allowed or root.drive.upper() not in allowed):
        raise ValueError(f'Research data must live on {", ".join(allowed)} under local Windows policy; requested {selected}, resolved {root}. {HINT}')
    if root.exists() and not root.is_dir():raise NotADirectoryError(f'External data root is not a directory: {root}')
    if require_exists and not root.is_dir():raise FileNotFoundError(f'External data root is unavailable: {root}. No data will be created in the repository. {HINT}')
    if for_write and not Path(root.anchor).is_dir():raise FileNotFoundError(f'External data drive/share is unavailable: {root.anchor}. No repository fallback is allowed. {HINT}')
    return root

def checked_data_destination(path, data_root=None):
    """Validate an absolute helper destination before any mkdir, write or network I/O."""
    root=resolve_data_root(data_root,for_write=True)
    destination=Path(path).resolve()
    if not destination.is_relative_to(root) or destination==root:
        raise ValueError(f'Data path escapes the configured external root: {destination}')
    if destination==REPO_ROOT or destination.is_relative_to(REPO_ROOT):
        raise ValueError(f'Research data cannot be written inside the repository: {destination}')
    return destination

def data_path(relative, data_root=None, *, must_exist=False, for_write=False):
    root=resolve_data_root(data_root,require_exists=must_exist,for_write=for_write)
    rel=Path(str(relative).replace('\\','/'))
    if rel.is_absolute() or rel.drive or PureWindowsPath(str(relative)).drive or '..' in rel.parts:
        raise ValueError(f'Manifest data path must stay relative to the external root: {relative}')
    path=(root/rel).resolve()
    if not path.is_relative_to(root) or path==root:
        raise ValueError(f'Data path escapes the configured external root: {relative}')
    if path==REPO_ROOT or path.is_relative_to(REPO_ROOT):
        raise ValueError(f'Data path resolves inside the code repository: {relative}')
    if must_exist and not path.exists():raise FileNotFoundError(f'Required external data is missing: {path}. {HINT}')
    return path

def configure_data_root_cli(description=None, *, for_write=False):
    """Shared standalone-script argument handling; never creates a directory."""
    parser=argparse.ArgumentParser(description=description)
    parser.add_argument('--data-root',type=Path,default=None,help=f'Absolute external source/cache root; overrides {ENV_NAME} and external user policy')
    args=parser.parse_args()
    root=resolve_data_root(args.data_root,require_exists=not for_write,for_write=for_write)
    os.environ[ENV_NAME]=str(root)
    return root


def evidence_path(name, repo_root=None, data_root=None, *, must_exist=False, for_write=False):
    """Resolve a logical evidence-relative path; never create or fall back.

    Only the explicit original-bulk names map to external original_results/.
    Compact results remain in repo/evidence/. Legacy receipts retain their keys;
    evidence/original_bulk_manifest.json documents the physical relocation.
    """
    rel=Path(str(name).replace('\\','/'))
    if rel.is_absolute() or rel.drive or PureWindowsPath(str(name)).drive or '..' in rel.parts or str(rel)=='.':
        raise ValueError(f'Evidence name must stay relative to evidence/: {name}')
    repo=Path(repo_root or REPO_ROOT).resolve()
    if rel.as_posix() in ORIGINAL_BULK_NAMES:
        root=resolve_data_root(data_root,require_exists=must_exist,for_write=for_write,repo_root=repo)
        return data_path('original_results/'+rel.as_posix(),root,must_exist=must_exist,for_write=for_write)
    base=repo/'evidence';path=(base/rel).resolve()
    if not path.is_relative_to(base) or path==base:
        raise ValueError(f'Evidence path escapes the repository evidence directory: {name}')
    if must_exist and not path.exists():raise FileNotFoundError(f'Required compact evidence is missing: {path}')
    return path


def logical_evidence_key(path, repo_root=None, data_root=None):
    """Return a legacy logical evidence/ key, not a claim of physical location."""
    repo=Path(repo_root or REPO_ROOT).resolve();path=Path(path).resolve()
    if path.is_relative_to(repo/'evidence'):
        return path.relative_to(repo).as_posix()
    if path.name in ORIGINAL_BULK_NAMES and path==evidence_path(path.name,repo_root=repo,data_root=data_root):
        return 'evidence/'+path.name
    raise ValueError(f'Path is not registered evidence: {path}')


def prepare_original_results(data_root=None):
    """Explicit writer setup after validating external storage; no import side effect."""
    directory=data_path('original_results',data_root,for_write=True)
    directory.mkdir(parents=True,exist_ok=True)
    return directory


def verify_original_bulk_manifest(repo_root=None, data_root=None, *, allow_missing=False, allow_reconstructed_gzip=False):
    """Verify retained original bytes; absent results may be explicitly rebuilt."""
    import gzip, hashlib
    repo=Path(repo_root or REPO_ROOT).resolve()
    manifest=repo/'evidence/original_bulk_manifest.json'
    if not manifest.is_file():raise FileNotFoundError(f'Original bulk relocation manifest is missing: {manifest}')
    payload=json.loads(manifest.read_text(encoding='utf-8-sig'));rows=[]
    expected={item['logical_path'].removeprefix('evidence/') for item in payload['files']}
    if expected!=ORIGINAL_BULK_NAMES:raise ValueError('Original bulk manifest does not match the registered file set')
    for item in payload['files']:
        name=item['logical_path'].removeprefix('evidence/')
        if item['data_relative_path']!='original_results/'+name:raise ValueError('Unexpected original bulk storage mapping')
        path=evidence_path(name,repo_root=repo,data_root=data_root)
        if not path.is_file():
            if not allow_missing:raise FileNotFoundError(f'Required original bulk result is missing: {path}')
            rows.append({'logical_path':item['logical_path'],'status':'absent_will_rebuild'});continue
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        exact=path.stat().st_size==item['bytes'] and digest==item['sha256']
        if not exact:
            decoded_ok=(allow_reconstructed_gzip and path.suffix=='.gz' and item.get('decoded_sha256')
                        and hashlib.sha256(gzip.decompress(path.read_bytes())).hexdigest()==item['decoded_sha256'])
            if not decoded_ok:raise RuntimeError(f'Original bulk result checksum mismatch: {path}')
        rows.append({'logical_path':item['logical_path'],'data_relative_path':item['data_relative_path'],
                     'status':'verified_original_bytes' if exact else 'verified_decoded_gzip_wrapper_differs',
                     'sha256':digest,'historical_sha256':item['sha256'],'original_bytes_identical':exact})
    return rows
