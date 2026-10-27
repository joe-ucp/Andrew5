"""Record or verify the public scientific files without copying them."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'PUBLIC_FILES.json'
ROOT_NAMES = {
    'README.md', 'CHANGELOG.md', 'VERSION', 'REUSE_TERMS.md',
    '.gitignore', '.gitattributes', 'storage.example.json', 'environment.json',
    'source_manifest.tsv', 'sample_manifest.tsv', 'cache_manifest.tsv',
    'primary_sample_manifest.tsv', 'primary_source_manifest.tsv',
    'literature_source_manifest.tsv', 'revision_source_manifest.tsv',
    'reference_donor_mapping.tsv', 'same_system_compatibility_ledger.tsv',
}
PUBLIC_TREES = {'code', 'scripts_primary', 'tests', 'docs', 'config',
                'evidence', 'manuscript', 'diagnostics', '.github'}
PUBLIC_SUFFIXES = {'.py', '.txt', '.md', '.tsv', '.csv', '.json', '.tex',
                   '.pdf', '.docx', '.png', '.html', '.yml', '.yaml'}


def approved_paths():
    result = []
    candidates = [ROOT / name for name in ROOT_NAMES if (ROOT / name).exists()]
    for name in sorted(PUBLIC_TREES):
        tree = ROOT / name
        if tree.exists():
            if tree.is_symlink() or (hasattr(tree, 'is_junction') and tree.is_junction()):
                raise ValueError('Public trees must not be filesystem links: ' + name)
            candidates.extend(tree.rglob('*'))
    for path in candidates:
        relative = path.relative_to(ROOT)
        if relative.parts[0] not in PUBLIC_TREES and relative.as_posix() not in ROOT_NAMES:
            continue
        if any(part in {'__pycache__', '.pytest_cache'} for part in relative.parts):
            continue
        if path.is_symlink() or (hasattr(path, 'is_junction') and path.is_junction()):
            raise ValueError('Public files must not link outside their tree: ' + relative.as_posix())
        if not path.is_file():
            continue
        if path.suffix not in PUBLIC_SUFFIXES and relative.as_posix() not in ROOT_NAMES:
            continue
        if path.stat().st_size > 20 * 1024 * 1024:
            raise ValueError('Unexpectedly large public file: ' + relative.as_posix())
        result.append(relative.as_posix())
    return sorted(result)


def record(relative):
    path = ROOT / relative
    return {'path': relative, 'bytes': path.stat().st_size,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def verify_manifest():
    payload = json.loads(MANIFEST.read_text(encoding='utf-8'))
    paths = approved_paths()
    if [r['path'] for r in payload['files']] != paths:
        raise ValueError('Public file set differs from the published manifest')
    for item in payload['files']:
        if record(item['path']) != item:
            raise ValueError('Public checksum mismatch: ' + item['path'])
    return {'status': 'PASS', 'listed_files': len(payload['files']),
            'listed_bytes': sum(r['bytes'] for r in payload['files']),
            'manifest': MANIFEST.name, 'files_modified': 0}


def check_git_candidates():
    """Reject tracked or unignored files outside the declared public file set."""
    output = subprocess.run(['git', 'ls-files', '--cached', '--others',
                             '--exclude-standard', '-z'], cwd=ROOT,
                            check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    candidates = {p for p in output.stdout.decode('utf-8').split('\0') if p}
    expected = set(approved_paths()) | {MANIFEST.name}
    extra = sorted(candidates - expected)
    missing = sorted(expected - candidates)
    if extra or missing:
        raise ValueError('Git/public file mismatch: ' + json.dumps({'unexpected': extra, 'ignored_or_missing_public': missing}))
    return {'status': 'PASS', 'candidate_files': len(candidates),
            'unexpected_files': 0, 'git_index_modified': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true', help='Record a reviewed public file set')
    mode.add_argument('--check', action='store_true', help='Read-only public file and checksum check')
    parser.add_argument('--git-check', action='store_true', help='Also inspect tracked and unignored Git candidates')
    args = parser.parse_args()
    if args.write:
        payload = {'schema_version': 2,
                   'version': (ROOT / 'VERSION').read_text(encoding='utf-8').strip(),
                   'scope': 'Scientific code, methods, compact results and manuscript files.',
                   'files': [record(p) for p in approved_paths()]}
        MANIFEST.write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')
    report = verify_manifest()
    if args.git_check:
        report['git_candidates'] = check_git_candidates()
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
