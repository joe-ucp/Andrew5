"""Published scientific configuration and external output locations."""
from pathlib import Path
import json
from revision_benchmark import sha
from data_paths import data_path
from public_runtime import protocol, frozen_score

ROOT = Path(__file__).resolve().parents[1]
EV = ROOT / 'evidence/revision_02'
OLD = ROOT / 'evidence/revision_01'


def load_protocol():
    return protocol('revision_02'), frozen_score()


def external(relative):
    path = data_path('revision_02/' + relative, for_write=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def write_json(path, value):
    path = Path(path)
    if path.exists():
        raise FileExistsError('No implicit evidence overwrite: ' + str(path))
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
