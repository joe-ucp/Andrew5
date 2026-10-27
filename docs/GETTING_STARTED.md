# Getting started

Run commands from the repository folder. Use Python 3.12 and the pinned [requirements](../code/requirements.txt). Reading the supplied [PDF](../manuscript/paper.pdf) or [Word document](../manuscript/paper.docx) needs no Python environment.

## Set up an external environment

Keep your environment, package cache, temporary files and generated outputs outside the repository. Replace the example paths with absolute paths on your computer.

PowerShell:

```powershell
$ResearchRoot = 'D:\research-data\signature-study'
$ResearchEnvironment = Join-Path $ResearchRoot 'environment'
$ResearchTemporary = Join-Path $ResearchRoot 'tmp'
New-Item -ItemType Directory -Force -Path $ResearchRoot,$ResearchTemporary | Out-Null
$env:TEMP = $ResearchTemporary
$env:TMP = $ResearchTemporary
$env:STRIKE007_TEST_TMPDIR = $ResearchTemporary
$env:MPLCONFIGDIR = Join-Path $ResearchRoot 'matplotlib-cache'
py -3.12 -m venv $ResearchEnvironment
$ResearchPython = Join-Path $ResearchEnvironment 'Scripts\python.exe'
& $ResearchPython -m pip install --cache-dir (Join-Path $ResearchRoot 'pip-cache') -r code/requirements.txt
```

Linux or another POSIX shell:

```sh
RESEARCH_ROOT=/absolute/external/research-data
RESEARCH_ENVIRONMENT="$RESEARCH_ROOT/environment"
mkdir -p "$RESEARCH_ROOT/tmp" "$RESEARCH_ROOT/pip-cache"
export TMPDIR="$RESEARCH_ROOT/tmp"
export STRIKE007_TEST_TMPDIR="$RESEARCH_ROOT/tmp"
export MPLCONFIGDIR="$RESEARCH_ROOT/matplotlib-cache"
python3.12 -m venv "$RESEARCH_ENVIRONMENT"
RESEARCH_PYTHON="$RESEARCH_ENVIRONMENT/bin/python"
"$RESEARCH_PYTHON" -m pip install --cache-dir "$RESEARCH_ROOT/pip-cache" -r code/requirements.txt
```

Use the environment's Python executable in subsequent commands. The shortened `python` commands below assume it is active.

## Check the files and run synthetic tests

```sh
python code/reproduce.py check
python code/run_tests.py
```

The file check verifies all code, documentation and artifacts listed in `PUBLIC_FILES.json`, plus registered scientific inputs. It requires no external research dataset and writes no scientific outputs. The default test suite uses small synthetic fixtures and reports platform-specific skips. No dataset is loaded during test collection.

## Run a synthetic example

```sh
python code/reproduce.py demo --output-dir /absolute/external/research-data/demo
```

Choose a new output folder. The example creates 96 synthetic cells and writes an HTML report, a JSON report and per-cell values outside the repository. Inspect the group detection means, occupancy and rank remainder, cell eligibility and state comparisons. These generated values illustrate the interface; they are not research results.

For your own count matrix, see [DIAGNOSTIC.md](DIAGNOSTIC.md). For retained scientific results, see [RESULTS.md](RESULTS.md).

## Reproduce a scientific computation

```sh
python code/reproduce.py stages
python code/reproduce.py influence --output-dir /absolute/external/research-data/influence
```

The influence calculation uses compact gene contributions already included in the repository. It is not a reconstruction from raw counts. Other stages require the external raw files or caches declared in the stage listing. [REPRODUCIBILITY.md](REPRODUCIBILITY.md) describes these prerequisites and output checks.

Missing requested inputs cause an error. Storage is never silently redirected into the repository. An optional local storage policy can restrict allowed external locations; see the configuration section of the reproducibility guide.

Document authoring has separate dependencies in [requirements-artifacts.txt](../code/requirements-artifacts.txt) and requires document-rendering tools. They are unnecessary for the data-free checks and synthetic diagnostic.
