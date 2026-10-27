# Andrew5

Andrew5 connects an exact mathematical explanation of a gene-expression score to measured consequences in human islet data and a diagnostic you can run and inspect. The scientific application is a fixed 50-up/50-down alpha-cell reprogramming signature transferred across experimental settings.

The accompanying paper is **Detection imbalance shapes a transferred alpha cell reprogramming signature**.

[Read the paper](manuscript/paper.pdf) · [Word document](manuscript/paper.docx) · [Results guide](docs/RESULTS.md) · [Reproducibility](docs/REPRODUCIBILITY.md)

## What Andrew5 provides

**An exact mathematical explanation.** For the stated signed average-rank score in a fixed gene universe, including tied zeros, the score separates exactly into detection imbalance O and a centered positive-rank remainder R:

$$
s=O+R,\qquad O=\frac{1}{2}\left(\frac{k_U}{m}-\frac{k_D}{n}\right)=\frac{k_U-k_D}{100}\quad\text{when }m=n=50.
$$

Here kU and kD count detected genes in the up/down blocks, whose sizes are m and n. The decomposition accounts for the common zero-rank baseline and shows why a positive score change can include loss of down-program detection. The remainder is retained; sparsity alone does not justify dropping it. [Methods and definitions](docs/METHODS.md) give the score definition, component interpretation and bounds. This is an algebraic explanation; biological causes require separate evidence.

**Measured consequences in human islet data.** The study evaluates four drug-treated donors, reference preparations and independent labelled alpha/beta cells. Observed treatment contrasts and simulated count-thinning comparisons are reported separately. Strong state ranking, a higher mean score and a stable fixed-threshold decision are tested as distinct properties, with donor results, eligible-cell counts and unsuccessful controls retained.

**A runnable, inspectable diagnostic.** Supply a sparse count matrix, an explicit gene universe and up/down program, and optional cell metadata and thresholds. The diagnostic exposes detection, score components, eligibility and strict-threshold behavior in HTML, JSON and per-cell output. It describes the supplied score without fitting a replacement or assigning biological conversion. Start with the [synthetic example below](#quickstart), inspect the [worked reports](diagnostics/index.html), or use [your own inputs](docs/DIAGNOSTIC.md).

## Main findings

| Question | Finding |
|---|---|
| Does native reference ranking ensure useful detection? | Reference AUROC was **0.859**, but only **48.5%** of PM-preparation cells exceeded the strict alpha-control threshold of zero; the achieved alpha false-positive fraction was **1.90%**. |
| What contributes to positive treatment contrasts? | Reduced down-program detection contributes positively in **all four donors**. Up-program detection also falls in two. The score has an exact occupancy-plus-rank-remainder decomposition. |
| Does the score discriminate independent cell states? | Native Baron alpha-versus-beta AUROCs were **0.9947–0.9998**. These are mature-state labels, not conversion truth. |
| What happens to fixed decisions after thinning to 2,000 counts? | The source threshold keeps alpha-positive fractions below **0.8%**, while beta-state detection falls **6.2–11.7 percentage points** on the same eligible cells. Holding negative target-native thresholds fixed instead raises alpha-positive fractions to **12.9–20.9%**. |
| Is broader signature specificity established? | **No.** The matched-program control failed its declared candidate-support requirement. No specificity percentile was produced. |

The [results guide](docs/RESULTS.md) links these statements to exact tables, figures and code.

## Quickstart

From the repository folder, with Python 3.12:

```sh
python -m pip install -r code/requirements.txt
python code/reproduce.py check
python code/run_tests.py
```

The file check verifies all code, documentation and artifacts listed in `PUBLIC_FILES.json`, plus registered scientific inputs, without downloading research data. The tests use synthetic fixtures; platform-specific skips are reported explicitly. See [setup instructions](docs/GETTING_STARTED.md) to keep the environment, caches and temporary files outside the repository.

Run a small synthetic diagnostic into a new external folder:

```sh
python code/reproduce.py demo --output-dir /absolute/external/folder/demo
```

Use a platform-native absolute path on Windows. The example produces HTML/JSON descriptions and per-cell values; its synthetic numbers have no biological interpretation.

## Explore and reproduce

| Material | Location |
|---|---|
| Paper and authoritative narrative source | [manuscript/](manuscript/) |
| Publication figures | [manuscript/figures/](manuscript/figures/) |
| Exact numerical results | [evidence/](evidence/), indexed in [docs/RESULTS.md](docs/RESULTS.md) |
| Scientific definitions and source qualifications | [docs/METHODS.md](docs/METHODS.md), [docs/PROVENANCE.md](docs/PROVENANCE.md) |
| Scientific analysis code and loaders | [code/](code/), [scripts_primary/](scripts_primary/) |
| Your own supplied-signature diagnostic | [docs/DIAGNOSTIC.md](docs/DIAGNOSTIC.md) |
| Worked diagnostic reports | [diagnostics/](diagnostics/) |
| Source URLs, checksums and sample mappings | [source_manifest.tsv](source_manifest.tsv), [sample_manifest.tsv](sample_manifest.tsv), [primary_sample_manifest.tsv](primary_sample_manifest.tsv) |

`python code/reproduce.py stages` lists supported computations and their prerequisites. Reproducing a stage from retained caches is different from rebuilding its inputs from raw data; [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) explains the distinction.

Bulk inputs, sparse caches and generated runs live outside the repository. Source accessions include GSE117454 and GSE123844 for the reference, GSE266427 for the drug-exposed islets, GSE84133 for independent state controls, and GSE217837 for the observational SMOC1 comparison. Retrieval manifests retain the exact URLs and hashes.

## Scope

The algebra applies to the stated tied-rank score; empirical findings are limited to the programs, datasets and populations examined. Four donors provide four biological contrasts. Thinning draws do not add biological replication. The motivating publication's exact curated barcode labels and corrected expression layer were not recovered, so the present population is not an exact replication of its candidate-cell analysis.

Expression similarity and mature-state discrimination do not establish alpha ancestry, regulated secretion, immune protection, durability or clinical efficacy.

## Attribution and reuse

**UCP Technology LLC.** Version: [VERSION](VERSION). No general reuse grant is supplied; obtain permission as described in [REUSE_TERMS.md](REUSE_TERMS.md). Third-party software and data retain their own terms.
