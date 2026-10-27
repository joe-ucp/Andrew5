# Supplied-signature diagnostic

`code/signature_diagnostic.py` describes detection, rank baselines, occupancy, rank remainder and supplied-threshold behavior. It does not fit a signature, choose a threshold, produce a corrected score, assign pass/fail, infer ancestry or assess clinical efficacy.

Run commands from the repository folder after installing [the pinned dependencies](../code/requirements.txt). The small synthetic example is:

```text
python code/signature_diagnostic.py --example --data-root /absolute/external/root --output-prefix diagnostics/example --threshold 0
```

Use a platform-native absolute path on Windows. The example creates 96 synthetic cells, inputs, a JSON report, an escaped static HTML report and a per-cell TSV outside the repository. Reusing an output prefix is rejected. Example numbers are generated for interface illustration and have no biological interpretation.

For your data:

```text
python code/signature_diagnostic.py --data-root /absolute/external/root --matrix /absolute/external/root/input/counts.npz --genes /absolute/external/root/input/genes.txt --cells /absolute/external/root/input/cells.txt --metadata /absolute/external/root/input/metadata.tsv --signature signature.json --threshold 0 --baseline-treatment vehicle --comparison-treatment exposure --state-positive beta --state-negative alpha --output-prefix diagnostics/my_report
```

Inputs:

- `counts.npz`: SciPy sparse **genes × cells**, finite nonnegative integer counts. It may contain genes outside the score universe. Normalized expression is not accepted as count input when fractional.
- `genes.txt` and `cells.txt`: one unique identifier per line, no header, ordered to match matrix rows and columns. Duplicate IDs are rejected.
- `metadata.tsv` (optional): unique `cell_id` matching the cell list exactly. Its rows are reordered by that ID, and the report records this. Optional `group`, `donor`, `treatment` and `state` columns define descriptive strata. Optional integer `total_counts` preserves whole-library mass when the supplied matrix is already a feature subset; it must be at least each column sum. Without it, matrix column sums define totals.
- `signature.json`: object with explicit `up`, `down` and `universe` gene lists. Each list must be unique; blocks must be nonempty, disjoint and contained in the universe. Missing input genes are reported and zero-filled while retaining the supplied universe and block denominators. This is an explicit target definition, not evidence of cross-platform equivalence.
- `--threshold`: supplied strict **signed-rank** threshold. Alternatively `--thresholds FILE.json` supplies a method-to-number object using `signed_rank`, `signed_log1p_CPT10K`, `ucell_form_1500`, or `unclipped_ucell_difference`. Thresholds are method-specific and are never estimated here. Optional `--provenance FILE.json` carries caller-provided source/calibration context into the escaped report.

All raw inputs must be under the configured external data root. Small signature, threshold and provenance JSON files may be elsewhere. Output paths are relative to the external root; repository overlap, path escapes and overwrites are rejected. The storage configuration is described in [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

The JSON/HTML reports group means, cancelling zero-rank baselines and centered block terms, up/down detection, exact occupancy and remainder, and four declared scoring methods. UCell-form methods are non-evaluable below the fixed 1,500-gene universe requirement; they are explanatory implementations with a frozen universe, not claims of package equivalence. Zero-total cells are excluded with counts shown. The default minimum group size of 20 governs displayed rates and directional flags; raw threshold counts and descriptive means remain available. `--minimum-group-size` explicitly changes that reporting rule.

Treatment flags require the two supplied labels and describe score increases accompanying up/down detection losses within matching other metadata strata. They are arithmetic descriptions, not significance tests. State AUROCs require supplied positive/negative labels, compare those labels within other strata, and are suppressed below the group minimum. State discrimination never establishes conversion, function or immune safety.

For an explicitly paired comparison, optional `pair_id` metadata identifies the same cell or observation across the two supplied conditions. IDs must be unique within each condition and stratum. Strict-threshold transition counts are shown only when the sets match exactly; otherwise unmatched counts are reported without silently selecting an intersection. Label simulated thinning conditions explicitly so that these comparisons cannot be mistaken for biological treatment effects.

Source hashes, supplied provenance, code hashes and eligibility details accompany each report. This tool neither corrects background contamination nor diagnoses doublets. These and the source lineage of the labels require separate measurements and review.

The [worked examples](../diagnostics/index.html) demonstrate primary donor D4 and one retained paired Baron count draw. [RESULTS.md](RESULTS.md) distinguishes these examples from the five-draw scientific summaries.
