# Worked signature diagnostics

These reports demonstrate the reusable [diagnostic CLI](../docs/DIAGNOSTIC.md) on retained research inputs. Per-cell files and count inputs remain outside this repository. The JSON and HTML reports are compact views; exact comparison checks are in `evidence/revision_02/diagnostic_worked_examples.json`.

| Report | Population and comparison | Threshold context |
|---|---|---|
| [Primary donor D4](primary_D4_source_thresholds.html) | Existing broad-alpha proxy: 1,957 vehicle and 4,938 harmine/native-GLP-1 cells | Frozen source-reference thresholds; signed-rank threshold 0 |
| [Baron human3 native](baron_human3_native_target_thresholds.html) | All 1,130 alpha and 787 beta author-labelled cells | Previously frozen native-alpha q95; signed-rank threshold −0.0204642 |
| [Baron paired 2,000 counts: source thresholds](baron_human3_paired2000_source_thresholds.html) | Same 1,081 alpha and 758 beta cells before and after simulated count thinning | Source-reference thresholds unchanged; signed-rank threshold 0 |
| [Baron paired 2,000 counts: target-native thresholds](baron_human3_paired2000_target_native_thresholds.html) | The identical paired count draw and cells | All-native-alpha thresholds unchanged; no refitting to eligible/thinned cells |

The primary D4 mean signed score increases while both up- and down-block detection fall. In the Baron worked draw, the source threshold keeps alpha positive fractions below 1% while beta detection falls from 88.9% to 83.2%. The target-native threshold instead yields alpha positive fractions of 4.44% before and 12.86% after thinning on the same eligible cells. These are different operating points, not contradictory measurements. Neither establishes cell conversion.

The Baron example is **one already-retained draw**, seed 7800 with `SeedSequence([7800, 3, 2000])`, reproduced exactly. It is not an additional replicate or the five-draw summary used in the paper. All four method-score arrays agree exactly with the corresponding retained arrays. Pair IDs identify the same cells, and the simulated observation condition is explicitly distinguished from biological treatment.

Source state labels are not ancestry truth. The primary proxy is not the publication's unavailable curated candidate-cell population. Already harmonized frozen-universe count inputs retain upstream zero-filling and original whole-library totals; the report provenance records unavailable source features. The tool does not fit thresholds, correct scores, judge mechanism or issue pass/fail conclusions.
