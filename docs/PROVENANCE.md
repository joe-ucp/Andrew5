# Data and source provenance

Public source files are identified by URL, access date, size and SHA-256 in the [source manifest](../source_manifest.tsv) and the [primary](../primary_source_manifest.tsv), [literature](../literature_source_manifest.tsv) and [additional-source](../revision_source_manifest.tsv) manifests. The [execution input registry](../config/scientific_inputs.json) distinguishes compact scientific inputs, downloadable sources and derived files used by the reproduction entry point. The documentary manifests distinguish numerical inputs, metadata, snapshots and unsuccessful retrievals. Hashes authenticate acquired bytes; they do not resolve biological metadata discrepancies.

## Dataset roles

| Accession | Role | Important qualification |
|---|---|---|
| [GSE117454](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117454) | Bulk alpha-GFP and alpha-PDX1/MAFA reference for program selection | Four donor pairs are verified from the supplement; matching suffixes alone produces a false pair. |
| [GSE123844](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123844) | Preparation-level reference calibration | 158 alpha-control, 241 beta-control, 860 PM-preparation and 77 graft-derived cells. Donor labels are not mapped to barcodes. |
| [GSE266427](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266427) | Four drug-exposed human islet donors, four exposures each, 96 hours | Deposited integer matrices and numerical QC do not reconstruct the published curated candidate-cell population. |
| [GSE84133](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84133) | Independent labelled human islet state/count control | Four human files contain 8,569 labelled cells, including 2,326 alpha and 2,525 beta cells; labels are mature states, not conversion ancestry. |
| [GSE217837](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217837) | Observational SMOC1-expression comparison | Released matrix is noninteger expression; no SMOC1 perturbation or adverse beta-lineage label is established by these samples. |

## Reference pairing and independence

Furuyama supplementary Tables 8–9 identify the four numbered control/PM donor pairs as IIDP1697, IIDP1738, IIDP1770 and IIDP1776. Numbered pair 3 is invalid: control and PM samples have different donors. Letter-coded libraries without corroborated mappings are excluded from paired feature selection. [reference_donor_mapping.tsv](../reference_donor_mapping.tsv) retains the explicit pairing.

The supplementary single-cell row lists donors HP1282 and IIDP1908, distinct from the listed bulk donors, but does not assign each released barcode to a donor. Preparation-level discrimination cannot be promoted into mapped donor replication. The available graft single-cell preparation was author-flagged for limited quality/quantity and is descriptive only.

## Primary sample identity and expression layer

The [primary sample manifest](../primary_sample_manifest.tsv) preserves explicit GEO donor and treatment fields. Donor/file aliases are D1/R407, D2/HP21226, D3/HP21227 and D4/SAMN1433. Barcodes are library-specific; bare barcodes across samples are not merged.

The 16 matrices contain 253,462 submitted barcodes and 36,601 feature rows each, with 36,591 unique symbols. Submitted values are nonnegative integers. Numerical filtering retains 177,149 barcodes, compared with 109,881 cells in the publication's final curated population. The exact curated barcode/state labels and identified corrected expression layer were not recovered from the inspected public records and supplements.

Deposited CellRanger header versions vary across donors, while the article describes later processing including CellBender. Integer values and inherited headers alone do not prove correction was absent. The analysis identifies its input as the deposited integer-count layer rather than assuming equivalence to the final published object.

The primary comparison uses harmine 10 μM plus native GLP-1 5 nM versus vehicle. The article Results names GLP-1(7–36), whereas its Key Resources table names fragment 7–37; both descriptions differ from exendin-4. LY364947 is 3 μM in Results and 10 μM in GEO. A donor alias also differs between one supplement entry and GEO. Explicit GEO sample blocking is retained without guessing a typo correction. These unresolved details remain source constraints.

## Independent state and observational controls

Baron's four deposited human files total 8,569 cells, whereas the paper reports 8,629 filtered cells. The cohort includes a 17-year-old donor and a donor with type 2 diabetes, not four healthy adults. Original labels and cell identifiers are retained. Missing source features are explicitly zero-filled under the original-universe definition; the shared-feature sensitivity is reported separately.

The GSE217837 scRNA matrix has 26,671 genes and 10,732 barcodes, with extensive noninteger values despite its count-matrix filename. Sample records describe normal non-diabetic tissue, without donor-to-barcode mapping or a declared SMOC1 perturbation. Barcode suffixes are not assumed to identify donors. The matrix is used for descriptive expression association, not UMI-depth resampling or a validated dedifferentiation challenge.

## Interpreting source evidence

The motivating studies use different interventions, preparations and endpoints. Direct functional or immune findings from another study are not inherited by the cells analyzed here. The paper distinguishes inspected full methods/data from abstract-level contextual sources, including the 2026 preprint, and notes reuse of GSE266427 where it is not an independent transcriptomic replication.

The study uses four primary biological donor contrasts. Large cell counts, source-state labels and repeated numerical draws do not add donor replication or establish absolute tissue abundance. Source-layer equivalence, permanent lineage, regulated secretion and immune durability remain separate evidentiary requirements.
