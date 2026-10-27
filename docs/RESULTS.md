# Results guide

The [paper](../manuscript/paper.pdf) gives the scientific interpretation. This index connects its main results to full-precision tables and implementations. File names containing `revision_01`, `revision_02` or `revision_03` identify retained analysis output groups; use the substantive descriptions below to navigate.

## Main results

| Result | Numerical evidence | Figure and implementation |
|---|---|---|
| Native reference AUROC 0.859; PM-preparation detection 48.49% at strict threshold zero, alpha-positive fraction 1.90% | [Source operating points](../evidence/revision_02/source_operating_points.tsv), [calibration operating points](../evidence/calibration_operating_points.tsv) | [Figure 1](../manuscript/figures/rev1_figure1_reference.png); [reference analysis](../code/reference_analysis.py) |
| All four native donor contrasts are positive, with reduced down-program detection contributing positively in every donor | [Observed component contrasts](../evidence/revision_01/observed_component_contrasts.tsv), [zero-baseline contrasts](../evidence/revision_02/zero_baseline_contrasts.tsv) | [Mechanism figure](../manuscript/figures/rev2_mechanism.png); [score components](../code/revision_score.py) |
| True thinning and depth eligibility both affect the contrast | [Thinning donor contrasts](../evidence/revision_01/thinning_donor_contrasts.tsv), [fixed-cohort summary](../evidence/revision_02/fixed_cohort_summary.tsv) | [Population/depth figure](../manuscript/figures/rev1_figure2_depth_population.png), [baseline/floor supplement](../manuscript/figures/rev2_supplement_baseline_floor.png) |
| Native independent alpha/beta state discrimination remains strong: AUROC 0.9947–0.9998 | [Independent state benchmark](../evidence/revision_01/independent_state_benchmark.tsv), [independent thinning](../evidence/revision_01/independent_thinning.tsv) | [Independent state figure](../manuscript/figures/rev1_figure5_independent.png); [independent analysis](../code/revision_independent.py) |
| At 2,000 counts the fixed source threshold retains low alpha-positive fractions while beta detection falls; fixed negative native-target thresholds instead increase alpha-positive fractions | [Threshold summary](../evidence/revision_02/threshold_summary.tsv), [target-native thresholds](../evidence/revision_02/target_native_thresholds.tsv) | [Threshold transport figure](../manuscript/figures/rev2_threshold_transport.png) |
| Gene contributions include both supporting and opposing terms | [Gene/donor contrasts](../evidence/revision_02/gene_donor_contrasts.tsv), [contribution summary](../evidence/revision_02/gene_contribution_summary.tsv), [concentration](../evidence/revision_02/gene_concentration.tsv) | [Gene contributions figure](../manuscript/figures/rev2_gene_drivers.png); [gene arithmetic](../code/revision02_gene_math.py) |
| Matched-program specificity is untested because candidate support failed | [Candidate support](../evidence/revision_02/matched_candidate_support.tsv), [comparison table](../evidence/revision_02/matched_program_comparison.tsv) | [Matching analysis](../code/revision02_gene_specificity.py) |
| Exploratory top-five contribution subtraction reverses D1 and slightly D2, but no single subtraction reverses a donor | [All 400 single subtractions](../evidence/revision_03/influence_single.tsv), [all cumulative subtractions](../evidence/revision_03/influence_topk.tsv), [summary](../evidence/revision_03/influence_summary.tsv) | [Influence arithmetic](../code/revision03_influence.py) |

The signed score is exactly occupancy plus a positive-rank remainder for the declared fixed universe and tied-average-rank definition. Sparse detection alone does not justify discarding the remainder. A rise in the raw up-block contribution can reflect a larger common zero-rank baseline even while detected up genes decline; the complete signed score cancels that common baseline.

## How to read the comparisons

- **Reference preparation detection** classifies preparation membership. It does not label individually verified conversion events.
- **Native drug contrasts** compare treated and vehicle mean cell scores within each donor. Four donors are four biological units.
- **Thinning** removes sampled molecules without replacement from observed cells. Each comparison uses its stated eligible cells; the fixed high-depth cohort is a separate sensitivity analysis.
- **Baron labels** provide independent mature alpha/beta state controls. They are not permanent lineage or functional labels.
- **Threshold rules** are distinct: the source-zero rule, fixed negative target-native rule and within-depth recalibration have different interpretations. Recalibration on the evaluated cells is resubstitution.
- **Gene subtractions** are exploratory arithmetic on retained additive contributions, not a recomputed smaller signature, causal intervention or significance test.

## Conditional checks and limitations

Finite-library occupancy expectations are in [detection_expectation_samples.tsv](../evidence/revision_01/detection_expectation_samples.tsv). They predict occupancy conditional on observed molecules, not the full rank remainder or biological response.

The [nuisance figure](../manuscript/figures/rev1_supplement_nuisance.png) retains favorable D3/D4 exceedances and finite-pool/coupling sensitivity. Imposed RNA mixtures, doublet fractions and cycle compositions are conditional simulations, not measured contamination bounds or convergence guarantees.

The matched comparison produced no programs or specificity percentiles after five down-gene slots failed support. No relaxed replacement analysis is presented. Broader empirical specificity therefore remains unresolved.

The motivating publication's curated candidate-cell labels and identified corrected layer were not recovered. Expression proxies and deposited counts do not reconstruct those unavailable annotations. Observational SMOC1 expression does not provide beta-origin dedifferentiation truth. Ancestry, regulated secretion, immune control and durability require separate direct evidence.

For interactive inspection of supplied scores and operating rules, open the [worked diagnostic reports](../diagnostics/index.html) or use the [diagnostic CLI](DIAGNOSTIC.md).
