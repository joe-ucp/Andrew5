# Scientific methods

This page summarizes the fixed definitions used by the [paper](../manuscript/paper.pdf). Machine-readable [scientific parameters](../config/scientific_parameters.json) retain the seeds, matching rules, count depths and decision definitions. Full-precision outputs are indexed in [RESULTS.md](RESULTS.md); dataset qualification and unresolved source discrepancies are described in [PROVENANCE.md](PROVENANCE.md).

## Study units and primary comparison

The primary estimand is the treated-minus-vehicle difference in the mean cell score within each of four human islet donors at 96 hours. The treatment is harmine plus native GLP-1. Each donor is reported separately. Captured cells and repeated numerical draws are not additional biological replicates, and unequal captured-cell yields do not measure absolute tissue abundance.

GSE117454 supplies the bulk reference used to select genes. GSE123844 supplies preparation-level single-cell calibration. GSE266427 supplies drug-exposed islets. GSE84133 supplies independent author-labelled mature alpha/beta states. GSE217837 supplies an observational expression comparison, not a labelled perturbation or lineage challenge.

The original program was fixed before primary treatment scores were calculated. Later analyses were specified with earlier findings already known; they are not public preregistration or an independent biological holdout. The qualitative influence-subtraction outcome was already known before the local arithmetic calculation. Its retained results are therefore explicitly exploratory, not prospective confirmation.

## Gene selection and score

Four bulk donor pairs are verified from the source supplement: IIDP1697, IIDP1738, IIDP1770 and IIDP1776. Matching numbered sample suffixes alone is insufficient. Duplicate exact gene symbols are summed. The fixed universe has 12,051 genes shared across the reference and primary assay and expressed at CPM ≥1 in at least four training libraries, after the declared exclusions.

Selection excludes hormones, lineage/transgene readouts, annotation markers, mitochondrial/ribosomal features and the explicit cell-cycle list in [reference_analysis.py](../code/reference_analysis.py). Within each verified donor pair, log2(CPM+1) differences compare alpha-PDX1/MAFA with alpha-GFP. Genes must have the same strict sign in all four pairs. The largest positive and negative median changes select 50 genes per block, with gene-name tie breaking. The [frozen program](../evidence/frozen_score.json) is retained without refitting.

For average ranks over the fixed universe, including tied zeros,

$$s(x)=\frac{1}{m}\sum_{g\in U}\frac{r_g-1/2}{G}-\frac{1}{n}\sum_{g\in D}\frac{r_g-1/2}{G},\qquad m=n=50,\ G=12051.$$

The score is dimensionless. Exact half-rank integer arithmetic preserves mathematical zero before floating conversion.

## Detection and rank decomposition

Let K be the number of detected genes in the universe, and kU and kD the detected genes in the two blocks. The score is exactly $s=O+R$, where

$$O=\frac{1}{2}\left(\frac{k_U}{m}-\frac{k_D}{n}\right)=\frac{k_U-k_D}{100}.$$

R contains the centered ranks among detected genes. The common zero-rank baseline is $b=(G-K)/(2G)$. The separately signed block contributions are $A_U=b+k_U/(2m)+R_U$ and $A_D=-b-k_D/(2n)+R_D$, with $R_U+R_D=R$. This explains why a raw up-block contribution can increase even as up-gene detection falls. The complete signed score cancels b.

For K≥1, $|R|\leq (K-1)/(2G)\,(k_U/m+k_D/n)$; at K=0, R=0. Approximating s by O requires a small remainder in the setting examined. Sparsity alone does not establish that condition. Components are arithmetic allocations, not causal or lineage assignments.

## Qualification and cell populations

Primary numerical QC requires at least 1,000 whole-library counts, at least 400 detected original feature rows, and mitochondrial fraction below 0.20. Duplicate symbols are aggregated after QC. Across all 16 deposited libraries, 177,149 of 253,462 barcodes pass this numerical filter. These are not the publication's 109,881 curated cells.

Broad alpha cells require GCG≥3 counts and strict dominance over INS, SST, PPY and GHRL. Strict alpha additionally requires GCG≥10 and detected ARX or IRX2. Cycle-high alpha uses the donor vehicle 75th percentile of positive cell-cycle fractions, held fixed across treatment arms. These are current-expression proxies, not published cycling-alpha annotations or permanent ancestry. Other endocrine populations provide separate sensitivity comparisons.

## Calibration, sampling and fixed decisions

Native reference comparison distinguishes AUROC from detection at a threshold. The alpha-control 95th percentile uses linear interpolation and a strictly greater-than comparison, so tied zero scores affect the achieved false-positive fraction. PM-preparation membership is not individually verified conversion truth.

Original depth simulations use multinomial draws from empirical profiles and may exceed observed counts. They are observation-model simulations, not recovered molecules. The declared 80% preparation-detection use criterion is analyst-chosen, not a biological requirement or clinical standard; the reference result does not meet it. Treatment scores therefore remain descriptive measurements.

True thinning samples d molecules without replacement from cells with observed whole-library counts N≥d, preserving outside-universe mass. Four depths (1,000, 2,000, 5,000 and 10,000) and five fixed draws per depth are retained. Every depth-specific contrast is compared with native scores on the same eligible cells. A separate N≥10,000 cohort keeps membership fixed across depths. Eligibility changes are not silently treated as a depth correction.

For a gene with observed count c, the exact conditional detection probability is $1-\binom{N-c}{d}/\binom{N}{d}$. Averaging these probabilities within blocks predicts occupancy by linearity; it does not predict the rank remainder. Mean simulation-minus-expectation differences assess conditional numerical agreement, not a biological residual.

The source-reference threshold is held fixed during independent state comparisons; [source operating points](../evidence/revision_02/source_operating_points.tsv) retain the thresholds and achieved rates. A separate rule holds each [target donor's native alpha threshold](../evidence/revision_02/target_native_thresholds.tsv) fixed. Within-depth recalibration is explicitly resubstitution on the evaluated cells. Counts and ties are retained; displayed state rates require at least 20 cells, and AUROC requires at least 20 of each labelled state.

## Feature, scoring and gene-level checks

Four leave-one-training-donor-out programs use the same selection rule on three reference pairs, retaining the preselected universe. Together with the original they give five programs. Fractional ranks, signed mean log1p(counts per 10,000), and a fixed UCell-form score are all reported. The UCell-form calculation uses a rank cap of 1,500 and clips the signed combination at zero; an unclipped diagnostic is reported separately. These are specified formulas, not a claim of equivalence to every package version and option.

Gene-level signed occupancy contributions are +Δdetection/100 for up genes and −Δdetection/100 for down genes. Positive-rank contributions complete the score. Conditional-positive expression compares potentially different positive-cell subsets and does not identify a within-cell biological mechanism.

The matched-program design proposed 500 programs with 50 genes per block, excluding the original genes. Candidates had to lie within 0.025 in equal-donor mean detection and within 0.5 in mean log2(all-cell CP10K+0.01); at most 200 nearest candidates were retained per slot. Every slot needed at least ten candidates before sampling. Planned donor/block balance tolerances were 0.025 for mean detection and 0.25 for mean log2 abundance, with unique genes within each program. Five down-gene slots failed initial support, so no random programs, later balance checks or specificity percentiles were produced. No widened calipers or partial successful subset replaces this result.

The exploratory influence calculation sums each gene's signed occupancy and centered positive-rank contributions, then subtracts single terms or each donor's largest positive terms. All 400 single-term and 20 cumulative results are retained. The original universe and block denominators stay fixed, with no reranking or refitting. Selection is outcome-informed; this is not gene knockout, mediation, causal necessity or statistical robustness.

## Nuisance models and biological scope

RNA-mixture, doublet and cycle-composition scenarios impose declared fractions and empirical sampling models. They do not estimate contamination prevalence. Their envelopes and finite-pool/coupling sensitivity are conditional numerical results; favorable D3/D4 exceedances remain visible. They neither exhaust biological alternatives nor provide universal confidence bounds.

The observational SMOC1 comparison uses a noninteger expression layer and is never thinned as counts. Mature-state classification and expression association do not establish ancestry or regulated function. The paper's separate source-to-function condition is

$$Q_{\mathrm{obs},L}>Q_{\mathrm{seed},U}+Q_{\mathrm{other},U},\qquad Q_{\mathrm{seed},U}=N_{\beta,0,U}F_{\mathrm{exp},U}q_{\beta,U}.$$

It requires independently measured, jointly valid bounds on initial beta seeds, survival/expansion, per-cell regulated output and other sources in the same preparation, challenge interval and units. Exceedance rejects that bounded seed-plus-other account; it does not assign alpha ancestry. Permanent lineage information and immune/durability measurements remain separate requirements. None of these output bounds is supplied by this signature.
