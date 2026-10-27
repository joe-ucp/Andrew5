# What cross-sectional expression can and cannot identify

This appendix supplies a constructive check of a known identification limitation. It is **not the run's main novel result**, not a fitted biological model, and not a substitute for the external calibration and four donor contrasts.

## Observed experiment and retained information

GSE266427 contains four donor-matched treatment conditions at one 96-hour sampling time. Donors are D1-D4; the primary comparison is harmine plus native GLP-1(7-36) versus vehicle. Each deposited sample has captured expression profiles and library identifiers, but no persistent ancestry label, longitudinal cell pairing, absolute starting endocrine-cell count, or secretion measurement for each candidate cell. The drug study's cells are not the cells in the PDX1/MAFA reference. The latter can constrain expression resemblance without assigning the former's ancestry.

Let j index donor and treatment, S=1 denote capture plus the declared quality/annotation selection, and Y denote all observed expression and technical covariates. Let Z be historical lineage and R be a historical transition. The observed law is

    p_j(y | S=1) = sum_(z,r) p_j(y | Z=z,R=r,S=1) w_j(z,r | S=1).

Neither Z nor R is measured. The weights are proportions among captured/selected cells; they are not tissue fractions, absolute counts, or conversion rates. Differential capture, survival, and selection may affect them.

## Constructive observational equivalence

Partition the observed profiles into score/annotation strata k with observed counts n_jk, total n_j, and empirical probabilities p_jk=n_jk/n_j. Every retained nonempty stratum has p_jk>0. Choose any 0<=m_jk<=n_jk within a candidate stratum. Compare two latent allocations:

| Observed stratum k | Recruitment history | Loss-of-identity history |
|---|---:|---:|
| Stable alpha-related captured cells | n_jk-m_jk | n_jk-m_jk |
| Alpha-origin cells with a transition | m_jk | 0 |
| Beta-origin cells with identity loss | 0 | m_jk |
| Observed total | n_jk | n_jk |

Both allocations preserve every observed profile, score, donor, treatment and selection decision. Emission distributions can be chosen to equal the conditional empirical distribution of the assigned profiles. Thus both assign the same strictly positive observed-data likelihood:

    L_R = L_D = product_j [n_j! / product_k n_jk!] product_k p_jk^(n_jk).

The construction can retain full expression vectors, not merely scores: swapping unmeasured ancestry labels leaves each observed vector unchanged. Score strata make the accounting auditable without implying that a score is a sufficient statistic for biological history.

`code/verify_identifiability.py` constructs a CSV from the actual cell-score table. Fixed quarter-width bins span [-1,1]. In each nonnegative-score alpha-related bin it reassigns ceil(n_jk/2) latent labels; that fraction is an arbitrary witness choice, **not an estimated biological rate**. Other observations retain their category and count. The verifier reopens the CSV, checks both row-sum reconstructions and every sample total, and calculates identical conditional log-likelihoods. The result explicitly records whether the supplied table covers all quality-controlled cells, only alpha-related cells, or the selected scored alpha/beta subset; no unobserved capture proportion is invented. The current primary pipeline scores the selected alpha/beta subset, so its correct execution flag is `--input-scope selected_scored_cells`.

## Assumptions and scope

Equivalence holds in the nonparametric model without externally fixed lineage-specific expression laws, initial source counts, transition-rate bounds, or a valid persistent ancestry label. It demonstrates insufficient information, not equal biological probability. A beta-origin explanation could be biologically implausible under additional measured constraints; these observations alone do not supply those constraints. An unrestricted latent model always has more freedom than a specified mechanistic model, so the equal likelihood is not a successful empirical prediction by the beta-loss mechanism.

Donor pairing and a common treatment time improve the response contrast but do not track individual cells over time. The two histories can be embedded in donor aliquots with a common unmeasured source reservoir and different transition, survival or capture probabilities. The construction does not estimate such reservoirs or probabilities. No measured absolute cell number is violated because none is used; if absolute starting/source counts were supplied, they would constrain and might eliminate some witnesses.

A positive frozen score therefore remains a molecular state observation. Independent source restrictions can break equivalence: validated persistent lineage labels, informative longitudinal perturbations with known lineage specificity, or quantitative source/output bounds that a rival history cannot satisfy. Function and durable immune control still require their own measurements.

## Seed-contamination sensitivity without invented biology

For donor d and a fixed assay/time window, define:

    Q_seed,d^U = N_beta0,d^U * G_beta,d^U * q_beta,d^U.

N_beta0,d^U bounds the absolute starting beta-cell contamination, G_beta,d^U bounds viable descendants per starting beta cell under the exact exposure and interval, and q_beta,d^U bounds output per surviving descendant under the same glucose assay. This is a conservative source-capacity envelope. It is not a prediction obtained by multiplying favorable published means.

The measurements required are distinct:

1. **Initial seed:** absolute input cell count and validated beta-cell purity/sensitivity, including rare contaminant detection limits. A captured beta-cell RNA fraction is insufficient without a capture model and input count.
2. **Expansion/survival:** a lineage-resolved count of descendants under the same treatment, tissue support, time and immune conditions. A proliferation marker, unrelated graft growth factor, or an assumed doubling rate is not this bound.
3. **Per-cell output:** secretion over the declared glucose challenge, viable cell counts and a conservative bound covering heterogeneity and assay uncertainty. Intracellular insulin or an expression score cannot supply it.
4. **Attribution errors:** false-positive alpha labels on starting beta cells, label transfer, label loss, assay background and output from other endocrine sources. If a beta cell can acquire the purported alpha-lineage label, its maximum contribution must enter a separate false-attribution envelope.

The bounds need simultaneous coverage. Statistical independence is unnecessary if each bound has stated one-sided error alpha_N, alpha_G, alpha_q: a union bound gives joint coverage at least 1-(alpha_N+alpha_G+alpha_q). Shared assays and donor dependence must still be retained in estimation. A joint uncertainty model can be tighter, but independent multiplication of point estimates is invalid.

For total regulated output, choose q_beta^U to bound the same endpoint (for example, high-minus-low glucose C-peptide released over a specified interval) and measure low-glucose shutoff separately. Let Q_obs^L be a lower bound on observed regulated output and Q_other^U cover all other unassigned sources. A sufficient source-capacity inequality is

    Q_obs^L > Q_seed^U + Q_other^U.

Exceeding it rejects the specified existing-seed explanation under its bounds; it does not by itself prove alpha ancestry. For lineage-attributed output use the analogous false-label envelope and validate the persistent alpha label. When a bound is absent or unbounded, the inequality cannot be evaluated and ancestry remains unresolved.

All three seed factors have unit elasticity: a twofold relaxation of one upper bound doubles Q_seed^U. Their units must multiply to the measured output unit. This symbolic sensitivity identifies where better measurement matters without fabricating contamination levels, growth rates, secretory capacities or donor sample sizes. The future experiment must preserve glucose regulation, surviving alpha-cell counterregulation, withdrawal durability, and immune-challenge function in the same system.
