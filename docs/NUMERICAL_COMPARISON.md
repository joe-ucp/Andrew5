# Numerical comparison across platforms

The reproduction runner preserves the original scores and scientific tables. Its default numerical comparison uses relative tolerance `1e-12` and absolute tolerance `2e-14`, with byte equality reported separately. One narrowly defined comparison handles rank statistics of mathematically zero signed-log scores. It does not round scores, refit signatures or substitute a new scientific result.

## Observed cause

Paired Windows and Linux CPU runs used the same pinned source files, gene programs, counts and normalization. The Windows native scores reproduced the retained values exactly. Native signed-log scores differed between platforms by at most `2.7755575615628914e-16`; the original signed-rank and UCell-form scores were bit-identical. Five entries in the native independent-state benchmark exceeded the ordinary tolerance: one donor-1 beta score/depth Spearman correlation, and two AUROCs plus two beta score/depth correlations for donor 3. The number of affected entries is specific to this platform pair.

The isolation experiment found differences of at most `8.881784197001252e-16` in selected positive-count `log1p` terms. Exchanging identical logged operands between platforms produced bit-identical dot products. Thus the observed difference entered through the evaluated log terms and subsequent cancellation, rather than a different dot-product result for identical operands. Every changed rank involved a score whose exact cancellation was independently established from integer counts. Replacing only those proven zeros by zero in a diagnostic calculation eliminated the cross-platform rank differences; that diagnostic replacement is not applied to saved scores or tables.

For example, a donor-1 beta cell had positive up-block counts `[1, 1, 2, 5]` and the same down-block counts, with 50 genes in each block. Its mathematical signed-log score is zero. The computed residual was approximately `-1.7347e-18` on Windows and `5.2042e-18` on Linux. The donor's beta score/depth Spearman correlation was `0.16593763116267823` versus `0.1662586070888252`, a difference of `0.0003209759261469658`. Rank statistics can respond to a changed tie or ordering even when the score difference is tiny.

Floating-point addition is order-sensitive, and numerical implementations can differ in evaluated elementary functions and reductions. The [NumPy summation documentation](https://numpy.org/doc/stable/reference/generated/numpy.sum.html) describes precision differences between reduction strategies; the [IEEE floating-point background from NVIDIA](https://docs.nvidia.com/cuda/floating-point/index.html) explains how finite precision and operation ordering affect results. These are general background references. The cause above was measured in this pipeline, which used CPUs and no GPU.

## Integer-count certificate

For a cell with whole-library total `N`, define the original signed-log score

$$
L=\frac{1}{n_u}\sum_{g\in U}\log(1+10^4c_g/N)
-\frac{1}{n_d}\sum_{g\in D}\log(1+10^4c_g/N).
$$

Here `n_u` and `n_d` are the represented up- and down-block sizes for the specified gene-universe mode. They are the same denominators used by the scorer. For every positive integer count `c`, let `h_U(c)` and `h_D(c)` count genes with that count in the two blocks. The exact certificate is

$$
n_d h_U(c)=n_u h_D(c)\quad\text{for every positive count }c.
$$

This integer equality proves `L = 0`: terms with the same count have the same normalization and cancel by count class. It also handles unequal block sizes. The implementation checks nonnegative integer counts, unique gene identifiers, disjoint signature blocks, valid totals, nonempty represented blocks and equality of the evaluated log terms for equal counts within the cell. A completely undetected signature receives zero radius, because its score is already computed from zeros. A cell without the certificate also receives zero radius, regardless of how small its score is.

## Float64 radius around a certified zero

For a nonempty certified cell, let `a` be the number of detected signed genes, `w_g` be `1/n_u` or `-1/n_d`, and `ell_g` the evaluated normalized log term. Set

$$
u=\operatorname{eps}_{64}/2,\qquad
\gamma_{a+4}=\frac{(a+4)u}{1-(a+4)u},\qquad
A=\sum_g|w_g\ell_g|.
$$

The implemented radius is the next representable float above

$$
r=\frac{\gamma_{a+4}A}{1-\gamma_{a+4}}.
$$

The extra factor and outward step provide slack for evaluating the nonnegative absolute-term sum. The budget covers float64 weight, product and accumulation rounding under ordinary IEEE arithmetic without overflow or underflow. It does not require a general error bound for `log1p`: the count-class certificate cancels common evaluated log values as well as their mathematical counterparts. Unequal evaluated values for identical counts fail the common-term check. See [the implementation](../code/rank_metric_comparison.py) for the exact support checks.

Certified cells receive intervals `[-r, r]`, centered at the proven true zero. Other cells remain fixed at their freshly computed scores. Each fresh certified score must lie inside its own interval. The integration also rejects nonfinite, negative or excessively large radii. These restrictions prevent a numerical exception from accepting arbitrary score changes.

## Bounds on rank statistics

Sort the score intervals by their lower endpoints and form connected overlapping clusters. Nonoverlapping clusters have a fixed relative order. Canonical ranks assign the average rank to every cell in a cluster. A cluster containing a positive-radius cell can admit different internal orders or ties; the bound conservatively permits full internal rank splitting, even where that relaxes constraints imposed by fixed points. A cluster containing only fixed points retains its observed ties.

For beta-versus-alpha AUROC, a cross-state pair in an uncertain cluster can move from canonical tie credit `1/2` to any credit between zero and one. If `P` is the number of such pairs and `n_beta`, `n_alpha` are the state counts, the envelope is

$$
\operatorname{AUROC}_{\rm canonical}
\ \pm\ \frac{P}{2n_{\rm beta}n_{\rm alpha}},
$$

clipped to `[0, 1]`.

For Spearman correlation, let `c` be the centered canonical score ranks and `y` the centered average ranks of depth. A permitted score rank vector can be written `c + e`; within each cluster the entries of `e` sum to zero, so `c` and `e` are orthogonal. Write `V0 = c·c`, `Vy = y·y` and `T = c·y`. For each uncertain cluster of size `m`, the fully split centered ranks are `d_j = j - (m+1)/2`. Rearrangement gives a covariance-adjustment bound

$$
B=\sum_{\text{uncertain clusters}}\sum_{j=1}^{m}d_j y_{(j)},
\qquad
V_{\max}=\sum_{\text{uncertain clusters}}\frac{m(m^2-1)}{12},
$$

where `y_(j)` is sorted within its cluster. The covariance lies in `[T-B, T+B]`, while the denominator lies between `sqrt(V0 Vy)` and `sqrt((V0+Vmax) Vy)`. For the lower endpoint, divide `T-B` by the larger denominator when nonnegative and the smaller denominator when negative. For the upper endpoint, divide `T+B` by the smaller denominator when nonnegative and the larger denominator when negative. Clip the result to `[-1, 1]`. Degenerate rank variance fails certification. These are numerical rank envelopes, not confidence intervals or biological uncertainty estimates.

## Exact scope and audit record

The exception applies only to `signed_log1p_CPT10K` rows in the native [independent state benchmark](../evidence/revision_01/independent_state_benchmark.tsv), and only to these three columns:

- `beta_vs_alpha_AUROC`;
- `alpha_score_depth_r`;
- `beta_score_depth_r`.

It covers the declared programs and both declared universe modes. Thresholds, detection fractions, counts, means, thinning results, other scoring methods and all other tables retain the original `rtol=1e-12`, `atol=2e-14` comparison. The elapsed-seconds field of `population_audit.tsv` remains the sole non-scientific column excluded from comparison.

An eligible mismatch must be recomputed from the **fresh run's** native NPZ scores, state labels and library totals, with matching program and method indices and cell counts. Both the fresh and retained metric must lie inside the derived envelope, allowing only the ordinary final numerical comparison tolerance. Failure to recompute, missing generated arrays, an invalid radius or an out-of-envelope metric causes failure. The comparator never overwrites a saved table; an internal comparison copy is used only after a mismatch has passed the certificate.

The stage record stores `exact_zero_rank_metric_comparisons` with the donor, mode, program, method, metric, retained value, fresh value, recomputed value, lower and upper bounds, canonical value, cluster diagnostics and acceptance reason. Byte equality remains a separate field. [Synthetic tests](../tests/test_rank_metric_comparison.py) cover cancellation certificates, rank bounds and rejection cases. This policy explains a demonstrated numerical degeneracy; it is not a guarantee for every future library or platform implementation.

## Run and inspect

```sh
python code/reproduce.py stage independent --data-root /absolute/external/baron-data --output-dir /absolute/external/independent-run --fetch
```

Use a new external output folder. This standalone command fetches the four registered Baron count files when absent, uses the included fixed program prerequisites, and executes the complete independent native and thinning stage. Inspect the six generated comparison tables, `receipts/independent.json`, and the native NPZ files under `data/revision_01/independent/`. [Reproducibility](REPRODUCIBILITY.md) explains the broader prerequisite chain. A standalone independent-stage run does not establish that the complete 18-stage reconstruction has been rerun.
