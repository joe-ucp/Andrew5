"""Revision 01 representative nuisance precision and profile/depth coupling.

This extends conditional numerical diagnostics, never estimates real nuisance
rates or a biological confidence interval. Original full-family results remain
unchanged. Requires the exact approved revision plan before execution.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os
import numpy as np
import pandas as pd
from scipy import sparse
from data_paths import data_path, resolve_data_root, ENV_NAME, evidence_path, logical_evidence_key
from primary_analysis import map_universe, prob_profiles, score_sparse, full_gene
from reference_analysis import CELL_CYCLE
from revision_mechanism import read_raw_sample, load_protocol

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "evidence/revision_01"
from public_runtime import parameters_sha256
MODELS = ["beta_RNA_0", "beta_RNA_0.1", "beta_RNA_0.2",
          "doublet_fraction_0.1", "cycle_high_vs_low"]
MODES = ["independent_profile_depth", "paired_profile_own_depth"]
SIZES = [400, 1600, 6400]
DRAWS = [200, 1000]
BATCHES = 3


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda: f.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def sample_profile_pool(profiles, own_depths, global_depths, weights, rng,
                        n, paired, beta=None, fraction=0.):
    """Multinomial observation model with explicit profile/depth coupling.

    Sparse positive supports plus one exact outside category give the same
    multinomial law as the full vector, avoiding arithmetic on zero probabilities.
    This is a nuisance model, not without-replacement thinning.
    """
    profiles = profiles.tocsc()
    if beta is not None:
        beta = beta.tocsc()
    choices = rng.integers(0, profiles.shape[1], n)
    depths = np.asarray(own_depths)[choices] if paired else rng.choice(global_depths, n)
    beta_choices = rng.integers(0, beta.shape[1], n) if beta is not None else np.full(n, -1)
    indices, values, pointers = [], [], [0]
    G = len(weights)
    for j, cell in enumerate(choices):
        a, b = profiles.indptr[cell:cell+2]
        ids, probs = profiles.indices[a:b], profiles.data[a:b]
        if beta is not None and fraction != 0:
            c, d = beta.indptr[beta_choices[j]:beta_choices[j]+2]
            b_ids, b_probs = beta.indices[c:d], beta.data[c:d]
            support = np.union1d(ids, b_ids)
            mixed = np.zeros(len(support))
            mixed[np.searchsorted(support, ids)] += (1-fraction) * probs
            mixed[np.searchsorted(support, b_ids)] += fraction * b_probs
            ids, probs = support, mixed
        positive = probs > 0
        ids, probs = ids[positive], probs[positive]
        if not len(probs) or not np.isclose(probs.sum(), 1., rtol=0, atol=1e-12):
            raise ValueError("Empirical profile does not preserve whole-library probability mass")
        probs = probs / probs.sum()
        counts = rng.multinomial(int(depths[j]), probs)
        keep = (ids < G) & (counts > 0)
        indices.append(ids[keep])
        values.append(counts[keep])
        pointers.append(pointers[-1] + int(keep.sum()))
    matrix = sparse.csc_matrix((np.concatenate(values), np.concatenate(indices),
                                np.asarray(pointers)), shape=(G, n))
    return dict(score=score_sparse(matrix, weights), profile_index=choices,
                beta_profile_index=beta_choices, observed_depth=depths)


def nested_doublet_mask(n, fraction, rng, block_size=400):
    """Exact declared fraction in each block, hence every planned pool prefix."""
    if n % block_size or fraction * block_size != int(fraction * block_size):
        raise ValueError("Nested doublet mask requires integral block counts")
    selected = np.zeros(n, dtype=bool)
    for start in range(0, n, block_size):
        selected[start + rng.choice(block_size, int(fraction*block_size), replace=False)] = True
    return selected


def group_draws(scenario, baseline, n_treated, n_control, seed, count=1000):
    rng = np.random.default_rng(seed)
    result = np.empty(count)
    for begin in range(0, count, 25):
        end = min(begin+25, count)
        result[begin:end] = (rng.choice(scenario, (end-begin, n_treated)).mean(axis=1) -
                             rng.choice(baseline, (end-begin, n_control)).mean(axis=1))
    return result


def main():
    plan, cfg = load_protocol()
    design = plan["nuisance_precision"]
    if design["seed_base"] != 7900 or design["batches"] != BATCHES or design["pool_sizes"] != SIZES or design["group_draw_counts"] != DRAWS:
        raise ValueError("Approved fixed nuisance design differs from this implementation")
    names = ["nuisance_precision_batches.tsv", "nuisance_precision_summary.tsv",
             "nuisance_pool_drift.tsv", "nuisance_coupling_differences.tsv",
             "nuisance_population_audit.tsv", "nuisance_precision_status.json"]
    bulk = data_path("revision_01/nuisance", for_write=True)
    for path in [OUT/name for name in names] + [bulk/"all_group_draws.tsv.gz"]:
        if path.exists():
            raise FileExistsError("Refusing to overwrite revision output: "+str(path))
    manifest = pd.read_csv(ROOT/"primary_sample_manifest.tsv", sep="\t")
    vehicle = manifest[manifest.treatment_code == "D"].sort_values("donor")
    historical_paths = [evidence_path(name, must_exist=True) for name in ["primary_cell_scores.tsv.gz", "primary_results.tsv",
                         "countermodel_results.tsv", "countermodel_stability_summary.tsv", "frozen_score.json"]]
    original_hashes = {logical_evidence_key(p): sha(p) for p in historical_paths}
    cells = pd.read_csv(historical_paths[0], sep="\t", float_precision="round_trip")
    results = pd.read_csv(historical_paths[1], sep="\t", float_precision="round_trip")
    targets = results[(results.contrast == "H_GLP-D") & (results.annotation == "alpha_related") &
                      (results.sensitivity == "primary")].set_index("donor")
    sources = pd.read_csv(ROOT/"source_manifest.tsv", sep="\t").set_index("local_path")
    input_hashes = {}
    for row in vehicle.to_dict("records"):
        for field in ["matrix_path", "features_path", "barcodes_path"]:
            relative = row[field].replace("\\", "/")
            actual = sha(data_path(relative, must_exist=True))
            if actual != sources.loc[relative, "sha256"]:
                raise ValueError("Raw source hash mismatch: "+relative)
            input_hashes[relative] = actual
        for batch in range(BATCHES):
            for mode in MODES:
                path = bulk/f"{row['donor']}_batch{batch}_{mode}_pools.npz"
                if path.exists():
                    raise FileExistsError("Refusing to overwrite pool cache: "+str(path))
    OUT.mkdir(exist_ok=True)
    bulk.mkdir(parents=True, exist_ok=True)
    universe = cfg["universe"]
    idx = {g: i for i, g in enumerate(universe)}
    weights = np.zeros(len(universe))
    weights[[idx[g] for g in cfg["up"]]] = 1/len(cfg["up"])
    weights[[idx[g] for g in cfg["down"]]] = -1/len(cfg["down"])
    records, all_draws, audits, caches = [], [], [], []
    start_time = datetime.datetime.now(datetime.timezone.utc).isoformat()
    for row in vehicle.to_dict("records"):
        donor = row["donor"]
        number = int(donor[1:])
        x, features, barcodes = read_raw_sample(row)
        saved = cells[cells.accession == row["accession"]].set_index("barcode")
        selected = np.flatnonzero(barcodes.barcode.isin(saved.index).to_numpy())
        aligned = saved.loc[barcodes.barcode.iloc[selected]]
        if len(aligned) != len(saved):
            raise ValueError("Retained vehicle barcode mismatch")
        totals = np.asarray(x[:, selected].sum(axis=0)).ravel()
        if not np.array_equal(totals, aligned.library_counts):
            raise ValueError("Retained vehicle count totals do not match raw counts")
        cycle = np.asarray(x[np.flatnonzero(features.symbol.isin(CELL_CYCLE)), :][:, selected].sum(axis=0)).ravel() / np.maximum(totals, 1)
        if not np.allclose(cycle, aligned.cycle_fraction, atol=1e-15, rtol=0):
            raise ValueError("Retained vehicle cycle fractions do not match raw counts")
        mapped = map_universe(x[:, selected], features, universe)
        profiles = prob_profiles(mapped, totals)
        am = aligned.annotation.eq("alpha_related").to_numpy()
        bm = aligned.annotation.eq("beta_like").to_numpy()
        alpha, beta = profiles[:, am], profiles[:, bm]
        alpha_depths, alpha_cycle = totals[am], cycle[am]
        q25, q75 = np.quantile(alpha_cycle, [.25, .75])
        low = alpha_cycle <= q25
        high = (alpha_cycle >= q75) & (alpha_cycle > 0) & ~low
        cycle_valid = (alpha_cycle > 0).sum() >= 20 and low.sum() >= 20 and high.sum() >= 20
        target = targets.loc[donor]
        n0, n1 = int(target.n_control), int(target.n_treated)
        if n0 != am.sum():
            raise ValueError("Actual control group size changed")
        audits.append(dict(donor=donor,control_alpha_n=n0,treated_alpha_n=n1,beta_n=int(bm.sum()),
                           cycle_q25=q25,cycle_q75=q75,cycle_low_n=int(low.sum()),cycle_high_n=int(high.sum()),
                           cycle_positive_n=int((alpha_cycle>0).sum()),cycle_evaluable=bool(cycle_valid)))
        del x, mapped, profiles
        for batch in range(BATCHES):
            for mode_index, mode in enumerate(MODES):
                paired = mode_index == 1
                def pool(stream, prof=alpha, own=alpha_depths, beta_profiles=None, fraction=0., n=6400):
                    seed = np.random.SeedSequence([7900, number, batch, mode_index, stream])
                    return sample_profile_pool(prof, own, alpha_depths, weights, np.random.default_rng(seed),
                                               n, paired, beta_profiles, fraction)
                base = pool(0)
                p10 = pool(1, beta_profiles=beta, fraction=.1)
                p20 = pool(2, beta_profiles=beta, fraction=.2)
                hit = nested_doublet_mask(6400, .1, np.random.default_rng(np.random.SeedSequence([7900, number, batch, mode_index, 3])))
                mixture = pool(4, beta_profiles=beta, fraction=.5, n=int(hit.sum()))
                doublet = base["score"].copy()
                doublet[hit] = mixture["score"]
                pools = {"baseline":base["score"], "beta_RNA_0":base["score"].copy(),
                         "beta_RNA_0.1":p10["score"], "beta_RNA_0.2":p20["score"],
                         "doublet_fraction_0.1":doublet}
                metadata = {f"{label}_{key}":value for label,p in [("baseline",base),("mix10",p10),("mix20",p20),("doublet_replacements",mixture)]
                            for key,value in p.items() if key != "score"}
                metadata["doublet_mask"] = hit
                if cycle_valid:
                    lo = pool(5, prof=alpha[:,low], own=alpha_depths[low])
                    hi = pool(6, prof=alpha[:,high], own=alpha_depths[high])
                    pools["cycle_low"] = lo["score"]
                    pools["cycle_high_vs_low"] = hi["score"]
                    metadata.update({f"{label}_{key}":value for label,p in [("cycle_low",lo),("cycle_high",hi)]
                                     for key,value in p.items() if key != "score"})
                cache = bulk/f"{donor}_batch{batch}_{mode}_pools.npz"
                np.savez_compressed(cache, **pools, **metadata)
                caches.append(dict(path=str(cache),sha256=sha(cache)))
                for model_index, model in enumerate(MODELS):
                    if model not in pools:
                        continue
                    baseline = pools["cycle_low"] if model == "cycle_high_vs_low" else pools["baseline"]
                    for size in SIZES:
                        scenario, baseline_prefix = pools[model][:size], baseline[:size]
                        seed = np.random.SeedSequence([7900, number, batch, mode_index, model_index, 99])
                        draw = group_draws(scenario, baseline_prefix, n1, n0, seed, 1000)
                        effect = float(scenario.mean()-baseline_prefix.mean())
                        if model == "beta_RNA_0" and effect != 0:
                            raise AssertionError("Zero model must use identical baseline")
                        for j, value in enumerate(draw):
                            all_draws.append(dict(donor=donor,batch=batch,coupling=mode,model=model,pool_n=size,
                                                  draw=j,group_mean_delta=float(value)))
                        for count in DRAWS:
                            use = draw[:count]
                            records.append(dict(donor=donor,batch=batch,coupling=mode,model=model,pool_n=size,
                                                group_draws=count,n_control=n0,n_treated=n1,
                                                conditional_pool_mean_effect=effect,
                                                conditional_q95=float(np.quantile(use,.95)),
                                                draw_mean=float(use.mean()),draw_SD=float(use.std(ddof=1)),
                                                draw_min=float(use.min()),draw_max=float(use.max()),
                                                observed_native_delta=float(target.delta)))
                print(f"{donor} batch {batch+1}/3 {mode}: nested 400/1600/6400 pools complete",flush=True)
        del alpha, beta
    batch_table = pd.DataFrame(records)
    draw_table = pd.DataFrame(all_draws)
    summary = []
    for keys, part in batch_table.groupby(["donor","coupling","model","pool_n","group_draws"],sort=True):
        d,m,model,size,count = keys
        ds = draw_table[(draw_table.donor==d)&(draw_table.coupling==m)&(draw_table.model==model)&
                        (draw_table.pool_n==size)&(draw_table.draw<count)].group_mean_delta
        summary.append(dict(donor=d,coupling=m,model=model,pool_n=size,group_draws_per_batch=count,
                            batches=len(part),total_draws=len(ds),pooled_q95=float(ds.quantile(.95)),
                            pooled_draw_min=float(ds.min()),pooled_draw_max=float(ds.max()),
                            batch_q95_min=float(part.conditional_q95.min()),batch_q95_max=float(part.conditional_q95.max()),
                            conditional_effect_mean=float(part.conditional_pool_mean_effect.mean()),
                            conditional_effect_min=float(part.conditional_pool_mean_effect.min()),
                            conditional_effect_max=float(part.conditional_pool_mean_effect.max()),
                            observed_native_delta=float(part.observed_native_delta.iloc[0]),
                            scope="conditional numerical sensitivity; three batches do not prove convergence or coverage"))
    summary = pd.DataFrame(summary)
    drift = []
    for keys, group in summary.groupby(["donor","coupling","model"]):
        reference = group[(group.pool_n==6400)&(group.group_draws_per_batch==1000)].iloc[0]
        for _, row in group.iterrows():
            drift.append(dict(donor=keys[0],coupling=keys[1],model=keys[2],pool_n=row.pool_n,
                              group_draws_per_batch=row.group_draws_per_batch,pooled_q95=row.pooled_q95,
                              reference_q95_6400_1000=reference.pooled_q95,
                              absolute_q95_difference=abs(row.pooled_q95-reference.pooled_q95),
                              signed_q95_difference=row.pooled_q95-reference.pooled_q95))
    keys = ["donor","model","pool_n","group_draws_per_batch"]
    ind = summary[summary.coupling==MODES[0]].set_index(keys)
    paired = summary[summary.coupling==MODES[1]].set_index(keys)
    coupling = ind.index.to_frame(index=False)
    coupling["independent_q95"] = ind.pooled_q95.to_numpy()
    coupling["paired_q95"] = paired.loc[ind.index].pooled_q95.to_numpy()
    coupling["paired_minus_independent_q95"] = coupling.paired_q95-coupling.independent_q95
    for name, table in zip(names[:5], [batch_table,summary,pd.DataFrame(drift),coupling,pd.DataFrame(audits)]):
        table.to_csv(OUT/name,sep="\t",index=False)
    draw_table.to_csv(bulk/"all_group_draws.tsv.gz",sep="\t",index=False,compression={"method":"gzip","mtime":0})
    for path in historical_paths:
        if sha(path)!=original_hashes[logical_evidence_key(path)]:
            raise AssertionError("Original evidence changed: "+str(path))
    status = dict(status="COMPLETE_REPRESENTATIVE_NUISANCE_PRECISION_DIAGNOSTICS",
                  started_utc=start_time,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  scientific_parameters_sha256=parameters_sha256(),
                  code_sha256=sha(Path(__file__)),original_evidence_unchanged=True,
                  original_evidence_sha256=original_hashes,raw_sources_sha256=input_hashes,
                  pool_sizes=SIZES,group_draw_counts=DRAWS,batches=BATCHES,models=MODELS,coupling=MODES,
                  batch_rows=len(batch_table),summary_rows=len(summary),group_draw_rows=len(draw_table),
                  pools=caches,full_draw_path=str(bulk/"all_group_draws.tsv.gz"),
                  conventions=["SeedSequence([7900,numeric_donor,batch,coupling_index,stream])",
                               "Group stream SeedSequence([7900,donor,batch,coupling,model,99]); same seed across pool-size prefixes",
                               "First200 of1000 group draws; first400/1600 of6400 synthetic cells",
                               "Exactly40 doublet replacements per400-cell block; 50 percent beta RNA in replaced profiles",
                               "Cycle baseline bottom vehicle quartile; high top quartile with positive cycle and no overlap; at least20 per extreme and20positive",
                               "Independent mode uses whole vehicle-alpha depth distribution even for cycle strata",
                               "Paired mode uses each selected vehicle alpha profile's own full-library depth",
                               "Zero condition shares identical baseline pool; group draws may still differ"],
                  inference_limits=["Selected scenarios do not replace original full-family diagnostics",
                                    "Three independent batches and finite pools do not prove quantile convergence",
                                    "Nuisance rates are assumed, not measured; no biological confidence bound",
                                    "Cycle high-versus-low baseline differs from original high-versus-whole-vehicle comparison"])
    (OUT/names[-1]).write_text(json.dumps(status,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in status.items() if k not in ["pools","raw_sources_sha256"]},indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    from public_runtime import execute_stage
    execute_stage('nuisance', args.data_root, args.output_dir)
