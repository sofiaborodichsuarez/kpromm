"""
Extends p_simulation_full_comparison.py (kproto, PAM-Gower, HC-Gower,
k-PROMMSA) with two more established baselines: tandem FAMD+k-means and
tandem FAMD+Student-t-mixture (Peel & McLachlan 2000 t-distribution EM
applied to FAMD factor scores, mirroring the "tandem analysis" design used
in Jimeno, Roy & Tortora 2021's KAMILA/k-prototypes benchmark study).
kamila, clustMD, and DIBmix remain R-only (no Python port); reported via
the companion R script.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from p_simulation_factors import base_data
import methods_lib as ml

t0 = time.time()
p_grid = (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
N, K, NNUM, NCAT, CATLV, SEP = 220, 3, 4, 4, 4, 1.5
n_init, max_iter = 6, 10
reps = 10

factors = {
    "baseline":            dict(),
    "contam_mild":         dict(contam_frac=0.03, contam_scale=8.0),
    "contam_severe":       dict(contam_frac=0.10, contam_scale=8.0),
    "heavytail_t2":        dict(dist="t2"),
    "skew_lognormal":      dict(dist="lognormal"),
    "unequal_var_mild":    dict(cluster_sds=np.array([0.7, 1.0, 1.6])),
    "unequal_var_severe":  dict(cluster_sds=np.array([0.4, 1.0, 3.0])),
    "correlation_mild":    dict(corr=0.4),
    "correlation_high":    dict(corr=0.8),
    "cat_imbalance_mild":  dict(cat_imbalance=0.6),
    "cat_imbalance_severe":dict(cat_imbalance=0.9),
    "unequal_size_mild":   dict(priors=[0.6, 0.25, 0.15]),
    "unequal_size_severe": dict(priors=[0.8, 0.15, 0.05]),
}

results = {}
for name, kwargs in factors.items():
    acc = {m: [] for m in ["kproto", "pam", "hc", "famd_km", "famd_t", "p2", "bestp"]}
    tim = {m: [] for m in acc}
    for rep in range(reps):
        rng = np.random.default_rng(hash((name, rep, "fullcomp_v2")) % (2**32))
        X_num, X_cat, labels = base_data(N, K, NNUM, NCAT, CATLV, SEP, rng, **kwargs)

        a, t = ml.kproto_run(X_num, X_cat, labels, K, n_init, rng)
        acc["kproto"].append(a); tim["kproto"].append(t)

        a, t = ml.pam_gower_run(X_num, X_cat, labels, K, n_init, rng)
        acc["pam"].append(a); tim["pam"].append(t)

        a, t = ml.hc_gower_run(X_num, X_cat, labels, K)
        acc["hc"].append(a); tim["hc"].append(t)

        a, t = ml.famd_kmeans_run(X_num, X_cat, labels, K, n_init, rng)
        acc["famd_km"].append(a); tim["famd_km"].append(t)

        a, t = ml.famd_tmixture_run(X_num, X_cat, labels, K, n_init, rng)
        acc["famd_t"].append(a); tim["famd_t"].append(t)

        a, t = ml.kmpmnorm_p_run(X_num, X_cat, labels, K, 2.0, n_init, max_iter, seed=rep)
        acc["p2"].append(a); tim["p2"].append(t)

        a, t, _ = ml.kmpmnorm_bestp_run(X_num, X_cat, labels, K, p_grid, n_init, max_iter, seed=rep)
        acc["bestp"].append(a); tim["bestp"].append(t)

    row = {}
    for m in acc:
        row[f"{m}_acc"] = ml.ci95(acc[m])
        row[f"{m}_time"] = float(np.nanmean(tim[m]))
    results[name] = row
    print(f"{name:20s} " + "  ".join(f"{m}={row[f'{m}_acc'][0]:.3f}" for m in acc))

print(f"\nelapsed: {time.time()-t0:.1f}s")
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "p_simulation_full_comparison_v2_results.json")
with open(out_path, "w") as f:
    json.dump(results, f, indent=2, default=str)
print(f"wrote {out_path}")
