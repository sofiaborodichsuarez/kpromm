"""
Shared library of clustering-method runners used by both the synthetic
factor-sweep simulation and the real-labeled-dataset benchmark. Every
runner has the same signature: (X_num, X_cat, true_labels, k, n_init, rng)
-> (accuracy, elapsed_seconds), so a single driver loop can call any of
them interchangeably.

Methods implemented here (all run natively in Python):
  - k-prototypes (Huang 1998), via kmodes.KPrototypes
  - PAM on Gower's dissimilarity (Kaufman & Rousseeuw 1990; Gower 1971)
  - Hierarchical clustering (average linkage) on Gower's dissimilarity
  - Tandem FAMD + k-means (Pages 2004 FAMD; standard tandem-analysis design)
  - Tandem FAMD + Student-t mixture (Peel & McLachlan 2000 t-mixture EM)
  - k-PROMM / k-PROMMSA (Suarez-Alvarez et al. 2010/2012), via kmpmnorm_core

NOT implemented here (R-only):
  - kamila (Foss et al. 2016)
  - clustMD (McParland & Gormley 2017)
  - DIBmix (Costa, Papatsouma & Markos 2024)
These R-only methods are handled by the companion R scripts
simulation_comparison.R and real_data_benchmark.R.
"""
import time
import numpy as np
import pandas as pd

np.infty = np.inf  # NumPy 2.0 compat shim required by the `smm` package

import gower
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.optimize import linear_sum_assignment
from kmodes.kprototypes import KPrototypes
from sklearn.cluster import KMeans
import prince
from smm import SMM

import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from kmpmnorm_core import kpromm as _kpromm_fit


def assignment_accuracy(true_labels, pred_labels):
    labs_t = np.unique(true_labels)
    labs_p = np.unique(pred_labels)
    kk = max(len(labs_t), len(labs_p))
    conf = np.zeros((kk, kk))
    map_t = {v: i for i, v in enumerate(labs_t)}
    map_p = {v: i for i, v in enumerate(labs_p)}
    for t, pr in zip(true_labels, pred_labels):
        conf[map_p[pr], map_t[t]] += 1
    row, col = linear_sum_assignment(-conf)
    return conf[row, col].sum() / len(true_labels)


def _mixed_frame(X_num, X_cat):
    df = pd.DataFrame(X_num, columns=[f"num{i}" for i in range(X_num.shape[1])])
    for i in range(X_cat.shape[1]):
        df[f"cat{i}"] = X_cat[:, i].astype(str)
    return df


def gower_dist_matrix(X_num, X_cat):
    df = np.hstack([X_num.astype(object), X_cat.astype(object)])
    cat_features = np.array([False] * X_num.shape[1] + [True] * X_cat.shape[1])
    return gower.gower_matrix(df, cat_features=cat_features).astype(np.float64)


def kproto_run(X_num, X_cat, true_labels, k, n_init, rng):
    X = np.hstack([X_num.astype(object), X_cat.astype(object)])
    cat_idx = list(range(X_num.shape[1], X_num.shape[1] + X_cat.shape[1]))
    best_acc, best_cost = None, np.inf
    t0 = time.perf_counter()
    for _ in range(n_init):
        try:
            kp = KPrototypes(n_clusters=k, n_init=1, max_iter=100,
                              random_state=int(rng.integers(0, 1_000_000)), verbose=0)
            labels = kp.fit_predict(X, categorical=cat_idx)
        except Exception:
            continue
        if kp.cost_ < best_cost:
            best_cost = kp.cost_
            best_acc = assignment_accuracy(true_labels, labels)
    return best_acc, time.perf_counter() - t0


def pam_gower_run(X_num, X_cat, true_labels, k, n_init, rng):
    t0 = time.perf_counter()
    D = gower_dist_matrix(X_num, X_cat)
    N = D.shape[0]
    best_acc, best_cost = None, np.inf
    for _ in range(n_init):
        medoids = rng.choice(N, size=k, replace=False)
        for _ in range(100):
            assign = np.argmin(D[:, medoids], axis=1)
            new_medoids = medoids.copy()
            changed = False
            for c in range(k):
                members = np.where(assign == c)[0]
                if len(members) == 0:
                    continue
                sub = D[np.ix_(members, members)]
                new_med = members[np.argmin(sub.sum(axis=1))]
                if new_med != medoids[c]:
                    changed = True
                new_medoids[c] = new_med
            medoids = new_medoids
            if not changed:
                break
        assign = np.argmin(D[:, medoids], axis=1)
        cost = D[np.arange(N), medoids[assign]].sum()
        if cost < best_cost:
            best_cost = cost
            best_acc = assignment_accuracy(true_labels, assign)
    return best_acc, time.perf_counter() - t0


def hc_gower_run(X_num, X_cat, true_labels, k):
    t0 = time.perf_counter()
    D = gower_dist_matrix(X_num, X_cat)
    Z = linkage(squareform(D, checks=False), method="average")
    assign = fcluster(Z, t=k, criterion="maxclust")
    acc = assignment_accuracy(true_labels, assign)
    return acc, time.perf_counter() - t0


def _famd_scores(X_num, X_cat, n_components, rng):
    df = _mixed_frame(X_num, X_cat)
    n_comp = min(n_components, df.shape[1] - 1, df.shape[0] - 1)
    n_comp = max(n_comp, 2)
    famd = prince.FAMD(n_components=n_comp, random_state=int(rng.integers(0, 1_000_000)))
    famd = famd.fit(df)
    return famd.transform(df).values


def famd_kmeans_run(X_num, X_cat, true_labels, k, n_init, rng, n_components=5):
    t0 = time.perf_counter()
    scores = _famd_scores(X_num, X_cat, n_components, rng)
    km = KMeans(n_clusters=k, n_init=n_init, random_state=int(rng.integers(0, 1_000_000))).fit(scores)
    acc = assignment_accuracy(true_labels, km.labels_)
    return acc, time.perf_counter() - t0


def famd_tmixture_run(X_num, X_cat, true_labels, k, n_init, rng, n_components=5):
    t0 = time.perf_counter()
    scores = _famd_scores(X_num, X_cat, n_components, rng)
    try:
        mix = SMM(n_components=k, n_init=max(n_init, 5), n_iter=200,
                  random_state=int(rng.integers(0, 1_000_000)))
        mix.fit(scores)
        pred = mix.predict(scores)
        acc = assignment_accuracy(true_labels, pred)
    except Exception:
        acc = np.nan
    return acc, time.perf_counter() - t0


def kmpmnorm_p_run(X_num, X_cat, true_labels, k, p, n_init, max_iter, seed, norm_subsample=None):
    t0 = time.perf_counter()
    fit = _kpromm_fit(X_num, X_cat, k, p, normalize_num=True, normalize_cat=True,
                       n_init=n_init, max_iter=max_iter, seed=seed, norm_subsample=norm_subsample)
    acc = assignment_accuracy(true_labels, fit["assign"])
    return acc, time.perf_counter() - t0


def kmpmnorm_bestp_run(X_num, X_cat, true_labels, k, p_grid, n_init, max_iter, seed, norm_subsample=None):
    t0 = time.perf_counter()
    accs = {}
    for p in p_grid:
        fit = _kpromm_fit(X_num, X_cat, k, p, normalize_num=True, normalize_cat=True,
                           n_init=n_init, max_iter=max_iter, seed=seed, norm_subsample=norm_subsample)
        accs[p] = assignment_accuracy(true_labels, fit["assign"])
    return max(accs.values()), time.perf_counter() - t0, accs


def ci95(vals):
    vals = np.asarray([v for v in vals if v is not None and not np.isnan(v)])
    if len(vals) == 0:
        return (np.nan, np.nan, np.nan)
    m = vals.mean()
    se = vals.std(ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else 0.0
    return m, m - 1.96 * se, m + 1.96 * se
