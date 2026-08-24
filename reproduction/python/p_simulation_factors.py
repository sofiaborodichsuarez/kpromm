"""
Descriptive simulation study (project brief Section 3): how does the
accuracy-optimal Minkowski exponent p shift with data characteristics not
yet tested in Phase 3/4 -- contamination/outliers, heavy tails, skewness,
unequal cluster variances, correlation among numeric attributes, categorical
imbalance, and unequal cluster sizes? This does NOT test a selection
heuristic; it only characterizes, using ground-truth labels for evaluation,
how the best-performing p moves as each factor is dialed up. No real data
required.
"""
import numpy as np
from scipy.optimize import linear_sum_assignment
import json, time

def bisect_minimize(v, p, tol=1e-5, max_iter=50):
    if len(v) == 0:
        return 0.0
    a, b = v.min(), v.max()
    if b - a < tol:
        return (a + b) / 2
    for _ in range(max_iter):
        c = (a + b) / 2
        diff = c - v
        w = np.abs(v - c) ** (p - 1) if p > 1 else np.ones_like(v)
        g = np.sum(w * np.sign(diff))
        if abs(g) < 1e-9 or (b - a) < tol:
            return c
        if g > 0:
            b = c
        else:
            a = c
    return (a + b) / 2

def normalize_weights(X_num, X_cat, p):
    N = X_num.shape[0]
    alpha = np.zeros(X_num.shape[1])
    for j in range(X_num.shape[1]):
        v = X_num[:, j]
        s = (np.abs(v[:, None] - v[None, :]) ** p).sum() / (N * (N - 1))
        alpha[j] = 1.0 / s if s > 1e-12 else 0.0
    beta = np.zeros(X_cat.shape[1])
    for j in range(X_cat.shape[1]):
        v = X_cat[:, j]
        s = (v[:, None] != v[None, :]).astype(float).sum() / (N * (N - 1))
        beta[j] = 1.0 / s if s > 1e-12 else 0.0
    return alpha, beta

def kmpmnorm_fit(X_num, X_cat, alpha, beta, k, p, max_iter, rng):
    N, n_num = X_num.shape
    n_cat = X_cat.shape[1]
    init_idx = rng.choice(N, size=k, replace=False)
    Q_num = X_num[init_idx].copy()
    Q_cat = X_cat[init_idx].copy()
    labels = np.zeros(N, dtype=int)
    for it in range(max_iter):
        d = np.zeros((N, k))
        for m in range(k):
            d[:, m] = (alpha[None, :] * np.abs(X_num - Q_num[m][None, :]) ** p).sum(axis=1)
            d[:, m] += (beta[None, :] * (X_cat != Q_cat[m][None, :]).astype(float)).sum(axis=1)
        new_labels = d.argmin(axis=1)
        if it > 0 and np.array_equal(new_labels, labels):
            labels = new_labels
            break
        labels = new_labels
        for m in range(k):
            members = labels == m
            if members.sum() == 0:
                continue
            for j in range(n_num):
                Q_num[m, j] = bisect_minimize(X_num[members, j], p)
            for j in range(n_cat):
                vals, counts = np.unique(X_cat[members, j], return_counts=True)
                Q_cat[m, j] = vals[np.argmax(counts)]
    d = np.zeros((N, k))
    for m in range(k):
        d[:, m] = (alpha[None, :] * np.abs(X_num - Q_num[m][None, :]) ** p).sum(axis=1)
        d[:, m] += (beta[None, :] * (X_cat != Q_cat[m][None, :]).astype(float)).sum(axis=1)
    obj = d[np.arange(N), labels].sum()
    return labels, obj

def assignment_accuracy(true_labels, pred_labels, k):
    conf = np.zeros((k, k), dtype=int)
    for t, p_ in zip(true_labels, pred_labels):
        conf[p_, t] += 1
    row, col = linear_sum_assignment(-conf)
    return conf[row, col].sum() / len(true_labels)

def best_of_accuracy(X_num, X_cat, true_labels, k, p, n_init, max_iter, rng):
    alpha, beta = normalize_weights(X_num, X_cat, p)
    best_obj, best_acc = np.inf, None
    for _ in range(n_init):
        lab, obj = kmpmnorm_fit(X_num, X_cat, alpha, beta, k, p, max_iter, rng)
        if obj < best_obj:
            best_obj = obj
            best_acc = assignment_accuracy(true_labels, lab, k)
    return best_acc

# ---------- data generators, one per factor ----------
def base_data(n, k, n_num, n_cat, cat_levels, sep, rng, priors=None, corr=0.0,
              cluster_sds=None, dist="gaussian", contam_frac=0.0, contam_scale=10.0,
              cat_imbalance=0.0):
    if priors is None:
        labels = rng.integers(0, k, size=n)
    else:
        labels = rng.choice(k, size=n, p=priors)
    centers = sep * rng.standard_normal((k, n_num))

    if cluster_sds is None:
        cluster_sds = np.ones(k)

    X_num = np.zeros((n, n_num))
    if corr > 0 and n_num > 1:
        cov = np.full((n_num, n_num), corr)
        np.fill_diagonal(cov, 1.0)
        L = np.linalg.cholesky(cov)
    for i in range(n):
        m = labels[i]
        sd = cluster_sds[m]
        if dist == "gaussian":
            z = rng.standard_normal(n_num)
        elif dist == "t2":
            z = rng.standard_t(df=2, size=n_num) / np.sqrt(2)  # roughly unit-scaled heavy tails
        elif dist == "lognormal":
            z = rng.lognormal(mean=0, sigma=0.6, size=n_num) - np.exp(0.6**2/2)  # centred skewed
        else:
            z = rng.standard_normal(n_num)
        if corr > 0 and n_num > 1:
            z = L @ z
        X_num[i] = centers[m] + sd * z

    if contam_frac > 0:
        n_contam = int(n * contam_frac)
        idx = rng.choice(n, size=n_contam, replace=False)
        signs = rng.choice([-1, 1], size=(n_contam, n_num))
        X_num[idx] += signs * contam_scale

    dom = rng.integers(0, cat_levels, size=(k, n_cat))
    X_cat = np.zeros((n, n_cat), dtype=int)
    cat_purity = 0.7
    for j in range(n_cat):
        matches = rng.random(n) < cat_purity
        if cat_imbalance > 0:
            # skew the "non-match" draw toward a single alternate category instead of uniform
            probs = np.full(cat_levels, (1 - cat_imbalance) / (cat_levels - 1))
            alt_cat = (dom[:, j] + 1) % cat_levels
            probs_full = np.tile(probs, (n, 1))
            for i in range(n):
                probs_full[i, alt_cat[labels[i]]] = cat_imbalance
            draws = np.array([rng.choice(cat_levels, p=probs_full[i]) for i in range(n)])
        else:
            draws = rng.integers(0, cat_levels, size=n)
        X_cat[:, j] = np.where(matches, dom[labels, j], draws)
    return X_num, X_cat, labels

def sweep_p(X_num, X_cat, labels, k, p_grid, n_init, max_iter, rng):
    return {p: best_of_accuracy(X_num, X_cat, labels, k, p, n_init, max_iter, rng) for p in p_grid}

if __name__ == "__main__":
    t0 = time.time()
    p_grid = (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
    N, K, NNUM, NCAT, CATLV, SEP = 220, 3, 4, 4, 4, 1.5
    n_init, max_iter = 6, 10
    reps = 3

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
        accs_reps = []
        for rep in range(reps):
            rng = np.random.default_rng(hash((name, rep)) % (2**32))
            X_num, X_cat, labels = base_data(N, K, NNUM, NCAT, CATLV, SEP, rng, **kwargs)
            accs = sweep_p(X_num, X_cat, labels, K, p_grid, n_init, max_iter, rng)
            accs_reps.append(accs)
        mean_accs = {p: float(np.mean([r[p] for r in accs_reps])) for p in p_grid}
        p_star = max(mean_accs, key=mean_accs.get)
        results[name] = dict(mean_accs=mean_accs, p_star=p_star)
        print(f"{name:20s} p*={p_star:.1f}  " + "  ".join(f"p={p}:{mean_accs[p]:.3f}" for p in p_grid))

    print(f"\nelapsed: {time.time()-t0:.1f}s")
    import os
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "p_simulation_factors_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"wrote {out_path}")
