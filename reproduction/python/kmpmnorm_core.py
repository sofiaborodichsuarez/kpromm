"""
Python implementation of k-PROMM / k-PROMMSA used for the
simulation study accompanying the manuscript.

Implements the mixed Minkowski objective and the unbiased
pairwise statistical-normalization estimators used in the paper.
"""
import numpy as np
from scipy.optimize import linear_sum_assignment


def mixed_dist_pow_p(X_num, X_cat, Q_num, Q_cat, p, alpha, beta):
    """Vectorized: distance^p from every record in X to prototype Q (1 x attrs)."""
    d = np.zeros(X_num.shape[0])
    if X_num.shape[1] > 0:
        d += (alpha[None, :] * np.abs(X_num - Q_num[None, :]) ** p).sum(axis=1)
    if X_cat.shape[1] > 0:
        d += (beta[None, :] * (X_cat != Q_cat[None, :]).astype(float)).sum(axis=1)
    return d


def est_alpha_num(col, p):
    """Unbiased pairwise estimator of E|X1-X2|^p, eq. (3.3)/(3.16)-analogue."""
    n = len(col)
    if n < 2:
        return 1.0
    # O(n^2) exact pairwise (fine for the normalization subsample size used)
    diffs = np.abs(col[:, None] - col[None, :]) ** p
    iu = np.triu_indices(n, k=1)
    m = diffs[iu].mean()
    if m <= 1e-12:
        return 0.0
    return 1.0 / m


def est_beta_cat(col):
    """Unbiased pairwise estimator, eq. (3.16): 1 - sum(n_r(n_r-1))/(N(N-1))."""
    n = len(col)
    if n < 2:
        return 1.0
    _, counts = np.unique(col, return_counts=True)
    m = 1.0 - (counts * (counts - 1)).sum() / (n * (n - 1))
    if m <= 1e-12:
        return 0.0
    return 1.0 / m


def compute_normalization(X_num, X_cat, p, normalize_num=True, normalize_cat=True,
                           norm_subsample=None, rng=None):
    """
    norm_subsample: if set, estimate alpha/beta from a random subsample of this
    size (unbiased estimator remains valid on a subsample; only variance changes).
    Needed here because the O(n^2) pairwise estimate is too slow/memory-heavy at
    n ~ 4-40k in a pure-Python/NumPy implementation; the R code, run natively, does
    not need this restriction for datasets in the low thousands.
    """
    m = X_num.shape[1]
    l = X_cat.shape[1]
    alpha = np.ones(m)
    beta = np.ones(l)

    idx = np.arange(X_num.shape[0] if m > 0 else X_cat.shape[0])
    if norm_subsample is not None and len(idx) > norm_subsample:
        rng = rng or np.random.default_rng(0)
        idx = rng.choice(idx, size=norm_subsample, replace=False)

    if normalize_num and m > 0:
        for j in range(m):
            alpha[j] = est_alpha_num(X_num[idx, j], p)
    if normalize_cat and l > 0:
        for j in range(l):
            beta[j] = est_beta_cat(X_cat[idx, j])
    return alpha, beta


def minimize_phi_j(x, p, tol=1e-8, max_iter=200):
    n = len(x)
    if n == 0:
        return np.nan
    if n == 1:
        return x[0]
    if abs(p - 2) < 1e-12:
        return x.mean()
    if abs(p - 1) < 1e-12:
        return np.median(x)
    a0, b0 = x.min(), x.max()
    if a0 == b0:
        return a0
    for _ in range(max_iter):
        c0 = (a0 + b0) / 2
        d = x - c0
        s = np.where(d > 0, -1.0, np.where(d < 0, 1.0, 0.0))
        u = np.sum(p * np.abs(d) ** (p - 1) * np.where(d == 0, -1.0, s))
        v = np.sum(p * np.abs(d) ** (p - 1) * np.where(d == 0, 1.0, s))
        if (b0 - a0) < tol:
            break
        if u > 0 and v > 0:
            b0 = c0
        elif u < 0 and v < 0:
            a0 = c0
        else:
            break
    return (a0 + b0) / 2


def weighted_mode(x):
    vals, counts = np.unique(x, return_counts=True)
    return vals[np.argmax(counts)]


def run_once(X_num, X_cat, k, p, alpha, beta, max_iter, rng):
    N = X_num.shape[0] if X_num.shape[1] > 0 else X_cat.shape[0]
    m, l = X_num.shape[1], X_cat.shape[1]
    init_idx = rng.choice(N, size=k, replace=False)
    Q_num = X_num[init_idx].copy() if m > 0 else np.zeros((k, 0))
    Q_cat = X_cat[init_idx].copy() if l > 0 else np.zeros((k, 0), dtype=object)

    assign = np.full(N, -1)
    obj_prev = np.inf
    it_used = 0
    for it in range(max_iter):
        it_used = it + 1
        D = np.zeros((N, k))
        for s in range(k):
            D[:, s] = mixed_dist_pow_p(X_num, X_cat, Q_num[s], Q_cat[s], p, alpha, beta)
        new_assign = D.argmin(axis=1)
        obj = D[np.arange(N), new_assign].sum()

        if np.array_equal(new_assign, assign) or obj >= obj_prev - 1e-10:
            assign = new_assign
            obj_prev = obj
            break
        assign = new_assign
        obj_prev = obj

        for s in range(k):
            members = assign == s
            if members.sum() == 0:
                continue
            if m > 0:
                for j in range(m):
                    Q_num[s, j] = minimize_phi_j(X_num[members, j], p)
            if l > 0:
                for j in range(l):
                    Q_cat[s, j] = weighted_mode(X_cat[members, j])
    return dict(assign=assign, objective=obj_prev, Q_num=Q_num, Q_cat=Q_cat, iterations=it_used)


def kpromm(X_num, X_cat, k, p, normalize_num=True, normalize_cat=True,
           n_init=8, max_iter=25, seed=1, norm_subsample=None):
    rng_norm = np.random.default_rng(seed)
    alpha, beta = compute_normalization(X_num, X_cat, p, normalize_num, normalize_cat,
                                         norm_subsample=norm_subsample, rng=rng_norm)
    best = None
    for attempt in range(n_init):
        rng = np.random.default_rng(seed * 10_000 + attempt)
        run = run_once(X_num, X_cat, k, p, alpha, beta, max_iter, rng)
        if best is None or run["objective"] < best["objective"]:
            best = run
    best["alpha"] = alpha
    best["beta"] = beta
    best["p"] = p
    best["k"] = k
    return best


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
