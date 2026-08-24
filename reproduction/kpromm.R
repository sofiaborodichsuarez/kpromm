
# kpromm.R
#
# Standalone R implementation of k-PROMM / k-PROMMSA mixed
# Minkowski-metric clustering, accompanying:
#   Borodich Suarez, S. Making Mixed-Type Clustering Decision-Ready:
#   k-PROMMSA, Software, and Practical Guidance.
#
# Based on the mathematical formulation in Suarez-Alvarez et al. (2012).


if (!requireNamespace("clue", quietly = TRUE)) {
  stop("Package 'clue' is required (Hungarian algorithm for accuracy). ",
       "Install with install.packages('clue').")
}


# 1. Mixed Minkowski p-metric (raised to the p-th power), eq. (19)-(20)


#' p-th power of the mixed Minkowski distance between one record and one prototype
#'
#' @param x_num numeric vector (numeric attributes of the record)
#' @param x_cat character/factor vector (categorical attributes of the record)
#' @param q_num numeric vector (numeric part of the prototype)
#' @param q_cat character vector (categorical part of the prototype, i.e. modes)
#' @param p Minkowski power (p >= 1)
#' @param alpha per-attribute numeric normalization weights (length = length(x_num)); 1 if unnormalized
#' @param beta per-attribute categorical normalization weights (length = length(x_cat)); 1 if unnormalized
mixed_dist_pow_p <- function(x_num, x_cat, q_num, q_cat, p, alpha, beta) {
  d_num <- if (length(x_num) > 0) sum(alpha * abs(x_num - q_num)^p) else 0
  d_cat <- if (length(x_cat) > 0) sum(beta * (x_cat != q_cat)) else 0  # omega^p == omega for p>=1
  d_num + d_cat
}


# 2. Statistical normalization weights, eq. (23)-(27)
#    alpha_j = 1 / E|X1j - X2j|^p   (numeric),  unbiased pairwise estimator
#    beta_j  = 1 / E[omega(Y1j,Y2j)] (categorical), unbiased pairwise estimator


#' Unbiased pairwise estimate of E|X1j - X2j|^p for one numeric column
#'
#' @param col numeric column
#' @param p Minkowski power
#' @param sample_n if not NULL and \code{length(col) > sample_n}, the O(N^2)
#'   pairwise estimator is computed on a random subsample of this size instead
#'   of the full column. eq.(25)/(3.3) is itself a sample mean over pairs, so
#'   this subsampling is just a smaller, still-unbiased, still-consistent
#'   estimate of the same population quantity, it changes the Monte Carlo
#'   noise of the estimate, not what it targets. Needed because the direct
#'   O(N^2) \code{dist()} call is a ~455M-pair, several-GB computation at
#'   N=30,000+ (e.g. the full-scale Adult dataset), repeated on every call to
#'   \code{kpromm()}; NULL (the default, full-data estimator) is unchanged
#'   from the original implementation and is what every other result in the
#'   paper was computed with.
.est_alpha_num <- function(col, p, sample_n = NULL) {
  n <- length(col)
  if (n < 2) return(1)
  if (!is.null(sample_n) && n > sample_n) {
    col <- col[sample.int(n, sample_n)]
  }
  # 2/(N(N-1)) * sum_{r<s} |x_r - x_s|^p   -- O(N^2), fine for benchmark-size data
  d <- as.numeric(dist(col, method = "minkowski", p = p))^p
  m <- mean(d)  # mean of the N(N-1)/2 unique pairs == the eq.(25) estimator
  if (m == 0) {
    warning("Zero-variance numeric attribute dropped from normalization (alpha undefined).")
    return(0)  # weight 0 -> attribute effectively removed, per the paper's rule
  }
  1 / m
}

#' Unbiased pairwise estimate of E[omega(Y1j,Y2j)] for one categorical column
#'
#' NOTE (fixed): the RSPA 2012 paper gives two estimators for this quantity:
#' eq. (3.15), the "sampling mean" plug-in E-hat = (1/N^2) sum_{r,s=1}^N u(y_r,y_s),
#' which is BIASED and equals the closed form 1 - sum(p_r^2); and eq. (3.16),
#' E-hat = 2/(N(N-1)) sum_{r<s} u(y_r,y_s), which the paper states is unbiased and
#' recommends "for small databases" (p.2642). The numeric weight .est_alpha_num
#' below already uses the analogous unbiased pairwise form (eq. 3.3). This function
#' previously implemented the biased plug-in (3.15) instead -- inconsistent with
#' both the numeric weight and the paper's own recommendation. Corrected to use
#' the unbiased pairwise closed form, 1 - sum(n_r*(n_r-1)) / (N*(N-1)), which is
#' algebraically identical to 2/(N(N-1)) sum_{r<s} u(y_r,y_s) but O(distinct
#' levels) instead of O(N^2). The two estimators coincide asymptotically; they can
#' differ by a few percent at the sample sizes typical of a mixed-data clustering
#' application, so the choice is not cosmetic.
.est_beta_cat <- function(col) {
  n <- length(col)
  if (n < 2) return(1)
  tab <- table(col)
  # unbiased pairwise estimator, eq. (3.16): 1 - sum(n_r(n_r-1)) / (N(N-1))
  m <- 1 - sum(tab * (tab - 1)) / (n * (n - 1))
  if (m == 0) {
    warning("Constant categorical attribute dropped from normalization (beta undefined).")
    return(0)
  }
  1 / m
}

#' Compute normalization weights for a mixed dataset
#'
#' @param X_num data.frame/matrix of numeric attributes (may have 0 columns)
#' @param X_cat data.frame/matrix of categorical attributes (may have 0 columns)
#' @param p Minkowski power
#' @param normalize_num logical: apply statistical normalization to numeric attrs
#' @param normalize_cat logical: apply statistical normalization to categorical attrs
compute_normalization <- function(X_num, X_cat, p, normalize_num = TRUE, normalize_cat = TRUE,
                                   norm_sample_n = NULL) {
  m <- ncol(X_num); l <- ncol(X_cat)
  alpha <- rep(1, m)
  beta  <- rep(1, l)
  if (normalize_num && m > 0) {
    alpha <- vapply(seq_len(m), function(j) .est_alpha_num(X_num[[j]], p, sample_n = norm_sample_n), numeric(1))
  }
  if (normalize_cat && l > 0) {
    # .est_beta_cat is already O(N + distinct levels) via table(), not O(N^2);
    # no subsampling needed here even at large N.
    beta <- vapply(seq_len(l), function(j) .est_beta_cat(X_cat[[j]]), numeric(1))
  }
  list(alpha = alpha, beta = beta)
}


# 3. Prototype recalculation
#    p = 2: closed-form (weighted) mean
#    p != 2: subgradient bisection on each attribute independently (paper Sec.4)
#    categorical: weighted mode (as in standard k-prototypes)


#' Minimizer of sum_i alpha * |x_i - t|^p over t, via subgradient bisection
.minimize_phi_j <- function(x, p, tol = 1e-8, max_iter = 200) {
  n <- length(x)
  if (n == 0) return(NA_real_)
  if (n == 1) return(x[1])
  if (abs(p - 2) < 1e-12) return(mean(x))
  if (abs(p - 1) < 1e-12) return(stats::median(x))

  a0 <- min(x); b0 <- max(x)
  if (a0 == b0) return(a0)

  subgrad <- function(t) {
    d <- x - t
    s <- ifelse(d > 0, -1, ifelse(d < 0, 1, 0))
    u <- sum(p * abs(d)^(p - 1) * ifelse(d == 0, -1, s))
    v <- sum(p * abs(d)^(p - 1) * ifelse(d == 0,  1, s))
    c(u, v)
  }

  for (iter in seq_len(max_iter)) {
    c0 <- (a0 + b0) / 2
    uv <- subgrad(c0)
    if ((b0 - a0) < tol) break
    if (uv[1] > 0 && uv[2] > 0) {
      b0 <- c0
    } else if (uv[1] < 0 && uv[2] < 0) {
      a0 <- c0
    } else {
      break  # 0 is in the subgradient interval -> c0 is the minimizer
    }
  }
  (a0 + b0) / 2
}

#' Weighted mode of a categorical vector (ties broken by first occurrence)
.weighted_mode <- function(x) {
  tab <- table(x)
  names(tab)[which.max(tab)]
}


# 4. Main clustering routine: one run of k-PROMMSA from one random init


.run_once <- function(X_num, X_cat, k, p, alpha, beta, max_iter = 100, seed = NULL,
                       init_Q_num = NULL, init_Q_cat = NULL) {
  if (!is.null(seed)) set.seed(seed)
  N_num <- nrow(X_num)
  N_cat <- nrow(X_cat)
  N <- max(N_num, N_cat)
  if (N < 1L) {
    stop("At least one numeric or categorical clustering attribute is required.")
  }
  if (k > N) {
    stop("k cannot exceed the number of observations.")
  }
  m <- ncol(X_num)
  l <- ncol(X_cat)

  # init_Q_num/init_Q_cat let a caller supply a starting prototype set
  # directly instead of drawing a fresh uniform-random init here. Used by
  # bee_kpromm() to run a short local refinement from a specific candidate
  # solution (a "bee"). Default (NULL) reproduces the original behavior
  # exactly. Every existing, already-verified result in this paper used
  # only the random-init path, unaffected by this addition.
  if (!is.null(init_Q_num) || !is.null(init_Q_cat)) {
    Q_num <- if (m > 0) init_Q_num else matrix(nrow = k, ncol = 0)
    Q_cat <- if (l > 0) init_Q_cat else matrix(nrow = k, ncol = 0)
  } else {
    init_idx <- sample.int(N, k)
    Q_num <- if (m > 0) as.matrix(X_num[init_idx, , drop = FALSE]) else matrix(nrow = k, ncol = 0)
    Q_cat <- if (l > 0) as.matrix(X_cat[init_idx, , drop = FALSE]) else matrix(nrow = k, ncol = 0)
  }

  assign_vec <- rep(NA_integer_, N)
  obj_prev <- Inf

  # Vectorized versions of X_num/X_cat, built once per call (not once per
  # iteration/prototype): identical arithmetic to the original per-record,
  # per-prototype loop over mixed_dist_pow_p(), just computed a whole column
  # at a time instead of via N*k individual R-level function calls. This
  # matters at large N (e.g. Adult, N=30,000+), where the original nested
  # for-loop's per-row data.frame subsetting overhead dominates runtime and
  # makes it impractical to afford enough restarts (n_init) to reliably find
  # the global optimum.
  X_num_mat <- if (m > 0) as.matrix(X_num) else matrix(nrow = N, ncol = 0)
  X_cat_mat <- if (l > 0) as.matrix(X_cat) else matrix(nrow = N, ncol = 0)

  for (it in seq_len(max_iter)) {
    # --- assignment step ---
    dist_mat <- matrix(0, nrow = N, ncol = k)
    for (s in seq_len(k)) {
      d_num <- if (m > 0) {
        diffs <- sweep(X_num_mat, 2, as.numeric(Q_num[s, ]), "-")
        rowSums(sweep(abs(diffs)^p, 2, alpha, "*"))
      } else rep(0, N)
      d_cat <- if (l > 0) {
        mism <- sweep(X_cat_mat, 2, as.character(Q_cat[s, ]), FUN = `!=`)
        storage.mode(mism) <- "double"
        rowSums(sweep(mism, 2, beta, "*"))
      } else rep(0, N)
      dist_mat[, s] <- d_num + d_cat
    }
    new_assign <- apply(dist_mat, 1, which.min)
    obj <- sum(dist_mat[cbind(seq_len(N), new_assign)])

    if (identical(new_assign, assign_vec) || obj >= obj_prev - 1e-10) {
      assign_vec <- new_assign
      obj_prev <- obj
      break
    }
    assign_vec <- new_assign
    obj_prev <- obj

    # --- update step ---
    for (s in seq_len(k)) {
      members <- which(assign_vec == s)
      if (length(members) == 0) next  # keep old prototype for empty cluster
      if (m > 0) {
        for (j in seq_len(m)) {
          Q_num[s, j] <- .minimize_phi_j(X_num[members, j], p)
        }
      }
      if (l > 0) {
        for (j in seq_len(l)) {
          Q_cat[s, j] <- .weighted_mode(X_cat[members, j])
        }
      }
    }
  }

  list(assign = assign_vec, objective = obj_prev, Q_num = Q_num, Q_cat = Q_cat, iterations = it)
}


# 5. Clustering accuracy via the assignment problem, eq. (33)-(35)


#' Best-matching accuracy between a cluster assignment and true labels
clustering_accuracy <- function(assign, true_labels) {
  assign <- as.integer(factor(assign))
  true_labels <- factor(true_labels)
  k <- max(assign)
  klab <- nlevels(true_labels)
  kk <- max(k, klab)
  eff <- matrix(0, kk, kk)
  for (i in seq_len(k)) {
    for (j in seq_len(klab)) {
      eff[i, j] <- sum(assign == i & true_labels == levels(true_labels)[j])
    }
  }
  # Hungarian algorithm maximizes weight by minimizing (max(eff) - eff)
  perm <- clue::solve_LSAP(max(eff) - eff, maximum = FALSE)
  acc <- sum(eff[cbind(seq_len(kk), perm)]) / length(assign)
  acc
}


# 6. Driver: mimics the old ini-file workflow


#' Fit k-PROMM / k-PROMMSA to a mixed dataset
#'
#' @param data data.frame containing all attributes (including id/class columns to be ignored)
#' @param format character vector, one entry per column of `data`: "N" (numeric),
#'   "S" (categorical/string), or "0" (ignore): mirrors the old FORMAT= string.
#' @param k number of clusters (CLUSTER_COUNT)
#' @param p Minkowski power (MINKOVSKI_METRIC_PARAM)
#' @param normalize_num,normalize_cat statistical normalization on/off
#'   (NUMERICAL_NORMALIZATION / CATEGORICAL_NORMALIZATION)
#' @param test_attribute optional column name/index holding the true class label,
#'   used only for accuracy evaluation (TEST_ATTRIBUTE), never for clustering
#' @param n_init number of random restarts (100 in the old tool)
#' @param seed base RNG seed
#' @param norm_sample_n optional subsample size for the numeric normalization
#'   weight estimator (see \code{.est_alpha_num}); NULL (default) reproduces
#'   the original full-data O(N^2) estimator used throughout the paper. Only
#'   needed for datasets where N is large enough (tens of thousands of rows)
#'   that the full pairwise computation is impractically slow.
kpromm <- function(data, format, k, p,
                    normalize_num = TRUE, normalize_cat = TRUE,
                    test_attribute = NULL, n_init = 100, seed = 1,
                    norm_sample_n = NULL) {
  stopifnot(length(format) == ncol(data))
  format <- toupper(format)

  cat_cols <- which(format == "S")
  num_cols <- which(format == "N")

  if (!is.null(test_attribute)) {
    ta <- if (is.character(test_attribute)) which(names(data) == test_attribute) else test_attribute
    true_labels <- data[[ta]]
    num_cols <- setdiff(num_cols, ta)
    cat_cols <- setdiff(cat_cols, ta)
  } else {
    true_labels <- NULL
  }

  X_num <- if (length(num_cols) > 0) as.data.frame(lapply(data[num_cols], as.numeric)) else data.frame()
  X_cat <- if (length(cat_cols) > 0) as.data.frame(lapply(data[cat_cols], as.character), stringsAsFactors = FALSE) else data.frame()

  norm <- compute_normalization(X_num, X_cat, p, normalize_num, normalize_cat,
                                 norm_sample_n = norm_sample_n)

  best <- NULL
  results <- data.frame(attempt = integer(0), objective = numeric(0), accuracy = numeric(0))

  for (attempt in seq_len(n_init)) {
    run <- .run_once(X_num, X_cat, k, p, norm$alpha, norm$beta, seed = seed + attempt)
    acc <- if (!is.null(true_labels)) clustering_accuracy(run$assign, true_labels) else NA_real_
    results <- rbind(results, data.frame(attempt = attempt, objective = run$objective, accuracy = acc))
    if (is.null(best) || run$objective < best$objective) best <- run
  }

  best_acc <- if (!is.null(true_labels)) clustering_accuracy(best$assign, true_labels) else NA_real_

  structure(list(
    call = match.call(),
    p = p, k = k,
    normalize_num = normalize_num, normalize_cat = normalize_cat,
    alpha = norm$alpha, beta = norm$beta,
    best = best,
    best_objective = best$objective,
    best_accuracy = best_acc,
    all_attempts = results,
    num_cols = names(data)[num_cols],
    cat_cols = names(data)[cat_cols]
  ), class = "kpromm")
}

print.kpromm <- function(x, ...) {
  cat(sprintf("k-PROMM%s fit: k=%d, p=%.2f\n",
              if (x$normalize_num || x$normalize_cat) "SA" else "", x$k, x$p))
  cat(sprintf("Best objective: %.4f\n", x$best_objective))
  if (!is.na(x$best_accuracy)) cat(sprintf("Best-objective accuracy: %.4f\n", x$best_accuracy))
  invisible(x)
}

#' Fit k-PROMMSA (statistically normalized k-PROMM)
#'
#' Convenience wrapper around \code{kpromm()} with normalization forced on
#' (\code{normalize_num = TRUE, normalize_cat = TRUE}), matching the
#' statistically normalized variant (k-PROMMSA) used throughout the paper.
#' Equivalent to calling \code{kpromm(..., normalize_num = TRUE,
#' normalize_cat = TRUE)} directly; provided as a named entry point so the
#' normalized and unnormalized variants (k-PROMMSA vs.\ k-PROMM) are called by
#' distinct, self-documenting function names rather than by remembering to
#' set two logical arguments correctly. Not exported from the R package
#' (only \code{kpromm()} and \code{screen_p()} are, matching the manuscript's
#' Software section exactly); kept here for convenience in this standalone
#' reproduction script.
#'
#' @inheritParams kpromm
kprommsa <- function(data, format, k, p,
                      test_attribute = NULL, n_init = 100, seed = 1) {
  kpromm(data, format, k, p,
         normalize_num = TRUE, normalize_cat = TRUE,
         test_attribute = test_attribute, n_init = n_init, seed = seed)
}


# 3. Labels-free ensemble screening procedure for p (Section 4.3 of the paper)


#' Labels-free screening procedure for shortlisting candidate values of p
#'
#' Implements the manuscript's ensemble screening procedure, in two variants
#' selected via \code{pool_restarts}.
#'
#' With \code{pool_restarts = TRUE} (the default): for each p in
#' \code{p_grid}, fit k-PROMMSA with \code{n_init} restarts and pool every
#' restart's partition across the full grid; find the pool member with
#' highest average pairwise similarity (Adjusted Rand Index) to the rest:
#' the "central partition"; score each p by the mean ARI of all of its
#' restart partitions (not just its best) to the central partition. This is
#' the variant actually validated in the manuscript (the 24-scenario ground-
#' truth check below) and used for the empirical application.
#'
#' With \code{pool_restarts = FALSE} (the "basic" variant): for each p,
#' retain only its single best-objective partition, identify the central
#' partition among those best-objective partitions alone, and score each p
#' by its own partition's similarity to that center.
#'
#' IMPORTANT (honesty note, consistent with the paper): this function
#' NARROWS the search over p using only the data and the clustering
#' objective, it does not RELIABLY SELECT the accuracy-optimal p. In the
#' paper's own validation (24 simulated scenarios with known ground truth,
#' using \code{pool_restarts = TRUE}), the top-1 shortlisted p matched the
#' true accuracy-optimal p only 25% of the time (vs. ~14% chance for a
#' 7-point grid), though the expected accuracy cost of trusting the top
#' choice was usually small. This function therefore does NOT auto-select a
#' single p and silently fit only that one; it always returns the full
#' ranking (and the fitted candidate for every p), so the analyst can
#' inspect how much the ranking would have changed their answer, and
#' combine it with domain knowledge (Section~4 recommendations) rather than
#' treat it as a black box.
#'
#' @param data,format,k,test_attribute,n_init,seed as in \code{kpromm()}
#' @param p_grid candidate values of p to screen (default: the paper's 7-point grid)
#' @param normalize_num,normalize_cat passed through to \code{kpromm()} (default: k-PROMMSA)
#' @param pool_restarts logical; if \code{TRUE} (default), pool all \code{n_init}
#'   restart partitions per p (the cross-restart variant); if \code{FALSE},
#'   use only each p's single best-objective partition (the basic variant).
#' @return a list with \code{ranking} (data.frame of p, ari_to_central, rank,
#'   ordered best-first), \code{fits} (named list of the full \code{kpromm}
#'   fit object for every p in \code{p_grid}, so the shortlisted fit(s) can be
#'   used directly without refitting), \code{central_p} (the p whose
#'   partition is most representative of the ensemble), and
#'   \code{pool_restarts} (echoes the argument)
screen_p <- function(data, format, k,
                             p_grid = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0),
                             normalize_num = TRUE, normalize_cat = TRUE,
                             test_attribute = NULL, n_init = 20, seed = 1,
                             pool_restarts = TRUE) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required (Adjusted Rand Index). ",
         "Install with install.packages('mclust').")
  }
  fits <- lapply(p_grid, function(p) {
    kpromm(data, format, k, p, normalize_num = normalize_num, normalize_cat = normalize_cat,
           test_attribute = test_attribute, n_init = n_init, seed = seed)
  })
  names(fits) <- as.character(p_grid)

  if (!pool_restarts) {
    # Basic variant: one best-objective partition per p.
    assigns <- lapply(fits, function(f) f$best$assign)
    np <- length(p_grid)
    sim <- matrix(0, np, np)
    for (i in seq_len(np)) for (j in seq_len(np)) {
      sim[i, j] <- if (i == j) 1 else mclust::adjustedRandIndex(assigns[[i]], assigns[[j]])
    }
    avg_sim <- rowMeans(sim)  # avg similarity to all others (diagonal = 1 by construction)
    central_idx <- which.max(avg_sim)

    ranking <- data.frame(
      p = p_grid,
      ari_to_central = sim[, central_idx],
      stringsAsFactors = FALSE
    )
    ranking <- ranking[order(-ranking$ari_to_central), ]
    ranking$rank <- seq_len(nrow(ranking))

    return(list(ranking = ranking, fits = fits, central_p = p_grid[central_idx],
                pool_restarts = FALSE))
  }

  # Cross-restart variant (default; the procedure actually validated in the
  # manuscript): pool every restart's partition across the full p grid,
  # identify the pool member with highest average similarity to the whole
  # pool as the central partition, and score each p by the mean ARI of all
  # of its n_init restart partitions (not just its best) to that partition.
  format_u <- toupper(format)
  cat_cols <- which(format_u == "S"); num_cols <- which(format_u == "N")
  if (!is.null(test_attribute)) {
    ta <- if (is.character(test_attribute)) which(names(data) == test_attribute) else test_attribute
    num_cols <- setdiff(num_cols, ta); cat_cols <- setdiff(cat_cols, ta)
  }
  X_num <- if (length(num_cols) > 0) as.data.frame(lapply(data[num_cols], as.numeric)) else data.frame()
  X_cat <- if (length(cat_cols) > 0) as.data.frame(lapply(data[cat_cols], as.character), stringsAsFactors = FALSE) else data.frame()

  pool_labels <- list(); pool_p <- c()
  for (p in p_grid) {
    norm <- compute_normalization(X_num, X_cat, p, normalize_num, normalize_cat)
    for (attempt in seq_len(n_init)) {
      run <- .run_once(X_num, X_cat, k, p, norm$alpha, norm$beta, seed = seed + attempt)
      pool_labels[[length(pool_labels) + 1]] <- run$assign
      pool_p <- c(pool_p, p)
    }
  }
  M <- length(pool_labels)
  sim <- matrix(0, M, M)
  for (i in seq_len(M)) for (j in seq_len(M)) {
    sim[i, j] <- if (i == j) 1 else mclust::adjustedRandIndex(pool_labels[[i]], pool_labels[[j]])
  }
  central_idx <- which.max(rowMeans(sim))
  central_partition <- pool_labels[[central_idx]]

  ens_score <- sapply(p_grid, function(p) {
    idx <- which(pool_p == p)
    mean(sapply(idx, function(i) mclust::adjustedRandIndex(pool_labels[[i]], central_partition)))
  })
  ranking <- data.frame(p = p_grid, ari_to_central = ens_score, stringsAsFactors = FALSE)
  ranking <- ranking[order(-ranking$ari_to_central), ]
  ranking$rank <- seq_len(nrow(ranking))

  list(ranking = ranking, fits = fits, central_p = pool_p[central_idx], pool_restarts = TRUE)
}


# NOTE: an earlier working version also included an exploratory bee-colony
# global-search variant (bee_kpromm()) developed during investigations on the
# Adult dataset. Adult and this variant are not part of the final manuscript
# or its reproducible analyses and are therefore not included in this release.

