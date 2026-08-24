# R/kpromm.R
#
# Package source for kpromm: k-PROMM / k-PROMMSA mixed Minkowski-metric
# clustering, as described in:
#   Suarez-Alvarez, M.M. (2010). Design and analysis of clustering algorithms
#   for numerical, categorical and mixed data. PhD thesis, Cardiff University.
#   Suarez-Alvarez, M.M., Pham, D-T., Prostov, M.Y. & Prostov, Y.I. (2012).
#   Statistical approach to normalization of feature vectors and clustering
#   of mixed data sets. Proceedings of the Royal Society A, 468(2145),
#   2630-2651.
#   Borodich Suarez, S. Making Mixed-Type Clustering Decision-Ready:
#   k-PROMMSA, Software, and Practical Guidance (manuscript).
#
# Public API corresponding to the manuscript's two documented functions:
#   kpromm(): fit k-PROMM/k-PROMMSA at one value of p
#   screen_p(): labels-free ensemble screening over a grid of p
# Everything else in this file is internal (not exported).
#
# NOTE ON PROVENANCE: this reimplements the logic of a legacy Borland C++
# tool (ZCLU_BC.exe, ini-driven: FORMAT, WEIGHTS, CLUSTER_COUNT,
# MINKOVSKI_METRIC_PARAM, NUMERICAL_NORMALIZATION, CATEGORICAL_NORMALIZATION,
# TEST_ATTRIBUTE) purely from the mathematical specification in the papers
# above, not from the executable itself.
#
# NOTE ON SCOPE: an earlier working version also included an exploratory
# bee-colony global-search variant (bee_kpromm()) developed during
# investigations on the Adult dataset. Adult and this variant are not part
# of the manuscript or its reproducible analyses and are therefore
# not included in this release.


# 1. Mixed Minkowski p-metric (raised to the p-th power)

#' p-th power of the mixed Minkowski distance between one record and one prototype
#'
#' Internal. Not exported.
#'
#' @param x_num numeric vector (numeric attributes of the record)
#' @param x_cat character/factor vector (categorical attributes of the record)
#' @param q_num numeric vector (numeric part of the prototype)
#' @param q_cat character vector (categorical part of the prototype, i.e. modes)
#' @param p Minkowski power (p >= 1)
#' @param alpha per-attribute numeric normalization weights (length = length(x_num)); 1 if unnormalized
#' @param beta per-attribute categorical normalization weights (length = length(x_cat)); 1 if unnormalized
#' @keywords internal
#' @noRd
mixed_dist_pow_p <- function(x_num, x_cat, q_num, q_cat, p, alpha, beta) {
  d_num <- if (length(x_num) > 0) sum(alpha * abs(x_num - q_num)^p) else 0
  d_cat <- if (length(x_cat) > 0) sum(beta * (x_cat != q_cat)) else 0  # omega^p == omega for p>=1
  d_num + d_cat
}


# 2. Statistical normalization weights
#    alpha_j = 1 / E|X1j - X2j|^p   (numeric),  unbiased pairwise estimator
#    beta_j  = 1 / E[omega(Y1j,Y2j)] (categorical), unbiased pairwise estimator


#' Unbiased pairwise estimate of E|X1j - X2j|^p for one numeric column
#'
#' Internal. Not exported.
#'
#' @param col numeric column
#' @param p Minkowski power
#' @param sample_n if not NULL and \code{length(col) > sample_n}, the O(N^2)
#'   pairwise estimator is computed on a random subsample of this size instead
#'   of the full column (still unbiased, still consistent -- changes the
#'   Monte Carlo noise of the estimate, not what it targets). NULL (default)
#'   reproduces the full-data estimator used throughout the manuscript.
#' @keywords internal
#' @noRd
.est_alpha_num <- function(col, p, sample_n = NULL) {
  n <- length(col)
  if (n < 2) return(1)
  if (!is.null(sample_n) && n > sample_n) {
    col <- col[sample.int(n, sample_n)]
  }
  # 2/(N(N-1)) * sum_{r<s} |x_r - x_s|^p   is O(N^2), fine for benchmark-size data
  d <- as.numeric(stats::dist(col, method = "minkowski", p = p))^p
  m <- mean(d)  # mean of the N(N-1)/2 unique pairs
  if (m == 0) {
    warning("Zero-variance numeric attribute dropped from normalization (alpha undefined).")
    return(0)  # weight 0 -> attribute effectively removed, per the paper's rule
  }
  1 / m
}

#' Unbiased pairwise estimate of E[omega(Y1j,Y2j)] for one categorical column
#'
#' Internal. Not exported. Uses the unbiased pairwise closed form,
#' 1 - sum(n_r*(n_r-1)) / (N*(N-1)), algebraically identical to the paper's
#' eq. (3.16) but O(distinct levels) instead of O(N^2).
#'
#' @param col a factor/character column
#' @keywords internal
#' @noRd
.est_beta_cat <- function(col) {
  n <- length(col)
  if (n < 2) return(1)
  tab <- table(col)
  m <- 1 - sum(tab * (tab - 1)) / (n * (n - 1))
  if (m == 0) {
    warning("Constant categorical attribute dropped from normalization (beta undefined).")
    return(0)
  }
  1 / m
}

#' Compute normalization weights for a mixed dataset
#'
#' Internal. Not exported.
#'
#' @param X_num data.frame/matrix of numeric attributes (may have 0 columns)
#' @param X_cat data.frame/matrix of categorical attributes (may have 0 columns)
#' @param p Minkowski power
#' @param normalize_num logical: apply statistical normalization to numeric attrs
#' @param normalize_cat logical: apply statistical normalization to categorical attrs
#' @param norm_sample_n optional subsample size, see \code{.est_alpha_num}
#' @keywords internal
#' @noRd
compute_normalization <- function(X_num, X_cat, p, normalize_num = TRUE, normalize_cat = TRUE,
                                   norm_sample_n = NULL) {
  m <- ncol(X_num); l <- ncol(X_cat)
  alpha <- rep(1, m)
  beta  <- rep(1, l)
  if (normalize_num && m > 0) {
    alpha <- vapply(seq_len(m), function(j) .est_alpha_num(X_num[[j]], p, sample_n = norm_sample_n), numeric(1))
  }
  if (normalize_cat && l > 0) {
    beta <- vapply(seq_len(l), function(j) .est_beta_cat(X_cat[[j]]), numeric(1))
  }
  list(alpha = alpha, beta = beta)
}


# 3. Prototype recalculation
#    p = 2: closed-form (weighted) mean
#    p != 2: subgradient bisection on each attribute independently
#    categorical: weighted mode (as in standard k-prototypes)


#' Minimizer of sum_i alpha * |x_i - t|^p over t, via subgradient bisection
#' @keywords internal
#' @noRd
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
#' @keywords internal
#' @noRd
.weighted_mode <- function(x) {
  tab <- table(x)
  names(tab)[which.max(tab)]
}


# 4. Main clustering routine: one run of k-PROMM(SA) from one random init


#' @keywords internal
#' @noRd
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

  X_num_mat <- if (m > 0) as.matrix(X_num) else matrix(nrow = N, ncol = 0)
  X_cat_mat <- if (l > 0) as.matrix(X_cat) else matrix(nrow = N, ncol = 0)

  for (it in seq_len(max_iter)) {
    # assignment step
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

    # update step
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


# 5. Clustering accuracy via the linear assignment problem


#' Best-matching accuracy between a cluster assignment and true labels
#'
#' Internal. Not exported (kept as an internal helper used by \code{kpromm()}
#' when \code{test_attribute} is supplied; not part of the manuscript's
#' documented two-function public API).
#'
#' @param assign integer/factor cluster assignment
#' @param true_labels ground-truth labels, same length as \code{assign}
#' @keywords internal
#' @noRd
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
  sum(eff[cbind(seq_len(kk), perm)]) / length(assign)
}


# 6. Public API


#' Fit k-PROMM / k-PROMMSA to a mixed-type dataset
#'
#' Fits the k-PROMM (unnormalized) or k-PROMMSA (statistically normalized)
#' clustering algorithm to a data frame containing numeric and/or categorical
#' attributes, using the mixed Minkowski p-metric.
#'
#' @param data data.frame containing all attributes (including any id/class
#'   columns to be ignored via \code{format = "0"} or excluded via
#'   \code{test_attribute})
#' @param format character vector, one entry per column of \code{data}:
#'   \code{"N"} (numeric), \code{"S"} (categorical/string), or \code{"0"}
#'   (ignore)
#' @param k number of clusters
#' @param p Minkowski exponent (p >= 1)
#' @param normalize_num,normalize_cat logical, independently toggle
#'   statistical normalization of the numeric and categorical halves of the
#'   distance. Both default to \code{TRUE} (k-PROMMSA); setting both to
#'   \code{FALSE} recovers unnormalized k-PROMM.
#' @param test_attribute optional column name/index holding the true class
#'   label, used only to report accuracy, never to fit the clustering
#' @param n_init number of random restarts; the best-objective partition
#'   across restarts is returned
#' @param seed base RNG seed; restart \code{i} uses \code{seed + i}
#' @param norm_sample_n optional subsample size for the numeric
#'   normalization-weight estimator; \code{NULL} (default) uses the full-data
#'   estimator. Only useful for datasets large enough (tens of thousands of
#'   rows) that the full O(N^2) pairwise computation is impractically slow.
#'
#' @return An object of class \code{"kpromm"}: a list with the fitted
#'   normalization weights (\code{alpha}, \code{beta}), the best-objective
#'   partition (\code{best}), its objective value (\code{best_objective}) and
#'   accuracy against \code{test_attribute} if supplied (\code{best_accuracy}),
#'   and the full per-restart objective/accuracy trace (\code{all_attempts}).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 60
#' df <- data.frame(
#'   x1 = c(rnorm(n / 2, 0), rnorm(n / 2, 5)),
#'   x2 = c(rnorm(n / 2, 0), rnorm(n / 2, 5)),
#'   cat1 = sample(c("a", "b"), n, replace = TRUE),
#'   true_label = rep(1:2, each = n / 2)
#' )
#' fit <- kpromm(df, format = c("N", "N", "S", "0"), k = 2, p = 2,
#'               test_attribute = "true_label", n_init = 10)
#' print(fit)
#' }
#'
#' @export
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

#' Print method for \code{kpromm} fit objects
#' @param x an object of class \code{"kpromm"}
#' @param ... unused, for S3 consistency
#' @export
print.kpromm <- function(x, ...) {
  cat(sprintf("k-PROMM%s fit: k=%d, p=%.2f\n",
              if (x$normalize_num || x$normalize_cat) "SA" else "", x$k, x$p))
  cat(sprintf("Best objective: %.4f\n", x$best_objective))
  if (!is.na(x$best_accuracy)) cat(sprintf("Best-objective accuracy: %.4f\n", x$best_accuracy))
  invisible(x)
}

#' Fit k-PROMMSA (convenience wrapper, not part of the manuscript's documented API)
#'
#' Convenience wrapper around \code{kpromm()} with normalization forced on
#' (\code{normalize_num = TRUE, normalize_cat = TRUE}). Not exported: the
#' manuscript's Software section documents exactly two public functions,
#' \code{kpromm()} and \code{screen_p()}; this wrapper is kept
#' internal so the package's exported surface matches the manuscript exactly.
#' Equivalent to \code{kpromm(..., normalize_num = TRUE, normalize_cat = TRUE)}.
#'
#' @inheritParams kpromm
#' @keywords internal
#' @noRd
kprommsa <- function(data, format, k, p,
                      test_attribute = NULL, n_init = 100, seed = 1) {
  kpromm(data, format, k, p,
         normalize_num = TRUE, normalize_cat = TRUE,
         test_attribute = test_attribute, n_init = n_init, seed = seed)
}


# 7. Labels-free ensemble screening procedure for p


#' Labels-free screening procedure for shortlisting candidate values of p
#'
#' Implements the manuscript's ensemble screening procedure in the two
#' variants defined there ("A labels-free screening procedure for p"),
#' selected via \code{pool_restarts}.
#'
#' With \code{pool_restarts = TRUE} (the default): for each \code{p} in
#' \code{p_grid}, fit k-PROMMSA with \code{n_init} restarts and pool
#' \emph{every} restart's partition (not only the best-objective one)
#' across the full grid; find the pool member with highest average
#' pairwise similarity (Adjusted Rand Index) to the entire pool -- the
#' "central partition"; score each \code{p} by the \emph{mean} ARI of all
#' of its own restart partitions to the central partition. This is the
#' variant actually used for the manuscript's ground-truth validation and
#' its real-data application, and is now the default so that this function
#' reproduces the paper's reported screening evidence out of the box.
#'
#' With \code{pool_restarts = FALSE} (the "basic" variant): for each
#' \code{p}, keep only the single best-objective partition across
#' \code{n_init} restarts; find the central partition among these
#' \code{p}-indexed partitions only; score each \code{p} by its own
#' best-objective partition's ARI to that central partition. This is
#' cheaper (its pairwise-ARI cost scales with the number of \code{p}
#' values rather than with \code{n_init} times the number of \code{p}
#' values) but was not the variant validated in the manuscript.
#'
#' \strong{Honesty note (consistent with the manuscript):} this function
#' narrows the search over \code{p} using only the data and the clustering
#' objective; it does not reliably select the accuracy-optimal \code{p}. In
#' the manuscript's own validation (24 simulated scenarios with known ground
#' truth, using \code{pool_restarts = TRUE}), the top-1 shortlisted \code{p}
#' matched the true accuracy-optimal \code{p} only 25\% of the time (versus
#' ~14\% chance for a 7-point grid), though the expected accuracy cost of
#' trusting the top choice was usually small. This function therefore does
#' not auto-select a single \code{p} and silently fit only that one; it
#' always returns the full ranking (and the fitted candidate for every
#' \code{p}), so the analyst can inspect how much the ranking would have
#' changed their answer and combine it with domain knowledge rather than
#' treat it as a black box.
#'
#' @param data,format,k,test_attribute,n_init,seed as in \code{kpromm()}
#' @param p_grid candidate values of p to screen (default: the manuscript's
#'   7-point grid)
#' @param normalize_num,normalize_cat passed through to \code{kpromm()}
#'   (default: k-PROMMSA)
#' @param pool_restarts logical; if \code{TRUE} (default), pool all
#'   \code{n_init} restart partitions per \code{p} into the similarity
#'   comparison (the variant validated in the manuscript); if \code{FALSE},
#'   use only each \code{p}'s best-objective partition (cheaper, but not
#'   the validated variant). Cost with pooling is quadratic in
#'   \code{n_init * length(p_grid)}; for large \code{n_init} on a wide
#'   grid, allow extra runtime.
#'
#' @return a list with \code{ranking} (data.frame of \code{p},
#'   \code{ari_to_central}, \code{rank}, ordered best-first), \code{fits}
#'   (named list of the full \code{kpromm} fit object for every \code{p} in
#'   \code{p_grid}, so the shortlisted fit(s) can be used directly without
#'   refitting), \code{central_p} (the \code{p} whose partition is most
#'   representative of the ensemble), and \code{pool_restarts} (which
#'   variant was used)
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 60
#' df <- data.frame(
#'   x1 = c(rnorm(n / 2, 0), rnorm(n / 2, 5)),
#'   x2 = c(rnorm(n / 2, 0), rnorm(n / 2, 5)),
#'   cat1 = sample(c("a", "b"), n, replace = TRUE)
#' )
#' scr <- screen_p(df, format = c("N", "N", "S"), k = 2,
#'                         p_grid = c(1, 2, 3), n_init = 5)
#' scr$ranking
#' }
#'
#' @export
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
