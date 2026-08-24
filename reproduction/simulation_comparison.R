
# simulation_comparison.R
#
# Runs the paper's Section 5/6 factor-sweep simulation NATIVELY in R, comparing
# k-PROMMSA (kpromm.R) against seven established mixed-type clustering methods:
# clustMixType::kproto, kamila::kamila, PAM on Gower's dissimilarity (PAM-
# Gower, via cluster::pam + cluster::daisy(metric="gower")), average-
# linkage hierarchical clustering on the same Gower dissimilarity (HC-Gower,
# via stats::hclust), and two tandem-analysis comparators: FactoMineR::FAMD
# followed by stats::kmeans (FAMD+k-means), and FactoMineR::FAMD followed by
# a multivariate Student-t mixture via teigen::teigen (FAMD+t-mixture). kproto
# and kamila are the two field-standard, actively maintained R packages for
# mixed-type clustering (per Jimeno, Roy & Tortora 2021, "Clustering
# Mixed-Type Data: A Benchmark Study on KAMILA and K-Prototypes", who use
# exactly this pairing); PAM-Gower and HC-Gower are added because two
# comparators are not enough to support a general claim of outperformance,
# and are, alongside k-prototypes, among the most commonly used mixed-type
# clustering baselines in the applied literature; the two FAMD tandem
# comparators match the benchmark designs of Jimeno et al. (2021) and Costa,
# Papatsouma & Markos (2024, DIBmix). clustMD is deliberately NOT included
# here (see real_data_benchmark.R), its Monte Carlo EM fit was not stable
# enough across this many synthetic conditions/replicates within a reasonable
# runtime; it is run instead on the smaller, fixed real datasets.
#
# What this produces, and where it goes in the manuscript:
#   - `results_table`  -> Table \ref{tab:speed-ci} / \ref{tab:speed-time} in
#                         Section 5 (Simulation comparison): accuracy (mean,
#                         95% CI) and mean wall-clock time per method, across
#                         all 13 conditions.
#   - Drop the printed/saved CSV numbers straight into the LaTeX table; no
#     further processing needed.
#
# Requires: install.packages(c("clustMixType", "kamila", "clue", "cluster",
#           "MASS", "FactoMineR", "teigen"))
# Put this file and kpromm.R in the same folder, then: Rscript simulation_comparison.R
# Runtime: with reps = 10 and n_init = 6, expect ~30-45 minutes depending on
# machine (the FAMD+t-mixture comparator is the slowest addition); reduce
# `reps` for a quick check first.


source("kpromm.R")
library(clustMixType)
library(kamila)
library(clue)
library(cluster)
library(FactoMineR)
library(teigen)

set.seed(1)

# 1. Data-generating process (ports p_simulation_factors.py::base_data)
# Mirrors the Python simulation exactly: K=3 Gaussian-mixture clusters (or
# heavy-tailed / skewed / contaminated variants) over 4 numeric attributes,
# 4 categorical attributes (4 levels each, 0.7 within-cluster purity, with an
# optional imbalance skew), N=220 by default.
generate_data <- function(n, k, n_num, n_cat, cat_levels, sep,
                           priors = NULL, corr = 0, cluster_sds = NULL,
                           dist = "gaussian", contam_frac = 0, contam_scale = 10,
                           cat_imbalance = 0) {
  if (is.null(priors)) {
    labels <- sample.int(k, n, replace = TRUE) - 1L
  } else {
    labels <- sample.int(k, n, replace = TRUE, prob = priors) - 1L
  }
  centers <- matrix(rnorm(k * n_num, sd = sep), nrow = k)
  if (is.null(cluster_sds)) cluster_sds <- rep(1, k)

  X_num <- matrix(0, n, n_num)
  L <- NULL
  if (corr > 0 && n_num > 1) {
    cov <- matrix(corr, n_num, n_num); diag(cov) <- 1
    L <- t(chol(cov))
  }
  for (i in seq_len(n)) {
    m <- labels[i] + 1L
    sd <- cluster_sds[m]
    if (dist == "gaussian") {
      z <- rnorm(n_num)
    } else if (dist == "t2") {
      z <- rt(n_num, df = 2) / sqrt(2)
    } else if (dist == "lognormal") {
      z <- rlnorm(n_num, meanlog = 0, sdlog = 0.6) - exp(0.6^2 / 2)
    } else {
      z <- rnorm(n_num)
    }
    if (!is.null(L)) z <- as.vector(L %*% z)
    X_num[i, ] <- centers[m, ] + sd * z
  }

  if (contam_frac > 0) {
    n_contam <- floor(n * contam_frac)
    idx <- sample.int(n, n_contam)
    signs <- matrix(sample(c(-1, 1), n_contam * n_num, replace = TRUE), nrow = n_contam)
    X_num[idx, ] <- X_num[idx, ] + signs * contam_scale
  }

  dom <- matrix(sample.int(cat_levels, k * n_cat, replace = TRUE) - 1L, nrow = k)
  X_cat <- matrix(0L, n, n_cat)
  cat_purity <- 0.7
  for (j in seq_len(n_cat)) {
    matches <- runif(n) < cat_purity
    if (cat_imbalance > 0) {
      alt_cat <- (dom[, j] + 1L) %% cat_levels
      draws <- integer(n)
      base_p <- (1 - cat_imbalance) / (cat_levels - 1)
      for (i in seq_len(n)) {
        probs <- rep(base_p, cat_levels)
        probs[alt_cat[labels[i] + 1L] + 1L] <- cat_imbalance
        draws[i] <- sample.int(cat_levels, 1, prob = probs) - 1L
      }
    } else {
      draws <- sample.int(cat_levels, n, replace = TRUE) - 1L
    }
    X_cat[, j] <- ifelse(matches, dom[labels + 1L, j], draws)
  }
  list(X_num = X_num, X_cat = X_cat, labels = labels)
}

# 2. Assignment-based accuracy (eq. AccNumb3 in the manuscript)
# Reuse kpromm.R's own clustering_accuracy(assign, true_labels) rather than
# reimplementing it: it re-factors both label vectors internally, so it is
# robust to any indexing/labelling convention (0- or 1-based, or arbitrary
# category codes) returned by kproto/kamila/kpromm alike.
assignment_accuracy <- function(true_labels, pred_labels, k) {
  clustering_accuracy(pred_labels, true_labels)
}

ci95 <- function(x) {
  m <- mean(x); se <- sd(x) / sqrt(length(x))
  c(mean = m, lo = m - 1.96 * se, hi = m + 1.96 * se)
}

# 3. The 13 factor-sweep conditions (identical to the Python version)
conditions <- list(
  baseline             = list(),
  contam_mild          = list(contam_frac = 0.03, contam_scale = 8),
  contam_severe        = list(contam_frac = 0.10, contam_scale = 8),
  heavytail_t2         = list(dist = "t2"),
  skew_lognormal       = list(dist = "lognormal"),
  unequal_var_mild     = list(cluster_sds = c(0.7, 1.0, 1.6)),
  unequal_var_severe   = list(cluster_sds = c(0.4, 1.0, 3.0)),
  correlation_mild     = list(corr = 0.4),
  correlation_high     = list(corr = 0.8),
  cat_imbalance_mild   = list(cat_imbalance = 0.6),
  cat_imbalance_severe = list(cat_imbalance = 0.9),
  unequal_size_mild    = list(priors = c(0.6, 0.25, 0.15)),
  unequal_size_severe  = list(priors = c(0.8, 0.15, 0.05))
)

N <- 220; K <- 3; NNUM <- 4; NCAT <- 4; CATLV <- 4; SEP <- 1.5
P_GRID <- c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)
n_init <- 6      # restarts per method per rep (matches Python; raise for final numbers)
max_iter <- 10
reps <- 10       # replicates per condition, for the 95% CI

results <- list()

for (cond_name in names(conditions)) {
  kp_acc <- c(); kp_time <- c()
  km_acc <- c(); km_time <- c()
  pam_acc <- c(); pam_time <- c()
  hc_acc <- c(); hc_time <- c()
  famdkm_acc <- c(); famdkm_time <- c()
  famdt_acc <- c(); famdt_time <- c()
  p2_acc <- c(); p2_time <- c()
  bp_acc <- c(); bp_time <- c()

  for (rep in seq_len(reps)) {
    set.seed(1000 * rep + which(names(conditions) == cond_name))
    d <- do.call(generate_data, c(list(n = N, k = K, n_num = NNUM, n_cat = NCAT,
                                        cat_levels = CATLV, sep = SEP), conditions[[cond_name]]))
    X_num <- d$X_num; X_cat <- d$X_cat; labels <- d$labels

    df_num <- as.data.frame(X_num)
    df_cat <- as.data.frame(lapply(seq_len(NCAT), function(j) factor(X_cat[, j])))
    names(df_cat) <- paste0("cat", seq_len(NCAT))
    df_mixed <- cbind(df_num, df_cat)

    # k-prototypes (clustMixType)
    t0 <- Sys.time()
    kp <- tryCatch(
      kproto(df_mixed, k = K, nstart = n_init, verbose = FALSE),
      error = function(e) NULL
    )
    kp_time <- c(kp_time, as.numeric(Sys.time() - t0, units = "secs"))
    kp_acc <- c(kp_acc, if (!is.null(kp)) assignment_accuracy(labels, kp$cluster - 1L, K) else NA)

    # kamila
    # NOTE: argument names below are for kamila package versions circa
    # 0.1.x (conVar, catFactor, numClust, numInit, maxIter, verbose). If your
    # installed version errors with "unused argument", run ?kamila and adjust,
    # the package has changed its interface across versions.
    t0 <- Sys.time()
    km <- tryCatch(
      kamila(conVar = as.data.frame(scale(X_num)), catFactor = df_cat,
             numClust = K, numInit = n_init, verbose = FALSE),
      error = function(e) NULL
    )
    km_time <- c(km_time, as.numeric(Sys.time() - t0, units = "secs"))
    km_acc <- c(km_acc, if (!is.null(km)) assignment_accuracy(labels, km$finalMemb - 1L, K) else NA)

    # Gower dissimilarity, computed once and shared by PAM-Gower and HC-Gower
    gower_d <- daisy(df_mixed, metric = "gower")

    # PAM-Gower (Kaufman & Rousseeuw 1990), via cluster::pam
    t0 <- Sys.time()
    pam_fit <- tryCatch(
      pam(gower_d, k = K, nstart = n_init),
      error = function(e) NULL
    )
    pam_time <- c(pam_time, as.numeric(Sys.time() - t0, units = "secs"))
    pam_acc <- c(pam_acc, if (!is.null(pam_fit)) assignment_accuracy(labels, pam_fit$clustering - 1L, K) else NA)

    # HC-Gower: average-linkage hierarchical clustering on Gower dissimilarity,
    # cut at K clusters. Deterministic (no restarts needed).
    t0 <- Sys.time()
    hc_fit <- tryCatch({
      hc <- hclust(gower_d, method = "average")
      cutree(hc, k = K)
    }, error = function(e) NULL)
    hc_time <- c(hc_time, as.numeric(Sys.time() - t0, units = "secs"))
    hc_acc <- c(hc_acc, if (!is.null(hc_fit)) assignment_accuracy(labels, hc_fit - 1L, K) else NA)

    # FAMD scores, shared by the two tandem comparators below
    famd_fit <- tryCatch(
      FAMD(df_mixed, ncp = min(5, ncol(df_mixed) - 1), graph = FALSE),
      error = function(e) NULL
    )
    famd_scores <- if (!is.null(famd_fit)) famd_fit$ind$coord else NULL

    # FAMD + k-means (tandem analysis; Jimeno et al. 2021; Costa et al. 2024)
    t0 <- Sys.time()
    famdkm_fit <- tryCatch(
      if (!is.null(famd_scores)) kmeans(famd_scores, centers = K, nstart = n_init) else NULL,
      error = function(e) NULL
    )
    famdkm_time <- c(famdkm_time, as.numeric(Sys.time() - t0, units = "secs"))
    famdkm_acc <- c(famdkm_acc, if (!is.null(famdkm_fit)) assignment_accuracy(labels, famdkm_fit$cluster - 1L, K) else NA)

    # FAMD + Student-t mixture (Peel & McLachlan 2000), via teigen
    t0 <- Sys.time()
    famdt_fit <- tryCatch(
      if (!is.null(famd_scores)) teigen(famd_scores, Gs = K, verbose = FALSE) else NULL,
      error = function(e) NULL
    )
    famdt_time <- c(famdt_time, as.numeric(Sys.time() - t0, units = "secs"))
    famdt_acc <- c(famdt_acc, if (!is.null(famdt_fit)) assignment_accuracy(labels, famdt_fit$classification - 1L, K) else NA)

    # k-PROMMSA at p = 2
    fmt <- c(rep("N", NNUM), rep("S", NCAT), "0")
    df_all <- cbind(df_num, df_cat, y = 0)
    t0 <- Sys.time()
    fit2 <- kpromm(df_all, fmt, k = K, p = 2.0, normalize_num = TRUE,
                    normalize_cat = TRUE, n_init = n_init, seed = rep)
    p2_time <- c(p2_time, as.numeric(Sys.time() - t0, units = "secs"))
    p2_acc <- c(p2_acc, assignment_accuracy(labels, fit2$best$assign - 1L, K))

    # k-PROMMSA, best of the p-grid
    t0 <- Sys.time()
    accs_p <- sapply(P_GRID, function(p) {
      f <- kpromm(df_all, fmt, k = K, p = p, normalize_num = TRUE,
                  normalize_cat = TRUE, n_init = n_init, seed = rep)
      assignment_accuracy(labels, f$best$assign - 1L, K)
    })
    bp_time <- c(bp_time, as.numeric(Sys.time() - t0, units = "secs"))
    bp_acc <- c(bp_acc, max(accs_p))

    cat(sprintf("  %s rep %d/%d done\n", cond_name, rep, reps))
  }

  results[[cond_name]] <- list(
    kproto_acc = ci95(kp_acc), kproto_time = mean(kp_time),
    kamila_acc = ci95(km_acc), kamila_time = mean(km_time),
    pam_acc = ci95(pam_acc), pam_time = mean(pam_time),
    hc_acc = ci95(hc_acc), hc_time = mean(hc_time),
    famdkm_acc = ci95(famdkm_acc[!is.na(famdkm_acc)]), famdkm_time = mean(famdkm_time),
    famdt_acc = ci95(famdt_acc[!is.na(famdt_acc)]), famdt_time = mean(famdt_time),
    p2_acc = ci95(p2_acc), p2_time = mean(p2_time),
    bestp_acc = ci95(bp_acc), bestp_time = mean(bp_time)
  )
  r <- results[[cond_name]]
  cat(sprintf("%s: kproto=%.3f kamila=%.3f pam=%.3f hc=%.3f famdkm=%.3f famdt=%.3f p2=%.3f bestp=%.3f\n",
              cond_name, r$kproto_acc["mean"], r$kamila_acc["mean"],
              r$pam_acc["mean"], r$hc_acc["mean"],
              r$famdkm_acc["mean"], r$famdt_acc["mean"],
              r$p2_acc["mean"], r$bestp_acc["mean"]))
}

# 4. Build the results table (drop straight into Table tab:speed-ci)
results_table <- do.call(rbind, lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(
    condition = nm,
    kproto_acc = r$kproto_acc["mean"], kproto_lo = r$kproto_acc["lo"], kproto_hi = r$kproto_acc["hi"],
    kproto_ms = r$kproto_time * 1000,
    kamila_acc = r$kamila_acc["mean"], kamila_lo = r$kamila_acc["lo"], kamila_hi = r$kamila_acc["hi"],
    kamila_ms = r$kamila_time * 1000,
    pam_acc = r$pam_acc["mean"], pam_lo = r$pam_acc["lo"], pam_hi = r$pam_acc["hi"],
    pam_ms = r$pam_time * 1000,
    hc_acc = r$hc_acc["mean"], hc_lo = r$hc_acc["lo"], hc_hi = r$hc_acc["hi"],
    hc_ms = r$hc_time * 1000,
    p2_acc = r$p2_acc["mean"], p2_lo = r$p2_acc["lo"], p2_hi = r$p2_acc["hi"],
    p2_ms = r$p2_time * 1000,
    bestp_acc = r$bestp_acc["mean"], bestp_lo = r$bestp_acc["lo"], bestp_hi = r$bestp_acc["hi"],
    bestp_ms = r$bestp_time * 1000
  )
}))
write.csv(results_table, "simulation_comparison_results.csv", row.names = FALSE)
saveRDS(results, "simulation_comparison_results.rds")
cat("\nSaved: simulation_comparison_results.csv / .rds\n")
cat("Use these columns directly for the manuscript's Table (Section 5): condition, ",
    "kproto_acc [kproto_lo, kproto_hi] kproto_ms, kamila_acc [...] kamila_ms, ",
    "pam_acc [...] pam_ms, hc_acc [...] hc_ms, ",
    "p2_acc [...] p2_ms, bestp_acc [...] bestp_ms.\n",
    "Expect HC-Gower (hc_acc) to outperform k-PROMMSA under 'unequal_size_severe' -> ",
    "this matches the Python-side finding (Section 5.3) and should be reported, not adjusted away.\n")
