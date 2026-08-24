
# kamila_simulation_ci.R


library(kamila)
library(clue)

set.seed(1)

# Accuracy metric (identical to kpromm.R's clustering_accuracy(),
# reproduced here standalone so this script needs no other project files)
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
  perm <- clue::solve_LSAP(max(eff) - eff, maximum = FALSE)
  sum(eff[cbind(seq_len(kk), perm)]) / length(assign)
}
assignment_accuracy <- function(true_labels, pred_labels) {
  clustering_accuracy(pred_labels, true_labels)
}
ci95 <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(c(mean = mean(x), lo = mean(x), hi = mean(x)))
  m <- mean(x); se <- sd(x) / sqrt(length(x))
  c(mean = m, lo = m - 1.96 * se, hi = m + 1.96 * se)
}

# Data-generating process, identical to simulation_comparison.R
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

# The 13 factor-sweep conditions, identical to simulation_comparison.R 
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
n_init <- 6      # matches simulation_comparison.R's kamila numInit
reps <- 10       # matches simulation_comparison.R's reps

results <- list()

for (cond_name in names(conditions)) {
  km_acc <- c()

  for (rep in seq_len(reps)) {
    # Identical seeding rule to simulation_comparison.R, so this draws the
    # same synthetic datasets that kamila's already-published point estimates
    # in Table \ref{tab:speed-ci} were computed on.
    set.seed(1000 * rep + which(names(conditions) == cond_name))
    d <- do.call(generate_data, c(list(n = N, k = K, n_num = NNUM, n_cat = NCAT,
                                        cat_levels = CATLV, sep = SEP), conditions[[cond_name]]))
    X_num <- d$X_num; X_cat <- d$X_cat; labels <- d$labels
    df_cat <- as.data.frame(lapply(seq_len(NCAT), function(j) factor(X_cat[, j])))

    km <- tryCatch(
      kamila(conVar = as.data.frame(scale(X_num)), catFactor = df_cat,
             numClust = K, numInit = n_init, verbose = FALSE),
      error = function(e) NULL
    )
    km_acc <- c(km_acc, if (!is.null(km)) assignment_accuracy(labels, km$finalMemb - 1L) else NA)

    cat(sprintf("  %s rep %d/%d done\n", cond_name, rep, reps))
  }

  results[[cond_name]] <- ci95(km_acc)
  r <- results[[cond_name]]
  cat(sprintf("%-22s kamila mean=%.4f  lo=%.4f  hi=%.4f  (halfwidth=%.4f)\n",
              cond_name, r["mean"], r["lo"], r["hi"], (r["hi"] - r["lo"]) / 2))
}

results_table <- do.call(rbind, lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(condition = nm, kamila_mean = r["mean"], kamila_lo = r["lo"], kamila_hi = r["hi"])
}))
write.csv(results_table, "kamila_simulation_ci_results.csv", row.names = FALSE)
cat("\nSaved kamila_simulation_ci_results.csv\n")
