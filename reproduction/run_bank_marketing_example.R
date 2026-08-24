
# run_bank_marketing_example.R
#
# Reproduces the R-side analysis for the Bank Marketing empirical
# application reported in the manuscript. The preprocessing and clustering
# configuration follow the specifications described in the paper.
#
# Requires: install.packages(c("clue", "clustMixType", "kamila", "mclust"))
# ("mclust" is only used for one Adjusted Rand Index number in Section 5; the
# rest of the script runs without it.)
# Put kpromm.R and bank-additional.csv in the same folder as this script,
# or edit the paths below.


source("kpromm.R")

# 1. Load & preprocess (mirrors the Python preprocessing exactly)
d <- read.csv("bank-additional.csv", sep = ";", stringsAsFactors = FALSE)

# duration is dropped: the dataset's own documentation states it leaks the
# outcome (duration = 0 implies y = "no") and is not known before a call is
# made, so it is not a legitimate clustering input for a decision-support use.
num_cols <- c("age", "campaign", "pdays", "previous",
              "emp.var.rate", "cons.price.idx", "cons.conf.idx",
              "euribor3m", "nr.employed")
cat_cols <- c("job", "marital", "education", "default", "housing", "loan",
              "contact", "month", "day_of_week", "poutcome")

# kpromm() expects one data.frame + a FORMAT vector; y and duration are
# marked "0" (ignored) and pdays/previous etc. are left as raw numeric values,
# including the pdays = 999 ("never previously contacted") sentinel -- this
# sentinel is deliberately NOT recoded, because part of the point of the
# application is to show what statistical normalization and low p do to a
# raw, uncorrected sentinel-coded variable, which is exactly the kind of
# artifact real analysts encounter and often do not think to fix by hand.
keep <- c(num_cols, cat_cols, "y")
d2 <- d[, keep]
fmt <- c(rep("N", length(num_cols)), rep("S", length(cat_cols)), "0")

K <- 4
P_GRID <- c(1.0, 1.5, 2.0, 2.5, 3.0, 4.0)

results <- list()
for (normalize in c(FALSE, TRUE)) {
  for (p in P_GRID) {
    fit <- kpromm(d2, fmt, k = K, p = p,
                  normalize_num = normalize, normalize_cat = normalize,
                  test_attribute = "y",  # used only to report accuracy-style
                                          # cross-tab diagnostics below, never
                                          # fed into the clustering itself
                  n_init = 20, seed = 1)
    # Twenty restarts are used because the raw pdays sentinel can induce
    # competing near-tied local optima; the best-objective restart is retained.
    key <- paste0("p", p, "_norm", normalize)
    results[[key]] <- fit
    sizes <- table(fit$best$assign)
    cat(sprintf("p=%.1f normalized=%s objective=%.1f sizes=%s\n",
                p, normalize, fit$best_objective, paste(sizes, collapse = ",")))
  }
}

# 2. Compare cluster sizes/objective against the Python run
# Expected pattern (from the Python port, seed-equivalent but not seed-identical
# since R's RNG stream differs from numpy's, compare PATTERNS, not exact
# cluster membership):
#   - normalize=FALSE: cluster sizes barely change across p >= 2, and one small
#     cluster (~160 records, ~3.9% of the sample) should closely track clients
#     with pdays != 999 (previously contacted) i.e. the raw pdays sentinel
#     dominates the clustering almost by itself.
#   - normalize=TRUE: cluster sizes should be much more balanced across all p,
#     and no single cluster should closely track pdays != 999 alone.
#
# 3. Profile a chosen solution against the held-out y
profile_solution <- function(fit, label) {
  a <- fit$best$assign
  cat("\n===", label, "===\n")
  print(data.frame(
    cluster = sort(unique(a)),
    n = as.integer(table(a)),
    subscribe_rate = tapply(d2$y == "yes", a, mean),
    prev_contacted_rate = tapply(d2$pdays != 999, a, mean)
  ))
}

profile_solution(results[["p2_normFALSE"]], "Unnormalized, p=2 (naive baseline)")
profile_solution(results[["p2_normTRUE"]],  "Normalized, p=2")

# 4. Save for comparison against the Python/manuscript results
saveRDS(results, "bank_marketing_kpromm_results.rds")
cat("\nSaved: bank_marketing_kpromm_results.rds\n")

# 5. Competitor comparison: clustMixType::kproto and kamila
# These are the two field-standard, actively maintained R packages for
# mixed-type clustering (the pairing used by Jimeno, Roy & Tortora 2021's
# benchmark study). This is what belongs in the manuscript's Section 6.3
# corroboration, run this block (install.packages(c("clustMixType","kamila"))
# first) and use kp_acc / km_acc / kp_ari / km_ari below directly; do not use
# the earlier Python/kmodes stand-in numbers in the final submission.
library(clustMixType)
library(kamila)
library(clue)

d3 <- d2
d3[cat_cols] <- lapply(d3[cat_cols], as.factor)

kp <- kproto(d3[, c(num_cols, cat_cols)], k = K, nstart = 8)
cat("\nclustMixType::kproto sizes:", table(kp$cluster), "\n")
cat("clustMixType::kproto lambda (auto gamma):", kp$lambda, "\n")
kp_acc <- clustering_accuracy(kp$cluster, d2$y)
kp_prevcontact_rate <- tapply(d2$pdays != 999, kp$cluster, mean)
cat("clustMixType::kproto accuracy vs. y:", round(kp_acc, 4), "\n")
cat("clustMixType::kproto prev.-contacted rate by cluster:\n"); print(round(kp_prevcontact_rate, 3))
# Expected, per the Python/kmodes run: sizes close to (243, 160, 2771, 945) in
# some order, and a ~160-record cluster with prev.-contacted rate close to
# 1.0 i.e. clustMixType's own default weighting should reproduce the same
# pdays-sentinel-driven pathology as unnormalized k-PROMM, not fix it.
kp_vs_unnorm_ari <- {
  a <- as.integer(factor(kp$cluster)); b <- as.integer(factor(results[["p2_normFALSE"]]$best$assign))
  # simple ARI via mclust if available, else report the confusion table only
  if (requireNamespace("mclust", quietly = TRUE)) mclust::adjustedRandIndex(a, b) else NA
}
cat("Adjusted Rand Index, kproto vs. unnormalized k-PROMM (p=2):", kp_vs_unnorm_ari,
    if (is.na(kp_vs_unnorm_ari)) "(install 'mclust' for this number)" else "", "\n")

conVars <- scale(d2[, num_cols])  # kamila expects continuous vars pre-scaled
catVarsFac <- as.data.frame(lapply(d2[, cat_cols], as.factor))
km <- kamila(conVar = as.data.frame(conVars), catFactor = catVarsFac,
             numClust = K, numInit = 8)
cat("\nkamila sizes:", table(km$finalMemb), "\n")
km_acc <- clustering_accuracy(km$finalMemb, d2$y)
km_prevcontact_rate <- tapply(d2$pdays != 999, km$finalMemb, mean)
cat("kamila accuracy vs. y:", round(km_acc, 4), "\n")
cat("kamila prev.-contacted rate by cluster:\n"); print(round(km_prevcontact_rate, 3))
# kamila pre-scales continuous variables by design (unlike kproto), so it may
# or may not reproduce the same failure mode, report whichever it actually
# does, do not assume either direction.

cat("\n=== Section 6.3 corroboration summary (drop straight into the manuscript) ===\n")
cat(sprintf("kproto max prev.-contacted rate in any cluster: %.3f\n", max(kp_prevcontact_rate)))
cat(sprintf("kamila max prev.-contacted rate in any cluster: %.3f\n", max(km_prevcontact_rate)))
cat("If either is close to 1.0 with a small cluster size, that package reproduces\n",
    "the pdays-artifact pathology; report the ARI against unnormalized k-PROMM (above)\n",
    "as the headline number, exactly as computed for kproto.\n")
