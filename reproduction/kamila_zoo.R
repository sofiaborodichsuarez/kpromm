
# kamila_zoo.R



library(kamila)
library(clue)

set.seed(1)

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

# Identical to load_zoo() in real_data_benchmark.R.
load_zoo <- function(path = "zoo.csv") {
  df <- read.csv(path, header = TRUE, stringsAsFactors = FALSE)
  num_cols <- c("legs")
  cat_cols <- c("hair","feathers","eggs","milk","airborne","aquatic","predator",
                "toothed","backbone","breathes","venomous","fins","tail",
                "domestic","catsize")
  labels <- df$class_type - 1L
  list(X_num = df[num_cols], X_cat = lapply(df[cat_cols], factor),
       labels = labels, K = 7, name = "Zoo")
}

N_INIT <- 10   # matches real_data_benchmark.R
REPS   <- 5    # matches real_data_benchmark.R

d <- load_zoo()
X_num <- as.data.frame(d$X_num)
X_cat <- as.data.frame(lapply(d$X_cat, factor))
labels <- d$labels; K <- d$K
cat(sprintf("Zoo: N=%d, K=%d\n", nrow(X_num), K))

km_acc <- c()
for (rep in seq_len(REPS)) {
  set.seed(2000 + rep)   # matches real_data_benchmark.R's per-rep seed
  cat(sprintf(" rep %d/%d\n", rep, REPS))
  km <- tryCatch(
    kamila(conVar = as.data.frame(scale(X_num)), catFactor = X_cat,
           numClust = K, numInit = N_INIT, verbose = FALSE),
    error = function(e) NULL
  )
  km_acc <- c(km_acc, if (!is.null(km)) assignment_accuracy(labels, km$finalMemb - 1L) else NA)
}

r <- ci95(km_acc)
cat(sprintf("\nkamila on Zoo -> mean=%.4f  lo=%.4f  hi=%.4f  (n=%d valid reps of %d)\n",
            r["mean"], r["lo"], r["hi"], sum(!is.na(km_acc)), REPS))

write.csv(data.frame(dataset = "zoo", kamila_mean = r["mean"], kamila_lo = r["lo"],
                      kamila_hi = r["hi"], n_reps = sum(!is.na(km_acc))),
          "kamila_zoo_results.csv", row.names = FALSE)
cat("Saved kamila_zoo_results.csv\n")
