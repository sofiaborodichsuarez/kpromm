# =============================================================================
# real_data_benchmark.R
#
# This final benchmark replaces an earlier development version that used CMC
# as the third dataset. CMC was dropped because no method (of eight tried)
# cleared its majority-class baseline by more than a marginal amount.
# Zoo is used instead:
#
#   Zoo (UCI id 111): N=101 animals, 16 predictor attributes (15 boolean:
#   hair, feathers, eggs, milk, airborne, aquatic, predator, toothed,
#   backbone, breathes, venomous, fins, tail, domestic, catsize -- plus 1
#   numeric, "legs", values in {0,2,4,5,6,8}), K=7 true taxonomic classes
#   (mammal/bird/reptile/fish/amphibian/insect/invertebrate). No missing
#   values. Majority-class baseline is 0.406 (41/101 mammals) -- much less
#   severe than CMC's near-uninformative regime, and on the Python side of
#   this benchmark (6 comparators run natively) every method clears it
#   comfortably (range 0.64-0.92), confirming real, recoverable label-aligned
#   structure. Source: archive.ics.uci.edu/dataset/111/zoo, file zoo.data
#   (or the equivalent zoo.csv used on the Python side of this benchmark --
#   same 101 rows, same columns, just CSV with a header row instead of
#   UCI's whitespace-free .data format).
#
# MODELING NOTE: "legs" is treated as numeric rather than categorical because
# it is a count-valued predictor. All 15 boolean attributes are treated as
# categorical (2-level factors), consistent with the treatment of binary
# attributes elsewhere in the benchmark.
#
# Methods run:
# k-prototypes (clustMixType::kproto), kamila::kamila, PAM on Gower's
# dissimilarity, average-linkage HC on Gower's dissimilarity, FAMD+k-means,
# FAMD+Student-t mixture (teigen), clustMD (attempted, expected to fail as
# on the other three datasets -- kept in the script for a complete record,
# not expected to be reportable), and k-PROMMSA at p=2 and best-of-grid
# (kpromm.R). The exploratory bee_kpromm() variant is not included because
# it is not part of the final manuscript or reproduction pipeline.
#
# Requires: install.packages(c("clustMixType", "kamila", "clue", "cluster",
#           "FactoMineR", "teigen", "clustMD"))
# Put this file, kpromm.R, and the four raw data files in the same folder:
#   processed.cleveland.data -- archive.ics.uci.edu/dataset/45/heart+disease
#   crx.data                 -- archive.ics.uci.edu/dataset/27/credit+approval
#   dermatology.data         -- archive.ics.uci.edu/dataset/33/dermatology
#   zoo.csv                  -- archive.ics.uci.edu/dataset/111/zoo (header row:
#     animal_name,hair,feathers,eggs,milk,airborne,aquatic,predator,toothed,
#     backbone,breathes,venomous,fins,legs,tail,domestic,catsize,class_type)
# dermatology.data is downloaded automatically if absent; the other listed
# data files should be placed alongside this script.
#   Rscript real_data_benchmark.R
# =============================================================================

source("kpromm.R")
library(clustMixType)
library(kamila)
library(clue)
library(cluster)
library(FactoMineR)
library(teigen)
library(clustMD)

GOWER_SUBSAMPLE_N <- 3000   # never triggers for these 4 datasets; kept for
                            # consistency with earlier scripts
N_INIT   <- 10
MAX_ITER <- 15
REPS     <- 5
P_GRID   <- c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)

assignment_accuracy <- function(true_labels, pred_labels) {
  clustering_accuracy(pred_labels, true_labels)
}
ci95 <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return(c(mean = mean(x), lo = mean(x), hi = mean(x)))
  m <- mean(x); se <- sd(x) / sqrt(length(x))
  c(mean = m, lo = m - 1.96 * se, hi = m + 1.96 * se)
}

# ---- Dataset loaders --------------------------------------------------------

# 1. Cleveland Heart Disease. Download: archive.ics.uci.edu/dataset/45/heart+disease
#    (file: processed.cleveland.data). Save alongside this script.
load_heart <- function(path = "processed.cleveland.data") {
  cols <- c("age","sex","cp","trestbps","chol","fbs","restecg","thalach",
            "exang","oldpeak","slope","ca","thal","num")
  df <- read.csv(path, header = FALSE, col.names = cols, na.strings = "?")
  df <- na.omit(df)
  num_cols <- c("age","trestbps","chol","thalach","oldpeak")
  cat_cols <- c("sex","cp","fbs","restecg","exang","slope","ca","thal")
  labels <- as.integer(df$num > 0)
  list(X_num = df[num_cols], X_cat = lapply(df[cat_cols], factor),
       labels = labels, K = 2, name = "Heart Disease (Cleveland)")
}

# 2. Credit Approval (Australian), crx.data. Download:
#    archive.ics.uci.edu/dataset/27/credit+approval (file: crx.data).
load_credit_approval <- function(path = "crx.data") {
  cols <- paste0("A", 1:16)
  df <- read.csv(path, header = FALSE, col.names = cols, na.strings = "?",
                  stringsAsFactors = FALSE)
  num_cols <- c("A2","A3","A8","A11","A14","A15")
  cat_cols <- c("A1","A4","A5","A6","A7","A9","A10","A12","A13")
  for (c in num_cols) {
    df[[c]] <- as.numeric(df[[c]])
    df[[c]][is.na(df[[c]])] <- median(df[[c]], na.rm = TRUE)
  }
  for (c in cat_cols) {
    mode_val <- names(sort(table(df[[c]]), decreasing = TRUE))[1]
    df[[c]][is.na(df[[c]])] <- mode_val
    df[[c]] <- factor(df[[c]])
  }
  labels <- as.integer(df$A16 == "+")
  list(X_num = df[num_cols], X_cat = df[cat_cols],
       labels = labels, K = 2, name = "Credit Approval (Australian)")
}

# 3. Zoo. Download: archive.ics.uci.edu/dataset/111/zoo (or use zoo.csv from
#    the Python side of this benchmark -- same data, header row included).
#    17 predictor columns (animal_name dropped) + class_type (1-7).
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

# 4. Dermatology. Source: UCI Machine Learning Repository, dataset 33
#    (https://archive.ics.uci.edu/dataset/33/dermatology), distributed under
#    CC BY 4.0. Ilter, N. & Guvenir, H. (1998). Dermatology [Dataset]. UCI
#    Machine Learning Repository. https://doi.org/10.24432/C5FK5P
#    34 feature columns (1-33 clinical/histopathological, 34 = Age) + class
#    (35, values 1-6). Age has 8 missing values coded "?" -- dropped, leaving
#    N=358. All severity-score attributes treated as categorical rather than
#    numeric; this follows the preprocessing specification used in the final
#    manuscript.
#
#    If dermatology.data is not present in this directory, it is downloaded
#    automatically from the UCI repository on first use (requires internet
#    access; ~25 KB). To fetch it manually instead, download the dataset zip
#    from the URL above and extract dermatology.data into this folder.
.download_dermatology <- function(path) {
  url <- "https://archive.ics.uci.edu/static/public/33/dermatology.zip"
  cat("dermatology.data not found -- downloading from UCI (", url, ")...\n", sep = "")
  tmp_zip <- tempfile(fileext = ".zip")
  ok <- tryCatch({
    download.file(url, tmp_zip, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)
  if (!ok || !file.exists(tmp_zip)) {
    stop("Automatic download of dermatology.data failed (no internet access, or ",
         "the UCI URL has changed). Download the dataset manually from ",
         "https://archive.ics.uci.edu/dataset/33/dermatology and place ",
         "dermatology.data in this directory, then re-run.")
  }
  extracted <- unzip(tmp_zip, files = "dermatology.data", exdir = dirname(normalizePath(path, mustWork = FALSE)))
  if (length(extracted) == 0 || !file.exists(path)) {
    stop("Downloaded dermatology.zip but could not find dermatology.data inside it. ",
         "Download and extract it manually from https://archive.ics.uci.edu/dataset/33/dermatology.")
  }
  cat("Saved ", path, "\n", sep = "")
  invisible(path)
}

load_dermatology <- function(path = "dermatology.data") {
  if (!file.exists(path)) .download_dermatology(path)
  cols <- c("erythema","scaling","definite_borders","itching","koebner_phenomenon",
            "polygonal_papules","follicular_papules","oral_mucosal_involvement",
            "knee_elbow_involvement","scalp_involvement","family_history",
            "melanin_incontinence","eosinophils_infiltrate","pnl_infiltrate",
            "fibrosis_papillary_dermis","exocytosis","acanthosis","hyperkeratosis",
            "parakeratosis","clubbing_rete_ridges","elongation_rete_ridges",
            "thinning_suprapapillary_epidermis","spongiform_pustule",
            "munro_microabcess","focal_hypergranulosis","disappearance_granular_layer",
            "vacuolisation_damage_basal_layer","spongiosis","saw_tooth_appearance_retes",
            "follicular_horn_plug","perifollicular_parakeratosis",
            "inflammatory_monoluclear_infiltrate","band_like_infiltrate",
            "age","class")
  df <- read.csv(path, header = FALSE, col.names = cols, na.strings = "?")
  df <- na.omit(df)
  num_cols <- c("age")
  cat_cols <- setdiff(cols, c("age", "class"))
  labels <- df$class - 1L
  list(X_num = df[num_cols], X_cat = lapply(df[cat_cols], factor),
       labels = labels, K = 6, name = "Dermatology (erythemato-squamous disease)")
}

DATASETS <- list(
  heart = load_heart,
  credit_approval = load_credit_approval,
  zoo = load_zoo,
  dermatology = load_dermatology
)

# ---- Runner: all eight methods on one dataset ------------------------------

run_one_dataset <- function(dataset_name, loader) {
  d <- loader()
  X_num <- as.data.frame(d$X_num)
  X_cat <- as.data.frame(lapply(d$X_cat, factor))
  labels <- d$labels; K <- d$K
  N <- nrow(X_num)
  cat(sprintf("\n=== %s: N=%d, K=%d ===\n", d$name, N, K))

  df_mixed <- cbind(X_num, X_cat)
  use_subsample <- N > GOWER_SUBSAMPLE_N

  acc <- list(kproto = c(), kamila = c(), pam = c(), hc = c(),
              famdkm = c(), famdt = c(), clustmd = c(), p2 = c(), bestp = c())

  for (rep in seq_len(REPS)) {
    set.seed(2000 + rep)
    cat(sprintf(" rep %d/%d\n", rep, REPS))

    kp <- tryCatch(kproto(df_mixed, k = K, nstart = N_INIT, verbose = FALSE), error = function(e) NULL)
    acc$kproto <- c(acc$kproto, if (!is.null(kp)) assignment_accuracy(labels, kp$cluster - 1L) else NA)

    km <- tryCatch(
      kamila(conVar = as.data.frame(scale(X_num)), catFactor = X_cat,
             numClust = K, numInit = N_INIT, verbose = FALSE),
      error = function(e) NULL)
    acc$kamila <- c(acc$kamila, if (!is.null(km)) assignment_accuracy(labels, km$finalMemb - 1L) else NA)

    if (use_subsample) {
      idx <- sample.int(N, GOWER_SUBSAMPLE_N)
    } else {
      idx <- seq_len(N)
    }
    df_sub <- df_mixed[idx, , drop = FALSE]
    labels_sub <- labels[idx]
    gower_d <- daisy(df_sub, metric = "gower")

    pam_fit <- tryCatch(pam(gower_d, k = K, nstart = N_INIT), error = function(e) NULL)
    acc$pam <- c(acc$pam, if (!is.null(pam_fit)) assignment_accuracy(labels_sub, pam_fit$clustering - 1L) else NA)

    hc_fit <- tryCatch(cutree(hclust(gower_d, method = "average"), k = K), error = function(e) NULL)
    acc$hc <- c(acc$hc, if (!is.null(hc_fit)) assignment_accuracy(labels_sub, hc_fit - 1L) else NA)

    famd_fit <- tryCatch(FAMD(df_mixed, ncp = min(5, ncol(df_mixed) - 1), graph = FALSE), error = function(e) NULL)
    famd_scores <- if (!is.null(famd_fit)) famd_fit$ind$coord else NULL

    famdkm_fit <- tryCatch(
      if (!is.null(famd_scores)) kmeans(famd_scores, centers = K, nstart = N_INIT) else NULL,
      error = function(e) NULL)
    acc$famdkm <- c(acc$famdkm, if (!is.null(famdkm_fit)) assignment_accuracy(labels, famdkm_fit$cluster - 1L) else NA)

    famdt_fit <- tryCatch(
      if (!is.null(famd_scores)) teigen(famd_scores, Gs = K, verbose = FALSE) else NULL,
      error = function(e) NULL)
    acc$famdt <- c(acc$famdt, if (!is.null(famdt_fit)) assignment_accuracy(labels, famdt_fit$classification - 1L) else NA)

    cmd_data <- as.matrix(cbind(
      X_num[idx, , drop = FALSE],
      as.data.frame(lapply(X_cat[idx, , drop = FALSE], function(col) as.integer(factor(col))))
    ))
    n_num <- ncol(X_num)
    cmd_fit <- tryCatch(
      clustMD(X = cmd_data, G = K, CnsIndx = n_num, OrdIndx = n_num,
              Nnorms = 100, MaxIter = 100, model = "EII",
              store.params = FALSE, scale = TRUE, startCL = "kmeans",
              autoStop = TRUE, ma.band = 30, stop.tol = 0.01),
      error = function(e) NULL)
    acc$clustmd <- c(acc$clustmd,
                      if (!is.null(cmd_fit)) assignment_accuracy(labels_sub, cmd_fit$cl - 1L) else NA)

    fmt <- c(rep("N", ncol(X_num)), rep("S", ncol(X_cat)), "0")
    df_all <- cbind(X_num, X_cat, y = 0)
    fit2 <- tryCatch(
      kpromm(df_all, fmt, k = K, p = 2.0, normalize_num = TRUE, normalize_cat = TRUE,
             n_init = N_INIT, seed = rep),
      error = function(e) NULL)
    acc$p2 <- c(acc$p2, if (!is.null(fit2)) assignment_accuracy(labels, fit2$best$assign - 1L) else NA)

    accs_p <- sapply(P_GRID, function(p) {
      f <- tryCatch(
        kpromm(df_all, fmt, k = K, p = p, normalize_num = TRUE, normalize_cat = TRUE,
               n_init = N_INIT, seed = rep),
        error = function(e) NULL)
      if (is.null(f)) NA else assignment_accuracy(labels, f$best$assign - 1L)
    })
    acc$bestp <- c(acc$bestp, max(accs_p, na.rm = TRUE))
  }

  result <- lapply(acc, ci95)
  cat(sprintf("%s -> kproto=%.3f kamila=%.3f pam=%.3f hc=%.3f famdkm=%.3f famdt=%.3f clustmd=%.3f p2=%.3f bestp=%.3f\n",
              dataset_name, result$kproto["mean"], result$kamila["mean"], result$pam["mean"],
              result$hc["mean"], result$famdkm["mean"], result$famdt["mean"],
              result$clustmd["mean"], result$p2["mean"], result$bestp["mean"]))
  list(N = N, K = K, name = d$name, subsampled_gower_clustmd = use_subsample,
       gower_subsample_n = if (use_subsample) GOWER_SUBSAMPLE_N else N, acc = result)
}

# ---- Run all four datasets and save ----------------------------------------

all_results <- list()
for (nm in names(DATASETS)) {
  all_results[[nm]] <- tryCatch(
    run_one_dataset(nm, DATASETS[[nm]]),
    error = function(e) { cat(sprintf("SKIPPED %s: %s\n", nm, conditionMessage(e))); NULL }
  )
}

saveRDS(all_results, "real_data_benchmark_results.rds")
results_table <- do.call(rbind, lapply(names(all_results), function(nm) {
  r <- all_results[[nm]]
  if (is.null(r)) return(NULL)
  data.frame(dataset = nm, N = r$N, K = r$K,
             t(sapply(r$acc, function(a) a["mean"])))
}))
write.csv(results_table, "real_data_benchmark_results.csv", row.names = FALSE)
cat("\nSaved real_data_benchmark_results.csv / .rds\n")
cat("This is the authoritative source for the manuscript's Table",
    "\"Real-data validation\" -- 4 rows: heart, credit_approval, zoo,\n",
    "dermatology.\n")
