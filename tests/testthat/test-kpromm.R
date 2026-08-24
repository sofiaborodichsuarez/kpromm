# Basic sanity tests for kpromm()/screen_p(). These are NOT a
# reproduction of the manuscript's simulation results (that requires the
# full simulation_comparison.R / real_data_benchmark.R in the accompanying
# reproduction/ folder, which take much longer to run), they check that
# the package's exported functions behave correctly on small, fast,
# well-understood synthetic inputs: correct output structure, sane values,
# and expected qualitative behaviour (e.g. recovering well-separated
# clusters with high accuracy).
#
# These tests are run by devtools::test() and R CMD check.
# For version 1.0.0, all 25 tests passed under R 4.6.1 on
# macOS Tahoe 26.1 (Apple Silicon).

make_toy_data <- function(seed = 1, n_per_cluster = 40) {
  set.seed(seed)
  n <- n_per_cluster
  data.frame(
    x1 = c(stats::rnorm(n, mean = 0, sd = 1), stats::rnorm(n, mean = 8, sd = 1)),
    x2 = c(stats::rnorm(n, mean = 0, sd = 1), stats::rnorm(n, mean = 8, sd = 1)),
    cat1 = c(sample(c("a", "b"), n, replace = TRUE, prob = c(0.9, 0.1)),
             sample(c("a", "b"), n, replace = TRUE, prob = c(0.1, 0.9))),
    true_label = rep(c(1L, 2L), each = n)
  )
}

test_that("kpromm() recovers well-separated clusters with high accuracy", {
  df <- make_toy_data()
  fit <- kpromm(df, format = c("N", "N", "S", "0"), k = 2, p = 2,
                test_attribute = "true_label", n_init = 10, seed = 1)

  expect_s3_class(fit, "kpromm")
  expect_true(fit$best_accuracy > 0.85)
  expect_length(fit$alpha, 2)   # x1, x2
  expect_length(fit$beta, 1)    # cat1
  expect_equal(nrow(fit$all_attempts), 10)
  expect_true(is.finite(fit$best_objective))
})

test_that("kpromm() with normalize_num/normalize_cat = FALSE gives unit weights", {
  df <- make_toy_data()
  fit <- kpromm(df, format = c("N", "N", "S", "0"), k = 2, p = 2,
                normalize_num = FALSE, normalize_cat = FALSE,
                test_attribute = "true_label", n_init = 5, seed = 1)

  expect_true(all(fit$alpha == 1))
  expect_true(all(fit$beta == 1))
})

test_that("kpromm() errors when format length doesn't match ncol(data)", {
  df <- make_toy_data()
  # format has 3 entries but df has 4 columns -- should hit the stopifnot()
  expect_error(
    kpromm(df, format = c("N", "N", "S"), k = 2, p = 2, n_init = 1)
  )
})

test_that("kpromm() works with p != 2 (bisection prototype update)", {
  df <- make_toy_data()
  fit <- kpromm(df, format = c("N", "N", "S", "0"), k = 2, p = 1.5,
                test_attribute = "true_label", n_init = 5, seed = 1)
  expect_s3_class(fit, "kpromm")
  expect_true(fit$best_accuracy > 0.8)
})

test_that("kpromm() works with numeric-only and categorical-only data", {
  df <- make_toy_data()

  fit_num_only <- kpromm(df[, c("x1", "x2", "true_label")],
                          format = c("N", "N", "0"), k = 2, p = 2,
                          test_attribute = "true_label", n_init = 5, seed = 1)
  expect_s3_class(fit_num_only, "kpromm")
  expect_length(fit_num_only$beta, 0)

  fit_cat_only <- kpromm(df[, c("cat1", "true_label")],
                          format = c("S", "0"), k = 2, p = 2,
                          test_attribute = "true_label", n_init = 5, seed = 1)
  expect_s3_class(fit_cat_only, "kpromm")
  expect_length(fit_cat_only$alpha, 0)
})

test_that("print.kpromm() runs without error and mentions k and p", {
  df <- make_toy_data()
  fit <- kpromm(df, format = c("N", "N", "S", "0"), k = 2, p = 2,
                test_attribute = "true_label", n_init = 3, seed = 1)
  out <- capture.output(print(fit))
  expect_true(any(grepl("k=2", out)))
  expect_true(any(grepl("p=2.00", out)))
})

test_that("screen_p() returns a ranking covering the full p_grid", {
  skip_if_not_installed("mclust")
  df <- make_toy_data()
  scr <- screen_p(df, format = c("N", "N", "S", "0"), k = 2,
                          p_grid = c(1, 2, 3), test_attribute = "true_label",
                          n_init = 5, seed = 1)

  expect_named(scr, c("ranking", "fits", "central_p", "pool_restarts"))
  expect_true(scr$pool_restarts)
  expect_equal(sort(scr$ranking$p), c(1, 2, 3))
  expect_equal(nrow(scr$ranking), 3)
  expect_true(scr$central_p %in% c(1, 2, 3))
  expect_length(scr$fits, 3)
  # ranking should be sorted best-first (rank 1 has the highest ari_to_central)
  expect_equal(scr$ranking$rank[which.max(scr$ranking$ari_to_central)], 1)
})

test_that("screen_p() errors clearly if mclust is unavailable", {
  # This test only makes a meaningful assertion in an environment without
  # mclust installed; it's a no-op (trivially passes) otherwise.
  if (!requireNamespace("mclust", quietly = TRUE)) {
    df <- make_toy_data()
    expect_error(
      screen_p(df, format = c("N", "N", "S", "0"), k = 2, n_init = 2),
      "mclust"
    )
  } else {
    succeed("mclust is installed; skipping the unavailable-package check")
  }
})
