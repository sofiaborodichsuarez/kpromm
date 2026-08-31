# kpromm
[![DOI](https://zenodo.org/badge/1345135977.svg)](https://doi.org/10.5281/zenodo.22084226)

An R package implementing **k-PROMM** and **k-PROMMSA**: clustering
algorithms for mixed numerical/categorical data that combine a tunable
Minkowski exponent *p* with statistical normalization of feature vectors.

This package accompanies the manuscript *"Robust Mixed-Type Clustering for
Operational Applications: k-PROMMSA, Parameter Choice and Practical
Guidance"* (Sofia Borodich Suarez, NEOMA Business School). It implements
exactly the two functions the manuscript's Software section documents:

```r
kpromm(data, format, k, p, normalize_num = TRUE, normalize_cat = TRUE,
       test_attribute = NULL, n_init = 100, seed = 1, norm_sample_n = NULL)

screen_p(data, format, k, p_grid = c(1, 1.5, 2, 2.5, 3, 3.5, 4),
                normalize_num = TRUE, normalize_cat = TRUE,
                test_attribute = NULL, n_init = 20, seed = 1,
                pool_restarts = TRUE)
```

## Installation

Not yet on CRAN. Install the development version from source:

```r
# from this folder
install.packages(".", repos = NULL, type = "source")

# or, once pushed to GitHub:
# install.packages("remotes")
# remotes::install_github("sofiaborodichsuarez/kpromm")
```

Dependencies: `clue` (required), `mclust` (only needed for
`screen_p()`).

## Quick example

```r
library(kpromm)

set.seed(1)
n <- 100
df <- data.frame(
  x1   = c(rnorm(n / 2, 0), rnorm(n / 2, 6)),
  x2   = c(rnorm(n / 2, 0), rnorm(n / 2, 6)),
  cat1 = c(sample(c("a", "b"), n / 2, replace = TRUE, prob = c(.9, .1)),
           sample(c("a", "b"), n / 2, replace = TRUE, prob = c(.1, .9))),
  true_label = rep(1:2, each = n / 2)
)

# k-PROMMSA (normalized) at p = 2
fit <- kpromm(df, format = c("N", "N", "S", "0"), k = 2, p = 2,
              test_attribute = "true_label", n_init = 20)
print(fit)

# Unnormalized k-PROMM, for comparison
fit_unnorm <- kpromm(df, format = c("N", "N", "S", "0"), k = 2, p = 2,
                      normalize_num = FALSE, normalize_cat = FALSE,
                      test_attribute = "true_label", n_init = 20)

# Labels-free screening over a grid of p (no ground truth used)
scr <- screen_p(df, format = c("N", "N", "S", "0"), k = 2,
                        n_init = 10)
scr$ranking
```

## What this package deliberately does NOT include

- **`bee_kpromm()`**, an exploratory bee-colony global-search variant
  developed during earlier investigations on the UCI Adult dataset. Adult
  was dropped from the manuscript's real-data section in favor of Dermatology,
  and this variant is neither part of the manuscript nor required
  for its reproduction. It is therefore not included in this release.
- **`kprommsa()`**, a one-line convenience wrapper that forces
  `normalize_num = normalize_cat = TRUE`. Kept in the source as an internal,
  non-exported function so the package's *exported* surface matches the
  manuscript's Software section exactly (two functions).
- **A single `normalize` argument.** An earlier draft of the manuscript
  described a single `normalize` toggle, but the actual implementation has
  always used two independent arguments, `normalize_num` and
  `normalize_cat`. The manuscript text has been corrected to match the code
  rather than the other way around.

## Package status

Version 1.0.0 was tested under R 4.6.1 on macOS Tahoe 26.1
(Apple Silicon). The `testthat` suite passes all 25 tests, and the
package passes `R CMD check` with 0 errors, 0 warnings, and 0 notes.

No CRAN submission has been made or is implied. This repository provides
the open-source implementation and reproduction materials accompanying
the manuscript.

## Citation

Version 1.0.0 is permanently archived on Zenodo.

DOI: [10.5281/zenodo.22084227](https://doi.org/10.5281/zenodo.22084227)

See `CITATION.cff` for citation metadata.

## License

MIT (see `LICENSE`).
