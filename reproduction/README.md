# Reproduction materials for "Making Mixed-Type Clustering Decision-Ready: k-PROMMSA, Software, and Practical Guidance"

This folder contains the reproduction materials referenced by the manuscript's
Software section. It is a cleaned release containing the scripts and data required
for the reproducible analyses reported in the manuscript; exploratory and superseded
development files are not included.

**Scope note.** Most of the manuscript's headline numbers (the 13-condition factor-sweep
simulation, all 7 non-`kamila` columns of Table `tab:speed-ci`) were computed with the
Python implementation in `python/`, which the
manuscript's Software section also explicitly documents providing ("I release an
open-source R implementation, together with the Python implementation used for the
simulations"). This R code is: (a) the reference implementation the Software section
describes and whose interface (`kpromm()`, `screen_p()`) the paper documents, and (b)
the actual source of every `kamila`-derived number in the paper (`kamila` is R-only, no
Python equivalent exists), including Table `tab:speed-ci`'s `kamila` column, Table
`tab:ci-full`'s `kamila` confidence intervals, and Table `tab:realdata-full`'s `kamila`
column across all four real datasets.

## Files (R)

| File | Reproduces | Requires |
|---|---|---|
| `kpromm.R` | Core `kpromm()`/`screen_p()` implementation described in the Software section. Trimmed from the original working file: the unused bee-colony global-search variant (`bee_kpromm()`) has been removed, see note at the end of the file. | `clue`; `mclust` only if using `screen_p()` |
| `simulation_comparison.R` | The 13-condition factor-sweep simulation (Table `tab:speed-ci`, `kamila` column and its Appendix C confidence interval) | `kamila`, `clue` |
| `real_data_benchmark.R` | Table `tab:realdata-full`'s `kamila` column, all four datasets (Heart Disease, Credit Approval, Zoo, Dermatology) | `kamila`, `clue`, plus data files below |
| `run_bank_marketing_example.R` | The empirical application (Section "Empirical application"), `kpromm`/`screen_p` side | `clue` |
| `kamila_simulation_ci.R` | A focused, faster re-run of just `kamila`'s simulation numbers with proper replicate-level 95% CIs (the exact source of the CI half-widths in Appendix C's `kamila` column) | `kamila`, `clue` |
| `kamila_zoo.R` | A focused, faster re-run of just `kamila` on Zoo with a 95% CI | `kamila`, `clue` |

## Files (Python)

The `python/` subfolder is the actual Python implementation used to produce the
non-`kamila` columns of Table `tab:speed-ci` (the 13-condition factor-sweep simulation).
It is self-contained four files, no dependency on anything outside this folder.

| File | Role |
|---|---|
| `python/kmpmnorm_core.py` | Core k-PROMM/k-PROMMSA fitting routine (the Python-side equivalent of `kpromm()` in `kpromm.R`) |
| `python/methods_lib.py` | Runner functions for every comparator method (k-prototypes, PAM-Gower, HC-Gower, FAMD+k-means, FAMD+Student-t-mixture) plus the k-PROMM/k-PROMMSA runners used by the driver script below |
| `python/p_simulation_factors.py` | Synthetic mixed-type data generator (`base_data()`) implementing the 13 factor-sweep conditions |
| `python/p_simulation_full_comparison_v2.py` | Driver script: runs all 7 non-`kamila` methods across all 13 conditions, 10 replicates each, and writes `p_simulation_full_comparison_v2_results.json`. This is the actual script that produced Table `tab:speed-ci`'s non-`kamila` columns. |

Requires: `numpy`, `pandas`, `scipy`, `scikit-learn`, `kmodes`, `gower`, `prince`, `smm`
(all pip-installable). Run from within `python/`:

```bash
cd python
pip install numpy pandas scipy scikit-learn kmodes gower prince smm
python p_simulation_full_comparison_v2.py
```

Takes roughly 15-20 minutes for the full 13 conditions x 10 replicates; writes results
next to the script.

## Data files included

- `processed.cleveland.data` (Heart Disease, UCI)
- `crx.data` (Credit Approval, UCI)
- `zoo.csv` (Zoo, UCI, CSV with header row)
- `bank-additional.csv` (Bank Marketing, UCI, semicolon-delimited)

## Data file downloaded automatically

- `dermatology.data` (Dermatology, UCI dataset 33, CC BY 4.0; Ilter, N. & Guvenir, H.
  (1998). Dermatology [Dataset]. UCI Machine Learning Repository.
  https://doi.org/10.24432/C5FK5P) is required by `real_data_benchmark.R` but is not
  bundled in this folder. `load_dermatology()` in `real_data_benchmark.R` downloads and
  extracts it automatically from https://archive.ics.uci.edu/static/public/33/dermatology.zip
  on first use if it is not already present, so no manual step is required as long as the
  machine running the script has internet access. To fetch it manually instead (e.g. on
  an offline machine), download that zip from
  https://archive.ics.uci.edu/dataset/33/dermatology and extract `dermatology.data` into
  this directory before running the script.

## How to run

All scripts expect to be run with this directory as the R working directory (they use
`source("kpromm.R")` and `read.csv("<file>")` with bare relative filenames, no
subfolders). From R:

```r
setwd("path/to/reproduction")
install.packages(c("clue", "mclust", "kamila"))  # cluster, clustMixType, teigen also
                                                    # needed by real_data_benchmark.R /
                                                    # simulation_comparison.R for the
                                                    # non-kamila comparators they also run
source("simulation_comparison.R")       # ~30-45 min, all 8 methods
source("real_data_benchmark.R")         # auto-downloads dermatology.data if missing
source("run_bank_marketing_example.R")
source("kamila_simulation_ci.R")        # faster: kamila only, ~a few minutes
source("kamila_zoo.R")                  # faster: kamila only, Zoo only, <1 min
```

If you only want `kamila`'s numbers (the R-only pieces the manuscript actually depends
on), running just `kamila_simulation_ci.R` and `kamila_zoo.R` plus the `kamila` portion
of `real_data_benchmark.R` is much faster than the full 8-method scripts.

