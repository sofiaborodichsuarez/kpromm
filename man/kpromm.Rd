
\name{kpromm}
\alias{kpromm}
\title{Fit k-PROMM / k-PROMMSA to a mixed-type dataset}
\usage{
kpromm(data, format, k, p,
       normalize_num = TRUE, normalize_cat = TRUE,
       test_attribute = NULL, n_init = 100, seed = 1,
       norm_sample_n = NULL)
}
\arguments{
\item{data}{data.frame containing all attributes (including any id/class
columns to be ignored via \code{format = "0"} or excluded via
\code{test_attribute})}

\item{format}{character vector, one entry per column of \code{data}:
\code{"N"} (numeric), \code{"S"} (categorical/string), or \code{"0"}
(ignore)}

\item{k}{number of clusters}

\item{p}{Minkowski exponent (p >= 1)}

\item{normalize_num, normalize_cat}{logical, independently toggle
statistical normalization of the numeric and categorical halves of the
distance. Both default to \code{TRUE} (k-PROMMSA); setting both to
\code{FALSE} recovers unnormalized k-PROMM.}

\item{test_attribute}{optional column name/index holding the true class
label, used only to report accuracy, never to fit the clustering}

\item{n_init}{number of random restarts; the best-objective partition
across restarts is returned}

\item{seed}{base RNG seed; restart \code{i} uses \code{seed + i}}

\item{norm_sample_n}{optional subsample size for the numeric
normalization-weight estimator; \code{NULL} (default) uses the full-data
estimator. Only useful for datasets large enough (tens of thousands of
rows) that the full O(N^2) pairwise computation is impractically slow.}
}
\value{
An object of class \code{"kpromm"}: a list with the fitted normalization
weights (\code{alpha}, \code{beta}), the best-objective partition
(\code{best}), its objective value (\code{best_objective}) and accuracy
against \code{test_attribute} if supplied (\code{best_accuracy}), and the
full per-restart objective/accuracy trace (\code{all_attempts}).
}
\description{
Fits the k-PROMM (unnormalized) or k-PROMMSA (statistically normalized)
clustering algorithm to a data frame containing numeric and/or categorical
attributes, using the mixed Minkowski p-metric.
}
\examples{
\donttest{
set.seed(1)
n <- 60
df <- data.frame(
  x1 = c(rnorm(n / 2, 0), rnorm(n / 2, 5)),
  x2 = c(rnorm(n / 2, 0), rnorm(n / 2, 5)),
  cat1 = sample(c("a", "b"), n, replace = TRUE),
  true_label = rep(1:2, each = n / 2)
)
fit <- kpromm(df, format = c("N", "N", "S", "0"), k = 2, p = 2,
              test_attribute = "true_label", n_init = 10)
print(fit)
}
}
