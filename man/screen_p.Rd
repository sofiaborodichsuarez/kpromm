
\name{screen_p}
\alias{screen_p}
\title{Labels-free screening procedure for shortlisting candidate values of p}
\usage{
screen_p(data, format, k,
                 p_grid = c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0),
                 normalize_num = TRUE, normalize_cat = TRUE,
                 test_attribute = NULL, n_init = 20, seed = 1,
                 pool_restarts = TRUE)
}
\arguments{
\item{data, format, k, test_attribute, n_init, seed}{as in \code{\link{kpromm}}}

\item{p_grid}{candidate values of p to screen (default: the manuscript's
7-point grid)}

\item{normalize_num, normalize_cat}{passed through to \code{\link{kpromm}}
(default: k-PROMMSA)}

\item{pool_restarts}{logical; if \code{TRUE} (default), use the
cross-restart variant: pool every restart's partition across the full
\code{p_grid} and score each \code{p} by the mean Adjusted Rand Index of
all of its restart partitions to the pool's central partition. If
\code{FALSE}, use the basic variant: retain only each \code{p}'s single
best-objective partition and score by that partition's similarity to a
central partition identified among the best-objective partitions alone.}
}
\value{
a list with \code{ranking} (data.frame of \code{p}, \code{ari_to_central},
\code{rank}, ordered best-first), \code{fits} (named list of the full
\code{kpromm} fit object for every \code{p} in \code{p_grid}, so the
shortlisted fit(s) can be used directly without refitting), \code{central_p}
(the \code{p} whose partition is most representative of the ensemble), and
\code{pool_restarts} (echoes the argument, so the variant used is recorded
alongside the result)
}
\description{
Implements the manuscript's labels-free ensemble screening procedure for
candidate values of \code{p}. With \code{pool_restarts = TRUE} (the
default), all restart partitions are pooled across the full \code{p_grid};
the partition with highest average pairwise Adjusted Rand Index to the
pool is identified as the central partition, and each \code{p} is scored
by the mean Adjusted Rand Index of its restart partitions to that central
partition. With \code{pool_restarts = FALSE}, only the best-objective
partition for each \code{p} is retained, providing the cheaper basic
variant.
}
\details{
\strong{Honesty note (consistent with the manuscript):} this function
narrows the search over \code{p} using only the data and the clustering
objective; it does not reliably select the accuracy-optimal \code{p}. In
the manuscript's own validation (24 simulated scenarios with known ground
truth), the top-1 shortlisted \code{p} matched the true accuracy-optimal
\code{p} only 25\% of the time (versus ~14\% chance for a 7-point grid),
though the expected accuracy cost of trusting the top choice was usually
small. This function therefore does not auto-select a single \code{p} and
silently fit only that one; it always returns the full ranking (and the
fitted candidate for every \code{p}), so the analyst can inspect how much
the ranking would have changed their answer and combine it with domain
knowledge rather than treat it as a black box.

Requires the \pkg{mclust} package (for the Adjusted Rand Index).
}
\examples{
\donttest{
set.seed(1)
n <- 60
df <- data.frame(
  x1 = c(rnorm(n / 2, 0), rnorm(n / 2, 5)),
  x2 = c(rnorm(n / 2, 0), rnorm(n / 2, 5)),
  cat1 = sample(c("a", "b"), n, replace = TRUE)
)
scr <- screen_p(df, format = c("N", "N", "S"), k = 2,
                        p_grid = c(1, 2, 3), n_init = 5)
scr$ranking
}
}
