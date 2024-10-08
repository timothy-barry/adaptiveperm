% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_adaptive_perm_test.R
\name{run_permutation_test}
\alias{run_permutation_test}
\title{Run permutation test}
\usage{
run_permutation_test(
  Y_list,
  x,
  side = "two_tailed",
  h = 15L,
  alpha = 0.1,
  test_statistic = "MW",
  adaptive = TRUE,
  B = NULL
)
}
\arguments{
\item{Y_list}{(required) a list of response vectors}

\item{x}{(required) a binary vector of "treatments"}

\item{side}{(optional; default \code{"two_tailed"}) the sidedness of the test, one of \code{"left"}, \code{"right"}, or \code{"two_tailed"}}

\item{h}{(optional; default \code{15}) the tuning parameter for the adaptive permutation test}

\item{alpha}{(optional; default \code{0.1}) the nominal false discovery rate}

\item{test_statistic}{(optional; default \code{"sum_over_treated_units"}) a string indicating the test statistic to use. Current options include \code{"sum_over_treated_units"} and \code{"MW"}}

\item{adaptive}{(optional; default \code{TRUE}) a logical indicating whether to run an adaptive permutation test (\code{TRUE}) or a classical permutation test (\code{FALSE})}

\item{B}{(optional; default \code{5 * length(Y)/alpha}) if \code{adaptive} is set to \code{FALSE}, the number of permutations to compute for each hypothesis}
}
\value{
the results data frame containing columns \code{p_value} and \code{rejected}
}
\description{
\code{run_permutation_test()} runs a permutation test to test for association between a binary "treatment" vector \eqn{x \in \{0, 1\}^n} and a collection of response vectors \eqn{y_1, \dots, y_n \in R^n}.
}
\details{
\code{Y_list} is a list of length \eqn{n} that contains the response vectors. Each vector should be of class \code{numeric}.

\code{x} is a binary "treatment" vector containing only \code{0}s and \code{1}s. \code{x} should be of class \code{numeric} or \code{integer}.

Users can select the test statistic to use via the \code{test_statistic} argument. Current options include \code{"MW"} (for the Mann-Whitney test statistic) and \code{"sum_over_treated_units"} (for the sum-over-treated-units test statistic). The sum-over-treated-units test statistic is defined as follows: \deqn{(1/\sqrt{s})\sum_{X_i = 1} Y_i,} where \deqn{s = \sum_{i=1}^n X_i} is the number of "treated" units within \eqn{X}.

\code{side} is the sidedness of the test; users can select between \code{"left"}, \code{"right"}, or \code{"two_tailed"} for left-, right-, and two-tailed tests, respectively.

\code{h} is the tuning parameter used within the adaptive permutation test. We stop computing permutations for a given hypothesis when the number of "losses" hits \code{h}, where, for a right-tailed test, a loss is defined as a null statistic exceeding the original statistic.

\code{alpha} is the nominal false discovery rate (FDR) of the test. \code{adaptive} is a flag indicating whether to run an adaptive permutation test (\code{TRUE}) or a classical permutation test (\code{FALSE}). Finally, \code{B} is the number of permutations to compute for each hypothesis if \code{adaptive} is set to \code{FALSE}. The default value of \code{B} is \eqn{5m/\alpha}, where \eqn{m} is the number of hypotheses to test.
}
\examples{
n <- 50L
m <- 500L
x <- rbinom(n = n, size = 1, prob = 0.5)
under_null <- sample(c(rep(TRUE, 0.9 * m), rep(FALSE, 0.1 * m)))
Y_list <- sapply(X = seq_len(m), FUN = function(i) {
 mus <- x * (if (under_null[i]) 0 else 1)
 rnorm(n = n, mean = mus)
}, simplify = FALSE)

# MW statistic
res <- run_permutation_test(Y_list, x, test_statistic = "MW")
evaluate_simulation_results(res, under_null)

# sum over treated units statistic
res <- run_permutation_test(Y_list, x, test_statistic = "sum_over_treated_units")
evaluate_simulation_results(res, under_null)
}
