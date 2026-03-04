#' Run permutation test
#'
#' `run_permutation_test()` runs a permutation test to test for association between a binary "treatment" vector \eqn{x \in \{0, 1\}^n} and a collection of response vectors \eqn{y_1, \dots, y_n \in R^n}.
#'
#'`Y_list` is a list of response vectors (one per hypothesis). Each vector should be of class `numeric`.
#'
#' `x` is a binary "treatment" vector containing only `0`s and `1`s. `x` should be of class `numeric` or `integer`.
#'
#' Users can select the test statistic to use via the `test_statistic` argument. Current options include `"MW"` (for the Mann-Whitney test statistic) and `"sum_over_treated_units"` (for the sum-over-treated-units test statistic). The sum-over-treated-units test statistic is defined as follows: \deqn{(1/\sqrt{s})\sum_{X_i = 1} Y_i,} where \deqn{s = \sum_{i=1}^n X_i} is the number of "treated" units within \eqn{X}.
#'
#' `side` is the sidedness of the test; users can select between `"left"`, `"right"`, or `"two_tailed"` for left-, right-, and two-tailed tests, respectively.
#'
#' `h` is the tuning parameter used within the adaptive permutation test. We stop computing permutations for a given hypothesis when the number of "losses" hits `h`, where, for a right-tailed test, a loss is defined as a null statistic exceeding the original statistic.
#'
#' `alpha` is the nominal false discovery rate (FDR) of the test. `method` controls which permutation procedure to run. Select `"classical"` for a fixed-`B` permutation test, `"anytime"` for the currently implemented adaptive test with early rejection and futility, or `"besag_clifford"` for adaptive stopping by futility only. Finally, `B` is the maximum number of permutations used in the `"classical"` and `"besag_clifford"` methods. The default value of `B` is \eqn{5m/\alpha}, where \eqn{m} is the number of hypotheses to test.
#'
#' @param Y_list (required) a list of response vectors
#' @param x (required) a binary vector of "treatments"
#' @param test_statistic (optional; default `"MW"`) a string indicating the test statistic to use. Current options include `"sum_over_treated_units"` and `"MW"`
#' @param side (optional; default `"two_tailed"`) the sidedness of the test, one of `"left"`, `"right"`, or `"two_tailed"`
#' @param h (optional; default `10`) the tuning parameter for the adaptive permutation test
#' @param alpha (optional; default `0.1`) the nominal false discovery rate
#' @param method (optional; default `"anytime"`) one of `"classical"`, `"anytime"`, or `"besag_clifford"`
#' @param B (optional; default `5 * length(Y)/alpha`) if `method` is `"classical"` or `"besag_clifford"`, the maximum number of permutations to compute for each hypothesis
#'
#' @return the results data frame containing columns `p_value` and `rejected`
#' @export
#'
#' @examples
#' n <- 50L
#' m <- 500L
#' x <- rbinom(n = n, size = 1, prob = 0.5)
#' under_null <- sample(c(rep(TRUE, 0.9 * m), rep(FALSE, 0.1 * m)))
#' Y_list <- sapply(X = seq_len(m), FUN = function(i) {
#'  mus <- x * (if (under_null[i]) 0 else 1)
#'  rnorm(n = n, mean = mus)
#' }, simplify = FALSE)
#'
#' # MW statistic
#' res <- run_permutation_test(Y_list, x, test_statistic = "MW")
#' evaluate_simulation_results(res, under_null)
#'
#' # sum over treated units statistic
#' res <- run_permutation_test(Y_list, x, test_statistic = "sum_over_treated_units")
#' evaluate_simulation_results(res, under_null)
run_permutation_test <- function(Y_list, x, side = "two_tailed", h = 10L, alpha = 0.1, test_statistic = "MW", method = "anytime", B = NULL) {
  # run checks
  if (!(test_statistic %in% c("MW", "sum_over_treated_units"))) {
    stop("`test_statistic` not recognized. Select between `MW` and `sum_over_treated_units`.")
  }
  if (!(method %in% c("classical", "anytime", "besag_clifford"))) {
    stop("`method` not recognized. Select between `classical`, `anytime`, and `besag_clifford`.")
  }
  # verify x is binary
  x <- as.integer(x)
  if (!all(x %in% c(0L, 1L))) {
    stop("x must be a binary vector")
  }
  # refactor x (if necessary)
  if (mean(x) > 0.5) x <- 1L - x
  # check that there is enough data for permutation testing
  if (2/choose(length(x), sum(x)) > 5e-4) {
    warning("Your sample size may be too small for permutation testing.")
  }
  side_code <- get_side_code(side)
  # set B (if not already set)
  if (is.null(B)) B <- as.integer(round(5 * length(Y_list)/alpha))
  B <- as.integer(B)
  # prepare to launch function
  funct_name <- switch(test_statistic,
                       "MW" = "run_permutation_test_mw",
                       "sum_over_treated_units" = "run_permutation_test_sum_over_treated_units")
  out <- do.call(what = funct_name, args = list(Y_list, x, side_code, h, alpha, method, B))
  return(out)
}


run_permutation_test_mw <- function(Y_list, x, side_code, h, alpha, method, B) {
  # iterate over hypotheses, performing precomputation
  precomp_list <- lapply(X = Y_list, FUN = function(y) {
    r <- rank(y)
    n_s1 <- as.double(sum(x == 1L))
    n_s2 <- as.double(length(x) - n_s1)
    n_ties <- as.numeric(table(r))
    sigma <- sqrt((n_s1 * n_s2/12) * ((n_s1 + n_s2 + 1) - sum(n_ties^3 - n_ties)/((n_s1 + n_s2) * (n_s1 + n_s2 - 1))))
    list(r = r, sigma = sigma, side_code = side_code)
  })
  # run the permutation test
  if (method == "anytime") {
    result <- run_adaptive_permutation_test_cpp(precomp_list, x, side_code, h, alpha, "compute_mw_test_statistic")
    p_values <- result$p_values; rejected <- result$rejected
  } else if (method == "besag_clifford") {
    p_values <- run_bc_permutation_test_cpp(precomp_list, x, side_code, h, B, "compute_mw_test_statistic")
    rejected <- stats::p.adjust(p_values, method = "BH") < alpha
  } else {
    p_values <- run_permutation_test_cpp(precomp_list, x, side_code, B, "compute_mw_test_statistic")
    rejected <- stats::p.adjust(p_values, method = "BH") < alpha
  }
  out <- data.frame(p_value = p_values, rejected = rejected)
  return(out)
}


run_permutation_test_sum_over_treated_units <- function(Y_list, x, side_code, h, alpha, method, B) {
  precomp_list <- lapply(Y_list, FUN = function(y) list(y = y))
  if (method == "anytime") {
    result <- run_adaptive_permutation_test_cpp(precomp_list, x, side_code, h, alpha, "compute_sum_over_treated_units")
    p_values <- result$p_values; rejected <- result$rejected
  } else if (method == "besag_clifford") {
    p_values <- run_bc_permutation_test_cpp(precomp_list, x, side_code, h, B, "compute_sum_over_treated_units")
    rejected <- stats::p.adjust(p_values, method = "BH") < alpha
  } else {
    p_values <- run_permutation_test_cpp(precomp_list, x, side_code, B, "compute_sum_over_treated_units")
    rejected <- stats::p.adjust(p_values, method = "BH") < alpha
  }
  out <- data.frame(p_value = p_values, rejected = rejected)
  return(out)
}
