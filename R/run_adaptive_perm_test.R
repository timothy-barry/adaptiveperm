#' Run adaptive permutation test
#'
#' @return the results data frame
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
run_permutation_test <- function(Y_list, x, side = "two_tailed", h = 15L, alpha = 0.1, test_statistic = "MW", adaptive = TRUE, B = NULL) {
  # run checks
  if (!(test_statistic %in% c("MW", "sum_over_treated_units"))) {
    stop("`test_statistic` not recognized. Select between `MW` and `sum_over_treated_units`.")
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
  if (is.null(B)) B <- round(5 * length(Y_list)/alpha)
  # prepare to launch function
  funct_name <- switch(test_statistic,
                       "MW" = "run_permutation_test_mw",
                       "sum_over_treated_units" = "run_permutation_test_sum_over_treated_units")
  out <- do.call(what = funct_name, args = list(Y_list, x, side_code, h, alpha, adaptive, B))
  return(out)
}


run_permutation_test_mw <- function(Y_list, x, side_code, h, alpha, adaptive, B) {
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
  if (adaptive) {
    result <- run_adaptive_permutation_test_cpp(precomp_list, x, side_code, h, alpha, "compute_mw_test_statistic")
    p_values <- result$p_values; rejected <- result$rejected
  } else {
    p_values <- run_permutation_test_cpp(precomp_list, x, side_code, B, "compute_mw_test_statistic")
    rejected <- stats::p.adjust(p_values, method = "BH") < alpha
  }
  out <- data.frame(p_value = p_values, rejected = rejected)
  return(out)
}


run_permutation_test_sum_over_treated_units <- function(Y_list, x, side_code, h, alpha, adaptive, B) {
  precomp_list <- lapply(Y_list, FUN = function(y) list(y = y))
  if (adaptive) {
    result <- run_adaptive_permutation_test_cpp(precomp_list, x, side_code, h, alpha, "compute_sum_over_treated_units")
    p_values <- result$p_values; rejected <- result$rejected
  } else {
    p_values <- run_permutation_test_cpp(precomp_list, x, side_code, B, "compute_sum_over_treated_units")
    rejected <- stats::p.adjust(p_values, method = "BH") < alpha
  }
  out <- data.frame(p_value = p_values, rejected = rejected)
  return(out)
}
