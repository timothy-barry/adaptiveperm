get_side_code <- function(side) {
  if (!(side %in% c("left", "right", "two_tailed"))) stop("`side` not recognized.")
  switch(side, left = -1L, right = 1L, two_tailed = 0L)
}


#' Evaluate simulation results
#'
#' `evaluate_simulation_results()` evaluates the results of `run_permutation_test()` when run on simulated data.
#'
#' @param result_df the results data frame outputted by `run_permutation_test()`
#' @param under_null a logical vector indicating whether a given hypothesis is under the null (`TRUE`) or under the alternative (`FALSE`)
#'
#' @return a vector containing the false discovery proportion (FDP) and true positive proportion (TPP)
#' @export
#'
#' @inherit run_permutation_test examples
evaluate_simulation_results <- function(result_df, under_null) {
  fdp <- mean(under_null[which(result_df$rejected)])
  tpp <- mean(result_df$rejected[which(!under_null)])
  c(fdp = fdp, tpp = tpp)
}
