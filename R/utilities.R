get_side_code <- function(side) {
  if (!(side %in% c("left", "right", "two_tailed"))) stop("`side` not recognized.")
  switch(side, left = -1L, right = 1L, two_tailed = 0L)
}


evaluate_simulation_results <- function(result_df, under_null) {
  fdp <- mean(under_null[which(result_df$rejected)])
  tpp <- mean(result_df$rejected[which(!under_null)])
  c(fdp = fdp, tpp = tpp)
}
