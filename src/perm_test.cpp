#include <Rcpp.h>
#include <boost/random.hpp>
#include <stdexcept>
#include <algorithm>
#include <string>
#include "utilities.h"
using namespace Rcpp;


// [[Rcpp::export]]
std::vector<double> run_permutation_test_cpp(List precomp_list, IntegerVector x, int side_code, int B, std::string test_stat_str) {
  // define variables and objects
  int n = x.length(), m = precomp_list.length(), n_trt;
  std::vector<double> original_statistics(m), p_values(m), n_right_losses(m, 0), n_left_losses(m, 0);
  std::vector<int> trt_idxs;
  double curr_test_stat, B_doub = static_cast<double>(B);

  // select the test statistic
  double (*funct)(List, const std::vector<int>&, int) = nullptr;
  if (test_stat_str == "compute_score_stat") {
    funct = compute_score_stat;
  } else if (test_stat_str == "compute_mean_over_treated_units") {
    funct = compute_mean_over_treated_units;
  } else if (test_stat_str == "compute_mw_test_statistic") {
    funct = compute_mw_test_statistic;
  } else {
    throw std::invalid_argument("Test statistic not recognized.");
  }

  // populate the trt_idxs vector
  for (int i = 0; i < x.length(); i++) if (x[i] == 1) trt_idxs.push_back(i);
  n_trt = trt_idxs.size();
  if (n_trt == 0) throw std::invalid_argument("Zero treatment units.");

  // compute the original test statistics
  for (int i = 0; i < m; i++) {
    original_statistics[i] = funct(precomp_list(i), trt_idxs, n_trt);
  }

  // define objects related to random permutations
  boost::random::mt19937 generator(4);
  boost::random::uniform_real_distribution<double> distribution(0, 1);
  double n_doub = static_cast<double>(n);
  std::vector<double> i_doub_array(n_trt);
  for (int i = 0; i < n_trt; i++) i_doub_array[i] = static_cast<double>(i);
  std::vector<int> random_samp(n);
  for (int i = 0; i < n; i++) random_samp[i] = i;

  // iterate through time, computing null test statistics and updating n_right_losses
  for (int k = 0; k < B; k++) {
    // generate a random permutation
    draw_wor_sample(generator, distribution, i_doub_array, random_samp, n_trt, n_doub);
    // iterate over hypotheses
    for (int i = 0; i < m; i ++) {
      // compute the test statistic
      curr_test_stat = funct(precomp_list(i), random_samp, n_trt);
      // determine whether we have a loss; if so, increment n_right_losses or n_left_losses
      n_right_losses[i] += (curr_test_stat >= original_statistics[i] ? 1.0 : 0.0);
      n_left_losses[i] += (curr_test_stat <= original_statistics[i] ? 1.0 : 0.0);
    }
  }

  // compute p-values
  if (side_code == 1) {
    for (int i = 0; i < m; i ++) p_values[i] = (1.0 + n_right_losses[i])/(1.0 + B_doub);
  } else if (side_code == -1) {
    for (int i = 0; i < m; i ++) p_values[i] = (1.0 + n_left_losses[i])/(1.0 + B_doub);
  } else {
    for (int i = 0; i < m; i ++) p_values[i] = std::min(1.0, 2.0 * std::min((1.0 + n_right_losses[i])/(1.0 + B_doub), (1.0 + n_left_losses[i])/(1.0 + B_doub)));
  }

  return p_values;
}
