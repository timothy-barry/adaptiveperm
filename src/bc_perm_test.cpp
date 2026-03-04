#include <Rcpp.h>
#include <boost/random.hpp>
#include <stdexcept>
#include <algorithm>
#include <string>
#include "utilities.h"
using namespace Rcpp;


// [[Rcpp::export]]
std::vector<double> run_bc_permutation_test_cpp(List precomp_list, IntegerVector x, int side_code, int h, int B, std::string test_stat_str) {
  // define variables and objects
  int n = x.length(), m = precomp_list.length(), n_trt;
  std::vector<bool> active_set(m, true), stopped_early(m, false);
  std::vector<double> stop_times(m, 0), original_statistics(m), p_values(m), n_right_losses(m, 0), n_left_losses(m, 0), n_losses(m, 0);
  std::vector<int> trt_idxs;
  double curr_test_stat, h_doub = static_cast<double>(h), B_doub = static_cast<double>(B), t, multiplier = (side_code == 0 ? 2.0 : 1.0);

  if (h <= 0) throw std::invalid_argument("h must be positive.");
  if (B <= 0) throw std::invalid_argument("B must be positive.");

  // select the test statistic
  double (*funct)(List, const std::vector<int>&, int) = nullptr;
  if (test_stat_str == "compute_sum_over_treated_units") {
    funct = compute_sum_over_treated_units;
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

  // iterate through time and stop only by futility
  for (int k = 0; k < B; k++) {
    t = static_cast<double>(k + 1);
    // generate a random permutation
    draw_wor_sample(generator, distribution, i_doub_array, random_samp, n_trt, n_doub);
    // iterate over hypotheses
    for (int i = 0; i < m; i ++) {
      if (active_set[i]) {
        // compute the test statistic
        curr_test_stat = funct(precomp_list(i), random_samp, n_trt);
        // determine whether we have a loss; if so, increment n_right_losses or n_left_losses
        n_right_losses[i] += (curr_test_stat >= original_statistics[i] ? 1.0 : 0.0);
        n_left_losses[i] += (curr_test_stat <= original_statistics[i] ? 1.0 : 0.0);
      }
    }

    // update n_losses
    if (side_code == -1) {
      n_losses = n_left_losses;
    } else if (side_code == 0) {
      for (int i = 0; i < m; i ++) n_losses[i] = std::min(n_left_losses[i], n_right_losses[i]);
    } else {
      n_losses = n_right_losses;
    }

    // move elements from active set into futility set
    for (int i = 0; i < m; i++) {
      if (active_set[i] && (n_losses[i] >= h_doub)) {
        active_set[i] = false;
        stopped_early[i] = true;
        stop_times[i] = t;
      }
    }

    // if no hypotheses remain active, stop iterating
    bool any_active = false;
    for (int i = 0; i < m; i++) {
      if (active_set[i]) {
        any_active = true;
        break;
      }
    }
    if (!any_active) break;
  }

  // compute the p-values using the Besag-Clifford rule
  for (int i = 0; i < m; i ++) {
    if (stopped_early[i]) {
      p_values[i] = std::min(1.0, multiplier * h_doub/stop_times[i]);
    } else {
      p_values[i] = std::min(1.0, multiplier * (n_losses[i] + 1.0)/(B_doub + 1.0));
    }
  }

  return p_values;
}
