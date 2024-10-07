#include <Rcpp.h>
#include <boost/random.hpp>
#include <stdexcept>
#include <algorithm>
#include <string>
#include "utilities.h"
using namespace Rcpp;


// [[Rcpp::export]]
List run_adaptive_permutation_test_cpp(List precomp_list, IntegerVector x, int side_code, int h, double alpha, std::string test_stat_str) {
  // define variables and objects
  int n = x.length(), m = precomp_list.length(), n_trt;
  std::vector<bool> active_set(m, true), futility_set(m, false), rejected_set(m, false);
  std::vector<double> stop_times(m), original_statistics(m), p_values(m), n_right_losses(m, 0), n_left_losses(m, 0), n_losses(m, 0);
  std::vector<int> trt_idxs;
  double curr_test_stat, t = 0, h_doub = static_cast<double>(h), m_doub = static_cast<double>(m), threshold, n_in_active_set, max_n_losses_active_set;

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

  // iterate through time
  for (int k = 0; k < 50000; k++) {
    // increment time
    t++;
    // generate a random permutation
    draw_wor_sample(generator, distribution, i_doub_array, random_samp, n_trt, n_doub);
    // iterate over hypotheses
    for (int i = 0; i < m; i ++) {
      // if in the active set, compute the test statistic and update gamma
      if (active_set[i]) {
        // compute the test statistic
        curr_test_stat = funct(precomp_list(i), random_samp, n_trt);
        // determine whether we have a loss; if so, increment n_right_losses
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
    max_n_losses_active_set = 0;
    n_in_active_set = 0;
    for (int i = 0; i < m; i++) {
      if (active_set[i]) {
        if (n_losses[i] >= h_doub) { // hit loss limit
          active_set[i] = false;
          futility_set[i] = true;
          stop_times[i] = t;
        } else {
          n_in_active_set ++;
          if (n_losses[i] > max_n_losses_active_set) max_n_losses_active_set = n_losses[i];
        }
      }
    }

    // check for rejection of active set hypotheses
    threshold = t + h_doub - ((side_code == 0 ? 2.0 : 1.0) * h_doub * m_doub)/(n_in_active_set * alpha);
    if (max_n_losses_active_set <= threshold) {
      for (int i = 0; i < m; i ++) {
        if (active_set[i]) {
          active_set[i] = false;
          rejected_set[i] = true;
          stop_times[i] = t;
        }
      }
      break;
    }
  }
  // compute the p-values
  for (int i = 0; i < m; i ++) {
    p_values[i] = std::min(1.0, (side_code == 0 ? 2.0 : 1.0) * h_doub/(stop_times[i] - static_cast<double>(n_losses[i]) + h_doub));
  }

  return List::create(Named("p_values") = p_values,
                      Named("rejected") = rejected_set);
}
