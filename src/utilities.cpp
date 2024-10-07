#include <boost/random.hpp>
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

/*
 * Functions for drawing WOR samples
 */
void draw_wor_sample(boost::random::mt19937& generator,
                     boost::random::uniform_real_distribution<double>& distribution,
                     const std::vector<double>& i_doub_array,
                     std::vector<int>& random_samp,
                     int n_trt, double n_doub) {
  double u;
  int pos;
  for (int i = 0; i < n_trt; i++) {
    u = distribution(generator);
    pos = floor((n_doub - i_doub_array[i]) * u);
    std::swap(random_samp[pos], random_samp[i]);
  }
  return;
}


// [[Rcpp::export]]
std::vector<std::vector<int>> generate_wor_sample_test(int n, int n_trt, int n_samples) {
  std::vector<std::vector<int>> out(n_samples);
  std::vector<int> random_samp(n);
  boost::random::mt19937 generator(4);
  boost::random::uniform_real_distribution<double> distribution(0, 1);
  double n_doub = static_cast<double>(n);
  std::vector<double> i_doub_array(n_trt);
  for (int i = 0; i < n_trt; i++) i_doub_array[i] = static_cast<double>(i);
  for (int i = 0; i < n; i++) random_samp[i] = i;
  for (int i = 0; i < n_samples; i ++) {
    draw_wor_sample(generator, distribution, i_doub_array, random_samp, n_trt, n_doub);
    out[i] = std::vector<int>(random_samp.begin(), random_samp.begin() + n_trt);
  }
  return out;
}

/*
 * Test statistics
 */
double compute_sum_over_treated_units(List precomp, const std::vector<int>& trt_idxs, int n_trt) {
  NumericVector r = precomp(0);
  double sum = 0, n_trt_doub = static_cast<double>(n_trt);
  for (int j = 0; j < n_trt; j ++) sum += r[trt_idxs[j]];
  return(sum/sqrt(n_trt_doub));
}


double compute_mw_test_statistic(List precomp, const std::vector<int>& trt_idxs, int n_trt) {
  NumericVector r = precomp(0);
  double sigma = precomp(1);
  int side_code = precomp(2);
  double n_trt_doub = static_cast<double>(n_trt);
  double n_cntrl_doub = static_cast<double>(r.size()) - n_trt_doub;
  double sum = 0;
  for (int i = 0; i < n_trt; i ++) sum += r[trt_idxs[i]];
  double statistic = sum - n_trt_doub * (n_trt_doub + 1.0)/2.0;
  double statistic_mean_0 = statistic - n_trt_doub * n_cntrl_doub/2.0;
  double correction;
  if (side_code == 0) {
    correction = (statistic_mean_0 == 0 ? 0 : (statistic_mean_0 > 0 ? 0.5 : -0.5));
  } else {
    correction = 0.5;
  }
  double z = (statistic_mean_0 - correction)/sigma;
  return z;
}
