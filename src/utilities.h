#include <Rcpp.h>
#include <boost/random.hpp>
using namespace Rcpp;

#ifndef DRAW_WOR_SAMPLE
#define DRAW_WOR_SAMPLE
void draw_wor_sample(boost::random::mt19937& generator, boost::random::uniform_real_distribution<double>& distribution, const std::vector<double>& i_doub_array, std::vector<int>& random_samp, int n_trt, double n_doub);
#endif

#ifndef COMPUTE_SCORE_STAT
#define COMPUTE_SCORE_STAT
double compute_score_stat(List precomp, const std::vector<int>& trt_idxs, int n_trt);
#endif

#ifndef COMPUTE_MEAN_OVER_TREATED_UNITS
#define COMPUTE_MEAN_OVER_TREATED_UNITS
double compute_mean_over_treated_units(List precomp, const std::vector<int>& trt_idxs, int n_trt);
#endif

#ifndef COMPUTE_MW_TEST_STATISTIC
#define COMPUTE_MW_TEST_STATISTIC
double compute_mw_test_statistic(List precomp, const std::vector<int>& trt_idxs, int n_trt);
#endif
