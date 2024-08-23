#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix findMatchesC(NumericVector feature_rts, double drt_max) {
  int n = feature_rts.size();
  std::deque<std::tuple<int, int, double>> matches;

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double drt = std::abs(feature_rts[i] - feature_rts[j]);
      if (drt < drt_max) {
        matches.push_back(std::make_tuple(i + 1, j + 1, drt)); // R is 1-indexed
      }
    }
  }

  int num_matches = matches.size();
  NumericMatrix rt_matches(num_matches, 3);

  for (int k = 0; k < num_matches; ++k) {
    rt_matches(k, 0) = std::get<0>(matches[k]);
    rt_matches(k, 1) = std::get<1>(matches[k]);
    rt_matches(k, 2) = std::get<2>(matches[k]);
  }

  return rt_matches;
}

