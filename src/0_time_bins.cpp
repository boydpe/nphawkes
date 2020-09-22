#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix get_time_bins( NumericMatrix time_mat, NumericVector time_breaks) {

  int n = pow(time_mat.size(), 0.5);
  int n_bins = time_breaks.size() - 1;
  NumericMatrix time_bins(n,n);

  for (int j = 0;j < n; j++) {
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n_bins; k++) {
        if (time_mat(i,j) > time_breaks[k] && time_mat(i,j) <= time_breaks[k + 1]) {
          time_bins(i,j) = k;
        }
      }
    }
  }
  return time_bins;
}
