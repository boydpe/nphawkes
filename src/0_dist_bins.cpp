#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix get_dist_bins( NumericMatrix dist_mat, NumericVector space_breaks) {

  int n = pow(dist_mat.size(), 0.5);
  int n_bins = space_breaks.size() - 1;
  NumericMatrix dist_bins(n,n);

  for (int j = 0;j < n; j++) {
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n_bins; k++) {
        if (dist_mat(i,j) > space_breaks[k] && dist_mat(i,j) <= space_breaks[k + 1]) {
          dist_bins(i,j) = k;
        }
      }
    }
  }
  return dist_bins;
}
