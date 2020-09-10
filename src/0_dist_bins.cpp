#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix get_dist_bins( NumericMatrix dist_mat, NumericVector h_bins) {
 
  int n = pow(dist_mat.size(), 0.5);
  int n_bins = h_bins.size() - 1;
  NumericMatrix dist_bins(n,n);
  
  for (int j = 0;j < n; j++) {
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n_bins; k++) {
        if (dist_mat(i,j) > h_bins[k] && dist_mat(i,j) <= h_bins[k + 1]) {
          dist_bins(i,j) = k;
        }
      }
    }
  }
  return dist_bins;
}