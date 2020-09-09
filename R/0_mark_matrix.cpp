#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector get_mark(NumericVector marks, NumericVector k_bins) {
  
  int n_i = marks.size();
  int n_j = k_bins.size() - 1;
  NumericVector mark_mat(n_i);

  for (int i = 0;i < n_i; i++) {
    for (int j = 0; j < n_j; j++) {
      if (marks[i] > k_bins[j] && marks[i] <= k_bins[j + 1]) {
        mark_mat[i] = j;
      }
    }
  } 
  return mark_mat;
}
