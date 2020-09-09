#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double check_p( NumericMatrix p0, NumericMatrix p) {
  
  int n = p0.ncol();
  NumericMatrix diffs(n,n);
  double max_diff;
  
  
  for ( int i = 0; i < n; i++){
    for ( int j = 0; j < i; j ++){
      diffs(i,j) = fabs(p(i,j) - p0(i,j));
    }
  }
  max_diff = max(diffs);
  return max_diff;
}
