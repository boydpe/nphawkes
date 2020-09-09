#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double calc_br( NumericMatrix p0, NumericVector times) {
  
  double br;
  double diag_sum;
  int n_i = pow(p0.size(), 0.5);
  int t = ceil(max(times));
  
  for (int i = 0; i < n_i; i++) {
      diag_sum += p0(i,i);
  }
  br = diag_sum/t;
  return br;
}
