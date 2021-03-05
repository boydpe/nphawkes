#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix calc_br_nonstat( NumericMatrix p0, NumericVector times, 
                       double z, NumericMatrix tau) {
  
  double diag_sum;
  int n_i = pow(p0.size(), 0.5);
  int m = tau.nrow();
  NumericMatrix br(m,3);
  
  for (int i = 0; i < n_i; i++) {
    diag_sum += p0(i,i);
  }
  
  for (int i = 0; i < m; i++){
    double val = diag_sum * tau(i,0)/ z; 
    br(i,0) = val;
    br(i,1) = tau(i,1);
    br(i,2) = tau(i,2); 
  }
  
  return br;
}
