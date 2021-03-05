#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double calc_z(NumericVector times, NumericVector x_pix,
                     NumericVector y_pix, NumericMatrix tau) {
  
  double difx = x_pix[1] - x_pix[0];
  double dify = y_pix[1] - y_pix[0];
  int n = tau.nrow();
  NumericVector zz(n);
  
  for (int i = 0; i < n; i++) {
    zz[i] = ceil(max(times)) * difx * dify * tau(i,0);
  }
  
  double z = sum(zz);
  return z; 
}
  
