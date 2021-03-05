#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix calc_tau(NumericVector x_pix, NumericVector y_pix,
                       NumericMatrix p0, NumericVector di,
                       NumericVector lon, NumericVector lat, 
                       NumericVector times) {
  
  int n = pow(p0.size(), 0.5);
  int m = x_pix.size();
  double pi = M_PI;
  NumericMatrix tau(pow(m,2),3);
  
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      double x = x_pix[i];
      double y = y_pix[j];
      double u = 0;
      
      for (int k = 0; k < n; k++) {
        u = u + (p0(k,k) * pow((2*pi*pow(di[k], 2)),-1) * 
          exp(-(pow((x - lon[k]), 2) + pow((y - lat[k]), 2)) / (2*pow(di[k], 2))));
      }
      tau(10*i + j, 0) = (1/ceil(max(times))) * u;
      tau(10*i + j, 1) = x;
      tau(10*i + j, 2) = y; 
    }
  }
  return tau;
}