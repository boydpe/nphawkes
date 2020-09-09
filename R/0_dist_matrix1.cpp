#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix get_dist(NumericVector lat, NumericVector lon) {
  int n_i = lat.size();
  NumericMatrix dist_mat(n_i, n_i);
  double R = 6371e3;
  double pi = M_PI;
  
  
  for (int i = 0; i < n_i; i++){
    for (int j =0; j < n_i; j++){
      double phi1 = lat[i]*pi/180; 
      double phi2 = lat[j]*pi/180;
      double latdiff = (lat[j] - lat[i])*pi/180;
      double londiff = (lon[j] - lon[i])*pi/180;
      
      double a = sin(latdiff/2)*sin(latdiff/2) + 
        cos(phi1)*cos(phi2)*sin(londiff/2)*sin(londiff/2);
      double b = 2*atan2(sqrt(a), sqrt(1-a));
      double d = R*b;
    
      dist_mat(i,j) = d;
    }
  }
  return dist_mat;
}
