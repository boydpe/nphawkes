#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix get_dist(NumericVector lat, NumericVector lon) {
  int n_i = lat.size();
  NumericMatrix dist_mat(n_i, n_i);
  
  for (int i = 0; i < n_i; i++){
    for (int j =0; j < n_i; j++){
      dist_mat(i,j) = pow(pow(lat[i] - lat[j], 2) + pow(lon[i] - lon[j], 2), 0.5);
    }
  }
  return dist_mat;
}
