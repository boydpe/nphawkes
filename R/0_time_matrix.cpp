#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix get_time(NumericVector times) {
  int n_i = times.size();
  NumericMatrix time_mat(n_i, n_i);
  
  for (int i = 0; i < n_i; i++){
    for (int j =0; j < n_i; j++){
      time_mat(i,j) = times[j] - times[i];
    }
  }
  return time_mat;
}
