#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix init_p0(NumericVector times) {
  int n_i = times.size();
  NumericMatrix p0(n_i, n_i);
  
  for (int i = 0; i < n_i; i++){
    for (int j = 0; j <= i; j++){
      p0(i,j) = (1.0/(i + 1));
    }
  } 
  return p0;
}

