#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector bins(NumericVector u,NumericVector v){
  int n_i = u.size();
  int n_j = v.size();
  NumericVector x(n_i);
  
  for (int i = 0; i <= n_i; i++) {
    for (int j = 0; j <= n_j; j++) {
      if(v[j] < u[i] && u[i] <= v[j + 1]){
        x[i] = j;
      }
    }
  } 
  return x;
}