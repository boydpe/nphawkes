#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_g(NumericMatrix p0, NumericVector g_bins, 
                       NumericMatrix  time_mat) {
  
  // numerator
  int n_l = g_bins.size() - 1;
  int n_i = pow(p0.size(), 0.5);
  NumericVector num_g(n_l);
  NumericVector g(n_l);
  NumericVector den_g(n_l);
  double diag_sum;
  
  for (int l = 0; l < n_l; l++) {
    for (int i = 0; i < n_i; i++){
      for (int j = i; j < n_i; j++){
        if(g_bins[l] < fabs(time_mat(i,j)) && fabs(time_mat(i,j)) <= g_bins[l + 1]){
          num_g[l] += p0(j,i);
        } else{
          num_g[l] += 0;
        }
      }
    }
  }
  
  // denominator
  
  for (int i = 0; i < n_l; i++){
    den_g[i] = (g_bins[i + 1] - g_bins[i]);
  }
  
  for (int i = 0; i < n_i; i++) {
    diag_sum += p0(i,i);
  }
  
  den_g = den_g * (sum(p0) - diag_sum);
  g = num_g/den_g;
  return g;
}

