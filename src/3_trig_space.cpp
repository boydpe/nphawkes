#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector get_h(NumericMatrix p0, NumericVector h_bins, 
                       NumericMatrix  dist_mat) {
  
  // numerator
  int n_l = h_bins.size() - 1;
  int n_i = pow(p0.size(), 0.5);
  NumericVector num_h(n_l);
  NumericVector h(n_l);
  NumericVector den_h(n_l);
  double diag_sum;
  
  if( sum(dist_mat) == 0){
    h[0] = 1;
  }
  
  else {
    
    for (int l = 0; l < n_l; l++) {
      for (int i = 0; i < n_i; i++){
        for (int j = i; j < n_i; j++){
          if(h_bins[l] < dist_mat(i,j) && dist_mat(i,j) <= h_bins[l+1]){
            num_h[l] += p0(j,i);
          } else{
            num_h[l] += 0;
          }
        }
      }
    }
    
    // denominator
    
    for (int i = 0; i < n_l; i++){
      den_h[i] = (h_bins[i+1] - h_bins[i]);
    }
    
    for (int i = 0; i < n_i; i++) {
      diag_sum += p0(i,i);
    }
    den_h = den_h * (sum(p0) - diag_sum);
    h = num_h / den_h;
    // return h;
  }
  return h;
}
