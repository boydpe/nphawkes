#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calc_d( NumericMatrix dist_mat2, int np) {
  
  int n = pow(dist_mat2.size(), 0.5);
  NumericVector d(n);
  
  for (int i = 0; i < n; i++){
    d[i] = dist_mat2(i,np) + 0.02;
  }

  // int r1 = r; 
  // int m = 0;
  // while (m < np) {
  //   m = 0; 
  //   for (int i = 1; i < n; i++) {
  //     if (dist_mat(i,0) < r1) {
  //       m = m + 1;
  //     }
  //   }
  //   r1 = r1 + dist_inc;
  // }
  // d(0,0) = r1; 
  // d(0,1) = m;
  // 
  // m = 0; 
  //   
  // for (int k = 1; k < n; k++){
  //   m = 0;
  //   int r2 = r; 
  //   
  //   while (m < np) {
  //     m = 0;
  //     for (int j = 0; j < k; j++) {
  //       if(dist_mat(k,j) < r2) {
  //         m = m + 1;
  //       } 
  //     }
  //     
  //     for (int i = k+1; i < n; i++) {
  //       if(dist_mat(i,k) < r2) {
  //         m = m + 1;
  //       }
  //     }
  //     r2 = r2 + dist_inc;
  //   }
  //   d(k,0) = r2;
  //   d(k,1) = m;
  // }
  
  return d; 
}
