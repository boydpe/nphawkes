#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix update_p(NumericMatrix p0, NumericMatrix time_mat, NumericMatrix dist_mat,
                       NumericVector mark_mat, NumericVector g, NumericVector h,
                       NumericVector k, NumericVector h_bins, NumericVector g_bins,
                       NumericVector k_bins, double br, NumericMatrix time_bins,
                       NumericMatrix dist_bins, NumericVector lat){

  int n_i = pow(p0.size(), 0.5);
  double num_p;
  double den_p ;
  NumericMatrix p(n_i, n_i);
  double y;
  // const double pi = 3.14159265358979323846;
  double pi = M_PI;
  double sum_mat;

  for (int i = 0; i < n_i; i++) {
    for ( int j = 0; j < n_i; j++){
      sum_mat += dist_mat(i,j);
    }
  }

  if ( sum(lat) != 0){
    for (int i = 0; i < n_i; i++){
      for (int j = 0; j < n_i; j++){
        if  ( j < i){
          num_p = (g[time_bins(j,i)]*h[dist_bins(i,j)] *
            k[mark_mat[j]]/(2*pi*dist_mat(i,j)));

          for (int  l = 0; l < i; l++){
            y += g[time_bins(l,i)]*h[dist_bins(i,l)]*k[mark_mat[l]]/
              (2*pi*dist_mat(i,l));
          }
          den_p = br + y;
          p(i,j) = num_p / den_p;
          y = 0;

        } else if (i == j) {
          num_p = br;

          for (int  l = 0; l < i; l++){
            y += g[time_bins(l,i)]*h[dist_bins(i,l)]*k[mark_mat[l]]/
              (2*pi*dist_mat(i,l));
          }
          den_p = br + y;
          p(i,j) = num_p / den_p;
          y = 0;
        } else {
          p(i,j) = 0;
        }
      }
    }
    p(0,0) = 1;
    return p;
  } else {

    for (int i = 0; i < n_i; i++){
      for (int j = 0; j < n_i; j++){
        if  ( j < i){
          num_p = g[time_bins(j,i)]*h[dist_bins(i,j)] *
            k[mark_mat[j]];

          for (int  l = 0; l < i; l++){
            y += g[time_bins(l,i)]*h[dist_bins(i,l)]*k[mark_mat[l]];
          }
          den_p = br +y;
          p(i,j) = num_p / den_p;
          y = 0;

        }else if (i == j) {
          num_p = br;

          for (int  l = 0; l < i; l++){
            y += g[time_bins(l,i)]*h[dist_bins(i,l)]*k[mark_mat[l]];
          }
          den_p = br + y;
          p(i,j) = num_p / den_p;
          y = 0;
        }

        else {
          p(i,j) = 0;//change to 0
        }
      }
    }
    p(0,0) = 1;
    return p;
  }
}
