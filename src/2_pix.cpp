#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector get_pix(NumericVector x_grid, NumericVector y_grid,
                      NumericVector lat, NumericVector lon,
                      NumericVector x_pix, NumericVector y_pix) {
  
  int n = lat.size();
  int m = x_grid.size();
  
  NumericVector x_loc(n);
  NumericVector y_loc(n);
  
  NumericMatrix locs(n,2);
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < (m-1); j++) {
      if(x_grid[j] < lon[i] && x_grid[j + 1] > lon[i]) {
        locs(i,0) = x_pix[j];
      }
    }
  }
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < (m-1); j++) {
      if(y_grid[j] < lat[i] && y_grid[j + 1] > lat[i]) {
        locs(i,1) = y_pix[j];
      }
    }
  }
  
  return locs;
}
  
  
