#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector get_mark(NumericVector marks, NumericVector mark_breaks) {

  int n_i = marks.size();
  int n_j = mark_breaks.size() - 1;
  NumericVector mark_mat(n_i);

  for (int i = 0;i < n_i; i++) {
    for (int j = 0; j < n_j; j++) {
      if (marks[i] > mark_breaks[j] && marks[i] <= mark_breaks[j + 1]) {
        mark_mat[i] = j;
      }
    }
  }
  return mark_mat;
}
