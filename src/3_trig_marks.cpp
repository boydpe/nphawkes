#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_k(NumericMatrix p0, NumericVector marks, NumericVector k_bins) {
  
  // PRELIMINARIES
  // Get the number of columns in p0
  // int n = p0.ncol(); // was .ncol()
  int n = pow(p0.size(), 0.5);
  // Get the number of bins
  int n_bins = k_bins.size() - 1; // was k_bins.size 
  // Create the vector for the numerator
  NumericVector k_num(n_bins);
  // Create the vector for the denominator
  NumericVector k_den(n_bins);
  // Create the vector for the estimates (output)
  NumericVector k(n_bins);
  // Create mark_bin as an integer valued object
  double mark_bin;
  
  // WORKING PORTION OF THE CODE
  // Start with column j and iterate
  
  // double mark_sum = sum(marks);
  // double mark_sum = std::accumulate(marks.begin(),
  //                    marks.end(), 0.0);

  if (sum(marks) == 0) {
    k[0] = 1;
  }
  else {
    for (int j = 0;j < n; j++) {

      // search for the bin the j-th mark is in
      for (int i = 0; i < n_bins; i++) {
        if (marks[j] > k_bins[i] && marks[j] <= k_bins[i + 1]) {
          mark_bin = i;
        }
      }

      // increment the denominator (count) for the bin by +1
      k_den[mark_bin] += 1;

      // sum the probabilities for each event i (row), occurring after event j
      for (int i = j + 1; i < n; i++) {
        k_num[mark_bin] += p0(i, j);
      }
    }
    // compute estimates of each bin
    for (int i = 0; i < n_bins; i++) {
      k[i] = k_num[i] / k_den[i];
    }
      //  return k; was here
  }
  return k;
}