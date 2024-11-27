#include <Rcpp.h>

// [[Rcpp::export]]
double fast_med(Rcpp::NumericVector xx) {
  Rcpp::NumericVector x = Rcpp::clone(Rcpp::na_omit(xx));
  std::size_t n = x.size() / 2;
  std::nth_element(x.begin(), x.begin() + n, x.end());

  if (x.size() % 2) return x[n]; 
  return (x[n] + *std::max_element(x.begin(), x.begin() + n)) / 2.;
}
