#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
double AUC_calc(NumericVector X, NumericVector W){
  int n = X.size();
  double numer = 0.0;
  double denom = 0.0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      numer += W[i] * (1 - W[j]) * ( (X[i] > X[j]) + 0.5 * (X[i] == X[j]) );
      denom += W[i] * (1 - W[j]);
    }
  }
  return numer / denom;
}
