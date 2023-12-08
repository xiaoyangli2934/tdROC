#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector AUC_calc2(NumericVector X, NumericVector W1, NumericVector W2){
  int n = X.size();
  double numer1 = 0.0;
  double denom1 = 0.0;
  double numer2 = 0.0;
  double denom2 = 0.0;
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      numer1 += W1[i] * (1 - W1[j]) * ( (X[i] > X[j]) + 0.5 * (X[i] == X[j]) );
      denom1 += W1[i] * (1 - W1[j]);
      numer2 += W1[i] * (1 - W1[j] - W2[j]) * ( (X[i] > X[j]) + 0.5 * (X[i] == X[j]) );
      denom2 += W1[i] * (1 - W1[j] - W2[j]);
    }
  }
  
  NumericVector result(2);
  
  result[0] =  numer1 / denom1;
  result[1] =  numer2 / denom2;
  
  return result;
}

