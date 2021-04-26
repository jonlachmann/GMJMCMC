/*#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int vector_in_matrix(NumericMatrix mat, NumericVector vec) {
  int row = 0;
  for (int i = 0; i < mat.nrow(); i++) {
    if (mat(i,_) == vec) {
      row = i;
    }
  }
  return row;
}*/