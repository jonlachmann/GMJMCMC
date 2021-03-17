#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::string set_alphas(std::string formula) {
  int alpha = 0;
  for (int i = 0; i < formula.length(); i++) {
    if (formula[i] == '?') {
      formula.replace(i, 1, "a");
      formula.insert(i+1, std::to_string(alpha));
      alpha++;
    }
  }
  return formula;
}