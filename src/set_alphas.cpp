#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List set_alphas(std::string formula) {
  int alpha = 0;
  for (int i = 0; i < formula.length(); i++) {
    if (formula[i] == '?') {
      alpha++;
      formula.replace(i, 1, "a");
      formula.insert(i+1, "["+std::to_string(alpha)+"]");
    }
  }
  List result;
  result["formula"] = formula;
  result["count"] = alpha;
  return result;
}