#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int vec_in_mat(NumericMatrix &mat, NumericVector &vec, int firstCol=0, int lastCol=0) {
    if (firstCol != 0) firstCol -= 1;
    if (lastCol == 0) lastCol = mat.ncol();
    int cols = lastCol - firstCol;
    int rows = mat.nrow();
    int equal = 0;
    for (int i = mat.nrow()-1; i > -1; i--) {
        // Reset equal counter
        equal = 0;
        for (int j = firstCol; j < lastCol; j++) {
            // Count how often the items are equal, avoiding if statements in the inner loop
            equal += (mat(i,j) == vec(j));
        }
        if (equal == cols) return i+1;
    }
    return 0;
}
