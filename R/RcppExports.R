# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

set_alphas <- function(formula) {
    .Call('_GMJMCMC_set_alphas', PACKAGE = 'GMJMCMC', formula)
}

vec_in_mat <- function(mat, vec, firstCol = 0L, lastCol = 0L) {
    .Call('_GMJMCMC_vec_in_mat', PACKAGE = 'GMJMCMC', mat, vec, firstCol, lastCol)
}

