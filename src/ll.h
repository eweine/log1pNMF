#ifndef LL_H
#define LL_H

#include <RcppArmadillo.h>
#include <Rcpp.h>


double get_loglik_quad_approx_full(
    const arma::mat U_T,
    const arma::mat V_T,
    const arma::vec U_cs,
    const arma::mat U_T_U,
    const std::vector<int> y_nz_vals,
    const std::vector<int> y_nz_rows_idx,
    const std::vector<int> y_nz_cols_idx,
    const double a1,
    const double a2
);

double get_loglik_quad_approx_sparse(
    const arma::mat U_T,
    const arma::mat V_T,
    const std::vector<int> y_nz_vals,
    const std::vector<int> y_nz_rows_idx,
    const std::vector<int> y_nz_cols_idx,
    const arma::vec s,
    const double a1,
    const double a2
);

double get_loglik_lin_approx_sparse(
    const arma::mat U_T,
    const arma::mat V_T,
    const arma::vec U_T_s,
    const std::vector<int> y_nz_vals,
    const std::vector<int> y_nz_rows_idx,
    const std::vector<int> y_nz_cols_idx,
    const arma::vec s,
    const double a
);

double get_loglik_exact(
    const arma::mat U_T,
    const arma::mat V_T,
    const std::vector<int> y_nz_vals,
    const std::vector<int> y_nz_rows_idx,
    const std::vector<int> y_nz_cols_idx,
    const arma::vec s,
    const int n,
    const int p
);

#endif
