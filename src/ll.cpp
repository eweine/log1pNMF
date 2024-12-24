#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>
#include "ll.h"

using namespace Rcpp;
using namespace arma;

double get_sparse_term_loglik_exact(
    const arma::mat U_T,
    const arma::mat V_T,
    const std::vector<int> nonzero_y,
    const std::vector<int> nonzero_y_i_idx,
    const std::vector<int> nonzero_y_j_idx,
    const arma::vec s,
    const int num_nonzero_y
) {

  double sum = 0.0;

  #pragma omp parallel for reduction(+:sum)
  for (int r = 0; r < num_nonzero_y; r++) {

    sum += nonzero_y[r] * log(
      exp(
        dot(
          U_T.col(nonzero_y_i_idx[r]), V_T.col(nonzero_y_j_idx[r])
        )
      ) - s[nonzero_y_i_idx[r]]
    );
  }

  return sum;
}

double get_sparse_term_loglik_exact_add_const(
    const arma::mat U_T,
    const arma::mat V_T,
    const std::vector<int> nonzero_y,
    const std::vector<int> nonzero_y_i_idx,
    const std::vector<int> nonzero_y_j_idx,
    const double c,
    const int num_nonzero_y
) {

  double sum = 0.0;

  #pragma omp parallel for reduction(+:sum)
  for (int r = 0; r < num_nonzero_y; r++) {

    sum += nonzero_y[r] * log(
      exp(
        dot(
          U_T.col(nonzero_y_i_idx[r]), V_T.col(nonzero_y_j_idx[r])
        )
      ) - c
    );
  }

  return sum;
}

double get_dense_term_loglik_exact(
    const arma::mat U_T,
    const arma::mat V_T,
    const int n,
    const int p
) {

  double sum = 0.0;

  #pragma omp parallel for reduction(+:sum)
  for (int i = 0; i < n; i++) {

    for (int j = 0; j < p; j++) {

      sum -= exp(dot(U_T.col(i), V_T.col(j)));

    }
  }

  return sum;
}

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
) {

  // Commenting out to avoid compilation error

  // double loglik = get_sparse_term_loglik_exact(
  //   U_T,
  //   V_T,
  //   y_nz_vals,
  //   y_nz_rows_idx,
  //   y_nz_cols_idx,
  //   y_nz_vals.size()
  // );

  double loglik = 0.0;

  double lin_term = a1 * arma::dot(U_cs, arma::sum(V_T, 1));

  arma::mat U_T_U_V_T = U_T_U * V_T;
  double quad_term = a2 * arma::accu(V_T % U_T_U_V_T);

  loglik = loglik - lin_term - quad_term;
  return(loglik);

}

double get_sparse_term_loglik_quad_sparse_approx(
    const arma::mat U_T,
    const arma::mat V_T,
    const std::vector<int> nonzero_y,
    const std::vector<int> nonzero_y_i_idx,
    const std::vector<int> nonzero_y_j_idx,
    const int num_nonzero_y,
    arma::vec s,
    const double a1,
    const double a2
) {

  double sp_term = 0.0;
  double lin_correction = 0.0;
  double quad_correction = 0;
  double cp;

  #pragma omp parallel for reduction(+:sp_term, lin_correction, quad_correction)
  for (int r = 0; r < num_nonzero_y; r++) {

    cp = dot(U_T.col(nonzero_y_i_idx[r]), V_T.col(nonzero_y_j_idx[r]));

    sp_term += nonzero_y[r] * log(exp(cp) - 1) -
      s[nonzero_y_i_idx[r]] * exp(cp);
    lin_correction += s[nonzero_y_i_idx[r]] * cp;
    quad_correction += s[nonzero_y_i_idx[r]] * cp * cp;

  }

  double ll = sp_term + a1 * lin_correction + a2 * quad_correction;

  return(ll);

}

double get_loglik_quad_approx_sparse(
    const arma::mat U_T,
    const arma::mat V_T,
    const std::vector<int> y_nz_vals,
    const std::vector<int> y_nz_rows_idx,
    const std::vector<int> y_nz_cols_idx,
    const arma::vec s,
    const double a1,
    const double a2
) {

  double loglik = get_sparse_term_loglik_quad_sparse_approx(
    U_T,
    V_T,
    y_nz_vals,
    y_nz_rows_idx,
    y_nz_cols_idx,
    y_nz_vals.size(),
    s,
    a1,
    a2
  );

  double lin_term = a1 * arma::dot(U_T * s, arma::sum(V_T, 1));

  arma::mat U_T_U_V_T = (U_T.each_row() % s.t()) * U_T.t() * V_T;
  double quad_term = a2 * arma::accu(V_T % U_T_U_V_T);

  loglik = loglik - lin_term - quad_term;
  return(loglik);

}


double get_sparse_term_loglik_lin_sparse_approx(
    const arma::mat U_T,
    const arma::mat V_T,
    const std::vector<int> nonzero_y,
    const std::vector<int> nonzero_y_i_idx,
    const std::vector<int> nonzero_y_j_idx,
    const arma::vec s,
    const int num_nonzero_y,
    const double a
) {

  double sp_term = 0.0;
  double lin_correction = 0.0;
  double cp;

  #pragma omp parallel for reduction(+:sp_term, lin_correction)
  for (int r = 0; r < num_nonzero_y; r++) {

    cp = dot(U_T.col(nonzero_y_i_idx[r]), V_T.col(nonzero_y_j_idx[r]));

    sp_term += nonzero_y[r] * log(exp(cp) - 1) -
      s[nonzero_y_i_idx[r]] * exp(cp);

    lin_correction += s[nonzero_y_i_idx[r]] * cp;

  }

  double ll = sp_term + a * lin_correction;

  return(ll);

}

double get_loglik_lin_approx_sparse(
    const arma::mat U_T,
    const arma::mat V_T,
    const arma::vec U_T_s,
    const std::vector<int> y_nz_vals,
    const std::vector<int> y_nz_rows_idx,
    const std::vector<int> y_nz_cols_idx,
    const arma::vec s,
    const double a
) {

  double loglik = get_sparse_term_loglik_lin_sparse_approx(
    U_T,
    V_T,
    y_nz_vals,
    y_nz_rows_idx,
    y_nz_cols_idx,
    s,
    y_nz_vals.size(),
    a
  );

  double lin_term = a * arma::dot(U_T_s, arma::sum(V_T, 1));

  loglik = loglik - lin_term;
  return(loglik);

}

double get_loglik_exact(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const std::vector<int>& y_nz_vals,
    const std::vector<int>& y_nz_rows_idx,
    const std::vector<int>& y_nz_cols_idx,
    const arma::vec& s,
    const int n,
    const int p
) {

  double loglik_sparse_term = get_sparse_term_loglik_exact(
    U_T,
    V_T,
    y_nz_vals,
    y_nz_rows_idx,
    y_nz_cols_idx,
    s,
    y_nz_vals.size()
  );

  double loglik_dense_term = get_dense_term_loglik_exact(
    U_T,
    V_T,
    n,
    p
  );

  double loglik = loglik_sparse_term + loglik_dense_term;

  return(loglik);

}

double get_loglik_exact_add_const(
    const arma::mat U_T,
    const arma::mat V_T,
    const std::vector<int> y_nz_vals,
    const std::vector<int> y_nz_rows_idx,
    const std::vector<int> y_nz_cols_idx,
    const double c,
    const int n,
    const int p
) {

  double loglik_sparse_term = get_sparse_term_loglik_exact_add_const(
    U_T,
    V_T,
    y_nz_vals,
    y_nz_rows_idx,
    y_nz_cols_idx,
    c,
    y_nz_vals.size()
  );

  double loglik_dense_term = get_dense_term_loglik_exact(
    U_T,
    V_T,
    n,
    p
  );

  double loglik = loglik_sparse_term + loglik_dense_term;

  return(loglik);

}

double get_loglik_nn_glmpca(
    const arma::mat U_T,
    const arma::mat V_T,
    const std::vector<int> y_nz_vals,
    const std::vector<int> y_nz_rows_idx,
    const std::vector<int> y_nz_cols_idx,
    const int n,
    const int p
) {

  double sp_term = 0;

  #pragma omp parallel for reduction(+:sp_term)
  for (int r = 0; r < y_nz_vals.size(); r++) {

    sp_term += y_nz_vals[r] *
        dot(
          U_T.col(y_nz_rows_idx[r]), V_T.col(y_nz_cols_idx[r])
        );

  }

  double ll = get_dense_term_loglik_exact(U_T, V_T, n, p) + sp_term;
  return ll;

}
