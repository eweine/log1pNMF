#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "ll.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

double get_sparse_term_loglik_exact(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const std::vector<int>& nonzero_y,
    const std::vector<int>& nonzero_y_i_idx,
    const std::vector<int>& nonzero_y_j_idx,
    const double cc_alpha,
    const arma::vec& s,
    const int num_nonzero_y
) {

  double sum = 0.0;

  #pragma omp parallel for reduction(+:sum)
  for (int r = 0; r < num_nonzero_y; r++) {

    sum += nonzero_y[r] * log(
      exp(
        dot(
          U_T.col(nonzero_y_i_idx[r]), V_T.col(nonzero_y_j_idx[r])
        ) / cc_alpha
      ) - s[nonzero_y_i_idx[r]]
    );
  }

  return sum;
}

double get_dense_term_loglik_exact(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const double cc_alpha,
    const int n,
    const int p
) {

  double sum = 0.0;

  #pragma omp parallel for reduction(+:sum)
  for (int i = 0; i < n; i++) {

    for (int j = 0; j < p; j++) {

      sum -= exp(dot(U_T.col(i), V_T.col(j)) / cc_alpha);

    }
  }

  return sum;
}

double get_sparse_term_loglik_quad_sparse_approx(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const std::vector<int>& nonzero_y,
    const std::vector<int>& nonzero_y_i_idx,
    const std::vector<int>& nonzero_y_j_idx,
    const int num_nonzero_y,
    const arma::vec& s,
    const double a1,
    const double a2
) {

  double ll = 0.0;
  double cp;

  #pragma omp parallel for reduction(+:ll) private(cp)
  for (int r = 0; r < num_nonzero_y; r++) {

    cp = dot(U_T.col(nonzero_y_i_idx[r]), V_T.col(nonzero_y_j_idx[r]));

    ll += nonzero_y[r] * log(expm1(cp)) -
      s[nonzero_y_i_idx[r]] * exp(cp) + a1 * s[nonzero_y_i_idx[r]] * cp +
      a2 * s[nonzero_y_i_idx[r]] * cp * cp;

  }

  return(ll);

}

double get_loglik_quad_approx_sparse(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const std::vector<int>& y_nz_vals,
    const std::vector<int>& y_nz_rows_idx,
    const std::vector<int>& y_nz_cols_idx,
    const arma::vec& s,
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


// [[Rcpp::export]]
double get_loglik_exact(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const std::vector<int>& y_nz_vals,
    const std::vector<int>& y_nz_rows_idx,
    const std::vector<int>& y_nz_cols_idx,
    const arma::vec& s,
    const double cc_alpha,
    const int n,
    const int p
) {

  double loglik_sparse_term = get_sparse_term_loglik_exact(
    U_T,
    V_T,
    y_nz_vals,
    y_nz_rows_idx,
    y_nz_cols_idx,
    cc_alpha,
    s,
    y_nz_vals.size()
  );

  double loglik_dense_term = get_dense_term_loglik_exact(
    U_T,
    V_T,
    cc_alpha,
    n,
    p
  );

  double loglik = loglik_sparse_term + loglik_dense_term;

  return(loglik);

}
