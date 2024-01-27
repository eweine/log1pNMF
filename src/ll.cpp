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
      ) - 1
    );
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

  double loglik = get_sparse_term_loglik_exact(
    U_T,
    V_T,
    y_nz_vals,
    y_nz_rows_idx,
    y_nz_cols_idx,
    y_nz_vals.size()
  );

  double lin_term = a1 * arma::dot(U_cs, arma::sum(V_T, 1));

  arma::mat U_T_U_V_T = U_T_U * V_T;
  double quad_term = a2 * arma::accu(V_T % U_T_U_V_T);

  loglik = loglik - lin_term - quad_term;
  return(loglik);

}
