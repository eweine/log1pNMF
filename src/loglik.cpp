#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;

// [[Rcpp::export]]
double get_sparse_term_loglik(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const std::vector<int> nonzero_y,
    const std::vector<int> nonzero_y_i_idx,
    const std::vector<int> nonzero_y_j_idx,
    const int num_nonzero_y
) {

  double sum = 0.0;
  double cp;
  int i;
  int j;

  for (int r = 0; r < num_nonzero_y; r++) {

    i = nonzero_y_i_idx[r];
    j = nonzero_y_j_idx[r];

    sum += nonzero_y[r] * log(exp(dot(U_T.col(i), V_T.col(j))) - 1);

  }

  return(sum);

}

// [[Rcpp::export]]
double get_loglik_exact(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const std::vector<int> nonzero_y,
    const std::vector<int> nonzero_y_i_idx,
    const std::vector<int> nonzero_y_j_idx,
    const int n,
    const int p
) {

  double sum = 0.0;
  double cp;

  int next_nz_i = nonzero_y_i_idx[0];
  int next_nz_j = nonzero_y_j_idx[0];
  unsigned int nz_idx = 0;

  for (int j = 0; j < p; j++) {

    for (int i = 0; i < n; i++) {

      cp = exp(dot(U_T.col(i), V_T.col(j)));

      if (i == next_nz_i && j == next_nz_j) {

        sum += nonzero_y[nz_idx] * log(cp - 1) - cp;
        nz_idx += 1;
        next_nz_i = nonzero_y_i_idx[nz_idx];
        next_nz_j = nonzero_y_j_idx[nz_idx];

      } else {

        sum -= cp;

      }

    }

  }

  return(sum);

}
