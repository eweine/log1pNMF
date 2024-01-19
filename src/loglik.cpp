#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
double get_loglik_exact(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const arma::vec& nonzero_y,
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

  for (int i = 0; i < n; i++) {

    for (int j = 0; j < n; j++) {

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
