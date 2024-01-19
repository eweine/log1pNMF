#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
arma::vec solve_pois_reg_log1p (
    const arma::mat X,
    const arma::vec y,
    const arma::uvec y_nz_idx,
    arma::vec b,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const double alpha,
    const double beta
) {

  double first_deriv;
  double second_deriv;
  double newton_dir;
  double newton_dec;
  vec eta = X * b;
  vec exp_eta = exp(eta);
  vec exp_eta_nz_m1 = exp_eta.elem(y_nz_idx) - 1;
  vec eta_proposed;
  vec exp_deriv_term;
  double t;
  double f_proposed;
  unsigned int i, j;
  double b_j_og;

  double current_lik = sum(exp_eta) - dot(
    y,
    log(exp_eta_nz_m1)
  );

  int num_indices = update_indices.size();

  for (int update_num = 1; update_num <= num_iter; update_num++) {

    for (i = 0; i < num_indices; i++) {
      j = update_indices[i];

      exp_deriv_term = exp_eta % X.col(j);

      first_deriv    = sum(exp_deriv_term) - dot(
        y,
        exp_deriv_term.elem(y_nz_idx) / exp_eta_nz_m1
      );

      second_deriv   = dot(exp_deriv_term, X.col(j)) + dot(
        y,
        (exp_deriv_term.elem(y_nz_idx) % exp_deriv_term.elem(y_nz_idx)) /
          square(exp_eta_nz_m1)
      );

      newton_dir     = first_deriv / second_deriv;

      // I need to handle the non-negativity constraint here
      if (newton_dir < 0) {

        t = 1.0;

      } else if (b[j] >= 1e-16) {

        t = std::min(b[j] / newton_dir, 1.0);

      } else {

        continue;

      }

      newton_dec    = alpha * first_deriv * newton_dir;
      b_j_og        = b[j];
      while (true) {
        b[j]             = b_j_og - t * newton_dir;
        eta_proposed     = eta + (b[j] - b_j_og) * X.col(j);
        exp_eta = exp(eta_proposed);
        exp_eta_nz_m1 = exp_eta.elem(y_nz_idx) - 1;
        f_proposed = sum(exp_eta) - dot(
            y,
            log(exp_eta_nz_m1)
        );
        if (f_proposed <= current_lik - t*newton_dec) {
          eta = eta_proposed;
          current_lik = f_proposed;
          break;
        } else {
          t *= beta;
        }
      }
    }
  }

  return(b);

}

// Y is an nxm matrix (each col is an n-dim data vec)
// X is an nxp matrix (each row is a p-dim covariate)
// B is a pxm matrix (each col is a p-dim reg coef)
// [[Rcpp::export]]
arma::mat regress_cols_of_Y_on_X_log1p_pois_exact(
  const arma::mat& X,
  Rcpp::List Y,
  Rcpp::List Y_nz_idx,
  arma::mat& B,
  const std::vector<int> update_indices,
  unsigned int num_iter,
  const double alpha,
  const double beta
) {

  for (int j = 0; j < B.n_cols; j++) {

    B.col(j) = solve_pois_reg_log1p (
      X,
      Y[j],
      Y_nz_idx[j],
      B.col(j),
      update_indices,
      num_iter,
      alpha,
      beta
    );

  }

  return(B);

}
