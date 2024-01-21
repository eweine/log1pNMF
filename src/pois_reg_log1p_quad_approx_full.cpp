#include <RcppArmadillo.h>

using namespace arma;


arma::vec solve_pois_reg_log1p_quad_approx_full (
    const arma::mat X_T,
    const arma::vec y,
    const arma::uvec y_nz_idx,
    const arma::vec X_cs_times_a1,
    const arma::mat X_T_X,
    const double a2,
    arma::vec b,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const double alpha,
    const double beta
) {

  const arma::mat X_nz = X_T.cols(y_nz_idx).t();

  double first_deriv;
  double second_deriv;
  double newton_dir;
  double newton_dec;
  vec eta_nz = X_nz * b;
  vec exp_eta_nz = exp(eta_nz);
  vec exp_eta_nz_m1 = exp_eta_nz - 1;
  vec eta_nz_proposed;
  vec exp_deriv_term_nz;
  vec quad_deriv_vec;
  double t;
  double f_proposed;
  unsigned int i, j;
  double b_j_og;
  double current_lik;
  double cross_term;
  double dot_b;

  double exact_lik = -dot(
    y,
    log(exp_eta_nz_m1)
  );

  int num_indices = update_indices.size();

  for (int update_num = 1; update_num <= num_iter; update_num++) {

    for (i = 0; i < num_indices; i++) {
      j = update_indices[i];

      dot_b = dot(X_T_X.col(j), b);

      cross_term = 2 * a2 * (dot_b - X_T_X(j, j) * b[j]);

      current_lik = exact_lik + b[j] * X_cs_times_a1[j] +
        a2 * X_T_X(j, j) * (b[j] * b[j]) +
        cross_term * b[j];

      exp_deriv_term_nz = exp_eta_nz % X_nz.col(j);

      first_deriv    = -dot(
        y,
        exp_deriv_term_nz / exp_eta_nz_m1
      ) + X_cs_times_a1[j] + 2 * a2 * dot_b;

      second_deriv   = dot(
        y,
        (exp_deriv_term_nz % exp_deriv_term_nz) /
          square(exp_eta_nz_m1)
      ) + 2 * a2 * X_T_X(j, j);

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
        eta_nz_proposed     = eta_nz + (b[j] - b_j_og) * X_nz.col(j);
        exp_eta_nz = exp(eta_nz_proposed);
        exp_eta_nz_m1 = exp_eta_nz - 1;

        exact_lik =  -dot(
          y,
          log(exp_eta_nz_m1)
        );

        f_proposed = exact_lik + b[j] * X_cs_times_a1[j] +
          a2 * X_T_X(j, j) * (b[j] * b[j]) +
          cross_term * b[j];

        if (f_proposed <= current_lik - t*newton_dec) {
          eta_nz = eta_nz_proposed;
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
arma::mat regress_cols_of_Y_on_X_log1p_quad_approx_full(
    const arma::mat X_T,
    Rcpp::List Y,
    Rcpp::List Y_nz_idx,
    const arma::vec X_cs_times_a1,
    const arma::mat X_T_X,
    arma::mat& B,
    const double a2,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const double alpha,
    const double beta
) {

  for (int j = 0; j < B.n_cols; j++) {

    B.col(j) = solve_pois_reg_log1p_quad_approx_full(
      X_T,
      Y[j],
      Y_nz_idx[j],
      X_cs_times_a1,
      X_T_X,
      a2,
      B.col(j),
      update_indices,
      num_iter,
      alpha,
      beta
    );

  }

  return(B);

}

