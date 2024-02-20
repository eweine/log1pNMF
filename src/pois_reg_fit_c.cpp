#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List solve_pois_reg_log1p_quad_approx_sparse_vec_s_fit_c (
    const arma::mat X_T,
    const arma::vec y,
    const arma::uvec y_nz_idx,
    const double sum_y,
    const arma::vec s_nz,
    const double sum_s,
    const double sum_s_0,
    const arma::vec X_T_s,
    const arma::mat X_T_diag_s_X,
    const double a0,
    const double a1,
    const double a2,
    arma::vec b,
    double c,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const double alpha,
    const double beta
) {

  const arma::mat X_T_nz = X_T.cols(y_nz_idx);
  const arma::mat X_nz = X_T_nz.t();
  const arma::vec X_T_0_s = X_T_s - X_T_nz * s_nz;
  const arma::vec X_T_0_s_times_a1 = a1 * X_T_0_s;
  const arma::mat X_0_T_diag_s_X_0 = X_T_diag_s_X -
    (X_T_nz.each_row() % s_nz.t()) * X_nz;

  double first_deriv;
  double second_deriv;
  double newton_dir;
  double newton_dec;
  vec eta_nz = X_nz * b;
  vec exp_eta_nz = exp(eta_nz);
  vec exp_eta_nz_m1 = exp_eta_nz - 1;
  vec eta_nz_proposed;
  vec exp_eta_nz_proposed;
  vec exp_eta_nz_m1_proposed;
  vec exp_deriv_term_nz;
  vec quad_deriv_vec;
  double s_nz_dot_exp_eta_nz = dot(s_nz, exp_eta_nz);
  double s_nz_dot_exp_eta_nz_proposed;
  double t;
  double f_proposed;
  unsigned int i, j;
  double b_j_og;
  double current_lik;
  double cross_term;
  double dot_b;
  double c_hat;

  double exact_lik;

  int num_indices = update_indices.size();

  for (int update_num = 1; update_num <= num_iter; update_num++) {

    exact_lik = c * s_nz_dot_exp_eta_nz - dot(
      y,
      log(exp_eta_nz_m1)
    );

    for (i = 0; i < num_indices; i++) {
      j = update_indices[i];

      dot_b = dot(X_0_T_diag_s_X_0.col(j), b);

      cross_term = 2 * a2 * (dot_b - X_0_T_diag_s_X_0(j, j) * b[j]);

      current_lik = exact_lik + b[j] * X_T_0_s_times_a1[j] * c +
        c * a2 * X_0_T_diag_s_X_0(j, j) * (b[j] * b[j]) +
        c * cross_term * b[j];

      exp_deriv_term_nz = exp_eta_nz % X_nz.col(j);

      first_deriv    = c * dot(s_nz, exp_deriv_term_nz) - dot(
        y,
        exp_deriv_term_nz / exp_eta_nz_m1
      ) + c * X_T_0_s_times_a1[j] + 2 * c * a2 * dot_b;

      second_deriv   = c * dot(s_nz, exp_deriv_term_nz % X_nz.col(j)) + dot(
        y,
        (exp_deriv_term_nz % exp_deriv_term_nz) /
          square(exp_eta_nz_m1)
      ) + 2 * c * a2 * X_0_T_diag_s_X_0(j, j);
      newton_dir     = first_deriv / second_deriv;

      // I need to handle the non-negativity constraint here
      if (newton_dir < 0) {

        t = 1.0;

      } else if (b[j] >= 1e-12) {

        t = std::min((b[j] - 1e-12) / newton_dir, 1.0);

      } else {

        continue;

      }

      newton_dec    = alpha * first_deriv * newton_dir;
      b_j_og        = b[j];
      while (true) {
        b[j]             = b_j_og - t * newton_dir;
        eta_nz_proposed     = eta_nz + (b[j] - b_j_og) * X_nz.col(j);
        exp_eta_nz_proposed = exp(eta_nz_proposed);
        exp_eta_nz_m1_proposed = exp_eta_nz_proposed - 1;
        s_nz_dot_exp_eta_nz_proposed = dot(s_nz, exp_eta_nz_proposed);

        exact_lik = c * s_nz_dot_exp_eta_nz_proposed - dot(
          y,
          log(exp_eta_nz_m1_proposed)
        );

        f_proposed = exact_lik + b[j] * X_T_0_s_times_a1[j] * c +
          c * a2 * X_0_T_diag_s_X_0(j, j) * (b[j] * b[j]) +
          c * cross_term * b[j];

        if (f_proposed <= current_lik - t*newton_dec) {
          eta_nz = eta_nz_proposed;
          exp_eta_nz = exp_eta_nz_proposed;
          exp_eta_nz_m1 = exp_eta_nz_m1_proposed;
          s_nz_dot_exp_eta_nz = s_nz_dot_exp_eta_nz_proposed;
          break;
        } else {
          t *= beta;

          if (t < 1e-12) {

            break;

          }

        }
      }
    }

    // This loop is called after each ccd iteration
    // Now I should be able to fit c

    double sqrd_term = sum(b.t() * X_0_T_diag_s_X_0 * b);

    c_hat = sum_y / (
      s_nz_dot_exp_eta_nz +
        a0 * sum_s_0 +
        a1 * dot(X_T_0_s, b) +
        a2 * sqrd_term -
        sum_s
    );

    c = std::max(1.0, c_hat);

  }

  List fit;

  fit["b"] = b;
  fit["c"] = c;

  return(fit);

}
