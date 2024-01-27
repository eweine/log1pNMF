#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>
#include "ll.h"
#include "utils.h"

using namespace Rcpp;
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
arma::mat regress_cols_of_Y_on_X_log1p_quad_approx_full_cpp(
    const arma::mat X_T,
    const std::vector<arma::vec> Y,
    const std::vector<arma::uvec> Y_nz_idx,
    const arma::vec X_cs_times_a1,
    const arma::mat X_T_X,
    arma::mat B,
    const double a2,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const double alpha,
    const double beta
) {

  #pragma omp parallel for
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


// [[Rcpp::export]]
List fit_factor_model_log1p_quad_approx_full_cpp_src(
  const std::vector<int> sc_x,
  const std::vector<int> sc_i,
  const std::vector<int> sc_j,
  const std::vector<int> sc_T_x,
  const std::vector<int> sc_T_i,
  const std::vector<int> sc_T_j,
  arma::mat U_T,
  arma::mat V_T,
  const double a1,
  const double a2,
  const int n,
  const int p,
  const int max_iter,
  const double alpha,
  const double beta,
  const int num_ccd_iter,
  const std::vector<int> update_indices
) {

  const std::vector<int> col_num_repeats = get_num_repeats_cpp(
      sc_j,
      p,
      sc_j.size()
  );

  const std::vector<arma::vec> y_cols_data = create_vals_vector_arma(
      p,
      col_num_repeats,
      sc_x
  );

  const std::vector<arma::uvec> y_cols_idx = create_vals_vector_uvec(
      p,
      col_num_repeats,
      sc_i
  );

  const std::vector<int> row_num_repeats = get_num_repeats_cpp(
      sc_T_j,
      n,
      sc_T_j.size()
  );

  const std::vector<arma::vec> y_rows_data = create_vals_vector_arma(
      n,
      row_num_repeats,
      sc_T_x
  );

  const std::vector<arma::uvec> y_rows_idx = create_vals_vector_uvec(
      n,
      row_num_repeats,
      sc_T_i
  );

  arma::vec U_cs = arma::sum(U_T, 1);
  arma::mat U_T_U = U_T * U_T.t();

  double loglik = get_loglik_quad_approx_full(
    U_T,
    V_T,
    U_cs,
    U_T_U,
    sc_x,
    sc_i,
    sc_j,
    a1,
    a2
  );

  std::vector<double> loglik_history;
  loglik_history.push_back(loglik);

  Rprintf("Fitting log1p factor model to %i x %i count matrix.\n",n,p);

  for (int iter = 0; iter < max_iter; iter++) {

    Rprintf("Iteration %i: objective = %+0.12e\n", iter, loglik);

    U_T = regress_cols_of_Y_on_X_log1p_quad_approx_full_cpp(
      V_T,
      y_rows_data,
      y_rows_idx,
      a1 * sum(V_T, 1),
      V_T * V_T.t(),
      U_T,
      a2,
      update_indices,
      num_ccd_iter,
      alpha,
      beta
    );

    U_cs = arma::sum(U_T, 1);
    U_T_U = U_T * U_T.t();

    V_T = regress_cols_of_Y_on_X_log1p_quad_approx_full_cpp(
      U_T,
      y_cols_data,
      y_cols_idx,
      a1 * U_cs,
      U_T_U,
      V_T,
      a2,
      update_indices,
      num_ccd_iter,
      alpha,
      beta
    );

    loglik = get_loglik_quad_approx_full(
      U_T,
      V_T,
      U_cs,
      U_T_U,
      sc_x,
      sc_i,
      sc_j,
      a1,
      a2
    );

    loglik_history.push_back(loglik);

  }

  List fit;

  fit["U"] = U_T.t();
  fit["V"] = V_T.t();
  fit["loglik"] = loglik_history;

  return(fit);

}


