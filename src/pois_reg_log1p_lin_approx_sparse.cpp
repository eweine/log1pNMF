#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>
#include "ll.h"
#include "utils.h"

using namespace Rcpp;
using namespace arma;


arma::vec solve_pois_reg_log1p_lin_approx_sparse (
    const arma::mat X_T,
    const arma::vec y,
    const arma::uvec y_nz_idx,
    const arma::vec s,
    const arma::vec X_cs_times_a,
    const double a,
    arma::vec b,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const double alpha,
    const double beta
) {

  const arma::mat X_T_nz = X_T.cols(y_nz_idx);
  const arma::mat X_nz = X_T_nz.t();
  const arma::vec X_0_cs_times_a = X_cs_times_a - a * arma::sum(X_T_nz, 1);

  double first_deriv;
  double second_deriv;
  double newton_dir;
  double newton_dec;
  vec eta_nz = X_nz * b;
  vec exp_eta_nz = exp(eta_nz);
  vec exp_eta_nz_m1 = exp_eta_nz - s;
  vec eta_nz_proposed;
  vec exp_deriv_term_nz;
  double t;
  double f_proposed;
  unsigned int i, j;
  double b_j_og;
  double current_lik;

  double exact_lik = sum(exp_eta_nz) - dot(
    y,
    log(exp_eta_nz_m1)
  );

  int num_indices = update_indices.size();

  for (int update_num = 1; update_num <= num_iter; update_num++) {

    for (i = 0; i < num_indices; i++) {
      j = update_indices[i];

      current_lik = exact_lik + b[j] * X_0_cs_times_a[j];

      exp_deriv_term_nz = exp_eta_nz % X_nz.col(j);

      first_deriv    = sum(exp_deriv_term_nz) - dot(
        y,
        exp_deriv_term_nz / exp_eta_nz_m1
      ) + X_0_cs_times_a[j];

      second_deriv   = dot(exp_deriv_term_nz, X_nz.col(j)) + dot(
        y,
        (exp_deriv_term_nz % exp_deriv_term_nz) /
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
        eta_nz_proposed     = eta_nz + (b[j] - b_j_og) * X_nz.col(j);
        exp_eta_nz = exp(eta_nz_proposed);
        exp_eta_nz_m1 = exp_eta_nz - s;

        exact_lik = sum(exp_eta_nz) - dot(
          y,
          log(exp_eta_nz_m1)
        );

        f_proposed = exact_lik + b[j] * X_0_cs_times_a[j];

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
arma::mat regress_cols_of_Y_on_X_log1p_lin_approx_sparse(
    const arma::mat X_T,
    const std::vector<arma::vec> Y,
    const std::vector<arma::uvec> Y_nz_idx,
    const arma::vec s,
    const bool common_size_factor,
    const arma::vec X_cs_times_a,
    const double a,
    arma::mat B,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const double alpha,
    const double beta
) {

  if (common_size_factor) {

    #pragma omp parallel for
    for (int j = 0; j < B.n_cols; j++) {

      arma::vec s_j(Y[j].n_elem);
      s_j.fill(s[j]);

      B.col(j) = solve_pois_reg_log1p_lin_approx_sparse(
        X_T,
        Y[j],
        Y_nz_idx[j],
        s_j,
        X_cs_times_a,
        a,
        B.col(j),
        update_indices,
        num_iter,
        alpha,
        beta
      );

    }

  } else {

    #pragma omp parallel for
    for (int j = 0; j < B.n_cols; j++) {

      B.col(j) = solve_pois_reg_log1p_lin_approx_sparse(
        X_T,
        Y[j],
        Y_nz_idx[j],
        s.elem(Y_nz_idx[j]),
        X_cs_times_a,
        a,
        B.col(j),
        update_indices,
        num_iter,
        alpha,
        beta
      );

    }

  }

  return(B);

}


// [[Rcpp::export]]
List fit_factor_model_log1p_lin_approx_sparse_cpp_src(
    const std::vector<int> sc_x,
    const std::vector<int> sc_i,
    const std::vector<int> sc_j,
    const std::vector<int> sc_T_x,
    const std::vector<int> sc_T_i,
    const std::vector<int> sc_T_j,
    const arma::vec s,
    arma::mat U_T,
    arma::mat V_T,
    const double a,
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

  double loglik = get_loglik_lin_approx_sparse(
    U_T,
    V_T,
    U_cs,
    sc_x,
    sc_i,
    sc_j,
    s,
    a
  );

  std::vector<double> loglik_history;
  loglik_history.push_back(loglik);

  Rprintf("Fitting log1p factor model to %i x %i count matrix.\n",n,p);

  for (int iter = 0; iter < max_iter; iter++) {

    Rprintf("Iteration %i: objective = %+0.12e\n", iter, loglik);

    U_T = regress_cols_of_Y_on_X_log1p_lin_approx_sparse(
      V_T,
      y_rows_data,
      y_rows_idx,
      s,
      true,
      a * sum(V_T, 1),
      a,
      U_T,
      update_indices,
      num_ccd_iter,
      alpha,
      beta
    );

    U_cs = arma::sum(U_T, 1);

    V_T = regress_cols_of_Y_on_X_log1p_lin_approx_sparse(
      U_T,
      y_cols_data,
      y_cols_idx,
      s,
      false,
      a * U_cs,
      a,
      V_T,
      update_indices,
      num_ccd_iter,
      alpha,
      beta
    );

    loglik = get_loglik_lin_approx_sparse(
      U_T,
      V_T,
      U_cs,
      sc_x,
      sc_i,
      sc_j,
      s,
      a
    );

    loglik_history.push_back(loglik);

  }

  List fit;

  fit["U"] = U_T.t();
  fit["V"] = V_T.t();
  fit["loglik"] = loglik_history;

  return(fit);

}
