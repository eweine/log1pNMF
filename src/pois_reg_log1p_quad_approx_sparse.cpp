#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "ll.h"
#include "utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;


arma::vec solve_pois_reg_log1p_quad_approx_sparse_scalar_s (
    const arma::mat& X_T,
    const arma::vec y,
    const arma::uvec y_nz_idx,
    const double s,
    const arma::vec& X_cs,
    const arma::mat& X_T_X,
    const double a1,
    const double a2,
    arma::vec b,
    const std::vector<int>& update_indices,
    int num_iter,
    const double alpha,
    const double beta
) {

  const arma::mat X_T_nz = X_T.cols(y_nz_idx);
  const arma::mat X_nz = X_T_nz.t();
  const arma::vec X_0_cs_times_sa1 = s * a1 * (X_cs - arma::sum(X_T_nz, 1));
  const arma::mat X_0_T_X_0 = X_T_X - X_T_nz * X_nz;

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
  double t;
  double f_proposed;
  int i, j;
  double b_j_og;
  double current_lik;
  double cross_term;
  double dot_b;

  double exact_lik = s * sum(exp_eta_nz) - dot(
    y,
    log(exp_eta_nz_m1)
  );

  int num_indices = update_indices.size();

  for (int update_num = 1; update_num <= num_iter; update_num++) {

    for (i = 0; i < num_indices; i++) {
      j = update_indices[i];

      dot_b = dot(X_0_T_X_0.col(j), b);

      cross_term = 2 * a2 * (dot_b - X_0_T_X_0(j, j) * b[j]);

      current_lik = exact_lik + b[j] * X_0_cs_times_sa1[j] +
        s * a2 * X_0_T_X_0(j, j) * (b[j] * b[j]) +
        s * cross_term * b[j];

      exp_deriv_term_nz = exp_eta_nz % X_nz.col(j);

      first_deriv    = s * sum(exp_deriv_term_nz) - dot(
        y,
        exp_deriv_term_nz / exp_eta_nz_m1
      ) + X_0_cs_times_sa1[j] + 2 * a2 * s * dot_b;

      second_deriv   = s * dot(exp_deriv_term_nz, X_nz.col(j)) + dot(
        y,
        (exp_deriv_term_nz % exp_deriv_term_nz) /
          square(exp_eta_nz_m1)
      ) + 2 * a2 * s * X_0_T_X_0(j, j);

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

        exact_lik = s * sum(exp_eta_nz_proposed) - dot(
          y,
          log(exp_eta_nz_m1_proposed)
        );

        f_proposed = exact_lik + b[j] * X_0_cs_times_sa1[j] +
          s * a2 * X_0_T_X_0(j, j) * (b[j] * b[j]) +
          s * cross_term * b[j];

        if (f_proposed <= current_lik - t*newton_dec) {
          eta_nz = eta_nz_proposed;
          eta_nz = eta_nz_proposed;
          exp_eta_nz = exp_eta_nz_proposed;
          exp_eta_nz_m1 = exp_eta_nz_m1_proposed;
          break;
        } else {
          t *= beta;
          if(t < 1e-12) {

            b[j] = b_j_og;
            break;

          }
        }
      }
    }
  }

  return(b);

}


arma::vec solve_pois_reg_log1p_quad_approx_sparse_vec_s (
    const arma::mat& X_T,
    const arma::vec y,
    const arma::uvec y_nz_idx,
    const arma::vec s_nz,
    const arma::vec& X_T_s,
    const arma::mat& X_T_diag_s_X,
    const double a1,
    const double a2,
    arma::vec b,
    const std::vector<int>& update_indices,
    int num_iter,
    const double alpha,
    const double beta
) {

  const arma::mat X_T_nz = X_T.cols(y_nz_idx);
  const arma::mat X_nz = X_T_nz.t();
  const arma::vec X_T_0_s_times_a1 = a1 * (X_T_s - X_T_nz * s_nz);
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
  double t;
  double f_proposed;
  int i, j;
  double b_j_og;
  double current_lik;
  double cross_term;
  double dot_b;

  double exact_lik = dot(s_nz, exp_eta_nz) - dot(
    y,
    log(exp_eta_nz_m1)
  );

  int num_indices = update_indices.size();

  for (int update_num = 1; update_num <= num_iter; update_num++) {

    for (i = 0; i < num_indices; i++) {
      j = update_indices[i];

      dot_b = dot(X_0_T_diag_s_X_0.col(j), b);

      cross_term = 2 * a2 * (dot_b - X_0_T_diag_s_X_0(j, j) * b[j]);

      current_lik = exact_lik + b[j] * X_T_0_s_times_a1[j] +
        a2 * X_0_T_diag_s_X_0(j, j) * (b[j] * b[j]) +
        cross_term * b[j];

      exp_deriv_term_nz = exp_eta_nz % X_nz.col(j);

      first_deriv    = dot(s_nz, exp_deriv_term_nz) - dot(
        y,
        exp_deriv_term_nz / exp_eta_nz_m1
      ) + X_T_0_s_times_a1[j] + 2 * a2 * dot_b;

      second_deriv   = dot(s_nz, exp_deriv_term_nz % X_nz.col(j)) + dot(
        y,
        (exp_deriv_term_nz % exp_deriv_term_nz) /
          square(exp_eta_nz_m1)
      ) + 2 * a2 * X_0_T_diag_s_X_0(j, j);
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

        exact_lik = dot(s_nz, exp_eta_nz_proposed) - dot(
          y,
          log(exp_eta_nz_m1_proposed)
        );

        f_proposed = exact_lik + b[j] * X_T_0_s_times_a1[j] +
          a2 * X_0_T_diag_s_X_0(j, j) * (b[j] * b[j]) +
          cross_term * b[j];

        if (f_proposed <= current_lik - t*newton_dec) {
          eta_nz = eta_nz_proposed;
          exp_eta_nz = exp_eta_nz_proposed;
          exp_eta_nz_m1 = exp_eta_nz_m1_proposed;
          break;
        } else {
          t *= beta;

          if (t < 1e-12) {

            b[j] = b_j_og;
            break;

          }

        }
      }
    }
  }

  return(b);

}

// Y is an nxm matrix (each col is an n-dim data vec)
// X is an nxp matrix (each row is a p-dim covariate)
// B is a pxm matrix (each col is a p-dim reg coef)
arma::mat regress_cols_of_Y_on_X_log1p_quad_approx_sparse_vec_s(
    const arma::mat& X_T,
    const std::vector<arma::vec>& Y,
    const std::vector<arma::uvec>& Y_nz_idx,
    const arma::vec& s,
    arma::mat& B,
    const double a1,
    const double a2,
    const std::vector<int>& update_indices,
    int num_iter,
    const double alpha,
    const double beta
) {

  const arma::mat X_T_s = X_T * s;
  const arma::mat X_T_diag_s_X = (X_T.each_row() % s.t()) * X_T.t();

  #pragma omp parallel for shared(B)
  for (int j = 0; j < static_cast<int>(B.n_cols); j++) {

    B.col(j) = solve_pois_reg_log1p_quad_approx_sparse_vec_s(
      X_T,
      Y[j],
      Y_nz_idx[j],
      s.elem(Y_nz_idx[j]),
      X_T_s,
      X_T_diag_s_X,
      a1,
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

arma::mat regress_cols_of_Y_on_X_log1p_quad_approx_sparse_scalar_s(
    const arma::mat& X_T,
    const std::vector<arma::vec>& Y,
    const std::vector<arma::uvec>& Y_nz_idx,
    const arma::vec& s,
    arma::mat& B,
    const double a1,
    const double a2,
    const std::vector<int>& update_indices,
    int num_iter,
    const double alpha,
    const double beta
) {

  const arma::mat X_cs = arma::sum(X_T, 1);
  const arma::mat X_T_X = X_T * X_T.t();

  // Commenting out parallelism for testing
  #pragma omp parallel for shared(B)
  for (int j = 0; j < static_cast<int>(B.n_cols); j++) {

    B.col(j) = solve_pois_reg_log1p_quad_approx_sparse_scalar_s(
      X_T,
      Y[j],
      Y_nz_idx[j],
      s(j),
      X_cs,
      X_T_X,
      a1,
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
List fit_factor_model_log1p_quad_approx_sparse_cpp_src(
    const std::vector<int>& sc_x,
    const std::vector<int>& sc_i,
    const std::vector<int>& sc_j,
    const std::vector<int>& sc_T_x,
    const std::vector<int>& sc_T_i,
    const std::vector<int>& sc_T_j,
    const arma::vec& s,
    arma::mat& U_T,
    arma::mat& V_T,
    const double a1,
    const double a2,
    const int n,
    const int p,
    const int max_iter,
    const double alpha,
    const double beta,
    const int num_ccd_iter,
    const std::vector<int>& update_indices,
    const bool verbose,
    const double tol
) {

  bool converged = false;

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

  double loglik = get_loglik_quad_approx_sparse(
    U_T,
    V_T,
    sc_x,
    sc_i,
    sc_j,
    s,
    a1,
    a2
  );

  double prev_lik = loglik;
  int K = U_T.n_rows;

  std::vector<double> loglik_history;
  loglik_history.reserve(max_iter);
  loglik_history.push_back(loglik);

  if (verbose) {

    Rprintf(
      "Optimizing rank %i log1p factor model to %i x %i count matrix.\n",K,n,p
    );

  }

  int barWidth = 50; // Adjust as you like

  for (int iter = 0; iter < max_iter; iter++) {

    if (verbose) {
      double fraction = (double)iter / (double)max_iter;
      int filled = (int)(barWidth * fraction);

      // Move to start of line
      Rprintf("\r[");

      // Fill progress
      for (int i = 0; i < barWidth; i++) {
        if (i < filled) {
          Rprintf("=");
        } else if (i == filled) {
          Rprintf(">");
        } else {
          Rprintf(" ");
        }
      }
      // Print percentage
      Rprintf("] %d%%", (int)(fraction * 100));

      // Flush to make sure it shows immediately
      R_FlushConsole();
    }
    // --------------------------------------

    arma::vec d = mean(U_T, 1) / mean(V_T, 1);

    U_T.each_col() %= arma::sqrt(1/d);
    V_T.each_col() %= arma::sqrt(d);

    V_T = regress_cols_of_Y_on_X_log1p_quad_approx_sparse_vec_s(
      U_T,
      y_cols_data,
      y_cols_idx,
      s,
      V_T,
      a1,
      a2,
      update_indices,
      num_ccd_iter,
      alpha,
      beta
    );

    U_T = regress_cols_of_Y_on_X_log1p_quad_approx_sparse_scalar_s(
      V_T,
      y_rows_data,
      y_rows_idx,
      s,
      U_T,
      a1,
      a2,
      update_indices,
      num_ccd_iter,
      alpha,
      beta
    );

    loglik = get_loglik_quad_approx_sparse(
      U_T,
      V_T,
      sc_x,
      sc_i,
      sc_j,
      s,
      a1,
      a2
    );

    loglik_history.push_back(loglik);

    if (loglik - prev_lik < tol) {

      converged = true;
      break;

    } else {

      prev_lik = loglik;

    }

  }

  if (verbose) {

    // Finalize the progress bar after the loop
    if (verbose) {
      // Move to start of line and overwrite
      Rprintf("\r[");
      for (int i = 0; i < barWidth; i++) {
        Rprintf("=");
      }
      Rprintf("] 100%%\n");
      R_FlushConsole();
    }

  }

  List fit = List::create(
    _["U"]         = U_T.t(),
    _["V"]         = V_T.t(),
    _["converged"] = converged,
    _["objective_trace"] = loglik_history
  );

  return(fit);

}
