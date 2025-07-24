#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "ll.h"
#include "utils.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;


arma::vec solve_pois_reg_log1p (
    const arma::mat& X,
    const arma::vec y,
    const arma::uvec y_nz_idx,
    const arma::vec& s,
    const double cc_alpha,
    arma::vec b,
    const std::vector<int>& update_indices,
    int num_iter,
    const double alpha,
    const double beta
) {

  double first_deriv;
  double second_deriv;
  double newton_dir;
  double newton_dec;
  vec eta = X * b * (1 / cc_alpha);
  vec exp_eta = exp(eta);
  vec exp_eta_nz_m1 = exp_eta.elem(y_nz_idx) - s;
  vec exp_eta_nz_m1_proposed;
  vec eta_proposed;
  vec exp_eta_proposed;
  vec exp_deriv_term;
  double t;
  double f_proposed;
  int i, j;
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

      // pick up here to re-define derivatives using chain rule
      first_deriv    = (1 / cc_alpha) * (sum(exp_deriv_term) - dot(
        y,
        exp_deriv_term.elem(y_nz_idx) / exp_eta_nz_m1
      ));

      second_deriv   = (1 / (cc_alpha * cc_alpha)) * (
        dot(exp_deriv_term, X.col(j)) + dot(
        y,
        (exp_deriv_term.elem(y_nz_idx) % exp_deriv_term.elem(y_nz_idx)) /
          square(exp_eta_nz_m1)
        ));

      newton_dir     = first_deriv / second_deriv;

      // I need to handle the non-negativity constraint here
      if (std::isfinite(newton_dir)) {
        
        if (newton_dir < 0) {
          
          t = 1.0;
          
        } else if (b[j] > 1e-12) {
          
          t = std::min((b[j] - 1e-12) / newton_dir, 1.0);
          
        } else {
          
          continue;
          
        }
        
      } else{
        
        continue;
        
      }

      newton_dec    = alpha * first_deriv * newton_dir;
      b_j_og        = b[j];
      while (true) {
        
        b[j]             = std::max(b_j_og - t * newton_dir, 1e-12);
        eta_proposed     = eta + (b[j] - b_j_og) * X.col(j) * (1 / cc_alpha);
        exp_eta_proposed = exp(eta_proposed);
        exp_eta_nz_m1_proposed = exp_eta_proposed.elem(y_nz_idx) - s;
        f_proposed = sum(exp_eta_proposed) - dot(
          y,
          log(exp_eta_nz_m1_proposed)
        );
        if (std::isfinite(f_proposed) && f_proposed <= current_lik - t*newton_dec) {
          eta = eta_proposed;
          exp_eta = exp_eta_proposed;
          exp_eta_nz_m1 = exp_eta_nz_m1_proposed;
          current_lik = f_proposed;
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

// Y is an nxm matrix (each col is an n-dim data vec)
// X is an nxp matrix (each row is a p-dim covariate)
// B is a pxm matrix (each col is a p-dim reg coef)
arma::mat regress_cols_of_Y_on_X_log1p_pois_exact(
    const arma::mat& X,
    const std::vector<arma::vec>& Y,
    const std::vector<arma::uvec>& Y_nz_idx,
    const arma::vec& s,
    const double cc_alpha,
    const bool common_size_factor,
    arma::mat& B,
    const std::vector<int>& update_indices,
    int num_iter,
    const double alpha,
    const double beta
) {

  if (common_size_factor) {

    #pragma omp parallel for shared(B)
    for (int j = 0; j < static_cast<int>(B.n_cols); j++) {

      arma::vec s_j(Y[j].n_elem);
      s_j.fill(s[j]);

      B.col(j) = solve_pois_reg_log1p (
        X,
        Y[j],
        Y_nz_idx[j],
        s_j,
        cc_alpha,
        B.col(j),
        update_indices,
        num_iter,
        alpha,
        beta
      );

    }

  } else {

    #pragma omp parallel for shared(B)
    for (int j = 0; j < static_cast<int>(B.n_cols); j++) {

      B.col(j) = solve_pois_reg_log1p (
        X,
        Y[j],
        Y_nz_idx[j],
        s.elem(Y_nz_idx[j]),
        cc_alpha,
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
List fit_factor_model_log1p_exact_cpp_src(
    const std::vector<int>& sc_x,
    const std::vector<int>& sc_i,
    const std::vector<int>& sc_j,
    const std::vector<int>& sc_T_x,
    const std::vector<int>& sc_T_i,
    const std::vector<int>& sc_T_j,
    const arma::vec& s,
    const double cc_alpha,
    arma::mat& U_T,
    arma::mat& V_T,
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

  double loglik = get_loglik_exact(
    U_T,
    V_T,
    sc_x,
    sc_i,
    sc_j,
    s,
    cc_alpha,
    n,
    p
  );

  double prev_lik = loglik;
  int K = U_T.n_rows - 1;

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

    U_T = regress_cols_of_Y_on_X_log1p_pois_exact(
      V_T.t(),
      y_rows_data,
      y_rows_idx,
      s,
      cc_alpha,
      true,
      U_T,
      update_indices,
      num_ccd_iter,
      alpha,
      beta
    );

    V_T = regress_cols_of_Y_on_X_log1p_pois_exact(
      U_T.t(),
      y_cols_data,
      y_cols_idx,
      s,
      cc_alpha,
      false,
      V_T,
      update_indices,
      num_ccd_iter,
      alpha,
      beta
    );

    arma::vec d = mean(U_T, 1) / mean(V_T, 1);
    d(0) = 1.0;

    U_T.each_col() %= arma::sqrt(1/d);
    V_T.each_col() %= arma::sqrt(d);

    loglik = get_loglik_exact(
      U_T,
      V_T,
      sc_x,
      sc_i,
      sc_j,
      s,
      cc_alpha,
      n,
      p
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

