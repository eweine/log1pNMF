#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <chrono>
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
    arma::vec b,
    const std::vector<int>& update_indices,
    int num_iter,
    const double alpha,
    const double beta,
    double& col_obj
) {

  double first_deriv;
  double second_deriv;
  double newton_dir;
  double newton_dec;
  vec eta = X * b;
  vec exp_eta = exp(eta);
  // LF is the linear predictor WITHOUT the fixed size-factor offset (col 0). Only
  // its nonzeros are needed (the sparse rate/log term), so compute lf_nz directly
  // from the nonzero ROWS of X with the offset coefficient zeroed: a matrix-VECTOR
  // product over nnz rows -- never the dense n*p factor product, and with no
  // cancellation (the offset column is multiplied by exactly 0, not subtracted).
  vec b_no_offset = b;
  b_no_offset[0] = 0.0;
  vec lf_nz = X.rows(y_nz_idx) * b_no_offset;
  // rate at the nonzeros = exp(eta) - s = s*expm1(LF): stable for small LF, with
  // no cancellation and no underflow to log(0) = -Inf
  vec rate_nz = s % expm1(lf_nz);
  vec rate_nz_proposed;
  vec lf_nz_proposed;
  vec Xj_nz;
  uvec jcol(1);
  vec eta_proposed;
  vec exp_eta_proposed;
  vec exp_deriv_term;
  double t;
  double f_proposed;
  int i, j;
  double b_j_og;

  double current_lik = sum(exp_eta) - dot(
    y,
    log(rate_nz)
  );

  int num_indices = update_indices.size();

  for (int update_num = 1; update_num <= num_iter; update_num++) {

    for (i = 0; i < num_indices; i++) {

      j = update_indices[i];

      // nonzero rows of X.col(j) (O(nnz), no full-column copy), reused across the
      // line-search trials to update lf_nz
      jcol[0] = j;
      Xj_nz = X.submat(y_nz_idx, jcol);

      exp_deriv_term = exp_eta % X.col(j);

      first_deriv    = sum(exp_deriv_term) - dot(
        y,
        exp_deriv_term.elem(y_nz_idx) / rate_nz
      );

      second_deriv   = dot(exp_deriv_term, X.col(j)) + dot(
        y,
        (exp_deriv_term.elem(y_nz_idx) % exp_deriv_term.elem(y_nz_idx)) /
          square(rate_nz)
      );

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
        eta_proposed     = eta + (b[j] - b_j_og) * X.col(j);
        exp_eta_proposed = exp(eta_proposed);
        // col 0 (the offset) is never updated, so a change in b[j] shifts LF at the
        // nonzeros by the same increment; update only those (O(nnz), not O(n))
        lf_nz_proposed   = lf_nz + (b[j] - b_j_og) * Xj_nz;
        rate_nz_proposed = s % expm1(lf_nz_proposed);
        f_proposed = sum(exp_eta_proposed) - dot(
          y,
          log(rate_nz_proposed)
        );
        if (std::isfinite(f_proposed) && f_proposed <= current_lik - t*newton_dec) {
          eta = eta_proposed;
          exp_eta = exp_eta_proposed;
          lf_nz = lf_nz_proposed;
          rate_nz = rate_nz_proposed;
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

  // current_lik is the full negative-loglik contribution of this column at the
  // final b (dense exp sum plus the sparse log term). It is committed only on
  // an accepted line-search step (see above), so it is always in sync with the
  // returned b and can be handed back without any recomputation. Summed across
  // columns and negated, this equals get_loglik_exact(U_T, V_T).
  col_obj = current_lik;

  return(b);

}

// Y is an nxm matrix (each col is an n-dim data vec)
// X is an nxp matrix (each row is a p-dim covariate)
// B is a pxm matrix (each col is a p-dim reg coef)
// As above, but also accumulates the model objective. Each column solve
// returns its full negative-loglik contribution; summed over all columns and
// negated (obj_out) this equals get_loglik_exact at the current (U_T, V_T),
// so a caller running this as the final pass of an iteration avoids a separate
// O(n*p*K) log-likelihood evaluation.
arma::mat regress_cols_of_Y_on_X_log1p_pois_exact(
    const arma::mat& X,
    const std::vector<arma::vec>& Y,
    const std::vector<arma::uvec>& Y_nz_idx,
    const arma::vec& s,
    const bool common_size_factor,
    arma::mat& B,
    const std::vector<int>& update_indices,
    int num_iter,
    const double alpha,
    const double beta,
    double& obj_out
) {

  double obj_sum = 0.0;

  if (common_size_factor) {

    #pragma omp parallel for shared(B) reduction(+:obj_sum)
    for (int j = 0; j < static_cast<int>(B.n_cols); j++) {

      arma::vec s_j(Y[j].n_elem);
      s_j.fill(s[j]);

      double col_obj = 0.0;

      B.col(j) = solve_pois_reg_log1p (
        X,
        Y[j],
        Y_nz_idx[j],
        s_j,
        B.col(j),
        update_indices,
        num_iter,
        alpha,
        beta,
        col_obj
      );

      obj_sum += col_obj;

    }

  } else {

    #pragma omp parallel for shared(B) reduction(+:obj_sum)
    for (int j = 0; j < static_cast<int>(B.n_cols); j++) {

      double col_obj = 0.0;

      B.col(j) = solve_pois_reg_log1p (
        X,
        Y[j],
        Y_nz_idx[j],
        s.elem(Y_nz_idx[j]),
        B.col(j),
        update_indices,
        num_iter,
        alpha,
        beta,
        col_obj
      );

      obj_sum += col_obj;

    }

  }

  // negate: solver returns the negative-loglik contribution per column
  obj_out = -obj_sum;

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
    const double tol,
    const bool track_time,
    const double max_time
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
    n,
    p
  );

  double prev_lik = loglik;
  int K = U_T.n_rows - 1;

  std::vector<double> loglik_history;
  loglik_history.reserve(max_iter);
  loglik_history.push_back(loglik);

  // optional cumulative wall-clock (seconds) per iteration for benchmarking,
  // excluding the progress bar; off by default. The exact objective is already
  // the model log-likelihood here, so no separate exact-loglik trace is needed.
  std::vector<double> time_history;
  if (track_time) time_history.reserve(max_iter);
  long long algo_ns = 0;

  // seed at the initial iterate (index 0) so time_trace aligns with
  // objective_trace, which also records the pre-loop value
  if (track_time) time_history.push_back(0.0);

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

    // always time (cheap) so max_time works even when track_time is off
    std::chrono::steady_clock::time_point iter_t0 =
      std::chrono::steady_clock::now();

    double u_obj_unused = 0.0;
    U_T = regress_cols_of_Y_on_X_log1p_pois_exact(
      V_T.t(),
      y_rows_data,
      y_rows_idx,
      s,
      true,
      U_T,
      update_indices,
      num_ccd_iter,
      alpha,
      beta,
      u_obj_unused
    );

    // the V-pass is the last block update of the iteration, so its accumulated
    // objective equals get_loglik_exact at the current (U_T, V_T). The rescale
    // below is objective-invariant (it leaves every U_T.col(i)'*V_T.col(j)
    // unchanged), so this value is also the objective at the rescaled iterate,
    // and no separate O(n*p*K) evaluation is needed.
    V_T = regress_cols_of_Y_on_X_log1p_pois_exact(
      U_T.t(),
      y_cols_data,
      y_cols_idx,
      s,
      false,
      V_T,
      update_indices,
      num_ccd_iter,
      alpha,
      beta,
      loglik
    );

    arma::vec d = mean(U_T, 1) / mean(V_T, 1);
    d(0) = 1.0;

    U_T.each_col() %= arma::sqrt(1/d);
    V_T.each_col() %= arma::sqrt(d);

    algo_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::steady_clock::now() - iter_t0
    ).count();
    if (track_time) time_history.push_back(static_cast<double>(algo_ns) * 1e-9);

    loglik_history.push_back(loglik);

    // stop once the accumulated optimization time reaches max_time
    if (static_cast<double>(algo_ns) * 1e-9 >= max_time) break;

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
    _["objective_trace"] = loglik_history,
    _["time_trace"] = time_history
  );

  return(fit);

}

