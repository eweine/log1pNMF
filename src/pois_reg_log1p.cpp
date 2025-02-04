#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppParallel.h>
#include "ll.h"
#include "utils.h"

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;


arma::vec solve_pois_reg_log1p (
    const arma::mat& X,
    const arma::vec y,
    const arma::uvec y_nz_idx,
    const arma::vec& s,
    arma::vec b,
    const std::vector<int>& update_indices,
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
  vec exp_eta_nz_m1 = exp_eta.elem(y_nz_idx) - s;
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

      } else if (b[j] >= 1e-12) {

        t = std::min((b[j] - 1e-12) / newton_dir, 1.0);

      } else {

        continue;

      }

      newton_dec    = alpha * first_deriv * newton_dir;
      b_j_og        = b[j];
      while (true) {
        b[j]             = b_j_og - t * newton_dir;
        eta_proposed     = eta + (b[j] - b_j_og) * X.col(j);
        exp_eta = exp(eta_proposed);
        exp_eta_nz_m1 = exp_eta.elem(y_nz_idx) - s;
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


struct RegressColsWorker : public Worker
{
  // Read-only inputs
  const arma::mat          &X;
  const std::vector<arma::vec> &Y;
  const std::vector<arma::uvec> &Y_nz_idx;
  const arma::vec          &s;
  const bool                common_size_factor;
  const std::vector<int>   &update_indices;
  const unsigned int        num_iter;
  const double              alpha;
  const double              beta;

  // The matrix B we need to modify in parallel (each thread touches distinct columns)
  arma::mat                &B;

  // Constructor
  RegressColsWorker(
    const arma::mat &X_,
    const std::vector<arma::vec> &Y_,
    const std::vector<arma::uvec> &Y_nz_idx_,
    const arma::vec &s_,
    const bool common_size_factor_,
    arma::mat &B_,
    const std::vector<int> &update_indices_,
    unsigned int num_iter_,
    const double alpha_,
    const double beta_
  )
    : X(X_), Y(Y_), Y_nz_idx(Y_nz_idx_), s(s_),
      common_size_factor(common_size_factor_),
      update_indices(update_indices_),
      num_iter(num_iter_),
      alpha(alpha_), beta(beta_),
      B(B_)
  {}

  // operator() processes columns j in [begin, end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++) {
      // If common_size_factor == true, build s_j as a constant vector
      // Else, s_j = s.elem(Y_nz_idx[j]).
      if (common_size_factor) {
        arma::vec s_j(Y[j].n_elem, arma::fill::value(s[j]));
        B.col(j) = solve_pois_reg_log1p(
          X,
          Y[j],
          Y_nz_idx[j],
          s_j,
          B.col(j),   // initial B.col(j)
          update_indices,
          num_iter,
          alpha,
          beta
        );
      } else {
        // s_j is s indexed by Y_nz_idx[j]
        arma::vec s_j = s.elem(Y_nz_idx[j]);
        B.col(j) = solve_pois_reg_log1p(
          X,
          Y[j],
          Y_nz_idx[j],
          s_j,
          B.col(j),  // initial B.col(j)
          update_indices,
          num_iter,
          alpha,
          beta
        );
      }
    }
  }

  // No join() needed since each j is distinct.
  // We do not accumulate any shared sum, only writing into distinct columns.
  void join(const RegressColsWorker & /*rhs*/) {}
};


arma::mat regress_cols_of_Y_on_X_log1p_pois_exact_parallel(
    const arma::mat& X,
    const std::vector<arma::vec>& Y,
    const std::vector<arma::uvec>& Y_nz_idx,
    const arma::vec& s,
    const bool common_size_factor,
    arma::mat& B,
    const std::vector<int>& update_indices,
    unsigned int num_iter,
    const double alpha,
    const double beta
) {
  // Build the worker
  RegressColsWorker worker(
      X, Y, Y_nz_idx, s,
      common_size_factor,
      B,               // reference to B
      update_indices,
      num_iter,
      alpha,
      beta
  );

  // Parallelize over columns j in [0, B.n_cols)
  RcppParallel::parallelFor(0, B.n_cols, worker);

  // Return the updated B
  return B;
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
    const std::vector<int>& update_indices
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

  Rprintf("Fitting log1p factor model to %i x %i count matrix.\n",n,p);

  for (int iter = 0; iter < max_iter; iter++) {

    Rprintf("Iteration %i: objective = %+0.12e\n", iter, loglik);

    U_T = regress_cols_of_Y_on_X_log1p_pois_exact_parallel(
      V_T.t(),
      y_rows_data,
      y_rows_idx,
      s,
      true,
      U_T,
      update_indices,
      num_ccd_iter,
      alpha,
      beta
    );

    V_T = regress_cols_of_Y_on_X_log1p_pois_exact_parallel(
      U_T.t(),
      y_cols_data,
      y_cols_idx,
      s,
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
      n,
      p
    );

    if (loglik - prev_lik < 1e-8) {

      break;

    } else {

      prev_lik = loglik;

    }

  }

  List fit;

  fit["U"] = U_T.t();
  fit["V"] = V_T.t();

  return(fit);

}

