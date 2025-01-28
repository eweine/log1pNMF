#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppParallel.h>
#include "ll.h"
#include "utils.h"

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

// I'm going to create two functions here for ease of computation.
// The first function will use a scalar size factor
// The second function will use a vector size factor


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
    unsigned int num_iter,
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
  unsigned int i, j;
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
    unsigned int num_iter,
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
  unsigned int i, j;
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


// A worker that processes columns [begin, end) of B in parallel:
struct RegressColsVecSWorker : public Worker
{
  // Inputs needed
  const arma::mat            &X_T;
  const std::vector<arma::vec>   &Y;
  const std::vector<arma::uvec>  &Y_nz_idx;
  const arma::vec            &s;
  const arma::mat            &X_T_s;
  const arma::mat            &X_T_diag_s_X;
  const double                a1;
  const double                a2;
  const std::vector<int>     &update_indices;
  const unsigned int          num_iter;
  const double                alpha;
  const double                beta;

  // The matrix B to be updated
  arma::mat &B;

  // Constructor
  RegressColsVecSWorker(const arma::mat &X_T_,
                        const std::vector<arma::vec> &Y_,
                        const std::vector<arma::uvec> &Y_nz_idx_,
                        const arma::vec &s_,
                        const arma::mat &X_T_s_,
                        const arma::mat &X_T_diag_s_X_,
                        double a1_,
                        double a2_,
                        const std::vector<int> &update_indices_,
                        unsigned int num_iter_,
                        double alpha_,
                        double beta_,
                        arma::mat &B_)
    : X_T(X_T_), Y(Y_), Y_nz_idx(Y_nz_idx_),
      s(s_),
      X_T_s(X_T_s_),
      X_T_diag_s_X(X_T_diag_s_X_),
      a1(a1_), a2(a2_),
      update_indices(update_indices_),
      num_iter(num_iter_),
      alpha(alpha_), beta(beta_),
      B(B_)
  {}

  // operator(): updates columns j in [begin, end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++) {
      // Solve and store into B.col(j)
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
  }
};


arma::mat regress_cols_of_Y_on_X_log1p_quad_approx_sparse_vec_s_parallel(
    const arma::mat &X_T,
    const std::vector<arma::vec> &Y,
    const std::vector<arma::uvec> &Y_nz_idx,
    const arma::vec &s,
    arma::mat &B,
    const double a1,
    const double a2,
    const std::vector<int> &update_indices,
    unsigned int num_iter,
    const double alpha,
    const double beta
) {
  // Precompute the same stuff you did before:
  arma::mat X_T_s = X_T * s;  // (p x 1) times s is (s is n-dim, so X_T is p x n ?)
  arma::mat X_T_diag_s_X = (X_T.each_row() % s.t()) * X_T.t();

  // Build the worker
  RegressColsVecSWorker worker(
      X_T, Y, Y_nz_idx, s,
      X_T_s, X_T_diag_s_X,
      a1, a2,
      update_indices,
      num_iter,
      alpha, beta,
      B
  );

  // Parallel loop over columns [0, B.n_cols)
  RcppParallel::parallelFor(0, B.n_cols, worker);

  return B;
}

struct RegressColsScalarSWorker : public RcppParallel::Worker
{
  // Inputs
  const arma::mat            &X_T;
  const std::vector<arma::vec>   &Y;
  const std::vector<arma::uvec>  &Y_nz_idx;
  const arma::vec            &s;
  const arma::mat            &X_cs;
  const arma::mat            &X_T_X;
  const double                a1;
  const double                a2;
  const std::vector<int>     &update_indices;
  const unsigned int          num_iter;
  const double                alpha;
  const double                beta;

  // Matrix B to update
  arma::mat &B;

  // Constructor
  RegressColsScalarSWorker(const arma::mat &X_T_,
                           const std::vector<arma::vec> &Y_,
                           const std::vector<arma::uvec> &Y_nz_idx_,
                           const arma::vec &s_,
                           const arma::mat &X_cs_,
                           const arma::mat &X_T_X_,
                           double a1_,
                           double a2_,
                           const std::vector<int> &update_indices_,
                           unsigned int num_iter_,
                           double alpha_,
                           double beta_,
                           arma::mat &B_)
    : X_T(X_T_), Y(Y_), Y_nz_idx(Y_nz_idx_),
      s(s_),
      X_cs(X_cs_), X_T_X(X_T_X_),
      a1(a1_), a2(a2_),
      update_indices(update_indices_),
      num_iter(num_iter_),
      alpha(alpha_), beta(beta_),
      B(B_)
  {}

  // operator() processes columns j in [begin, end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++) {
      B.col(j) = solve_pois_reg_log1p_quad_approx_sparse_scalar_s(
        X_T,
        Y[j],
        Y_nz_idx[j],
        s(j),       // single scalar for column j
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
  }
};


arma::mat regress_cols_of_Y_on_X_log1p_quad_approx_sparse_scalar_s_parallel(
    const arma::mat &X_T,
    const std::vector<arma::vec> &Y,
    const std::vector<arma::uvec> &Y_nz_idx,
    const arma::vec &s,
    arma::mat &B,
    const double a1,
    const double a2,
    const std::vector<int> &update_indices,
    unsigned int num_iter,
    const double alpha,
    const double beta
) {
  // Precompute as in your original code
  arma::mat X_cs  = arma::sum(X_T, 1);  // sum each row
  arma::mat X_T_X = X_T * X_T.t();

  // Build the worker
  RegressColsScalarSWorker worker(
      X_T, Y, Y_nz_idx, s,
      X_cs, X_T_X,
      a1, a2,
      update_indices,
      num_iter,
      alpha, beta,
      B
  );

  // Parallel loop over [0, B.n_cols)
  RcppParallel::parallelFor(0, B.n_cols, worker);

  return B;
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

  std::vector<double> loglik_history;
  loglik_history.push_back(loglik);

  Rprintf("Fitting log1p factor model to %i x %i count matrix.\n",n,p);

  for (int iter = 0; iter < max_iter; iter++) {

    Rprintf("Iteration %i: objective = %+0.12e\n", iter, loglik);

    arma::vec d = mean(U_T, 1) / mean(V_T, 1);

    U_T.each_col() %= arma::sqrt(1/d);
    V_T.each_col() %= arma::sqrt(d);

    V_T = regress_cols_of_Y_on_X_log1p_quad_approx_sparse_vec_s_parallel(
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

    U_T = regress_cols_of_Y_on_X_log1p_quad_approx_sparse_scalar_s_parallel(
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

    if (loglik - prev_lik < 1e-8) {

      Rprintf("Log-likelihood change is very small!!!!");
      //break;

    } else {

      prev_lik = loglik;

    }

  }

  List fit;

  fit["U"] = U_T.t();
  fit["V"] = V_T.t();
  fit["loglik"] = loglik_history;

  return(fit);

}
