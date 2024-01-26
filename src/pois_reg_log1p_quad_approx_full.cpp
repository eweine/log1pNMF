#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;
using namespace RcppParallel;


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


// Ideally, this is the code that could be made to run in parallel
// I'm not sure how hard this will be
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

struct RegressColsOfYOnXlog1pQuadFull : public Worker {
  const arma::mat X_T;
  const std::vector<arma::vec> Y;
  const std::vector<arma::uvec> Y_nz_idx;
  const arma::vec X_cs_times_a1;
  const arma::mat X_T_X;
  arma::mat B;
  arma::mat B_out;
  const double a2;
  const std::vector<int> update_indices;
  unsigned int num_iter;
  const double alpha;
  const double beta;

  RegressColsOfYOnXlog1pQuadFull(
    const arma::mat X_T,
    const std::vector<arma::vec> Y,
    const std::vector<arma::uvec> Y_nz_idx,
    const arma::vec X_cs_times_a1,
    const arma::mat X_T_X,
    arma::mat B,
    arma::mat B_out,
    const double a2,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const double alpha,
    const double beta) : X_T(X_T),
    Y(Y), Y_nz_idx(Y_nz_idx), X_cs_times_a1(X_cs_times_a1), X_T_X(X_T_X),
    B(B), B_out(B_out), a2(a2), update_indices(update_indices),
    num_iter(num_iter), alpha(alpha), beta(beta) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++) {
      B_out.col(j) = solve_pois_reg_log1p_quad_approx_full (
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
  }
};

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::export]]
arma::mat regress_cols_of_Y_on_X_log1p_quad_approx_full_parallel(
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

  arma::mat B_out(B.n_rows, B.n_cols);

  RegressColsOfYOnXlog1pQuadFull updater(
      X_T,
      Y,
      Y_nz_idx,
      X_cs_times_a1,
      X_T_X,
      B,
      B_out,
      a2,
      update_indices,
      num_iter,
      alpha,
      beta
  );

  parallelFor(0, B.n_cols, updater);

  // It's not clear to me if this will contain the updates
  return B_out;

}


double get_loglik_cpp(
    const arma::mat U_T,
    const arma::mat V_T,
    const std::vector<int> y_nz_vals,
    const std::vector<int> y_nz_rows_idx,
    const std::vector<int> y_nz_cols_idx,
    const double a1,
    const double a2
) {

  double sum = 0.0;
  double lin_term;
  double quad_term;
  int i;
  int j;

  for (int r = 0; r < y_nz_vals.size(); r++) {

    i = y_nz_rows_idx[r];
    j = y_nz_cols_idx[r];

    sum += y_nz_vals[r] * log(exp(dot(U_T.col(i), V_T.col(j))) - 1);

  }

  lin_term = a1 * arma::dot(arma::sum(U_T, 1), arma::sum(V_T, 1));

  arma::mat U_T_U = U_T * U_T.t();
  arma::mat U_T_U_V_T = U_T_U * V_T;
  quad_term = a2 * arma::accu(V_T % U_T_U_V_T);

  sum = sum - lin_term - quad_term;
  return(sum);

}

// Now, I want to add the code that will run the full matrix factorization
std::vector<int> get_num_repeats_cpp(
    const std::vector<int> idx,
    const int p,
    const int total_idx
) {

  std::vector<int> idx_cts(p, 1);

  int prev_val = idx[0];
  int ctr = 0;

  for (int i = 1; i < total_idx; i++) {

    if (idx[i] != prev_val) {

      ctr += 1;
      prev_val = idx[i];

    } else {

      idx_cts[ctr] += 1;

    }

  }

  return(idx_cts);

}

std::vector<arma::vec> create_vals_vector_arma(
    const int num_vectors,
    const std::vector<int> vector_sizes,
    const std::vector<int> values
) {
  // Preallocate memory for the vector
  std::vector<arma::vec> vectorOfVectors;
  vectorOfVectors.reserve(num_vectors);

  int values_ctr = 0;

  // Create vectors in a for loop and add them to the vector
  for (int i = 0; i < num_vectors; i++) {
    // Create an Armadillo vector
    arma::vec vec(vector_sizes[i]);

    // Populate the Armadillo vector with values (e.g., incrementing values)
    for (int j = 0; j < vector_sizes[i]; j++) {
      vec(j) = values[values_ctr];
      values_ctr++;
    }

    // Add the Armadillo vector to the vector
    vectorOfVectors.push_back(vec);
  }

  return vectorOfVectors;
}

std::vector<arma::uvec> create_vals_vector_uvec(
    const int num_vectors,
    const std::vector<int> vector_sizes,
    const std::vector<int> values
) {
  // Preallocate memory for the vector
  std::vector<arma::uvec> vectorOfUvecs;
  vectorOfUvecs.reserve(num_vectors);

  int values_ctr = 0;

  // Create vectors in a for loop and add them to the vector
  for (int i = 0; i < num_vectors; i++) {
    // Create an Armadillo unsigned integer vector
    arma::uvec uvec(vector_sizes[i]);

    // Populate the Armadillo unsigned integer vector with values (e.g., incrementing values)
    for (int j = 0; j < vector_sizes[i]; j++) {
      uvec(j) = static_cast<unsigned int>(values[values_ctr]);
      values_ctr++;
    }

    // Add the Armadillo unsigned integer vector to the vector
    vectorOfUvecs.push_back(uvec);
  }

  return vectorOfUvecs;
}

// Now, I want to define a function to do the full matrix
// factorization

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

  std::vector<int> col_num_repeats = get_num_repeats_cpp(
      sc_j,
      p,
      sc_j.size()
  );

  std::vector<arma::vec> y_cols_data = create_vals_vector_arma(
      p,
      col_num_repeats,
      sc_x
  );

  std::vector<arma::uvec> y_cols_idx = create_vals_vector_uvec(
      p,
      col_num_repeats,
      sc_i
  );

  std::vector<int> row_num_repeats = get_num_repeats_cpp(
      sc_T_j,
      n,
      sc_T_j.size()
  );

  std::vector<arma::vec> y_rows_data = create_vals_vector_arma(
      n,
      row_num_repeats,
      sc_T_x
  );

  std::vector<arma::uvec> y_rows_idx = create_vals_vector_uvec(
      n,
      row_num_repeats,
      sc_T_i
  );

  double loglik = get_loglik_cpp(
    U_T,
    V_T,
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

    U_T = regress_cols_of_Y_on_X_log1p_quad_approx_full_parallel(
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

    V_T = regress_cols_of_Y_on_X_log1p_quad_approx_full_parallel(
      U_T,
      y_cols_data,
      y_cols_idx,
      a1 * sum(U_T, 1),
      U_T * U_T.t(),
      V_T,
      a2,
      update_indices,
      num_ccd_iter,
      alpha,
      beta
    );

    loglik = get_loglik_cpp(
      U_T,
      V_T,
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


// NOTE: The issue with this implementation is that it will be difficult
// to iterate through a list in parallel. Instead, I should use a vector
// of arma::uvecs or a vector of arma::vecs. Once I do this, I should be
// able to plugin the code above relatively easily.
// I can take care of this tomorrow.

