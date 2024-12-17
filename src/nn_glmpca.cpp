#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>
#include "ll.h"
#include "utils.h"

using namespace Rcpp;
using namespace arma;


arma::vec solve_pois_reg_nn_glmpca (
    const arma::mat X,
    const arma::vec m,
    arma::vec b,
    const std::vector<int> update_indices,
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
  vec eta_proposed;
  vec exp_deriv_term;
  double t;
  double f_proposed;
  unsigned int i, j;
  double b_j_og;

  double current_lik;

  int num_indices = update_indices.size();

  for (int update_num = 1; update_num <= num_iter; update_num++) {

    for (i = 0; i < num_indices; i++) {

      j = update_indices[i];

      current_lik = sum(exp_eta) - b(j) * m(j);

      exp_deriv_term = exp_eta % X.col(j);
      first_deriv    = sum(exp_deriv_term) - m(j);
      second_deriv   = dot(exp_deriv_term, X.col(j));

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

        f_proposed = sum(exp_eta) - b(j) * m(j);

        if (f_proposed <= current_lik - t*newton_dec) {
          eta = eta_proposed;
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
arma::mat regress_cols_of_Y_on_X_nn_glmpca_pois_exact(
    const arma::mat X,
    const arma::sp_mat Y,
    arma::mat B,
    const std::vector<int> update_indices,
    unsigned int num_iter,
    const double alpha,
    const double beta
) {

  arma::mat M = X.t() * Y;

  #pragma omp parallel for shared(B, M)
  for (int j = 0; j < B.n_cols; j++) {

    B.col(j) = solve_pois_reg_nn_glmpca (
      X,
      M.col(j),
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
List fit_factor_model_nn_glmpca_cpp_src(
    const arma::sp_mat Y,
    const std::vector<int> sc_x,
    const std::vector<int> sc_i,
    const std::vector<int> sc_j,
    arma::mat U_T,
    arma::mat V_T,
    const int n,
    const int p,
    const int max_iter,
    const double alpha,
    const double beta,
    const int num_ccd_iter,
    const std::vector<int> update_indices
) {

  // FIX THIS FOR GLMPCA
  double loglik = get_loglik_nn_glmpca(
    U_T,
    V_T,
    sc_x,
    sc_i,
    sc_j,
    n,
    p
  );

  double prev_lik = loglik;
  arma::sp_mat Y_T = Y.t();

  std::vector<double> loglik_history;
  loglik_history.push_back(loglik);

  Rprintf("Fitting nn GLMPCA factor model to %i x %i count matrix.\n",n,p);

  for (int iter = 0; iter < max_iter; iter++) {

    Rprintf("Iteration %i: objective = %+0.12e\n", iter, loglik);

    U_T = regress_cols_of_Y_on_X_nn_glmpca_pois_exact(
      V_T.t(),
      Y_T,
      U_T,
      update_indices,
      num_ccd_iter,
      alpha,
      beta
    );

    V_T = regress_cols_of_Y_on_X_nn_glmpca_pois_exact(
      U_T.t(),
      Y,
      V_T,
      update_indices,
      num_ccd_iter,
      alpha,
      beta
    );

    arma::vec d = mean(U_T, 1) / mean(V_T, 1);

    U_T.each_col() %= arma::sqrt(1/d);
    V_T.each_col() %= arma::sqrt(d);

    loglik = get_loglik_nn_glmpca(
      U_T,
      V_T,
      sc_x,
      sc_i,
      sc_j,
      n,
      p
    );

    loglik_history.push_back(loglik);

    if (loglik - prev_lik < 1e-8) {

      break;

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
