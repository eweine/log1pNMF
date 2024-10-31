#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <omp.h>
#include "update_add_const.h"


double update_a(
  double a,
  const arma::mat U_T,
  const arma::mat V_T,
  const std::vector<int> nonzero_y,
  const std::vector<int> nonzero_y_i_idx,
  const std::vector<int> nonzero_y_j_idx,
  const int num_nonzero_y,
  const arma::vec s,
  const int n,
  const int p,
  double alpha,
  double beta
) {

  //Rprintf("Updating constant...\n");
  double exp_term = 0.0;

  #pragma omp parallel for reduction(+:exp_term)
  for (int i = 0; i < n; i++) {

    for (int j = 0; j < p; j++) {

      exp_term += exp(a + dot(U_T.col(i), V_T.col(j)));

    }
  }

  double lin_term_f = 0.0;
  double lin_term_g = 0.0;
  double lin_term_h = 0.0;

  #pragma omp parallel for reduction(+:lin_term_f, lin_term_g, lin_term_h)
  for (int r = 0; r < num_nonzero_y; r++) {

    double exp_eta = exp(
      a + dot(U_T.col(nonzero_y_i_idx[r]), V_T.col(nonzero_y_j_idx[r]))
    );

    lin_term_f -= nonzero_y[r] * log(exp_eta - s[nonzero_y_i_idx[r]]);
    lin_term_g -= nonzero_y[r] * (exp_eta / (exp_eta - s[nonzero_y_i_idx[r]]));
    lin_term_h -= nonzero_y[r] * (
      (exp_eta * (exp_eta - s[nonzero_y_i_idx[r]]) -
        exp_eta * exp_eta) / pow(exp_eta - s[nonzero_y_i_idx[r]], 2)
    );

  }

  double current_lik = lin_term_f + exp_term;
  double g = lin_term_g + exp_term;
  double h = lin_term_h + exp_term;

  // now, do a line search until the new value of alpha has improved
  double new_lik = std::numeric_limits<double>::infinity();

  double t;

  double newton_dir = g / h;
  if (newton_dir < 0) {

    t = 1.0;

  } else if (a >= 1e-12) {

    t = std::min((a - 1e-12) / newton_dir, 1.0);

  } else {

    return a;

  }

  double newton_dec = alpha * g * newton_dir;
  double a_new = a;

  bool return_a = false;

  while (new_lik > current_lik - t * newton_dec) {

    a_new = a - t * newton_dir;
    exp_term = 0.0;

    #pragma omp parallel for reduction(+:exp_term)
    for (int i = 0; i < n; i++) {

      for (int j = 0; j < p; j++) {

        exp_term += exp(a_new + dot(U_T.col(i), V_T.col(j)));

      }
    }

    lin_term_f = 0.0;

    #pragma omp parallel for reduction(+:lin_term_f)
    for (int r = 0; r < num_nonzero_y; r++) {

      double exp_eta = exp(
        a_new + dot(U_T.col(nonzero_y_i_idx[r]), V_T.col(nonzero_y_j_idx[r]))
      );

      lin_term_f -= nonzero_y[r] * log(exp_eta - 1);

    }

    new_lik = lin_term_f + exp_term;
    t *= beta;

    if (t <= 1e-12) {

      return_a = true;
      break;

    }

  }

  if (return_a) {

    return a;

  }

  return a_new;

}

