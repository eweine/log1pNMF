#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
arma::vec solve_pois_reg_log1p (
    const arma::mat X, 
    const arma::vec y,
    const arma::uvec y_nz_idx,
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
  vec exp_eta_nz_m1 = exp_eta.elem(y_nz_idx) - 1;
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
        
      } else if (b[j] >= 1e-16) {
        
        t = std::min(b[j] / newton_dir, 1.0);
        
      } else {
        
        continue;
        
      }
      
      newton_dec    = alpha * first_deriv * newton_dir;
      b_j_og        = b[j];
      while (true) {
        b[j]             = b_j_og - t * newton_dir;
        eta_proposed     = eta + (b[j] - b_j_og) * X.col(j);
        exp_eta = exp(eta_proposed);
        exp_eta_nz_m1 = exp_eta.elem(y_nz_idx) - 1;
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
        }
      }
    }
  }
  
  return(b);
  
}
  
/*** R
set.seed(1)
n <- 5000
p <- 10
X <- matrix(
  data = abs(rnorm(n * p, sd = .15)),
  nrow = n,
  ncol = p
)
b <- abs(rnorm(p, sd = 0.75))
b[c(1, 2, 3, 4, 7, 8, 10)] <- 0
lambda <- exp(X %*% b) - 1
y <- rpois(n, lambda)
init <- abs(rnorm(p, sd = 0.5))
b_out <- solve_pois_reg_log1p (
  X, 
  y[y != 0],
  as.integer(which(y != 0) - 1),
  init, 
  0:9,
  100,
  .01,
  .25
)
*/
