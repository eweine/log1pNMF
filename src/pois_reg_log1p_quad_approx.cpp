#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::export]]
arma::vec solve_pois_reg_log1p_quad_approx (
    const arma::mat X_nz, 
    const arma::vec y,
    const arma::vec X_0_cs_times_a1,
    const arma::mat X_0_T_X_0,
    const double a2,
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
  
  double exact_lik = sum(exp_eta_nz) - dot(
    y,
    log(exp_eta_nz_m1)
  );
  
  int num_indices = update_indices.size();
  
  for (int update_num = 1; update_num <= num_iter; update_num++) {
    
    for (i = 0; i < num_indices; i++) {
      j = update_indices[i];
      
      current_lik = exact_lik + b[j] * X_0_cs_times_a1[j] + 
        a2 * sum(b.t() * X_0_T_X_0 * b);
      
      exp_deriv_term_nz = exp_eta_nz % X_nz.col(j);
      
      quad_deriv_vec = X_0_T_X_0 * b;
      
      first_deriv    = sum(exp_deriv_term_nz) - dot(
        y,
        exp_deriv_term_nz / exp_eta_nz_m1
      ) + X_0_cs_times_a1[j] + 2 * a2 * quad_deriv_vec[j];
      
      second_deriv   = dot(exp_deriv_term_nz, X_nz.col(j)) + dot(
        y,
        (exp_deriv_term_nz % exp_deriv_term_nz) /
          square(exp_eta_nz_m1)
      ) + 2 * a2 * X_0_T_X_0(j, j);
      
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
        
        exact_lik = sum(exp_eta_nz) - dot(
          y,
          log(exp_eta_nz_m1)
        );
        
        f_proposed = exact_lik + b[j] * X_0_cs_times_a1[j] + 
          a2 * sum(b.t() * X_0_T_X_0 * b);
        
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

y_nz_idx <- which(y != 0)
y_z_idx <- which(y == 0)

approx_poly <- pracma::polyApprox(
  exp,
  a = 0,
  b = log(7.9),
  n = 2
)

a1 <- approx_poly$p[2] 
a2 <- approx_poly$p[1]

b_out <-  solve_pois_reg_log1p_quad_approx (
  X[y_nz_idx, ], 
  y[y != 0],
  a1 * colSums(X[y_z_idx, ]),
  crossprod(X[y_z_idx, ]),
  a2,
  init, 
  0:9,
  100,
  .01,
  .25
)
*/
