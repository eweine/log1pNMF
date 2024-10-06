data {
  real z; // single observation
  real<lower=0> s; // known standard deviation
  real<lower=0> pi0; // mixture component on point mass at 0
  real<lower=0> mu; // rate parameter for the exponential distribution
}

parameters {
  real<lower=0> m; // mean parameter, constrained to be non-negative
}

transformed parameters {
  real log_p0; // log posterior probability of the 0 component
  real log_p1; // log posterior probability of the exponential component
  real p0; // posterior responsibility for the 0 component

  log_p0 = log(pi0) + normal_lpdf(z | 0, s);
  log_p1 = log(1 - pi0) + normal_lpdf(z | m, s);

  p0 = exp(log_p0 - log_sum_exp(log_p0, log_p1));
}

model {
  m ~ exponential(mu); // prior for m

  target += log_sum_exp(log_p0, log_p1);
}
