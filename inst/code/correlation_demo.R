# here, I could think about generating data from a mixture of zeros and
# multivariate log-normal distributions
set.seed(1)
# Function to generate an AR(1) covariance matrix
ar1_cov <- function(n, rho, sigma2 = 1) {
  # Create the matrix of |i - j| distances
  d <- abs(outer(1:n, 1:n, "-"))
  # Covariance entries: sigma2 * rho^|i-j|
  sigma2 * rho^d
}

# Example usage
n <- 400
m <- 400
k     <- 4      # dimension
rho   <- 0.5    # AR(1) correlation
sigma2 <- (0.75)^2     # variance

Sigma <- ar1_cov(k, rho, sigma2)
FF <- MASS::mvrnorm(
  n = m,
  mu = rep(0, k),
  Sigma = Sigma
)
FF <- matrix(data = pmax(
  1e-2,
  FF
), nrow = nrow(FF), ncol = ncol(FF))

LL <- matrix(
  data = pmax(rnorm(n * k, sd = 0.75), 1e-2),
  nrow = n, ncol = k
)
B <- tcrossprod(LL,FF)
Lambda <- exp(B) - 1

Y <- matrix(
  data = rpois(n * m, as.vector(Lambda)),
  nrow = n,
  ncol = m
)

Y <- as(Y, "CsparseMatrix")
library(passPCA)

compute_loglik <- function (fit, X) {
  return(sum(dpois(X,fitted(fit),log = TRUE)))
}

compute_condition_number <- function (x) {
  e <- eigen(cor(x))$values
  n <- length(e)
  return(e[1]/e[n])
}

hoyer <- function (x) {
  n <- length(x)
  return((sqrt(n) - sum(abs(x))/sqrt(sum(x^2)))/(sqrt(n - 1)))
}

ccs <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000)

# pre-allocate a data.frame to hold results
results <- data.frame(
  cc        = ccs,
  loglik    = NA_real_,
  cond_num  = NA_real_,
  hoyerL     = NA_real_,
  hoyerF = NA_real_
)

for (i in seq_along(ccs)) {
  cc_val <- ccs[i]
  set.seed(1)
  fit <- fit_poisson_log1p_nmf(
    Y           = Y,
    K           = k,
    cc          = cc_val,
    loglik      = "exact",
    init_method = "random",
    control     = list(maxiter = 1000)
  )
  
  results$loglik[i] <- compute_loglik(fit, as.matrix(Y))
  
  Lmat <- fit$LL
  Fmat <- fit$FF
  
  results$cond_num[i] <- compute_condition_number(Fmat)
  
  hoyersL <- apply(Lmat, 2, hoyer)
  hoyersF <- apply(Fmat, 2, hoyer)
  results$hoyerL[i] <- mean(hoyersL)
  results$hoyerF[i] <- mean(hoyersF)
}

library(ggplot2)
ga <- ggplot(data = results, aes(x = cc, y = loglik)) +
  geom_point(color = "dodgerblue") +
  geom_line(color = "dodgerblue") +
  scale_x_log10() +
  cowplot::theme_cowplot() +
  geom_point(
    data   = subset(results, cc == 1),
    shape  = 21,         # circle with border + fill
    colour = "black",    # border color
    fill   = NA,         # no fill
    size   = 2,          # adjust outer diameter
    stroke = 1.5         # adjust border thickness
  )

gb <- ggplot(data = results, aes(x = cc, y = hoyerL)) +
  ylab("Sparsity of L") + 
  geom_point(color = "darkorchid") +
  geom_line(color = "darkorchid") +
  scale_x_log10() +
  cowplot::theme_cowplot() +
  geom_point(
    data   = subset(results, cc == 1),
    shape  = 21,         # circle with border + fill
    colour = "black",    # border color
    fill   = NA,         # no fill
    size   = 2,          # adjust outer diameter
    stroke = 1.5         # adjust border thickness
  )

gc <- ggplot(data = results, aes(x = cc, y = hoyerF)) +
  ylab("Sparsity of F") + 
  geom_point(color = "olivedrab") +
  geom_line(color = "olivedrab") +
  scale_x_log10() +
  cowplot::theme_cowplot() +
  geom_point(
    data   = subset(results, cc == 1),
    shape  = 21,         # circle with border + fill
    colour = "black",    # border color
    fill   = NA,         # no fill
    size   = 2,          # adjust outer diameter
    stroke = 1.5         # adjust border thickness
  )

gd <- ggplot(data = results, aes(x = cc, y = cond_num)) +
  ylab("F Correlation Cond #") + 
  geom_point(color = "darkorange") +
  geom_line(color = "darkorange") +
  scale_x_log10() +
  cowplot::theme_cowplot() +
  geom_point(
    data   = subset(results, cc == 1),
    shape  = 21,         # circle with border + fill
    colour = "black",    # border color
    fill   = NA,         # no fill
    size   = 2,          # adjust outer diameter
    stroke = 1.5         # adjust border thickness
  )

library(ggpubr)
ggarrange(ga, gb, gc, gd, nrow = 2, ncol = 2)
