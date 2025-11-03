# install_github('linxihui/NNLM')
library(NNLM)

fit_nmf_rank1_init <- function(
  X,
  K
) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  fit_rank1 <- nnmf(
    A = X,
    k = 1
  )
  
  init_W <- cbind(
    fit_rank1$W,
    matrix(
      data = runif(n * (K - 1), min = 0, max = 0.1),
      nrow = n,
      ncol = K - 1
    )
  )
  
  
  init_H <- rbind(
    fit_rank1$H,
    matrix(
      data = runif(p * (K - 1), min = 0, max = 0.1),
      nrow = K - 1,
      ncol = p
    )
  )
  
  fit_rank_K <- nnmf(
    A = X,
    init = list(W = init_W, H = init_H),
    k = K
  )
  
  return(fit_rank_K)
  
}

set.seed(1)
n <- 500
p <- 500
K <- 5

# generate rank 5 data
LL <- matrix(
  data = runif(n * K, min = 0, max = 5),
  nrow = n,
  ncol = K
)

FF <- matrix(
  data = runif(p * K, min = 0, max = 5),
  nrow = p,
  ncol = K
)

M <- tcrossprod(LL, FF)

# add log-normal noise
E <- exp(
  matrix(
    rnorm(n = n * p),
    nrow = n,
    ncol = p
  )
)
X <- M + E

nmf_fit <- fit_nmf_rank1_init(X, K)

