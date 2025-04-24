# Short script to illustrate the use of passPCA on a toy data set.
library(passPCA)
set.seed(1)

# Simulate a data set from a Poisson NMF model with the shifted log
# link.
n <- 100
m <- 400
k <- 4
L <- matrix(0,n,k)
F <- matrix(0,m,k)
shift_factor <- 1
for (j in 1:k) {
  i <- sample(m,40)
  F[i,j] <- abs(rnorm(40))
}
for (i in 1:n) {
  j <- sample(k,1)
  L[i,j] <- 1
}
B <- tcrossprod(L,F)
ginv <- function (x, s = 1)
  s * exp(x/max(1,s) - 1)
X <- rpois(n*m,ginv(B,s = shift_factor))
X <- matrix(X,n,m)
storage.mode(X) <- "double"

# Fit a Poisson NMF model with the shifted log link.
fit <- fit_poisson_log1p_nmf(X,k,s = rep(1,n),cc = shift_factor,
                             loglik = "exact",
                             control = list(maxiter = 500,threads = 1))

# Compare the estimates to the truth.
L_est <- fit$LL
F_est <- fit$FF
kset <- apply(cor(L,L_est),2,which.max)
L_est <- L_est[,kset]
F_est <- F_est[,kset]
s     <- apply(L_est,2,max)
L_est <- L_est %*% diag(1/s)
F_est <- F_est %*% diag(s)
print(cor(L,L_est))
print(cor(F,F_est))
plot(F[,kset],fit$FF,pch = 20,xlab = "truth",ylab = "estimate")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")
