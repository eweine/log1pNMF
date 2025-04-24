# Short script to illustrate the use of passPCA on a toy data set.
library(passPCA)
set.seed(1)

# Simulate a data set from a Poisson NMF model with the shifted log
# link.
n <- 200
m <- 400
k <- 4
L <- matrix(0,n,k)
F <- matrix(0,m,k)
shift_factor <- 1
for (j in 1:k) {
  i <- sample(m,40)
  F[i,j] <- abs(rnorm(40))
}
for (i in seq(1,n/2)) {
  j <- sample(k,1)
  L[i,j] <- 1
}
for (i in seq(n/2+1,n)) {
  j <- sample(k,2)
  L[i,j] <- 0.5
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
kset <- apply(cor(L_est,L),2,which.max)
L_est <- L_est[,kset]
F_est <- F_est[,kset]
s     <- apply(L_est,2,max)
L_est <- L_est %*% diag(1/s)
F_est <- F_est %*% diag(s)
print(cor(L,L_est))
print(cor(F,F_est))
plot(F,F_est,pch = 20,xlab = "truth",ylab = "estimate")
abline(a = 0,b = 1,col = "magenta",lty = "dotted")

# Now fit a Poisson NMF model with the shifted log like for different
# values of the shift factor.
shift_factors <- 10^seq(-1,1,length.out = 20)
fits <- vector("list",length(shift_factors))
L_mse <- rep(0,length(shift_factors))
F_mse <- rep(0,length(shift_factors))
names(fits)  <- as.character(shift_factors)
names(L_mse) <- as.character(shift_factors)
names(F_mse) <- as.character(shift_factors)
for (shift_factor in shift_factors) {
  cat("shift_factor =",shift_factor,"\n")

  # Fit the model.
  fit <- fit_poisson_log1p_nmf(X,k,s = rep(1,n),cc = shift_factor,
                               loglik = "exact",
                               control = list(maxiter = 200,threads = 1,
                                              verbose = FALSE))
  fits[[as.character(shift_factor)]] <- fit

  # Compare the estimates to the truth.
  L_est <- fit$LL
  F_est <- fit$FF
  kset <- apply(cor(L_est,L),2,which.max)
  L_est <- L_est[,kset]
  F_est <- F_est[,kset]
  s     <- apply(L_est,2,max)
  L_est <- L_est %*% diag(1/s)
  F_est <- F_est %*% diag(s)
  L_mse[as.character(shift_factor)] <- mean((L_est - L)^2)
  F_mse[as.character(shift_factor)] <- mean((F_est - F)^2)
}
plot(shift_factors,L_mse,pch = 20,log = "x",col = "tomato",
     xlab = "shift factor",ylab = "MSE",ylim = c(0,0.3))
lines(shift_factors,L_mse,col = "tomato")
i <- which.min(L_mse)
points(shift_factors[i],L_mse[i],pch = 1,cex = 1.2,col = "black")
points(shift_factors,F_mse,col = "royalblue",pch = 20)
lines(shift_factors,F_mse,col = "royalblue")
i <- which.min(F_mse)
points(shift_factors[i],F_mse[i],pch = 1,cex = 1.2,col = "black")
