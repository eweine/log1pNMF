# Short script to illustrate the use of passPCA on a toy data set.
library(passPCA)
set.seed(3)

# Simulate a data set from a Poisson NMF model with the shifted log
# link.
n <- 200
m <- 400
k <- 4
L <- matrix(0,n,k)
F <- matrix(0,m,k)
shift_factor_true <- 2
for (j in 1:k) {
  i <- sample(m,100)
  F[i,j] <- abs(rnorm(100))
}
for (i in seq(1,n/2)) {
  j <- sample(k,1)
  L[i,j] <- 1
}
for (i in seq(n/2+1,n)) {
  j <- sample(k,2)
  x <- runif(1)
  L[i,j] <- c(x,1 - x)
}
L <- cbind(runif(n),L)
F <- cbind(abs(rnorm(m)),F)
i <- sample(m,40)
F[i,1] <- abs(10*rnorm(40))
B <- tcrossprod(L,F)
ginv <- function (x, s = 1)
  s * exp(x/max(1,s) - 1)
X <- rpois(n*m,ginv(B,s = shift_factor_true))
X <- matrix(X,n,m)
storage.mode(X) <- "double"

# Function to compute the log-likelihood for a given Poisson log1p NMF
# model.
compute_loglik <- function (fit, X) {
  return(sum(dpois(X,fitted(fit)),log = TRUE))
}

# Fit a Poisson NMF model with the shifted log link.
k   <- k + 1
fit <- fit_poisson_log1p_nmf(X,k,s = rep(1,n),cc = shift_factor_true,
                             init_LL = L,init_FF = F,loglik = "exact",
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
# plot(F,F_est,pch = 20,xlab = "truth",ylab = "estimate")
# abline(a = 0,b = 1,col = "magenta",lty = "dotted")

# Now fit a Poisson NMF model with the shifted log like for different
# values of the shift factor.
shift_factors <- round(10^seq(-1.3,1,length.out = 20),digits = 4)
fits   <- vector("list",length(shift_factors))
loglik <- rep(0,length(shift_factors))
L_cor  <- rep(0,length(shift_factors))
F_cor  <- rep(0,length(shift_factors))
names(fits)   <- as.character(shift_factors)
names(loglik) <- as.character(shift_factors)
names(L_cor)  <- as.character(shift_factors)
names(F_cor)  <- as.character(shift_factors)
for (shift_factor in shift_factors) {
  cat("shift_factor =",shift_factor,"\n")

  # Fit the model.
  fit <- fit_poisson_log1p_nmf(X,k,s = rep(1,n),cc = shift_factor,
                               init_LL = L,init_FF = F,loglik = "exact",
                               control = list(maxiter = 200,threads = 1,
                                              verbose = FALSE))
  fits[[as.character(shift_factor)]] <- fit

  # Compute the log-likelihood:
  loglik[as.character(shift_factor)] <- compute_loglik(fit,X)
  
  # Compare the estimates to the truth.
  L_est <- fit$LL
  F_est <- fit$FF
  kset <- apply(cor(L_est,L),2,which.max)
  L_est <- L_est[,kset]
  F_est <- F_est[,kset]
  s     <- apply(L_est,2,max)
  L_est <- L_est %*% diag(1/s)
  F_est <- F_est %*% diag(s)
  # L_cor[as.character(shift_factor)] <- mean((L_est - L)^2)
  # F_cor[as.character(shift_factor)] <- mean((F_est - F)^2)
  L_cor[as.character(shift_factor)] <- mean(diag(cor(L_est,L)))
  F_cor[as.character(shift_factor)] <- mean(diag(cor(F_est,F)))
}

# Plot shift factor vs. MSE for L.
par(mfrow = c(2,3))
i <- which.min(abs(shift_factor_true - shift_factors))
plot(shift_factors,L_cor,pch = 20,log = "x",col = "darkorchid",
     xlab = "shift factor",ylab = "mean correlation, L")
lines(shift_factors,L_cor,col = "darkorchid")
points(shift_factors[i],L_cor[i],pch = 1,cex = 1.2,col = "black")

# Plot shift factor vs. MSE for F.
plot(shift_factors,F_cor,col = "royalblue",pch = 20,log = "x",
     xlab = "shift factor",ylab = "mean correlation, F")
lines(shift_factors,F_cor,col = "royalblue")
points(shift_factors[i],F_cor[i],pch = 1,cex = 1.2,col = "black")

# Plot shift factor vs. log-likelihood.
plot(shift_factors,loglik,pch = 20,log = "x",col = "dodgerblue",
     xlab = "shift factor")
lines(shift_factors,loglik,col = "dodgerblue")
points(shift_factors[i],loglik[i],pch = 1,cex = 1.2,col = "black")

# Plot the condition number of F.
compute_condition_number <- function (x) {
  e <- eigen(cor(x))$values
  n <- length(e)
  return(e[1]/e[n])
}
cond_F <- sapply(fits,function (x) compute_condition_number(x$FF))
plot(shift_factors,cond_F,pch = 20,log = "x",col = "olivedrab",
     xlab = "shift factor",ylab = "orthogonality of F")
lines(shift_factors,cond_F,col = "olivedrab")
points(shift_factors[i],cond_F[i],pch = 1,cex = 1.2,col = "black")

# Plot the sparsity of F.
hoyer <- function (x) {
  n <- length(x)
  return((sqrt(n) - sum(abs(x))/sqrt(sum(x^2)))/(sqrt(n - 1)))
}
sparsity_F <- sapply(fits,function (x) median(apply(x$FF,2,hoyer)))
plot(shift_factors,sparsity_F,pch = 20,log = "x",col = "darkorange",
     xlab = "shift factor",ylab = "sparsity in F")
lines(shift_factors,sparsity_F,col = "darkorange")
points(shift_factors[i],sparsity_F[i],pch = 1,cex = 1.2,col = "black")

# Plot the sparsity of L.
sparsity_L <- sapply(fits,function (x) median(apply(x$LL,2,hoyer)))
plot(shift_factors,sparsity_L,pch = 20,log = "x",col = "tomato",
     xlab = "shift factor",ylab = "sparsity in L")
lines(shift_factors,sparsity_L,col = "tomato")
points(shift_factors[i],sparsity_L[i],pch = 1,cex = 1.2,col = "black")

