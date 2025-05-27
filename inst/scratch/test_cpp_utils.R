# Here, I want to explore different ways of dealing with
# sparse matrices in R so that I can easily feed them to C++
library(Matrix)
load("~/Documents/data/fastglmpca/raw_data/droplet.RData")

sc <- summary(counts)

col_num_repeats <- log1pNMF:::get_num_repeats(
  sc$j,
  ncol(counts),
  length(sc$j)
)

y_cols_list <- log1pNMF:::create_vals_list(
  ncol(counts),
  col_num_repeats,
  sc$x
)

y_cols_idx <- log1pNMF:::create_vals_list(
  ncol(counts),
  col_num_repeats,
  sc$i - 1 # account for 0 indexing in C++
)


sc_t <- summary(Matrix::t(counts))

row_num_repeats <- log1pNMF:::get_num_repeats(
  sc_t$j,
  nrow(counts),
  length(sc_t$j)
)

y_rows_list <- log1pNMF:::create_vals_list(
  nrow(counts),
  row_num_repeats,
  sc_t$x
)

y_rows_idx <- log1pNMF:::create_vals_list(
  nrow(counts),
  row_num_repeats,
  sc_t$i - 1 # account for 0 indexing in C++
)

# Now, I would like to test if I can run a series of regressions
# I have to be careful about how exactly I call this
# but I don't think it should be too hard

set.seed(1)
n <- 1250
p <- 5

X <- matrix(
  data = 0, nrow = n, ncol = p
)

X[, 1] <- abs(rnorm(n, sd = .5))
X[, 2] <- abs(rnorm(n, sd = 1))
X[, 3] <- rexp(n, rate = 3)
X[, 4] <- rgamma(n, shape = 3, rate = 3)
X[1:460, 5] <- abs(rnorm(460, sd = .01))
X[461:500, 5] <- rexp(40, 1)

b <- c(0, 0, .1, .3, .6)

lambda <- exp(X %*% b) - 1
y1 <- rpois(n, lambda)
y2 <- rpois(n, lambda)

y_list <- list()
y_idx_list <- list()

y_list[[1]] <- as.integer(y1[y1 > 0])
y_list[[2]] <- as.integer(y2[y2 > 0])

y_idx_list[[1]] <- as.integer(which(y1 > 0) - 1)
y_idx_list[[2]] <- as.integer(which(y2 > 0) - 1)

reg_mat_out <- log1pNMF:::regress_cols_of_Y_on_X_log1p_pois_exact(
  X,
  y_list,
  y_idx_list,
  matrix(data = abs(rnorm(10)), nrow = 5, ncol = 2),
  0:4,
  150,
  .01,
  .25
)

