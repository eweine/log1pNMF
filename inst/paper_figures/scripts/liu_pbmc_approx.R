library(Matrix)
library(dplyr)
library(passPCA)

load("/home/ericw456/pbmc/liu_data.Rdata")

counts <- counts[,Matrix::colSums(counts) > 0]
# require that all used genes appear in at least 5 cells
s <- Matrix::rowSums(counts)
s <- s / mean(s)
genes_to_use <- which(Matrix::colSums(counts>0)>4)
counts <- counts[,genes_to_use]
K <- 25
cc_vec <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4)

n <- nrow(counts)
p <- ncol(counts)

for (cc in cc_vec) {

  print(cc)

  set.seed(1)
  log1p_k1 <- fit_factor_model_log1p_exact(
    Y = counts,
    K = 1,
    maxiter = 10,
    s = cc * s,
    init_method = "frob_nmf"
  )

  init_LL <- log1p_k1$U %>%
    cbind(
      matrix(
        data = rexp(
          n = n * (K - 1), rate = 15
        ),
        nrow = n,
        ncol = K - 1
      )
    )

  init_FF <- log1p_k1$V %>%
    cbind(
      matrix(
        data = rexp(
          n = p * (K - 1), rate = 15
        ),
        nrow = p,
        ncol = K - 1
      )
    )

  tictoc::tic()
  set.seed(1)
  fit <- fit_factor_model_log1p_quad_approx_sparse(
    Y = counts,
    K = K,
    init_U = init_LL,
    init_V = init_FF,
    maxiter = 1000,
    s = cc * s,
    approx_method = "taylor"
  )
  total_time <- tictoc::toc()

  fit[["total_time"]] <- total_time$toc
  rownames(fit$U) <- rownames(counts)
  rownames(fit$V) <- colnames(counts)

  readr::write_rds(
    fit, glue::glue("liu_pbmc_log1p_c{cc}_k25_approx_taylor_1K_iter.rds")
  )

}
