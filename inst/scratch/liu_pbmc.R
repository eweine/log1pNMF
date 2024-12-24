library(Matrix)
library(dplyr)
counts <- counts[,Matrix::colSums(counts) > 0]
# require that all used genes appear in at least 5 cells
s <- Matrix::rowSums(counts)
s <- s / mean(s)
genes_to_use <- which(Matrix::colSums(counts>0)>4)
counts <- counts[,genes_to_use]

cc_vec <- c(0.0001, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

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


  set.seed(1)
  fit <- fit_factor_model_log1p_exact(
    Y = counts,
    K = K,
    init_U = init_LL,
    init_V = init_FF,
    maxiter = 100,
    s = cc * s
  )

  rownames(fit$U) <- rownames(counts)
  rownames(fit$V) <- colnames(counts)

  readr::write_rds(
    fit, glue::glue("liu_pbmc_log1p_c{cc}_k25.rds")
  )

}
