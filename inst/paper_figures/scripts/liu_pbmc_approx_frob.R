library(Matrix)
library(dplyr)
library(passPCA)

#load("/home/ericw456/pbmc/liu_data.Rdata")
load("~/Documents/data/passPCA/liu_data.Rdata")


counts <- counts[,Matrix::colSums(counts) > 0]
# require that all used genes appear in at least 5 cells
s <- Matrix::rowSums(counts)
s <- s / mean(s)
genes_to_use <- which(Matrix::colSums(counts>0)>4)
counts <- counts[,genes_to_use]
K <- 25
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

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

  set.seed(1)
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

  set.seed(1)
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
  fit_full_iter <- fit_factor_model_log1p_exact(
    Y = counts,
    K = K,
    init_U = init_LL,
    init_V = init_FF,
    maxiter = 1,
    s = cc * s
  )

  counts_trans <- MatrixExtra::mapSparse(counts, function(x){log1p(x / cc)})
  Y_dense <- as.matrix(counts_trans)

  set.seed(1)
  fit <- NNLM::nnmf(
    A = Y_dense,
    k = K,
    init = list(
      W = fit_full_iter$U,
      H = t(fit_full_iter$V)
    ),
    n.threads = 7
  )

  readr::write_rds(
    fit, glue::glue(
      "~/Documents/data/passPCA/paper_figures_res/liu_pbmc/liu_pbmc_log1p_c{cc}_k25_approx_frob_100_iter.rds"
    )
  )

}
