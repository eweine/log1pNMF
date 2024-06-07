# here, I want to try a few models to understand the effects of different
# fitting and initializations
# first, I'm interested in fitting both log1p and Poisson NMF models
# both with and without an "intercept" factor initialization
library(glue)
library(Matrix)

data_dir <- "~/Documents/data/passPCA/experiment_results"

m1 <- as.matrix(Matrix::readMM(
  '~/Downloads/GSE128243_RAW/GSM3669244_NKT_HS_Unstim1_matrix.mtx'
))
genes1 <- readr::read_tsv("~/Downloads/GSE128243_RAW/GSM3669244_NKT_HS_Unstim1_genes.tsv",
                          col_names = c("ensembl", "name"))
rownames(m1) <- genes1$ensembl

m2 <- as.matrix(Matrix::readMM('~/Downloads/GSE128243_RAW/GSM3669245_NKT_HS_Unstim2_matrix.mtx'))
genes2 <- readr::read_tsv(
  "~/Downloads/GSE128243_RAW/GSM3669245_NKT_HS_Unstim2_genes.tsv",
  col_names = c("ensembl", "name")
)
rownames(m2) <- genes2$ensembl

m3 <- as.matrix(Matrix::readMM('~/Downloads/GSE128243_RAW/GSM3669246_NKT_HS_Unstim3_matrix.mtx'))
genes3 <- readr::read_tsv(
  "~/Downloads/GSE128243_RAW/GSM3669246_NKT_HS_Unstim3_genes.tsv",
  col_names = c("ensembl", "name")
)
rownames(m3) <- genes3$ensembl


m4 <- as.matrix(Matrix::readMM('~/Downloads/GSE128243_RAW/GSM3669247_NKT_HS_Stim1_matrix.mtx'))
genes4 <- readr::read_tsv(
  "~/Downloads/GSE128243_RAW/GSM3669247_NKT_HS_Stim1_genes.tsv",
  col_names = c("ensembl", "name")
)
rownames(m4) <- genes4$ensembl

m5 <- as.matrix(Matrix::readMM('~/Downloads/GSE128243_RAW/GSM3669248_NKT_HS_Stim2_matrix.mtx'))
genes5 <- readr::read_tsv("~/Downloads/GSE128243_RAW/GSM3669248_NKT_HS_Stim2_genes.tsv",
                          col_names = c("ensembl", "name")
)
rownames(m5) <- genes5$ensembl

m6 <- as.matrix(Matrix::readMM('~/Downloads/GSE128243_RAW/GSM3669249_NKT_HS_Stim3_matrix.mtx'))
genes6 <- readr::read_tsv("~/Downloads/GSE128243_RAW/GSM3669249_NKT_HS_Stim3_genes.tsv",
                          col_names = c("ensembl", "name")
)
rownames(m6) <- genes6$ensembl

m <- cbind(
  m1, m2, m3, m4, m5, m6
)

samples <- c(
  rep("Unstim Donor 1", ncol(m1)),
  rep("Unstim Donor 2", ncol(m2)),
  rep("Unstim Donor 3", ncol(m3)),
  rep("Stim Donor 1", ncol(m4)),
  rep("Stim Donor 2", ncol(m5)),
  rep("Stim Donor 3", ncol(m6))
)

rm(m1, m2, m3, m4, m5, m6)
m <- as(m, "sparseMatrix")

#m <- m[(rowSums(m) >= 10) | (apply(m, 1, max) >= 5), ]
m <- m[, Matrix::colSums(m) > 0]
m <- m[Matrix::rowSums(m) > 0, ]

m <- Matrix::t(m)

gc()

set.seed(1)
log1p_fit_rand_init <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = m,
  K = 10,
  approx_range = c(0, 1.25),
  maxiter = 100,
  s = rowSums(m) / mean(rowSums(m))
)

readr::write_rds(log1p_fit_rand_init, glue("{data_dir}/log1p_nkt_rand_init_k10.rds"))

set.seed(1)
log1p_mod_k1 <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = m,
  K = 1,
  approx_range = c(0, 1.25),
  maxiter = 10,
  s = rowSums(m) / mean(rowSums(m))
)

set.seed(1)
U_init <- cbind(
  log1p_mod_k1$U,
  matrix(
    data = rexp(n = nrow(log1p_mod_k1$U) * 9, rate = 15), nrow = nrow(log1p_mod_k1$U)
  )
)

V_init <- cbind(
  log1p_mod_k1$V,
  matrix(
    data = rexp(n = nrow(log1p_mod_k1$V) * 9, rate = 15), nrow = nrow(log1p_mod_k1$V)
  )
)

log1p_mod_k1_init <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = m,
  K = 10, s = rowSums(m) / mean(rowSums(m)),
  approx_range = c(0, 1.25), maxiter = 100,
  init_U = U_init, init_V = V_init
)

readr::write_rds(log1p_mod_k1_init, glue("{data_dir}/log1p_nkt_k1_init_k10.rds"))
library(fastTopics)
# now, I want to fit Poisson NMF models
set.seed(1)
nmf_pois_rand_init <- fit_poisson_nmf(
  m,
  k = 10,
  numiter = 100,
  control = list(nc = 7),
  init.method = "random"
)

readr::write_rds(nmf_pois_rand_init, glue("{data_dir}/pois_nmf_nkt_rand_init_k10.rds"))

nmf_pois_k1 <- fastTopics:::fit_pnmf_rank1(m)

set.seed(1)
L_init <- cbind(
  nmf_pois_k1$L,
  matrix(
    data = rexp(n = nrow(nmf_pois_k1$L) * 9, rate = 15), nrow = nrow(nmf_pois_k1$L)
  )
)

F_init <- cbind(
  nmf_pois_k1$F,
  matrix(
    data = rexp(n = nrow(nmf_pois_k1$F) * 9, rate = 15), nrow = nrow(nmf_pois_k1$F)
  )
)

rownames(F_init) <- colnames(m)

set.seed(1)

fit0_pois_nmf <- init_poisson_nmf(
  X = m,
  L = L_init,
  F = F_init
)

nmf_pois_k1_init <- fit_poisson_nmf(
  m,
  numiter = 100,
  control = list(nc = 7),
  fit0 = fit0_pois_nmf
)

readr::write_rds(nmf_pois_k1_init, glue("{data_dir}/pois_nmf_nkt_k1_init_k10.rds"))

library(Matrix)
library(RcppML)

Y <- m / (rowSums(m) / mean(rowSums(m)))
Y_tilde <- MatrixExtra::mapSparse(Y, log1p)
set.seed(1)
frob_nmf_fit <- nmf(
  A = Y_tilde, k = 10
)

readr::write_rds(frob_nmf_fit, glue("{data_dir}/frob_nmf_nkt_rand_init_k10.rds"))

frob_nmf_k1 <- NNLM::nnmf(
  A = as.matrix(Y_tilde), k = 1
)

W_init <- cbind(
  frob_nmf_k1$W,
  matrix(data = rexp(n = 9 * nrow(Y_tilde), rate = 15), nrow = nrow(Y_tilde))
)

H_init <- rbind(
  frob_nmf_k1$H,
  matrix(data = rexp(n = 9 * ncol(Y_tilde), rate = 15), ncol = ncol(Y_tilde))
)

frob_nmf_k1_init <- NNLM::nnmf(
  A = as.matrix(Y_tilde), k = 10,
  n.threads = 0, max.iter = 100,
  init = list(
    W = W_init,
    H = H_init
  )
)

readr::write_rds(frob_nmf_k1_init, glue("{data_dir}/frob_nmf_nkt_k1_init_k10.rds"))

flash_fit <- flashier::flash(
  data = Y_tilde, greedy_Kmax = 10, backfit = T,
  ebnm_fn = ebnm::ebnm_point_exponential
)

readr::write_rds(flash_fit, glue("{data_dir}/frob_nmf_flash_k10.rds"))

max_col <- apply(flash_fit$L_pm, 2, max)
flash_LL <- sweep(flash_fit$L_pm, 2, max_col, FUN = "/")
flash_FF <- sweep(flash_fit$F_pm, 2, max_col, FUN = "*")

structure_plot(flash_LL, grouping = samples)



