library(dplyr)
library(Matrix)
library(fastTopics)

m <- pbmc_facs$counts
m <- m[,Matrix::colSums(m) >= 100]

set.seed(1)
log1p_mod <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = m,
  K = 6,
  maxiter = 100,
  approx_range = c(0, 1.25)
)

readr::write_rds(
  log1p_mod, "~/Documents/data/passPCA/experiment_results/pbmc_facs_k6_mle.rds"
)

set.seed(1)
ff <- passPCA::run_flash_log1p_with_greedy_init(
  Y = m,
  var_type = 2
)

readr::write_rds(
  ff, "~/Documents/data/passPCA/experiment_results/pbmc_facs_flash_greedy_init.rds"
)

set.seed(1)
ff2 <- passPCA::run_flash_log1p_with_MLE_init(
  Y = m,
  K = 6
)

readr::write_rds(
  ff2, "~/Documents/data/passPCA/experiment_results/pbmc_facs_flash_mle_init.rds"
)

# I would also like to fit a model here using the previous flash pipeline
m_tilde <- MatrixExtra::mapSparse(m, log1p)
n  <- nrow(m)
x  <- rpois(1e7, 1/n)
s1 <- sd(log(x + 1))

ff3 <- flashier::flash(m_tilde,
             ebnm_fn = ebnm::ebnm_point_exponential,
             var_type = 2,
             greedy_Kmax = 10,
             S = s1,
             backfit = TRUE)

readr::write_rds(
  ff3, "~/Documents/data/passPCA/experiment_results/pbmc_facs_flash_old.rds"
)

m1 <- as.matrix(Matrix::readMM(
  '~/Downloads/GSE128243_RAW/GSM3669244_NKT_HS_Unstim1_matrix.mtx'
))
genes1 <- readr::read_tsv("~/Downloads/GSE128243_RAW/GSM3669244_NKT_HS_Unstim1_genes.tsv",
                          col_names = c("ensembl", "name"))
rownames(m1) <- genes1$ensembl

m4 <- as.matrix(Matrix::readMM('~/Downloads/GSE128243_RAW/GSM3669247_NKT_HS_Stim1_matrix.mtx'))
genes4 <- readr::read_tsv(
  "~/Downloads/GSE128243_RAW/GSM3669247_NKT_HS_Stim1_genes.tsv",
  col_names = c("ensembl", "name")
)
rownames(m4) <- genes4$ensembl

m <- cbind(
  m1, m4
)

samples <- c(
  rep("Unstim", ncol(m1)),
  rep("Stim", ncol(m4))
)

rm(m1, m4)
m <- as(m, "sparseMatrix")
m <- Matrix::t(m)
m <- m[,Matrix::colSums(m) >= 100]


set.seed(1)
log1p_mod <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = m,
  K = 5,
  maxiter = 100,
  approx_range = c(0, 1.25)
)

readr::write_rds(
  log1p_mod, "~/Documents/data/passPCA/experiment_results/nkt_k5_mle.rds"
)

set.seed(1)
ff <- passPCA::run_flash_log1p_with_greedy_init(
  Y = m,
  var_type = 2
)

readr::write_rds(
  ff, "~/Documents/data/passPCA/experiment_results/nkt_flash_greedy_init.rds"
)

m_tilde <- MatrixExtra::mapSparse(m, log1p)
n  <- nrow(m)
x  <- rpois(1e7, 1/n)
s1 <- sd(log(x + 1))

ff3 <- flashier::flash(m_tilde,
                       ebnm_fn = ebnm::ebnm_point_exponential,
                       var_type = 2,
                       greedy_Kmax = 10,
                       S = s1,
                       backfit = TRUE)


readr::write_rds(
  ff3, "~/Documents/data/passPCA/experiment_results/nkt_flash_old.rds"
)
