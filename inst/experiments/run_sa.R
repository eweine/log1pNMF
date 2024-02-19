sla <- readr::read_csv("~/Downloads/paperList.txt")

sla <- sla[!is.na(sla$abstract),]
sla$docnum = 1:nrow(sla)
datax = readRDS('~/Downloads/sla_full.rds')
dim(datax$data)

library(Matrix)
datax$data = Matrix(datax$data,sparse = TRUE)

doc_to_use = order(rowSums(datax$data),decreasing = T)[1:round(nrow(datax$data)*0.6)]
mat = datax$data[doc_to_use,]
sla = sla[doc_to_use,]
samples = datax$samples
samples = lapply(samples, function(z){z[doc_to_use]})

word_to_use = which(colSums(mat>0)>4)
mat = mat[,word_to_use]
mat = Matrix(mat,sparse=TRUE)

library(passPCA)

cc <- 1
log1p_fit <- fit_factor_model_log1p_quad_approx_sparse(
  Y = mat,
  K = 50,
  maxiter = 1000,
  approx_range = c(0, 1.25),
  s = cc * as.vector(Matrix::rowSums(mat) / mean(Matrix::rowSums(mat)))
)

readr::write_rds(
  log1p_fit,
  "~/Documents/passPCA/inst/experiments/results/sa_text_c1_log1p.rds"
)

cc <- 5
log1p_fit <- fit_factor_model_log1p_quad_approx_sparse(
  Y = mat,
  K = 50,
  maxiter = 1000,
  approx_range = c(0, 1.25),
  s = cc * as.vector(Matrix::rowSums(mat) / mean(Matrix::rowSums(mat)))
)

readr::write_rds(
  log1p_fit,
  "~/Documents/passPCA/inst/experiments/results/sa_text_c5_log1p.rds"
)

cc <- 10
log1p_fit <- fit_factor_model_log1p_quad_approx_sparse(
  Y = mat,
  K = 50,
  maxiter = 1000,
  approx_range = c(0, 1.25),
  s = cc * as.vector(Matrix::rowSums(mat) / mean(Matrix::rowSums(mat)))
)

readr::write_rds(
  log1p_fit,
  "~/Documents/passPCA/inst/experiments/results/sa_text_c10_log1p.rds"
)

cc <- 50
log1p_fit <- fit_factor_model_log1p_quad_approx_sparse(
  Y = mat,
  K = 50,
  maxiter = 1000,
  approx_range = c(0, 1.25),
  s = cc * as.vector(Matrix::rowSums(mat) / mean(Matrix::rowSums(mat)))
)

readr::write_rds(
  log1p_fit,
  "~/Documents/passPCA/inst/experiments/results/sa_text_c50_log1p.rds"
)

cc <- 100
log1p_fit <- fit_factor_model_log1p_quad_approx_sparse(
  Y = mat,
  K = 50,
  maxiter = 1000,
  approx_range = c(0, 1.25),
  s = cc * as.vector(Matrix::rowSums(mat) / mean(Matrix::rowSums(mat)))
)

readr::write_rds(
  log1p_fit,
  "~/Documents/passPCA/inst/experiments/results/sa_text_c100_log1p.rds"
)

library(fastTopics)

nmf_fit <- fit_poisson_nmf(X = mat, k = 50, numiter = 1000)
readr::write_rds(
  nmf_fit,
  "~/Documents/passPCA/inst/experiments/results/sa_text_nmf.rds"
)

cc <- .1
log1p_fit <- fit_factor_model_log1p_exact(
  Y = mat,
  K = 50,
  maxiter = 1000,
  s = cc * as.vector(Matrix::rowSums(mat) / mean(Matrix::rowSums(mat)))
)

readr::write_rds(
  log1p_fit,
  "~/Documents/passPCA/inst/experiments/results/sa_text_c.1_log1p.rds"
)

cc <- 1000
log1p_fit <- fit_factor_model_log1p_quad_approx_sparse(
  Y = mat,
  K = 50,
  maxiter = 1000,
  approx_range = c(0, 1.25),
  s = cc * as.vector(Matrix::rowSums(mat) / mean(Matrix::rowSums(mat)))
)

readr::write_rds(
  log1p_fit,
  "~/Documents/passPCA/inst/experiments/results/sa_text_c1000_log1p.rds"
)

cc <- 250
log1p_fit <- fit_factor_model_log1p_quad_approx_sparse(
  Y = mat,
  K = 50,
  maxiter = 1000,
  approx_range = c(0, 1.25),
  s = cc * as.vector(Matrix::rowSums(mat) / mean(Matrix::rowSums(mat)))
)

readr::write_rds(
  log1p_fit,
  "~/Documents/passPCA/inst/experiments/results/sa_text_c250_log1p.rds"
)

cc <- 150
log1p_fit <- fit_factor_model_log1p_quad_approx_sparse(
  Y = mat,
  K = 50,
  maxiter = 1000,
  approx_range = c(0, 1.25),
  s = cc * as.vector(Matrix::rowSums(mat) / mean(Matrix::rowSums(mat)))
)

readr::write_rds(
  log1p_fit,
  "~/Documents/passPCA/inst/experiments/results/sa_text_c150_log1p.rds"
)
