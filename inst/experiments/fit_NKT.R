library(Matrix)

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
nmf_pois <- fit_poisson_nmf(
  m,
  k = 10,
  numiter = 100,
  control = list(nc = 7)
)

#readr::write_rds(nmf_pois, "~/Documents/passPCA/inst/experiments/results/pois_nmf_NKT.rds")

set.seed(1)
log1p_fit <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = m,
  K = 10,
  approx_range = c(0, 1.25),
  maxiter = 100,
  s = rowSums(m) / mean(rowSums(m))
)

readr::write_rds(log1p_fit, "~/Documents/passPCA/inst/experiments/results/pois_log1p_NKT.rds")

log1p_fit <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/pois_log1p_NKT.rds")

g1 <- structure_plot(
  nmf_pois,
  grouping = as.factor(samples),
  gap = 10
) + ggplot2::ggtitle("Pois NMF")

log1p_fit$L <- log1p_fit$U
log1p_fit$FF <- log1p_fit$V

log1p_fit$Lnorm = t(t(log1p_fit$L)/apply(log1p_fit$L,2,max))
log1p_fit$Fnorm = t(t(log1p_fit$FF)*apply(log1p_fit$L,2,max))

nmf_pois$Lnorm = t(t(nmf_pois$L)/apply(nmf_pois$L,2,max))
nmf_pois$Fnorm = t(t(nmf_pois$F)*apply(nmf_pois$L,2,max))

g2 <- structure_plot(
  log1p_fit$Lnorm,
  grouping = as.factor(samples),
  gap = 10,
  n = nrow(log1p_fit$Lnorm)
) + ggplot2::ggtitle("Log1p")

library(ggpubr)
ggarrange(
  g2, g1, nrow = 1, ncol = 2
)

rownames(log1p_fit$Fnorm) <- colnames(m)

par(mfrow=c(1,2))

hist(
  log1p_fit$Lnorm,
  breaks = seq(0, 1, .025),
  ylim = c(0, 30000),
  main = "Log1p",
  xlab = "Lnorm"
)
hist(
  nmf_pois$Lnorm,
  breaks = seq(0, 1, .025),
  ylim = c(0, 30000),
  main = "Pois NMF",
  xlab = "Lnorm"
)

par(mfrow=c(1,2))

hist(
  log1p_fit$Fnorm,
  xlim = c(0, 0.5),
  breaks = seq(0, max(log1p_fit$Fnorm) + .025, .025),
  ylim = c(0, 50000),
  main = "Log1p",
  xlab = "Fnorm"
)

hist(
  nmf_pois$Fnorm,
  xlim = c(0, 0.5),
  breaks = seq(0, max(nmf_pois$Fnorm) + .025, .025),
  ylim = c(0, 50000),
  main = "Pois NMF",
  xlab = "Fnorm"
)
