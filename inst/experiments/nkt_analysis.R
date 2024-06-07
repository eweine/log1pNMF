# I think there is potentially a lot of interesting analysis to be
# done on this natural killer T cell dataset.

# for one, I think it would be useful for me to do what Yusha
# described, to see if I am able to capture similar multiplicative
# effects. Once I've done that, then I think it would be interesting
# to look at the differences in results of the factor model methods
library(dplyr)
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
  rep("Unstim", ncol(m1)),
  rep("Unstim", ncol(m2)),
  rep("Unstim", ncol(m3)),
  rep("Stim", ncol(m4)),
  rep("Stim", ncol(m5)),
  rep("Stim", ncol(m6))
)

rm(m1, m2, m3, m4, m5, m6)
m <- as(m, "sparseMatrix")

#m <- m[(rowSums(m) >= 10) | (apply(m, 1, max) >= 5), ]
m <- m[, Matrix::colSums(m) > 0]
m <- m[Matrix::rowSums(m) > 0, ]

m <- Matrix::t(m)

stim_m <- m[samples == "Stim",]
unstim_m <- m[samples == "Unstim",]

s_stim <- sum(stim_m)
s_unstim <- sum(unstim_m)

lambda_g <- Matrix::colSums(m) / (s_stim + s_unstim)

lambda_stim <- Matrix::colSums(stim_m) / s_stim
lambda_unstim <- Matrix::colSums(unstim_m) / s_unstim

lambda_df <- data.frame(
  lambda_g = lambda_g,
  lambda_stim = lambda_stim,
  lambda_unstim = lambda_unstim
)

# filter out genes that are 0 in one condition
lambda_df <- lambda_df %>%
  dplyr::filter(
    lambda_stim > 0 & lambda_unstim > 0
  )

lambda_df <- lambda_df %>%
  dplyr::mutate(
    add_diff = abs(lambda_stim - lambda_unstim),
    mult_diff = abs(log(lambda_stim) - log(lambda_unstim))
  )

library(ggplot2)

ggplot(data = lambda_df) +
  geom_point(aes(x = lambda_g, y = add_diff)) +
  geom_smooth(aes(x = lambda_g, y = add_diff))


ggplot(data = lambda_df) +
  geom_point(aes(x = lambda_g, y = mult_diff)) +
  geom_smooth(aes(x = lambda_g, y = mult_diff))
