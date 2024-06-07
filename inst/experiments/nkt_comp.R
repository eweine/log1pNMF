library(Matrix)
library(dplyr)
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

condition <- c(
  rep("Unstim", ncol(m1)),
  rep("Unstim", ncol(m2)),
  rep("Unstim", ncol(m3)),
  rep("Stim", ncol(m4)),
  rep("Stim", ncol(m5)),
  rep("Stim", ncol(m6))
)

sample <- c(
  rep("1", ncol(m1)),
  rep("2", ncol(m2)),
  rep("3", ncol(m3)),
  rep("1", ncol(m4)),
  rep("2", ncol(m5)),
  rep("3", ncol(m6))
)

rm(m1, m2, m3, m4, m5, m6, genes1, genes2, genes3, genes4, genes5, genes6)
m <- as(m, "sparseMatrix")
gc()

m <- Matrix::t(m)

m <- m[,Matrix::colSums(m) >= 100]

samp_df <- data.frame(condition = condition, sample = sample, cell_id = 1:length(condition))
rm(condition, sample)
gc()
samp_df$s <- Matrix::rowSums(m)
samp_df$condition <- factor(samp_df$condition, levels = c("Unstim", "Stim"))

stim_mean <- c()
unstim_mean <- c()
p_vec <- c()

i <- 1

for (gene in colnames(m)) {
  
  print(i)
  i <- i + 1
  samp_df$count <- m[,gene]
  sum_df <- samp_df %>% 
    dplyr::group_by(condition) %>%
    dplyr::summarize(
      gene_cts = sum(count),
      tot_cts = sum(s)
    )
  
  fish_test <- fisher.test(matrix(c(
      sum_df$gene_cts[1], 
      sum_df$tot_cts[1]-sum_df$gene_cts[1], 
      sum_df$gene_cts[2], 
      sum_df$tot_cts[2]-sum_df$gene_cts[2]
    ), ncol=2)
  )
  
  p_vec <- c(p_vec, fish_test$p.value)
  unstim_mean <- c(unstim_mean, sum_df$gene_cts[1] / sum_df$tot_cts[1])
  stim_mean <- c(stim_mean, sum_df$gene_cts[2] / sum_df$tot_cts[2])
  
}

df_out <- data.frame(stim_mean = stim_mean, unstim_mean = unstim_mean, p_value = p_vec)

df_out$gene <- colnames(m)

df_out <- df_out %>%
  dplyr::filter(p_value < 5e-8)

df_out <- df_out %>%
  dplyr::mutate(diff = stim_mean - unstim_mean)

df_out$diffpm <- df_out$diff * 10000

big_diff <- df_out %>% dplyr::filter(diffpm >= 1)
plot(big_diff$unstim_mean * 10000, big_diff$diffpm)

library(glmGamPoi)
# mod <- glm_gp(
#   Matrix::t(m),
#   design = ~ 1 + condition,
#   col_data = samp_df,
#   reference_level = "Unstim",
#   subsample = TRUE,
#   on_disk = FALSE,
#   verbose = T
# )

dds <- DESeqDataSetFromMatrix(
  countData = Matrix::t(m),
  design = ~ 1 + condition,
  colData = samp_df
)

dds <- DESeq(dds, sfType = "poscounts", parallel = TRUE)
