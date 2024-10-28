library(dplyr)
dat <- readr::read_rds("~/Documents/data/lung_and_airway.rds")

cells <- dat@meta.data
genes <- dat@assays$RNA@meta.features

counts <- Matrix::t(dat@assays$RNA@counts)

genes <- genes %>% dplyr::filter(feature_type == "protein_coding")

counts <- counts[, colnames(counts) %in% rownames(genes)]

counts <- counts[,Matrix::colSums(counts) > 0]
counts <- counts[Matrix::rowSums(counts) > 0, ]

cells_sum <- cells %>%
  dplyr::group_by(author_day, author_cell_type) %>%
  dplyr::summarise(
    tot = n()
  )


counts <- as(counts, "CsparseMatrix")

n <- nrow(counts)
p <- ncol(counts)
K <- 11

rs <- Matrix::rowSums(counts)
s <- rs / mean(rs)

set.seed(1)
log1p_k1 <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = 1,
  maxiter = 10,
  approx_range = c(0, 1.25),
  s = s
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
log1p_k11 <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = K,
  maxiter = 100,
  approx_range = c(0, 1.25),
  s = s,
  init_U = init_LL,
  init_V = init_FF
)

#readr::write_rds(log1p_k11, "~/Documents/data/passPCA/experiment_results/log1p_pois_k11_mouse_embryo_lung.rds")

log1p_k11 <- readr::read_rds(
  "~/Documents/data/passPCA/experiment_results/log1p_pois_k11_mouse_embryo_lung.rds"
  )

normalize_bars <- function(LL) {

  max_col <- apply(LL, 2, max)
  sweep(LL, 2, max_col, FUN = "/")

}

LL <- normalize_bars(log1p_k11$U)
rownames(LL) <- rownames(cells)

library(fastTopics)

structure_plot(LL, grouping = cells$author_cell_type, gap = 15)

day_str <- substr(cells$author_day, 2, 5)
day_str <- ifelse(day_str == "0000", "1900", day_str)
day <- as.numeric(day_str)

corr_vec <- c()

for (k in 1:11) {

  corr_vec <- c(
    corr_vec,
    cor(LL[,k], day)
    )

}

# it would also be interesting to look at the correlations by cell type
cells$day <- day

corr_vec2 <- c()
ct_vec <- c()
factor_vec <- c()

for (ct in unique(cells$author_cell_type)) {

  for (k in 1:11) {

    cells_ct <- cells %>% dplyr::filter(author_cell_type == ct)
    LL_ct <- LL[which(rownames(counts) %in% rownames(cells_ct)), ]
    corr_vec2 <- c(
      corr_vec2,
      cor(LL_ct[,k], cells_ct$day)
    )
    factor_vec <- c(factor_vec, k)
    ct_vec <- c(ct_vec, ct)

  }

}

df <- data.frame(
  corr = corr_vec2,
  factor = factor_vec,
  ct = ct_vec
)

# let's first look at the alveolar type 1 cells
cells_av1 <- cells %>%
  dplyr::filter(author_cell_type == "Alveolar Type 1 cells") %>%
  dplyr::arrange(day)

LL_av1 <- LL[rownames(cells_av1), ]
library(fastTopics)
structure_plot(LL_av1, grouping = as.factor(as.character(cells_av1$day)))

library(ggplot2)

cells_av1$loading_4 <- LL_av1[,4]

ggplot(data = cells_av1) +
  geom_jitter(aes(x = day, y = loading_4), alpha = 0.1)



# let's look at the alveolar type 2 cells
cells_av2 <- cells %>%
  dplyr::filter(author_cell_type == "Alveolar Type 2 cells") %>%
  dplyr::arrange(day)

LL_av2 <- LL[rownames(cells_av2), ]

cells_av2 <- cells_av2 %>%
  dplyr::mutate(
    time = ifelse(day > 1875, "post_birth", "pre_birth")
  )

structure_plot(LL_av2, grouping = as.character(cells_av2$day))

library(ggplot2)

cells_av2$loading_9 <- LL_av2[,9]

cs <- cells_av2 %>%
  dplyr::group_by(day) %>%
  dplyr::summarise(
    avg_loading_9 = mean(loading_9)
  )

ggplot(data = cells_av1) +
  geom_point(aes(x = day, y = loading_4))

# there does appear to be an interesting trend above
# I think it's worth thinking about how some of the other methods would perform
# here. In particular, I want to think about Frobenius NMF with a log1p transform
# and Poisson ID link NMF

counts_tilde <- Matrix::Diagonal(x = 1/s) %*% counts
counts_tilde <- MatrixExtra::mapSparse(counts_tilde, log1p)

library(RcppML)
library(Matrix)

frob_nmf <- RcppML::nmf(
  A = counts_tilde, k = 11
)

LL_frob <- frob_nmf$w %*% diag(sqrt(frob_nmf$d))
LL_frob <- normalize_bars(LL_frob)

structure_plot(LL_frob, grouping = cells$author_cell_type, gap = 15)

rownames(LL_frob) <- rownames(cells)


corr_vec2 <- c()
ct_vec <- c()
factor_vec <- c()

for (ct in unique(cells$author_cell_type)) {

  for (k in 1:11) {

    cells_ct <- cells %>% dplyr::filter(author_cell_type == ct)
    LL_ct <- LL_frob[which(rownames(counts) %in% rownames(cells_ct)), ]
    corr_vec2 <- c(
      corr_vec2,
      cor(LL_ct[,k], cells_ct$day)
    )
    factor_vec <- c(factor_vec, k)
    ct_vec <- c(ct_vec, ct)

  }

}


df2 <- data.frame(
  corr = corr_vec2,
  factor = factor_vec,
  ct = ct_vec
)

counts_sum <- Matrix::summary(counts)

readr::write_csv(counts_sum, "~/Documents/data/passPCA/mouse_dev_lung_sum.csv")
cells$bc <- rownames(cells)
rownames(cells) <- NULL
readr::write_csv(cells, "~/Documents/data/passPCA/mouse_dev_lung_cells.csv")
genes <- genes %>% dplyr::filter(rownames(genes) %in% colnames(counts))
genes$ensembl <- rownames(genes)
rownames(genes) <- NULL
readr::write_csv(genes, "~/Documents/data/passPCA/mouse_dev_lung_genes.csv")


cnmf_usage <- readr::read_csv("/Users/eweine/Documents/data/passPCA/mouse_lung_dev_cnmf_usage.csv")
cnmf_LL <- cnmf_usage %>%
  dplyr::select(-bc) %>%
  as.matrix()
rownames(cnmf_LL) <- cnmf_usage$bc


structure_plot(
  cnmf_LL, grouping = cells$author_cell_type
)


corr_vec2 <- c()
ct_vec <- c()
factor_vec <- c()

for (ct in unique(cells$author_cell_type)) {

  for (k in 1:11) {

    cells_ct <- cells %>% dplyr::filter(author_cell_type == ct)
    LL_ct <- cnmf_LL[which(rownames(counts) %in% cells_ct$bc), ]
    corr_vec2 <- c(
      corr_vec2,
      cor(LL_ct[,k], cells_ct$day)
    )
    factor_vec <- c(factor_vec, k)
    ct_vec <- c(ct_vec, ct)

  }

}


df3 <- data.frame(
  corr = corr_vec2,
  factor = factor_vec,
  ct = ct_vec
)


