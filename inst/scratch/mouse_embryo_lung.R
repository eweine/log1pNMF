library(dplyr)
dat <- readr::read_rds("~/Documents/data/lung_and_airway.rds")

cells <- dat@meta.data
genes <- dat@assays$RNA@meta.features

counts <- Matrix::t(dat@assays$RNA@counts)

genes <- genes %>% dplyr::filter(feature_type == "protein_coding")

counts <- counts[, colnames(counts) %in% rownames(genes)]

day_str <- substr(cells$author_day, 2, 5)
day_str <- ifelse(day_str == "0000", "1900", day_str)
day <- as.numeric(day_str)

cells$day <- day

cells <- cells %>% dplyr::filter(day >= 1725)
cells <- cells %>% dplyr::filter(author_cell_type != "Lung progenitor cells")

set.seed(1)
df_list <- list()
i <- 1
for (d in unique(cells$author_day)) {

  for (ct in unique(cells$author_cell_type)) {

    ctd_df <- cells %>% dplyr::filter(author_day == d & author_cell_type == ct)

    if (nrow(ctd_df) > 500) {

      ctd_df <- ctd_df %>% dplyr::sample_n(500)

    }

    if (nrow(ctd_df) >= 25) {

      df_list[[i]] <- ctd_df
      i <- i + 1

    }

  }

}

dfo <- do.call(rbind, df_list)

counts <- counts[rownames(dfo), ]

counts <- counts[, Matrix::colSums(counts) > 0]
genes_to_use <- which(Matrix::colSums(counts>0)>4)
counts <- counts[,genes_to_use]

counts <- as(counts, "CsparseMatrix")

n <- nrow(counts)
p <- ncol(counts)
K <- 12

cc_vec <- c(0.0001, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

nmf_k1 <- fastTopics:::fit_pnmf_rank1(counts)

nmf_LL <- nmf_k1$L %>%
  cbind(
    matrix(
      data = rexp(
        n = n * (K - 1), rate = 15
      ),
      nrow = n,
      ncol = K - 1
    )
  )
rownames(nmf_LL) <- rownames(counts)

nmf_FF <- nmf_k1$F %>%
  cbind(
    matrix(
      data = rexp(
        n = p * (K - 1), rate = 15
      ),
      nrow = p,
      ncol = K - 1
    )
  )

rownames(nmf_FF) <- colnames(counts)

set.seed(1)
nmf_fit0 <- fastTopics::init_poisson_nmf(
  X = counts,
  L = nmf_LL,
  F = nmf_FF
)

nmf_fit <- fastTopics::fit_poisson_nmf(
  X = counts,
  fit0 = nmf_fit0
)

readr::write_rds(
  nmf_fit, "~/Documents/data/log1pNMF/lung_embryo_results/nmf_k12.rds"
)

library(log1pNMF)
library(Matrix)
cc_vec <- c(0.0001, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)
rs <- Matrix::rowSums(counts)
s <- rs / mean(rs)


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
    fit, glue::glue("~/Documents/data/log1pNMF/lung_embryo_results/log1p_c{cc}_k12.rds")
  )

}

fit_list <- list()

fit_list[["nmf"]] <- readr::read_rds("~/Documents/data/log1pNMF/lung_embryo_results/nmf_k12.rds")

cc_vec <- c(0.0001, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

for (cc in cc_vec) {

  print(cc)

  fit_list[[as.character(cc)]] <- readr::read_rds(
    glue::glue("~/Documents/data/log1pNMF/lung_embryo_results/log1p_c{cc}_k12.rds")
  )

}

readr::write_rds(
  fit_list, glue::glue("~/Documents/data/log1pNMF/lung_embryo_results/k12_fit_list.rds")
)

############ Code to fit full models above


set.seed(1)
log1p_k1 <- log1pNMF::fit_factor_model_log1p_quad_approx_sparse(
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
log1p_k11 <- log1pNMF::fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = K,
  maxiter = 100,
  approx_range = c(0, 1.25),
  s = s,
  init_U = init_LL,
  init_V = init_FF
)

library(fastTopics)
set.seed(1)
nmf <- fit_poisson_nmf(
  X = counts,
  k = K,
  control = list(nc = 6)
)

readr::write_rds(nmf, "~/Documents/data/log1pNMF/experiment_results/nmf_pois_k11_mouse_embryo_lung.rds")

#readr::write_rds(log1p_k11, "~/Documents/data/log1pNMF/experiment_results/log1p_pois_k11_mouse_embryo_lung.rds")

log1p_k11 <- readr::read_rds(
  "~/Documents/data/log1pNMF/experiment_results/log1p_pois_k11_mouse_embryo_lung.rds"
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
av1d <- cells_av1$day
names(av1d) <- rownames(cells_av1)
structure_plot(LL_av1, grouping = av1d)

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

readr::write_csv(counts_sum, "~/Documents/data/log1pNMF/mouse_dev_lung_sum.csv")
cells$bc <- rownames(cells)
rownames(cells) <- NULL
readr::write_csv(cells, "~/Documents/data/log1pNMF/mouse_dev_lung_cells.csv")
genes <- genes %>% dplyr::filter(rownames(genes) %in% colnames(counts))
genes$ensembl <- rownames(genes)
rownames(genes) <- NULL
readr::write_csv(genes, "~/Documents/data/log1pNMF/mouse_dev_lung_genes.csv")


cnmf_usage <- readr::read_csv("/Users/eweine/Documents/data/log1pNMF/mouse_lung_dev_cnmf_usage.csv")
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


