load("~/Downloads/mouse_brain_stim.Rdata")

library(dplyr)

set.seed(1)
cells <- cells %>% dplyr::filter(!is.na(maintype))


counts <- counts[rownames(counts) %in% cells$...1, ]
#counts <- counts[, Matrix::colSums(counts) > 0]
counts <- counts[Matrix::rowSums(counts) > 0, ]

counts <- counts[sample(x = rownames(counts), size = 7500, replace = FALSE), ]
counts <- counts[, Matrix::colSums(counts) > 0]
counts <- as(counts, "CsparseMatrix")

n <- nrow(counts)
p <- ncol(counts)
K <- 10

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
log1p_k10 <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = K,
  maxiter = 100,
  approx_range = c(0, 1.25),
  s = s,
  init_U = init_LL,
  init_V = init_FF
)

readr::write_rds(log1p_k10, "~/Documents/data/passPCA/experiment_results/mouse_brain_k10.rds")

cells_sub <- cells %>%
  dplyr::filter(
    `...1` %in% rownames(counts)
  )

library(fastTopics)

normalize_bars <- function(LL) {

  max_col <- apply(LL, 2, max)
  sweep(LL, 2, max_col, FUN = "/")

}

LL <- log1p_k10$U

LL <- normalize_bars(LL)

ct <- cells_sub$maintype
names(ct) <- cells_sub$...1
ct <- ct[rownames(counts)]

structure_plot(LL, grouping = ct)

flash_mod <- readr::read_rds("~/Documents/data/passPCA/experiment_results/mouse_brain_flash_fit_out.rds")

LL2 <- normalize_bars(flash_mod$LL)

structure_plot(LL2, grouping = ct, gap = 25)


# Y_sum <- Matrix::summary(counts)
# readr::write_csv(Y_sum, "~/Documents/data/mouse_brain_ijx.csv")
#
# cells <- cells %>%
#   dplyr::rename(bc = `...1`)
#
# readr::write_csv(cells, "~/Documents/data/passPCA/mouse_brain_cells.csv")
# readr::write_csv(as.data.frame(colnames(counts)), "~/Documents/data/passPCA/mouse_brain_genes.csv")
log1p_k15 <- readr::read_rds("~/Documents/data/passPCA/experiment_results/mouse_brain_k15.rds")

normalize_bars <- function(LL) {

  max_col <- apply(LL, 2, max)
  sweep(LL, 2, max_col, FUN = "/")

}

LL <- normalize_bars(log1p_k15$U)


structure_plot(LL, grouping = cells$maintype)

cell_stim <- paste(cells$maintype, cells$stim)

structure_plot(LL, grouping = cell_stim, gap = 25)

usage_df <- readr::read_csv(
  "~/Documents/data/passPCA/mouse_brain_cnmf_usage.csv"
)

top_genes_df <- readr::read_csv(
  "~/Documents/data/passPCA/mouse_brain_cnmf_top_usage.csv"
)

cNMF_LL <- usage_df %>%
  dplyr::select(-bc) %>%
  as.matrix()

structure_plot(cNMF_LL, grouping = cell_stim, gap = 25)

rs <- Matrix::rowSums(counts)
s <- rs / mean(rs)
Diag_S <- Matrix::Diagonal(x=1/s)

counts <- Diag_S %*% counts
counts <- MatrixExtra::mapSparse(counts, log1p)

library(RcppML)
nmf_fit <- nmf(
  A = counts, k = 15
)

readr::write_rds(nmf_fit, "~/Documents/data/passPCA/experiment_results/rcppML_nmf_k15.rds")

nmf_LL <- nmf_fit$w %*% diag(x = nmf_fit$d)
nmf_LL <- normalize_bars(nmf_LL)

structure_plot(nmf_LL, grouping = cell_stim, gap = 25)

