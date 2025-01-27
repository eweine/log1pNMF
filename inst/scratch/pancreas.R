library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(passPCA)
library(dplyr)
library(fastTopics)

set.seed(1)

load("~/Downloads/pancreas.RData")

i           <- which(sample_info$tech == "celseq2")
sample_info <- sample_info[i,]
counts      <- counts[i,]
sample_info <- transform(sample_info,celltype = factor(celltype))

counts <- counts[, Matrix::colSums(counts) > 0]
genes_to_use <- which(Matrix::colSums(counts>0)>4)
counts <- counts[,genes_to_use]
s <- Matrix::rowSums(counts)
s <- s / mean(s)

K <- 9
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

  tictoc::tic()
  set.seed(1)
  fit <- fit_factor_model_log1p_exact(
    Y = counts,
    K = K,
    init_U = init_LL,
    init_V = init_FF,
    maxiter = 100,
    s = cc * s
  )
  total_time <- tictoc::toc()

  fit[["total_time"]] <- total_time$toc
  rownames(fit$U) <- rownames(counts)
  rownames(fit$V) <- colnames(counts)

  readr::write_rds(
    fit, glue::glue("~/Documents/data/passPCA/pancreas/pancreas_log1p_c{cc}_k9_exact_100_iter.rds")
  )

}

fit0_nmf <- fastTopics:::fit_pnmf_rank1(counts)

set.seed(1)
init_LL <- fit0_nmf$L %>%
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
init_FF <- fit0_nmf$F %>%
  cbind(
    matrix(
      data = rexp(
        n = p * (K - 1), rate = 15
      ),
      nrow = p,
      ncol = K - 1
    )
  )

rownames(init_LL) <- rownames(counts)
rownames(init_FF) <- colnames(counts)

fit0_K <- init_poisson_nmf(
  X = counts, F = init_FF, L = init_LL
)

fit_nmf <- fit_poisson_nmf(
  X = counts,
  fit0 = fit0_K,
  control = list(list(nc = 8))
)

readr::write_rds(
  fit_nmf, glue::glue("~/Documents/data/passPCA/pancreas/pancreas_pois_nmf_k9_exact_100_iter.rds")
)

library(fastTopics)
fit_list <- list()

for (cc in cc_vec) {

  fit_list[[as.character(cc)]] <- readr::read_rds(
    glue::glue(
      "~/Documents/data/passPCA/pancreas/pancreas_log1p_c{cc}_k9_exact_100_iter.rds"
    )
  )

}

fit_list[["Inf"]] <- readr::read_rds(
  glue::glue("~/Documents/data/passPCA/pancreas/pancreas_pois_nmf_k9_exact_100_iter.rds")
)

normalize_bars <- function(LL) {

  max_col <- apply(LL, 2, max)
  sweep(LL, 2, max_col, FUN = "/")

}

celltype <- sample_info$celltype
celltype <-
  factor(celltype,
         c("acinar","ductal","activated_stellate","quiescent_stellate",
           "endothelial","macrophage","mast","schwann","alpha","beta",
           "delta","gamma","epsilon"))


plot_list <- list()

for (cc in cc_vec) {

  plot_list[[as.character(cc)]] <- structure_plot(
    normalize_bars(fit_list[[as.character(cc)]]$U),
    grouping = celltype,gap = 20,perplexity = 70,n = Inf
    )

}

plot_list[["Inf"]] <- structure_plot(
  fit_list[["Inf"]],
  grouping = celltype,gap = 20,perplexity = 70,n = Inf
  )

ggarrange(
  plotlist = plot_list,
  ncol = 1,
  labels = paste0("K = ", names(plot_list))
)
