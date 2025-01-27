library(dplyr)
library(Matrix)
RhpcBLASctl::blas_set_num_threads(7)

cl_to_use <- c(
  "COV434_OVARY", "UMUC1_URINARY_TRACT", "SH10TC_STOMACH", "COLO680N_OESOPHAGUS",
  "SKMEL2_SKIN", "NCIH226_LUNG", "CCFSTTG1_CENTRAL_NERVOUS_SYSTEM"
)

Y_bor <- Matrix::readMM("~/Downloads/Bortezomib_24hr_expt1/matrix.mtx")
bcs_bor <- readr::read_csv("~/Downloads/Bortezomib_24hr_expt1/barcodes.tsv", col_names = c("bc"))
genes_bor <- readr::read_csv("~/Downloads/Bortezomib_24hr_expt1/genes.tsv", col_names = c("gene"))
colnames(Y_bor) <- bcs_bor$bc
rownames(Y_bor) <- genes_bor$gene
md_bor <- readr::read_csv("~/Downloads/Bortezomib_24hr_expt1/classifications.csv")
md_bor <- md_bor %>% dplyr::filter(
  cell_quality == "normal" &
    singlet_ID %in% cl_to_use
  )
Y_bor <- Y_bor[,md_bor$barcode]
md_bor$treatment <- "bortezomib"


library(dplyr)

Y_tra <- Matrix::readMM("~/Downloads/Trametinib_24hr_expt1/matrix.mtx")
bcs_tra <- readr::read_csv("~/Downloads/Trametinib_24hr_expt1/barcodes.tsv", col_names = c("bc"))
genes_tra <- readr::read_csv("~/Downloads/Trametinib_24hr_expt1/genes.tsv", col_names = c("gene"))
colnames(Y_tra) <- bcs_tra$bc
rownames(Y_tra) <- genes_tra$gene
md_tra <- readr::read_csv("~/Downloads/Trametinib_24hr_expt1/classifications.csv")
md_tra <- md_tra %>% dplyr::filter(
  cell_quality == "normal" &
    singlet_ID %in% cl_to_use
)
Y_tra <- Y_tra[,md_tra$barcode]
md_tra$treatment <- "trametinib"

Y_ida <- Matrix::readMM("~/Downloads/Idasanutlin_24hr_expt1/matrix.mtx")
bcs_ida <- readr::read_csv("~/Downloads/Idasanutlin_24hr_expt1/barcodes.tsv", col_names = c("bc"))
genes_ida <- readr::read_csv("~/Downloads/Idasanutlin_24hr_expt1/genes.tsv", col_names = c("gene"))
colnames(Y_ida) <- bcs_ida$bc
rownames(Y_ida) <- genes_ida$gene
md_ida <- readr::read_csv("~/Downloads/Idasanutlin_24hr_expt1/classifications.csv")
md_ida <- md_ida %>% dplyr::filter(
  cell_quality == "normal" &
    singlet_ID %in% cl_to_use
)
Y_ida <- Y_ida[,md_ida$barcode]
md_ida$treatment <- "idasanutlin"

Y_unt <- Matrix::readMM("~/Downloads/DMSO_24hr_expt1/matrix.mtx")
bcs_unt <- readr::read_csv("~/Downloads/DMSO_24hr_expt1/barcodes.tsv", col_names = c("bc"))
genes_unt <- readr::read_csv("~/Downloads/DMSO_24hr_expt1/genes.tsv", col_names = c("gene"))
colnames(Y_unt) <- bcs_unt$bc
rownames(Y_unt) <- genes_unt$gene
md_unt <- readr::read_csv("~/Downloads/DMSO_24hr_expt1/classifications.csv")
md_unt <- md_unt %>% dplyr::filter(
  cell_quality == "normal" &
    singlet_ID %in% cl_to_use
)
Y_unt <- Y_unt[,md_unt$barcode]
md_unt$treatment <- "dmso"

Y <- cbind(
  as.matrix(Y_unt), as.matrix(Y_tra), as.matrix(Y_bor), as.matrix(Y_ida)
)
Y <- as(Y, "CsparseMatrix")
Y <- Y[Matrix::rowSums(Y) > 0, ]

md <- rbind(md_unt, md_tra, md_bor, md_ida)
Y <- Matrix::t(Y)

genes_to_use <- which(Matrix::colSums(Y>0)>4)
Y <- Y[, genes_to_use]

library(fastTopics)
library(passPCA)
s <- Matrix::rowSums(Y)
s <- s / mean(s)

K <- 14
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

n <- nrow(Y)
p <- ncol(Y)

for (cc in cc_vec) {

  print(cc)

  set.seed(1)
  log1p_k1 <- fit_factor_model_log1p_exact(
    Y = Y,
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
    Y = Y,
    K = K,
    init_U = init_LL,
    init_V = init_FF,
    maxiter = 100,
    s = cc * s
  )
  total_time <- tictoc::toc()

  fit[["total_time"]] <- total_time$toc
  rownames(fit$U) <- rownames(Y)
  rownames(fit$V) <- colnames(Y)

  readr::write_rds(
    fit, glue::glue("~/Documents/data/passPCA/cancer/cancer_log1p_c{cc}_k14_exact_100_iter_intercept.rds")
  )

}

fit0_nmf <- fastTopics:::fit_pnmf_rank1(Y)

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

rownames(init_LL) <- rownames(Y)
rownames(init_FF) <- colnames(Y)

fit0_K <- init_poisson_nmf(
  X = Y, F = init_FF, L = init_LL
)

fit_nmf <- fit_poisson_nmf(
  X = Y,
  fit0 = fit0_K,
  control = list(list(nc = 6))
)

readr::write_rds(
  fit_nmf, glue::glue("~/Documents/data/passPCA/cancer/cancer_pois_nmf_k14_exact_100_iter_intercept.rds")
)




normalize_bars <- function(LL) {

  max_col <- apply(LL, 2, max)
  sweep(LL, 2, max_col, FUN = "/")

}

md <- md %>%
  dplyr::mutate(cancer = sub("^[^_]*_", "", singlet_ID))

structure_plot(fit_nmf, grouping = paste(md$cancer, md$treatment), gap = 20)

fit_c1 <- readr::read_rds("~/Documents/data/passPCA/cancer/cancer_log1p_c1_k14_exact_100_iter_intercept.rds")

structure_plot(
  normalize_bars(fit_c1$U),
  grouping = paste(md$cancer, md$treatment),gap = 20,
)

