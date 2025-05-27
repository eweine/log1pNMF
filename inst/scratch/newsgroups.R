load("~/Downloads/newsgroups.RData")
library(dplyr)
library(Matrix)
library(fastTopics)
library(log1pNMF)
rm(counts)

# I did the following to the original dataset
# counts <- counts[, !(colnames(counts) %in% stopwords("SMART"))]
# colnames(counts) <- wordStem(colnames(counts))

counts <- readr::read_rds(
  "~/Downloads/newsgroups_stemmed_stopwords.rds"
)

topics <- topics[Matrix::rowSums(counts) > 9]
counts <- counts[Matrix::rowSums(counts) > 9, ]
genes_to_use <- which(Matrix::colSums(counts>0)>4)
counts <- counts[,genes_to_use]
s <- Matrix::rowSums(counts)
s <- s / mean(s)


K <- 25
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

n <- nrow(counts)
p <- ncol(counts)

# turn off BLAS threads and increase rcpp Parallel threads to 8
RcppParallel::setThreadOptions(numThreads = 8)
RhpcBLASctl::blas_set_num_threads(1)

for (cc in cc_vec) {

  print(cc)

  set.seed(1)
  log1p_k1 <- fit_factor_model_log1p_exact(
    Y = counts,
    K = 1,
    maxiter = 5,
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

  if (cc < 10) {

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

    method <- "exact"

  } else {

    tictoc::tic()
    set.seed(1)
    fit <- fit_factor_model_log1p_quad_approx_sparse(
      Y = counts,
      K = K,
      init_U = init_LL,
      init_V = init_FF,
      maxiter = 100,
      s = cc * s,
      approx_method = "taylor"
    )
    total_time <- tictoc::toc()
    method <- "approx"

  }



  fit[["total_time"]] <- total_time$toc
  rownames(fit$U) <- rownames(counts)
  rownames(fit$V) <- colnames(counts)

  readr::write_rds(
    fit, glue::glue("~/Documents/data/news_log1p_c{cc}_k{K}_{method}_100_iter.rds")
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
  fit_nmf, glue::glue("~/Documents/data/news_pois_nmf_k{K}_exact_100_iter.rds")
)
