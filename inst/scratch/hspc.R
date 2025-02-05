library(dplyr)
library(Matrix)
library(fastTopics)
library(passPCA)

data_dir <- "/home/ericweine/hspc"
counts <- readr::read_rds(glue::glue("{data_dir}/hspcs.rds"))

K <- 10
cc_vec <- c(1e-3, 1)
s <- Matrix::rowSums(counts)
s <- s / mean(s)

n <- nrow(counts)
p <- ncol(counts)

# turn off BLAS threads and increase rcpp Parallel threads to 48
RcppParallel::setThreadOptions(numThreads = 48)
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
    fit, glue::glue("{data_dir}/hspc_log1p_c{cc}_k{K}_exact_100_iter.rds")
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
  control = list(list(nc = 48))
)

readr::write_rds(
  fit_nmf, glue::glue("{data_dir}/hspc_pois_nmf_k{K}_exact_100_iter.rds")
)
