library(dplyr)
dat <- readr::read_rds("~/Downloads/Datasets/Liu/pure_select.rds")

ct <- readr::read_csv("~/Downloads/Datasets/Liu/Liu_purified_celltype.csv")

dmats <- lapply(dat, as.matrix)
dmat <- do.call(rbind, dmats)

counts <- as(dmat, "CsparseMatrix")

rm(dmats, dmat, dat)
gc()

save(counts, ct, file = "~/Documents/data/passPCA/liu_data.Rdata")

# now, the task is just to run the algorithm across
# different values of c. I'm not sure how small I should let
# c get.

counts <- counts[,Matrix::colSums(counts) > 0]
# require that all used genes appear in at least 5 cells
s <- Matrix::rowSums(counts)
s <- s / mean(s)
genes_to_use <- which(Matrix::colSums(counts>0)>4)
counts <- counts[,genes_to_use]

fit_list <- list()
n <- nrow(counts)
p <- ncol(counts)

nmf_k1 <- fastTopics:::fit_pnmf_rank1(counts)

K <- 25

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
  nmf_fit, "~/Documents/data/passPCA/liu_results/nmf_k25.rds"
)

library(passPCA)
library(Matrix)
cc_vec <- c(0.0001, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

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
    fit, glue::glue("~/Documents/data/passPCA/liu_results/log1p_c{cc}_k25.rds")
  )

}

fit_list <- list()
fit_list[["nmf"]] <- readr::read_rds("~/Documents/data/passPCA/liu_results/nmf_k25.rds")

for (cc in cc_vec) {

  fit_list[[as.character(cc)]] <- readr::read_rds(
    glue::glue("~/Documents/data/passPCA/liu_results/log1p_c{cc}_k25.rds")
  )

}

readr::write_rds(fit_list, "~/Documents/data/passPCA/liu_results/fit_list.rds")

# there are a variety of analyses I should run here
# the first two I would like to look at are the sparsity and
# log-likelihood for each model

Y <- as.matrix(counts)

ll_vec <- c()

for (cc in cc_vec) {

  fit <- fit_list[[as.character(cc)]]
  B <- fit$U %*% t(fit$V)
  Lambda <- cc * (exp(B) - 1)
  Lambda <- as.matrix(Matrix::Diagonal(x = s) %*% Lambda)

  ll <- sum(
    dpois(
      x = as.vector(Y),
      lambda = as.vector(Lambda),
      log = TRUE
    )
  )

  ll_vec <- c(ll_vec, ll)

}

ll_nmf <- sum(
  dpois(
    x = as.vector(Y),
    lambda = as.vector(fit_list$nmf$L %*% t(fit_list$nmf$F)),
    log = TRUE
  )
)

df <- data.frame(
  cc = cc_vec,
  ll_diff = max(ll_vec) - ll_vec + 1
)

ll_list <- list(
  df = df,
  nmf_ll = ll_nmf
)

ll_list$ll_vec <- ll_vec

readr::write_rds(ll_list, "~/Documents/data/passPCA/liu_results/loglik.rds")

library(ggplot2)
ggplot(data = df) +
  geom_point(aes(x = cc, y = ll_diff)) +
  geom_line(aes(x = cc, y = ll_diff)) +
  xlab("c") +
  ylab("Distance from best log-likelihood") +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept = max(ll_vec) - ll_nmf + 1, color = "red", linetype = "dashed")


#now, I want to examine the sparsity here...

normalize_bars <- function(LL) {

  max_col <- apply(LL, 2, max)
  sweep(LL, 2, max_col, FUN = "/")

}

L_sparsity_vec <- c()

for (cc in cc_vec) {

  fit <- fit_list[[as.character(cc)]]
  LL <- normalize_bars(fit$U)
  L_sparsity <- mean(LL < 1e-5)
  L_sparsity_vec <- c(L_sparsity_vec, L_sparsity)


}

nmf_sparsity <- mean(normalize_bars(fit_list$nmf$L) < 1e-5)
df <- data.frame(
  cc = cc_vec,
  l_sparsity = L_sparsity_vec
)

library(ggplot2)
ggplot(data = df) +
  geom_point(aes(x = cc, y = l_sparsity)) +
  geom_line(aes(x = cc, y = l_sparsity)) +
  xlab("c") +
  ylab("Sparsity of Loadings") +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept = nmf_sparsity, color = "red", linetype = "dashed") +
  cowplot::theme_cowplot()

F_sparsity_vec <- c()

for (cc in cc_vec) {

  fit <- fit_list[[as.character(cc)]]
  FF <- normalize_bars(fit$V)
  F_sparsity <- mean(FF < 1e-5)
  F_sparsity_vec <- c(F_sparsity_vec, F_sparsity)


}

nmf_sparsity <- mean(normalize_bars(fit_list$nmf$F) < 1e-5)
df <- data.frame(
  cc = cc_vec,
  f_sparsity = F_sparsity_vec
)

library(ggplot2)
ggplot(data = df) +
  geom_point(aes(x = cc, y = f_sparsity)) +
  geom_line(aes(x = cc, y = f_sparsity)) +
  xlab("c") +
  ylab("Sparsity of Factors") +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept = nmf_sparsity, color = "red", linetype = "dashed") +
  cowplot::theme_cowplot()

# this is all quite interesting to me...
# I next need to make sure that I can wrap all of this in a notebook
# and then get the loadings and the GO terms



