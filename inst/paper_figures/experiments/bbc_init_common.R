# Shared setup for the BBC initialization-stability experiment (reviewer
# request): how stable are the K = 10 BBC fits to initialization?
#
# For each c in {0.001, 1, Inf}, fit the model
#   (a) from N_RANDOM random initializations (seeds 1..N_RANDOM), and
#   (b) from the rank-1 initialization used in the paper,
# every fit run to convergence under the EXACT log-likelihood (finite c;
# c = Inf is plain Poisson NMF via fastTopics, whose likelihood is exact).
#
# Task layout (one fit per Slurm array task), ordered so that the 8-core
# c = 0.001 tasks form one contiguous block and the 1-core rest another:
#
#   ids  0-30   c = 0.001  (seeds 1-30 random, then rank1)   8 cores each
#   ids 31-61   c = 1      (seeds 1-30 random, then rank1)   1 core each
#   ids 62-92   c = Inf    (seeds 1-30 random, then rank1)   1 core each
#
# The document-term matrix is rebuilt at run time from the raw CSV via the
# tm / SnowballC pipeline -- identical to fit_bbc_models.R and
# exact_vs_approx_bbc_common.R -- and the fits use the same size factors
# s = rowSums / mean(rowSums) as the paper's BBC fits.
#
# REQUIRES on the server: the raw CSV at CSV_PATH, and packages tm,
# SnowballC, readr, Matrix, log1pNMF, fastTopics.

library(Matrix)

## ---- configuration ---------------------------------------------------------
CSV_PATH <- Sys.getenv("BBC_INIT_CSV",
  "/rafalab/eweine/log1pNMF/inst/paper_figures/data/bbc_news_text_complexity_summarization.csv")
OUT_DIR  <- Sys.getenv("BBC_INIT_OUT",
                       "/rafalab/eweine/log1p_experiments/bbc_init_matched")

K_FIT       <- 10L
CC_GRID     <- c(0.001, 1, Inf)
N_RANDOM    <- 30L                 # random-init seeds per value of c
TOL         <- 1e-6                # absolute change in the fit's own objective
MAXITER     <- 1000000L            # safety cap, finite c (convergence is the stop)
MAXITER_INF <- 100000L             # safety cap, fastTopics SCD
RANK1_PAD   <- 1e-5                # padding for the rank-1 init at c = Inf
                                   # (matches fit_bbc_models.R)

## ---- task id <-> (cc, init, seed) ------------------------------------------
## init "random" carries a seed 1..N_RANDOM; init "rank1" has seed NA.
SPEC <- do.call(rbind, lapply(CC_GRID, function(cc) rbind(
  data.frame(cc = cc, init = "random", seed = seq_len(N_RANDOM)),
  data.frame(cc = cc, init = "rank1",  seed = NA_integer_))))
SPEC$id <- seq_len(nrow(SPEC)) - 1L
N_TASKS <- nrow(SPEC)              # 93

fmtc <- function(x) ifelse(is.infinite(x), "Inf",
  format(x, scientific = FALSE, trim = TRUE, drop0trailing = TRUE))

task_spec <- function(id) {
  stopifnot(id >= 0L, id < N_TASKS)
  SPEC[id + 1L, ]
}

tag_of <- function(spec) {
  if (spec$init == "random")
    sprintf("bbcinit_c%s_random_seed%02d", fmtc(spec$cc), spec$seed)
  else
    sprintf("bbcinit_c%s_rank1", fmtc(spec$cc))
}

## ---- build the document-term matrix from the raw CSV -----------------------
## Same pipeline as fit_bbc_models.R / exact_vs_approx_bbc_common.R.
build_bbc <- function(csv_path = CSV_PATH) {
  message("Building BBC document-term matrix from ", csv_path)
  dat <- readr::read_csv(csv_path, show_col_types = FALSE)
  suppressPackageStartupMessages({ library(tm); library(SnowballC) })
  my_corpus <- VCorpus(VectorSource(dat$text))
  addspace  <- content_transformer(function(x, pattern) gsub(pattern, " ", x))
  my_corpus <- tm_map(my_corpus, addspace, "-")
  my_corpus <- tm_map(my_corpus, removeNumbers)
  my_corpus <- tm_map(my_corpus, content_transformer(tolower))
  my_corpus <- tm_map(my_corpus, removeWords, stopwords("SMART"))
  my_corpus <- tm_map(my_corpus, removePunctuation)
  my_corpus <- tm_map(my_corpus, stripWhitespace)
  my_corpus <- tm_map(my_corpus, stemDocument)

  dtm  <- DocumentTermMatrix(my_corpus)
  dtm2 <- Matrix::sparseMatrix(i = dtm$i, j = dtm$j, x = dtm$v)
  colnames(dtm2) <- dtm$dimnames$Terms
  dtm2 <- dtm2[, which(Matrix::colSums(dtm2 > 0) > 4)]

  Y <- as(dtm2, "CsparseMatrix")
  Y <- Y[Matrix::rowSums(Y) > 0, , drop = FALSE]
  Y <- Y[, Matrix::colSums(Y) > 0, drop = FALSE]
  message(sprintf("Data: %d documents x %d terms, %s nonzeros",
                  nrow(Y), ncol(Y), format(length(Y@x), big.mark = ",")))
  s <- Matrix::rowSums(Y)
  s <- s / mean(s)
  list(Y = Y, s = s)
}

pois_loglik <- function(Y, lam) {
  lam <- pmax(lam, 1e-300)
  sum(Y * log(lam) - lam - lgamma(as.matrix(Y) + 1))
}

## ---- matched random initialization -----------------------------------------
## One random factorization per seed, mapped into each model's parameter
## space: L_raw and F_raw have iid Uniform(0, 1) entries and depend ONLY on
## the seed (identical matrices for every c), and both are scaled by
## sqrt(tau_c), with the scalar tau_c solved so that the model's implied
## MEAN RATE at the initial point equals the mean observed count. Every
## model therefore starts from the same random direction, at rates on the
## data's scale; across-c differences cannot come from the init
## distribution. (For finite c the implied rates are s_i * c *
## (exp(tau * m_ij / alpha_c) - 1); for c = Inf they are tau * m_ij.)
matched_random_init <- function(Y, s, cc, K, seed) {
  set.seed(seed)                       # seed alone determines the draws
  L_raw <- matrix(stats::runif(nrow(Y) * K), nrow(Y), K)
  F_raw <- matrix(stats::runif(ncol(Y) * K), ncol(Y), K)
  M     <- tcrossprod(L_raw, F_raw)
  ybar  <- sum(Y) / prod(dim(Y))
  if (is.infinite(cc)) {
    tau <- ybar / mean(M)
  } else {
    alpha <- max(1, cc)
    ## clamped so overflow cannot produce non-finite values; the clamp is
    ## far above any plausible root and preserves monotonicity
    f <- function(tau)
      mean(s * cc * expm1(pmin(tau * M / alpha, 700))) - ybar
    upper <- 1
    while (f(upper) < 0) upper <- 2 * upper
    tau <- stats::uniroot(f, c(0, upper), tol = 1e-10)$root
  }
  L0 <- sqrt(tau) * L_raw
  F0 <- sqrt(tau) * F_raw
  rownames(L0) <- rownames(Y); rownames(F0) <- colnames(Y)
  colnames(L0) <- colnames(F0) <- paste0("k", seq_len(K))
  list(L = L0, F = F0, tau = tau)
}
