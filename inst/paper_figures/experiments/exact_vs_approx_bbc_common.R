# Shared setup for the exact-vs-approx timing experiment on the BBC news corpus
# (a documents x terms document-term matrix; a topic-model-style example).
#
# Unlike pancreas / MOCA, BBC is not stored as a saved matrix: the document-term
# matrix is rebuilt at run time from the raw CSV via the tm / SnowballC pipeline
# (identical to inst/scratch/bbc_approx_exact_comp.R and fit_bbc_models.R). The
# resulting sparse matrix `dtm2` is already documents x terms = observations x
# features, so no transpose is needed. Sourced by the job runner, which calls
# run_job(scheme, method).
#
# REQUIRES on the server: the raw CSV at csv_path, and packages tm, SnowballC,
# readr, Matrix. Set csv_path to wherever the CSV lives.

library(Matrix)
library(log1pNMF)

## ---- configuration ----------------------------------------------------------
csv_path       <- "/rafalab/eweine/log1pNMF/inst/paper_figures/data/raw_data/bbc_news_text_complexity_summarization.csv"
K              <- 10        # rank of the factorization (matches the BBC analysis)
cc             <- 1         # link-function tuning parameter c
maxiter_approx <- 100000L   # high safety cap; max_time / tol are the real stops
maxiter_exact  <- 100000L   # high safety cap; max_time / tol are the real stops
init_maxiter   <- 5         # rank-1 warm-up iterations (package default)
seed           <- 1
n_threads      <- 64        # keep identical across all datasets for consistency
out_dir        <- "/home/ericweine/log1p_experiments/"   # where the .rds outputs go
out_tag        <- sprintf("exact_vs_approx_bbc_K%d_c%s", K, format(cc))
tol            <- 1e-8      # a fit that converges within the budget stops early
max_time       <- 10 * 3600 # stop each fit after 10 h of optimization time

## ---- build the document-term matrix from the raw CSV -------------------------
message("Building BBC document-term matrix from ", csv_path)
dat <- readr::read_csv(csv_path, show_col_types = FALSE)

suppressPackageStartupMessages({ library(tm); library(SnowballC) })
set.seed(seed)
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
## keep words appearing in > 4 documents (matches the BBC analysis)
dtm2 <- dtm2[, which(Matrix::colSums(dtm2 > 0) > 4)]

Y <- as(dtm2, "CsparseMatrix")
## drop any all-zero documents / terms the fitter rejects
Y <- Y[Matrix::rowSums(Y) > 0, , drop = FALSE]
Y <- Y[, Matrix::colSums(Y) > 0, drop = FALSE]
n <- nrow(Y)
p <- ncol(Y)
message(sprintf("Data: %d documents x %d terms, %s nonzeros",
                n, p, format(length(Y@x), big.mark = ",")))

config <- list(
  csv_path = csv_path, K = K, cc = cc,
  maxiter_approx = maxiter_approx, maxiter_exact = maxiter_exact,
  init_maxiter = init_maxiter, seed = seed,
  n_threads = n_threads, tol = tol, max_time = max_time, n = n, p = p
)

## shared fitting/init/trace machinery (ctrl, fit_one, run_job, ...)
source("exact_vs_approx_helpers.R")
