#dat <- readr::read_csv("/home/ericw456/bbc/bbc_news_text_complexity_summarization.csv")
dat <- readr::read_csv("~/Downloads/bbc_news_text_complexity_summarization.csv")
library(tm)
library(SnowballC)

my_corpus <- VCorpus(VectorSource(dat$text))

addspace <- content_transformer(function(x, pattern) {
  return(gsub(pattern, " ", x))
})

my_corpus <- tm_map(my_corpus, addspace, "-")

my_corpus <- tm_map(my_corpus, removeNumbers)

# Transform to lower case (need to wrap in content_transformer)
my_corpus <- tm_map(my_corpus,content_transformer(tolower))
my_corpus <- tm_map(my_corpus, removeWords, stopwords("SMART"))
my_corpus <- tm_map(my_corpus, removePunctuation)
my_corpus <- tm_map(my_corpus, stripWhitespace)

my_corpus <- tm_map(my_corpus, stemDocument)

dtm <- DocumentTermMatrix(my_corpus)
counts <- Matrix::sparseMatrix(
  i = dtm$i,
  j = dtm$j,
  x = dtm$v
)


colnames(counts) <- dtm$dimnames$Terms
words_to_use <- which(Matrix::colSums(counts>0)>4)

counts <- counts[,words_to_use]

s <- Matrix::rowSums(counts)
s <- s / mean(s)

library(log1pNMF)
library(Matrix)
library(dplyr)

K <- 25
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

  fit[["total_time"]] <- total_time$toc
  rownames(fit$U) <- rownames(counts)
  rownames(fit$V) <- colnames(counts)

  readr::write_rds(
    fit, glue::glue("~/Documents/data/log1pNMF/bbc_log1p_c{cc}_k25_approx_taylor_100_iter.rds")
  )

}
