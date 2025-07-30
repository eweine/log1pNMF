dat <- readr::read_csv("../data/raw_data/bbc_news_text_complexity_summarization.csv")

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
dtm2 <- Matrix::sparseMatrix(
  i = dtm$i,
  j = dtm$j,
  x = dtm$v
)


colnames(dtm2) <- dtm$dimnames$Terms
words_to_use <- which(Matrix::colSums(dtm2>0)>4)

dtm2 <- dtm2[,words_to_use]

s <- Matrix::rowSums(dtm2)
s <- s / mean(s)

library(passPCA)
library(Matrix)
library(dplyr)
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

n <- nrow(dtm2)
p <- ncol(dtm2)

create_bin_thin_mats <- function(M) {
  
  x_train <- rbinom(n = length(M@x), size = M@x, prob = 0.5)
  x_test <- M@x - x_train
  
  browser()
  M_train <- M
  M_train@x <- as.numeric(x_train)
  M_train <- drop0(M_train)
  
  M_test <- M
  M_test@x <- as.numeric(x_test)
  M_test <- drop0(M_test)
  
  cs_train <- Matrix::colSums(M_train)
  cs_test <- Matrix::colSums(M_test)
  
  nz_idx <- which(cs_train > 0 & cs_test > 0)
  M_train <- M_train[, nz_idx]
  M_test <- M_test[, nz_idx]
  
  return(
    list(
      train = M_train,
      test = M_test
    )
  )
  
}

n_sims <- 10

for (i in 1:n_sims) {
  
  set.seed(i)
  data_split <- create_bin_thin_mats(dtm2)
  
  r1_fit <- fastTopics:::fit_pnmf_rank1(dtm2)
  
  set.seed(1)
  init_LL <- cbind(
    r1_fit$L,
    matrix(
      data = 1e-5,
      nrow = nrow(dtm2),
      ncol = K - 1
    )
  )
  rownames(init_LL) <- rownames(dtm2)
  
  set.seed(1)
  init_FF <- cbind(
    r1_fit$F,
    matrix(
      data = 1e-5,
      nrow = ncol(dtm2),
      ncol = K - 1
    )
  )
  rownames(init_FF) <- colnames(dtm2)
  
  fit0 <- init_poisson_nmf(X = dtm2, L = init_LL, F = init_FF)
  
  nmf_fit <- fit_poisson_nmf(
    X = dtm2,
    fit0 = fit0,
    numiter = 250,
    control = list(nc = 7)
  )
  
}

data_split <- create_bin_thin_mats(dtm2)
