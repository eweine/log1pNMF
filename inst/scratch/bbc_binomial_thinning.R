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

library(log1pNMF)
library(Matrix)
library(dplyr)
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

n <- nrow(dtm2)
p <- ncol(dtm2)
K <- 10

create_bin_thin_mats <- function(M) {
  
  x_train <- rbinom(n = length(M@x), size = M@x, prob = 0.5)
  x_test <- M@x - x_train
  
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

n_sims <- 1
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

sim_res <- list()

for (i in 1:n_sims) {
  
  set.seed(i)
  data_split <- create_bin_thin_mats(dtm2)
  
  r1_fit <- fastTopics:::fit_pnmf_rank1(data_split$train)
  
  set.seed(1)
  init_LL <- cbind(
    r1_fit$L,
    matrix(
      data = 1e-5,
      nrow = nrow(data_split$train),
      ncol = K - 1
    )
  )
  rownames(init_LL) <- rownames(data_split$train)
  
  set.seed(1)
  init_FF <- cbind(
    r1_fit$F,
    matrix(
      data = 1e-5,
      nrow = ncol(data_split$train),
      ncol = K - 1
    )
  )
  rownames(init_FF) <- colnames(data_split$train)
  
  fit0 <- init_poisson_nmf(X = data_split$train, L = init_LL, F = init_FF)
  
  nmf_fit <- fit_poisson_nmf(
    X = data_split$train,
    fit0 = fit0,
    numiter = 250,
    control = list(nc = 7)
  )
  
  train_ll <- tail(nmf_fit$progress$loglik, 1)
  s_train <- Matrix::rowSums(data_split$train)
  s_train <- s_train / mean(s_train)
  s_test <- Matrix::rowSums(data_split$test)
  s_test <- s_test / mean(s_test)
  test_Lambda <- Matrix::Diagonal(x = s_test/s_train) %*% nmf_fit$L %*% t(nmf_fit$F)
  test_ll <- sum(
      dpois(
      as.vector(
        as.matrix(data_split$test)), 
      lambda = as.vector(as.matrix(test_Lambda)), 
      log = TRUE
    )
  )
  
  sim_res[[i]] <- list()
  sim_res[[i]]$tm_train_ll <- train_ll
  sim_res[[i]]$tm_test_ll <- test_ll
  
  for (cc in cc_vec) {
    
    sim_res[[i]][[as.character(cc)]] <- list()
    print(cc)
    set.seed(1)
    fit <- fit_poisson_log1p_nmf(
      Y = data_split$train,
      K = K,
      cc = cc,
      loglik = ifelse(cc >= 100, "approx", "exact"),
      init_method = "rank1",
      control = list(
        maxiter = 250
      )
    )
    
    sim_res[[i]][[as.character(cc)]]$train_ll <- logLik(fit, Y = data_split$train)
    sim_res[[i]][[as.character(cc)]]$test_ll <- logLik(fit, Y = data_split$test, s = s_test)
    
  }
  
}

