load("../data/experiment_results.Rdata")
dat <- readr::read_csv("../data/raw_data/bbc_news_text_complexity_summarization.csv")

set.seed(1)
library(fastTopics)
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
#cc_vec <- c(1e3)

n <- nrow(dtm2)
p <- ncol(dtm2)

K <- 10

hoyer_sparsity <- function(x) {
  
  n <- length(x)
  (1 / (sqrt(n) - 1)) * (sqrt(n) - (sum(x) / (sqrt(sum(x ^ 2)))))
  
}

fit_list <- list()
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

for (cc in cc_vec) {
  
  fit_list[[as.character(cc)]] <- res_list$bbc[[as.character(cc)]]
  
  colnames(fit_list[[as.character(cc)]]$LL) <- paste0("k", 1:K)
  
}

fit_list[["Inf"]] <- res_list$bbc$`Inf`

plot_list <- list()
cc_vec <- c(1e-3)

log1p_loadings_order <- c(1, 2, 5, 4, 3, 8, 7, 6, 9, 10)
tm_loadings_order <- c(1, 8, 3, 9, 5, 10, 7, 2, 4, 6)


colnames(fit_list[[as.character(1e-3)]]$LL) <- paste0(
  "k", log1p_loadings_order
)

fit_list[[as.character(1e-3)]]$LL <- fit_list[[as.character(1e-3)]]$LL[,paste0("k", 1:10)]

colnames(fit_list[[as.character(1e-3)]]$FF) <- paste0(
  "k", log1p_loadings_order
)

fit_list[[as.character(1e-3)]]$FF <- fit_list[[as.character(1e-3)]]$FF[,paste0("k", 1:10)]



colnames(fit_list[[as.character(Inf)]]$L) <- paste0(
  "k", tm_loadings_order
)

fit_list[[as.character(Inf)]]$L <- fit_list[[as.character(Inf)]]$L[,paste0("k", 1:10)]

colnames(fit_list[[as.character(Inf)]]$F) <- paste0(
  "k", tm_loadings_order
)

fit_list[[as.character(Inf)]]$F <- fit_list[[as.character(Inf)]]$F[,paste0("k", 1:10)]

dat$doc_id <- 1:nrow(dat)
dat$k7 <- fit_list[[as.character(1e-3)]]$LL[,"k7"]
dat$tm_k7 <- fit_list[[as.character(Inf)]]$L[,"k7"]
dat <- dat %>% dplyr::select(text, labels, k7, doc_id)
word_vec_full <- (fit_list[[as.character(1e-3)]]$LL[947,] %*% t(fit_list[[as.character(1e-3)]]$FF))[1,]
lambda_vec_full <- 1e-3 * (exp(word_vec_full) - 1)

l_vec_k3_0 <- fit_list[[as.character(1e-3)]]$LL[947,]

l_vec_k3_0["k3"] <- 0

word_vec_k3_0 <- (l_vec_k3_0 %*% t(fit_list[[as.character(1e-3)]]$FF))[1,]
lambda_vec_k3_0 <- 1e-3 * (exp(word_vec_k3_0) - 1)
used_words <- colnames(dtm2)[which(dtm2[947,] > 0)]

diff <- lambda_vec_full - lambda_vec_k3_0
diff_used <- lambda_vec_full[used_words] - lambda_vec_k3_0[used_words]

cc <- 1e-3

LL_explain <- matrix(data = 0, nrow = nrow(fit_list[[as.character(cc)]]$LL), ncol = K)

for (i in 1:nrow(dtm2)) {
  
  word_vec_full <- (
    fit_list[[as.character(cc)]]$LL[i,] %*% 
      t(fit_list[[as.character(cc)]]$FF)
    )[1,]
  lambda_vec_full <- cc * (exp(word_vec_full) - 1)
  rate_vec_full <- lambda_vec_full / sum(lambda_vec_full)
  
  l <- fit_list[[as.character(cc)]]$LL[i,]
  
  for (k in 1:K) {
    
    lk_og <- l[k]
    l[k] <- 0
    
    word_vec_k0 <- (
      l %*% 
        t(fit_list[[as.character(cc)]]$FF)
    )[1,]
    
    lambda_vec_k0 <- cc * (exp(word_vec_k0) - 1)
    rate_vec_k0 <- lambda_vec_k0 / sum(lambda_vec_k0)
    
    l1_diff <- mean(abs(rate_vec_full - rate_vec_k0))
    l[k] <- lk_og
    LL_explain[i, k] <- l1_diff
    
  }
  
}

LL_explain <- LL_explain / Matrix::rowSums(LL_explain)

structure_plot(LL_explain, grouping = dat$labels)

