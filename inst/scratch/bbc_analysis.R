dat <- readr::read_csv("~/Downloads/bbc_news_text_complexity_summarization.csv")

library(tm)
library(SnowballC)

my_corpus <- VCorpus(VectorSource(dat$text))

addspace <- content_transformer(function(x, pattern) {
  return(gsub(pattern, " ", x))
})

my_corpus <- tm_map(my_corpus, addspace, "-")

my_corpus <- tm_map(my_corpus, removePunctuation)
my_corpus <- tm_map(my_corpus, removeNumbers)
my_corpus <- tm_map(my_corpus, removeWords, stopwords("english"))

# Transform to lower case (need to wrap in content_transformer)
my_corpus <- tm_map(my_corpus,content_transformer(tolower))
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
cc_vec <- c(1e-4)
fit_list <- list()
n <- nrow(dtm2)
p <- ncol(dtm2)

K <- 50

for (cc in cc_vec) {

  print(cc)

  set.seed(1)
  log1p_k1 <- fit_factor_model_log1p_exact(
    Y = dtm2,
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
    Y = dtm2,
    K = K,
    init_U = init_LL,
    init_V = init_FF,
    maxiter = 100,
    s = cc * s
  )

}



library(fastTopics)
nmf_k1 <- fastTopics:::fit_pnmf_rank1(X = dtm2)
L0 <- nmf_k1$L %>%
  cbind(
    matrix(
      data = rexp(
        n = n * (K - 1), rate = 15
      ),
      nrow = n,
      ncol = K - 1
    )
  )
F0 <- nmf_k1$F %>%
  cbind(
    matrix(
      data = rexp(
        n = p * (K - 1), rate = 15
      ),
      nrow = p,
      ncol = K - 1
    )
  )

rownames(L0) <- rownames(dtm2)
rownames(F0) <- colnames(dtm2)

nmf_fit0 <- fastTopics::init_poisson_nmf(
  X = dtm2,
  F = F0,
  L = L0
)
nmf <- fit_poisson_nmf(X = dtm2, k = 50, control = list(nc = 4))

fit_list <- readr::read_rds("~/Documents/data/passPCA/bbc_articles.rds")
fit_list[["nmf"]] <- nmf
fit_list[[as.character(1e-4)]] <- fit
readr::write_rds(fit_list, "~/Documents/data/passPCA/bbc_articles.rds")
