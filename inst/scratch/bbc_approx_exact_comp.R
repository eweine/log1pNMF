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

library(log1pNMF)


fit_list <- list()
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

fit_list <- list()

for (cc in cc_vec) {
  
  for (ll in c("exact", "approx")) {
    
    if (ll == "approx") {
      
      if (cc >= 1) {
        
        set.seed(1)
        fit <- fit_poisson_log1p_nmf(
          Y = dtm2, K = 10, cc = cc, loglik = ll, init_method = "rank1",
          control = list(maxiter = 250)
        )
        
        fit$final_ll <- logLik(fit, dtm2)
        fit_list[[ll]][[as.character(cc)]] <- fit
        
      }
      
    } else {
      
      set.seed(1)
      fit <- fit_poisson_log1p_nmf(
        Y = dtm2, K = 10, cc = cc, loglik = ll, init_method = "rank1",
        control = list(maxiter = 250)
      )
      
      fit$final_ll <- logLik(fit, dtm2)
      fit_list[[ll]][[as.character(cc)]] <- fit
      
    }
    
  }
  
}

K <- 10
library(fastTopics)
r1_fit <- fastTopics:::fit_pnmf_rank1(dtm2)

init_LL <- cbind(
  r1_fit$L,
  matrix(
    data = 1e-5,
    nrow = nrow(dtm2),
    ncol = K - 1
  )
)
rownames(init_LL) <- rownames(dtm2)

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

final_fit_list <- fit_list$exact

# use approximate methods here to avoid numerical issues
final_fit_list[[as.character(100)]] <- fit_list$approx$`100`
final_fit_list[[as.character(1000)]] <- fit_list$approx$`1000`

final_fit_list[[as.character(Inf)]] <- nmf_fit

readr::write_rds(final_fit_list, "~/Documents/data/bbc_fit_list.rds")
