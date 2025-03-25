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

library(fastTopics)
library(passPCA)
library(Matrix)
library(dplyr)

K <- 10

fit <- readr::read_rds(
  glue::glue(
    "~/Documents/data/bbc_log1p_c0.001_k{K}_exact_1000_iter.rds"
  )
)

normalized_structure_plot(fit, grouping = dat$labels)

LL <- passPCA:::normalize_bars(fit$LL)
dat$doc_id <- 1:nrow(dat)
dat$k7 <- LL[,"k_7"]
dat <- dat %>% dplyr::select(text, labels, k7, doc_id)

