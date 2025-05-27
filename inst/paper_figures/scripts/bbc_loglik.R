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

library(log1pNMF)
library(Matrix)
library(dplyr)

fit_list <- readr::read_rds(
  "~/Documents/data/log1pNMF/bbc_articles.rds"
)

cc_vec <- c(0.0001, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4)

ll_vec <- c()
ll_vec_approx <- c()

for (cc in cc_vec) {

  fit <- fit_list[[as.character(cc)]]
  B <- fit$U %*% t(fit$V)
  Lambda <- cc * (exp(B) - 1)
  Lambda <- as.matrix(Matrix::Diagonal(x = s) %*% Lambda)

  ll <- sum(
    dpois(
      x = as.vector(as.matrix(dtm2)),
      lambda = as.vector(Lambda),
      log = TRUE
    )
  )

  ll_vec <- c(ll_vec, ll)

  fit <- fit_list$approx[[as.character(cc)]]
  B <- fit$U %*% t(fit$V)
  Lambda <- cc * (exp(B) - 1)
  Lambda <- as.matrix(Matrix::Diagonal(x = s) %*% Lambda)

  ll <- sum(
    dpois(
      x = as.vector(as.matrix(dtm2)),
      lambda = as.vector(Lambda),
      log = TRUE
    )
  )
  ll_vec_approx <- c(ll_vec_approx, ll)

}

ll_nmf <- sum(
  dpois(
    x = as.vector(as.matrix(dtm2)),
    lambda = as.vector(fit_list$nmf$L %*% t(fit_list$nmf$F)),
    log = TRUE
  )
)

df <- data.frame(
  cc = cc_vec,
  ll_diff = c(
    max(ll_vec) - ll_vec + 1
  ),
  model = c(
    rep("log1p Exact", length(cc_vec))
  )
)

library(ggplot2)
ga <- ggplot(data = df) +
  geom_point(aes(x = cc, y = ll_diff)) +
  geom_line(aes(x = cc, y = ll_diff)) +
  xlab("c") +
  ylab("Distance from best log-likelihood") +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept = max(ll_vec) - ll_nmf + 1, color = "red", linetype = "dashed") +
  cowplot::theme_cowplot()

df <- data.frame(
  cc = cc_vec,
  ll_diff = pmax(c(
    ll_vec - ll_vec_approx
  ), 1)
)

library(ggplot2)
gb <- ggplot(data = df) +
  geom_point(aes(x = cc, y = ll_diff)) +
  geom_line(aes(x = cc, y = ll_diff)) +
  xlab("c") +
  ylab("Approximation error (log-likelihood)") +
  scale_x_log10() +
  scale_y_log10() +
  cowplot::theme_cowplot()

library(ggpubr)

f2 <- ggarrange(
  ga, gb, labels = "AUTO"
)

ggsave(
  "/Users/eweine/Documents/log1pNMF/inst/paper_figures/pdfs/bbc_loglik.pdf",
  f2,
  device = "pdf",
  width = 8,
  height = 4
)
