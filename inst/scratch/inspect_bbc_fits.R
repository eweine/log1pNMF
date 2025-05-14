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

library(passPCA)
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
  
  fit_list[[as.character(cc)]]$l_sparsity <- apply(
    fit_list[[as.character(cc)]]$LL, 2, hoyer_sparsity
  )
  
  fit_list[[as.character(cc)]]$f_sparsity <- apply(
    fit_list[[as.character(cc)]]$FF, 2, hoyer_sparsity
  )
  
  fit_list[[as.character(cc)]]$cor_mat <- cor(fit_list[[as.character(cc)]]$FF, method = "spearman")
  
}

fit_list[["Inf"]] <- res_list$bbc$`Inf`

fit_list[["Inf"]]$Ls <- Matrix::Diagonal(x = 1/s) %*% fit_list[["Inf"]]$L

fit_list[["Inf"]]$l_sparsity <- apply(
  fit_list[["Inf"]]$Ls, 2, hoyer_sparsity
)

fit_list[["Inf"]]$f_sparsity <- apply(
  fit_list[["Inf"]]$F, 2, hoyer_sparsity
)

fit_list[["Inf"]]$cor_mat <- cor(fit_list[["Inf"]]$F, method = "spearman")

l_sparsity_vec <- unlist(lapply(fit_list, function(x) {median(x$l_sparsity)}))
f_sparsity_vec <- unlist(lapply(fit_list, function(x) {median(x$f_sparsity)}))
cor_vec <- unlist(
  lapply(fit_list, function(x) {median(abs(x$cor_mat[lower.tri(x$cor_mat)]))})
)

df_sparsity_l <- data.frame(
  cc = as.numeric(names(l_sparsity_vec)),
  sparsity = l_sparsity_vec
) %>% filter(is.finite(cc))

df_sparsity_f <- data.frame(
  cc = as.numeric(names(f_sparsity_vec)),
  sparsity = f_sparsity_vec
) %>% filter(is.finite(cc))

df_cor <- data.frame(
  cc = as.numeric(names(cor_vec)),
  correlation = cor_vec
) %>% filter(is.finite(cc))

library(ggpubr)
library(ggplot2)

g1 <- ggplot(data = df_cor, aes(x = cc, y = correlation)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-3, 1, 1e3), transform = "log10") +
  xlab("c (log10 scale)") +
  ylab("Median Abs. Factor Correlation") +
  geom_hline(yintercept = cor_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.004, y=cor_vec["Inf"] + 0.05, label="ID Link", color="red"
  )

g2 <- ggplot(data = df_sparsity_l, aes(x = cc, y = sparsity)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-3, 1, 1e3), transform = "log10") +
  xlab("c (log10 scale)") +
  ylab("Median Loading Sparsity") +
  geom_hline(yintercept = l_sparsity_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.004, y=l_sparsity_vec["Inf"] + 0.05, label="ID Link", color="red"
  )

g3 <- ggplot(data = df_sparsity_f, aes(x = cc, y = sparsity)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-3, 1, 1e3), transform = "log10") +
  xlab("c (log10 scale)") +
  ylab("Median Factor Sparsity") +
  geom_hline(yintercept = f_sparsity_vec["Inf"], color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.004, y=f_sparsity_vec["Inf"] + 0.05, label="ID Link", color="red"
  )


g <- ggarrange(g1,g3,g2, nrow = 1, labels = "AUTO")

g <- annotate_figure(g,
                     top = text_grob("BBC K = 10", size = 20, face = "bold"))

plot_list <- list()
cc_vec <- c(1e-3)

colnames(fit_list[[as.character(1e-3)]]$LL) <- paste0(
  "k", c(1, 2, 4, 6, 9, 5, 3, 8, 7, 10)
)

fit_list[[as.character(1e-3)]]$LL <- fit_list[[as.character(1e-3)]]$LL[,paste0("k", 1:10)]

colnames(fit_list[[as.character(1e-3)]]$FF) <- paste0(
  "k", c(1, 2, 4, 6, 9, 5, 3, 8, 7, 10)
)

fit_list[[as.character(1e-3)]]$FF <- fit_list[[as.character(1e-3)]]$FF[,paste0("k", 1:10)]

topic_order <- rev(paste0(
  "k",
  c(9, 1, 2, 8, 5, 10, 3, 4, 7, 6)
))

sp <- normalized_structure_plot(
  fit_list[[as.character(1e-3)]],
  grouping = dat$labels,gap = 25,perplexity = 70,n = Inf, font.size = 14,
  topics = topic_order
)

dat$doc_id <- 1:nrow(dat)
dat$k3 <- fit_list[[as.character(1e-3)]]$LL[,"k3"]
dat <- dat %>% dplyr::select(text, labels, k3, doc_id)
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

