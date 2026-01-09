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
dat <- dat %>%
  dplyr::mutate(
    labels = case_when(
      labels == "business" ~ "Business", 
      labels == "entertainment" ~ "Entertainment",
      labels == "politics" ~ "Politics",
      labels == "sport" ~ "Sports",
      labels == "tech" ~ "Tech"
    )
  )

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

plot_list <- list()

for (cc in cc_vec) {
  
  set.seed(1)
  plot_list[[as.character(cc)]] <- normalized_structure_plot(
    fit_list[[as.character(cc)]], 
    grouping = dat$labels,gap = 25,perplexity = 70,n = Inf, font.size = 14
  )  + ggtitle(glue::glue("log1p Model Loadings (c = {cc})")) +
    theme(
      plot.title = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(angle = 0,hjust = 0.5, size = 11)
    ) + ylab("Membership") +
    guides(fill=guide_legend(title="Factor", ncol = 1)) +
    guides(colour = "none")
  
}

L <- log1pNMF:::normalize_bars(diag(1 / fit_list[[as.character(1)]]$s) %*% fit_list[[as.character(Inf)]]$L)
set.seed(1)
plot_list[[as.character(Inf)]] <- structure_plot(
  L, 
  grouping = dat$labels,gap = 25,perplexity = 70,n = Inf, font.size = 12
)  + ggtitle("log1p Model Loadings (c = \u221E)") +
  theme(
    plot.title = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 11)
  ) + ylab("Membership") +
  guides(fill=guide_legend(title="Factor", ncol = 1)) +
  guides(colour = "none")
  
g <- ggarrange(
  plotlist = plot_list, 
  nrow = 8, 
  ncol = 1, 
  common.legend = TRUE, 
  legend = "right"
  )

ggsave(
  "../images/supp_bbc_structure.png",
  g,
  device = "png",
  width = 10,
  height = 14
)

  