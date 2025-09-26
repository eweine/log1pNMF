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



set.seed(1)
sp <- normalized_structure_plot(
  fit_list[[as.character(1e-3)]],
  grouping = dat$labels,gap = 25,perplexity = 70,n = Inf, font.size = 12
)
plot_list[[glue::glue("c = 0.001")]] <- sp + 
  ggtitle("log1p Model Loadings (c = 0.001)") +
  theme(
    plot.title = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 12.5),
    legend.key.size = unit(0.6, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 12.5),
    legend.title = element_text(size = 12.5)
  ) + ylab("Membership") +
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none")

set.seed(1)
tm_sp <- structure_plot(
  log1pNMF:::normalize_bars(diag(1 / s) %*% fit_list[["Inf"]]$L),
  grouping = dat$labels,gap = 25,perplexity = 70,font.size = 12
)

plot_list[["Topic Model"]] <- tm_sp + 
  ggtitle("Topic Model Loadings") +
  theme(
    plot.title = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 12.5),
    legend.key.size = unit(0.6, "cm"),   # make color boxes bigger
    legend.text = element_text(size = 12.5),
    legend.title = element_text(size = 12.5)
  )  + ylab("Membership") +
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none")


g <- ggarrange(
  plotlist = plot_list,
  ncol = 1,
  labels = "AUTO",
  common.legend = TRUE,
  legend = "right"
)

ggsave(
  "../images/bbc_structure.png",
  g,
  device = "png",
  width = 11,
  height = 5
)

# now, I want to look for top and distinctive words for each factor

get_top_words <- function(f, n_top = 20) {
  
  words <- names(sort(f, decreasing = TRUE))[1:n_top]
  pasted_words <- paste(words, collapse = ", ")
  return(pasted_words)
  
}

top_words <- apply(fit_list$`0.001`$FF, 2, get_top_words)
top_words <- as.data.frame(top_words)
colnames(top_words) <- c("top_words")
top_words$factor <- rownames(top_words)
rownames(top_words) <- NULL
top_words <- top_words %>% dplyr::select(c("factor", "top_words"))

colnames(fit_list$`Inf`$F) <- paste0("k", 1:10)
top_words_tm <- apply(fit_list$`Inf`$F, 2, get_top_words)
top_words_tm <- as.data.frame(top_words_tm)
colnames(top_words_tm) <- c("top_words")
top_words_tm$factor <- rownames(top_words_tm)
rownames(top_words_tm) <- NULL
top_words_tm <- top_words_tm %>% dplyr::select(c("factor", "top_words"))

colnames(top_words_tm) <- c("Factor", "Top Words - Topic Model")
colnames(top_words) <- c("Factor", "Top Words - log1p Model c = 0.001")

top_words <- top_words_tm %>%
  dplyr::inner_join(top_words)

library(readr)
library(dplyr)
library(stringr)
library(knitr)
library(kableExtra)

top_words <- top_words %>%
  mutate(
    `Top Words - Topic Model` = str_wrap(`Top Words - Topic Model`, width = 40),
    `Top Words - log1p Model c = 0.001` = str_wrap(`Top Words - log1p Model c = 0.001`, width = 40)
  )

df_checker <- top_words %>%
  mutate(
    Factor = Factor,
    # For odd rows, shade only TopWords1
    `Top Words - Topic Model` = if_else(row_number() %% 2 == 1,
                        cell_spec(`Top Words - Topic Model`, "latex", background = "gray!10"),
                        `Top Words - Topic Model`
    ),
    # For even rows, shade only TopWords2
    `Top Words - log1p Model c = 0.001` = if_else(row_number() %% 2 == 0,
                        cell_spec(`Top Words - log1p Model c = 0.001`, "latex", background = "gray!10"),
                        `Top Words - log1p Model c = 0.001`
    )
  )

# 3. Make the LaTeX table
# - 'booktabs = TRUE' for better looking horizontal rules
# - 'escape = FALSE' avoids escaping characters like underscores
# - column_spec can set the width for each column to force wrapping
kable(df_checker, format = "latex", booktabs = TRUE, escape = FALSE) %>%
  kable_styling(
    position = "center"
  ) %>%
  column_spec(2, width = "5cm") %>%
  column_spec(3, width = "5cm")


