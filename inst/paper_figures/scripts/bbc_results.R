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

readr::write_rds(g, "../data/bbc_sparsity_ggplot.rds")

g <- ggarrange(g1,g3,g2, nrow = 1, labels = c("A", "B", "C"))

ggplot2::ggsave(
  "../pdfs/bbc_sparsity.pdf",
  g,
  device = "pdf",
  width = 11.5,
  height = 4
)

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
plot_list[[glue::glue("c = 0.001")]] <- sp$plot + 
  ggtitle(glue::glue("Loadings of log1p Link Poisson NMF With c = 0.001")) +
  theme(
    plot.title = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )

tm_sp <- structure_plot(
  fit_list[["Inf"]],
  grouping = dat$labels,gap = 25,perplexity = 70,font.size = 14,
  topics = topic_order, loadings_order = sp$loadings_order
)

plot_list[["Topic Model"]] <- tm_sp$plot + 
  ggtitle("Loadings of Identity Link Poisson NMF / Topic Model") +
  theme(
    plot.title = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )


g <- ggarrange(
  plotlist = plot_list,
  ncol = 1,
  labels = "AUTO"
)

ggsave(
  "../pdfs/bbc_structure.png",
  g,
  device = "png",
  width = 11,
  height = 6
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

library(kableExtra)

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


