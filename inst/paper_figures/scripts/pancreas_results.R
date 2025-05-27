library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggpubr)
library(passPCA)

load("../data/experiment_results.Rdata")

set.seed(1)

K <- 9

jaccard_index <- function(vec1, vec2) {
  
  set1 <- unique(vec1)
  
  set2 <- unique(vec2)
  
  intersection <- length(intersect(set1, set2))
  
  union <- length(union(set1, set2))
  
  return(intersection / union)
  
}

set.seed(1)

load("../data/raw_data/pancreas.RData")

i           <- which(sample_info$tech == "celseq2")
sample_info <- sample_info[i,]
counts      <- counts[i,]
sample_info <- transform(sample_info,celltype = factor(celltype))

counts <- counts[, Matrix::colSums(counts) > 0]
genes_to_use <- which(Matrix::colSums(counts>0)>4)
counts <- counts[,genes_to_use]

s <- Matrix::rowSums(counts)
s <- s / mean(s)

n_top <- 20

hoyer_sparsity <- function(x) {
  
  n <- length(x)
  (1 / (sqrt(n) - 1)) * (sqrt(n) - (sum(x) / (sqrt(sum(x ^ 2)))))
  
}

fit_list <- list()
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

for (cc in cc_vec) {
  
  fit_list[[as.character(cc)]] <- res_list$pancreas[[as.character(cc)]]
  
  colnames(fit_list[[as.character(cc)]]$FF) <- paste0("k", 1:K)
  colnames(fit_list[[as.character(cc)]]$LL) <- paste0("k", 1:K)
  
  F_norm <- passPCA:::normalize_bars(fit_list[[as.character(cc)]]$FF)
  top_list <- list()
  distinct_list <- list()
  jaccard_vec <- c()
  
  for (k in 1:K) {
    
    top_list[[k]] <- names(sort(F_norm[,k], decreasing = TRUE))[1:n_top]
    distinct_list[[k]] <- fastTopics:::get_distinctive_features(
      effects_matrix = F_norm,
      k = glue::glue("k{k}"), n = n_top, feature_sign = "positive"
    )
    
    jaccard_vec <- c(jaccard_vec, jaccard_index(
      top_list[[k]],
      distinct_list[[k]]
    ))
    
  }
  
  fit_list[[as.character(cc)]]$top_genes <- top_list
  fit_list[[as.character(cc)]]$distinct_genes <- distinct_list
  fit_list[[as.character(cc)]]$jaccard <- jaccard_vec
  
  fit_list[[as.character(cc)]]$l_sparsity <- apply(
    fit_list[[as.character(cc)]]$LL, 2, hoyer_sparsity
  )
  
  fit_list[[as.character(cc)]]$f_sparsity <- apply(
    fit_list[[as.character(cc)]]$FF, 2, hoyer_sparsity
  )
  
  fit_list[[as.character(cc)]]$cor_mat <- cor(fit_list[[as.character(cc)]]$FF, method = "spearman")
  
}

fit_list[["Inf"]] <- res_list$pancreas$`Inf`

fit_list[["Inf"]]$Ls <- Matrix::Diagonal(x = 1/s) %*% fit_list[["Inf"]]$L

F_norm <- passPCA:::normalize_bars(fit_list[["Inf"]]$F)
colnames(F_norm) <- paste0("k",1:K)
top_list <- list()
distinct_list <- list()
jaccard_vec <- c()

for (k in 1:K) {
  
  top_list[[k]] <- names(sort(F_norm[,k], decreasing = TRUE))[1:n_top]
  distinct_list[[k]] <- fastTopics:::get_distinctive_features(
    effects_matrix = F_norm,
    k = glue::glue("k{k}"), n = n_top, feature_sign = "positive"
  )
  
  jaccard_vec <- c(jaccard_vec, jaccard_index(
    top_list[[k]],
    distinct_list[[k]]
  ))
  
}

fit_list[["Inf"]]$top_genes <- top_list
fit_list[["Inf"]]$distinct_genes <- distinct_list
fit_list[["Inf"]]$jaccard <- jaccard_vec

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

jaccard_vec <- unlist(
  lapply(fit_list, function(x) {median(x$jaccard)})
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
    geom="text", x=0.004, y=cor_vec["Inf"] + 0.05, label="ID Link", color="red",
    size = 5
  )  + theme(
    axis.title = element_text(size = 14, face = "bold")
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
    geom="text", x=0.004, y=l_sparsity_vec["Inf"] + 0.05, label="ID Link", color="red",
    size = 5
  )  + theme(
    axis.title = element_text(size = 14, face = "bold")
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
    geom="text", x=0.004, y=f_sparsity_vec["Inf"] + 0.05, label="ID Link", color="red",
    size = 5
  ) + theme(
    axis.title = element_text(size = 14, face = "bold")
  )


g <- ggarrange(g1,g3,g2, nrow = 1, labels = c("D", "E", "F"))

g <- annotate_figure(g,
                     top = text_grob(glue::glue("Pancreas K = {K}"), size = 20, face = "bold"))

readr::write_rds(g, "../data/pancreas_sparsity_ggplot.rds")

g <- ggarrange(g1,g3,g2, nrow = 1, labels = c("A", "B", "C"))

ggplot2::ggsave(
  "../pdfs/pancreas_sparsity.pdf",
  g,
  device = "pdf",
  width = 11.5,
  height = 4
)

celltype <- sample_info$celltype
celltype <-
  factor(celltype,
         c("acinar","ductal","activated_stellate","quiescent_stellate",
           "endothelial","macrophage","mast","schwann","alpha","beta",
           "delta","gamma","epsilon"))

levels(celltype)[levels(celltype) == "activated_stellate"] <- "activated PSC"
levels(celltype)[levels(celltype) == "quiescent_stellate"] <- "quiescent PSC"



gamma_info <- sample_info %>%
  dplyr::filter(celltype == "gamma")

delta_info <- sample_info %>%
  dplyr::filter(celltype == "delta")

gamma_counts <- counts[gamma_info$id, ]
delta_counts <- counts[delta_info$id, ]

gamma_means <- Matrix::colMeans(gamma_counts)
delta_means <- Matrix::colMeans(delta_counts)

means_df <- data.frame(
  gene = names(gamma_means),
  gamma_expr = unname(gamma_means),
  delta_expr = unname(delta_means)
)

library(ggplot2)

genes_to_label <- c("RBP4", "SST", "PPY")
library(ggrepel)
# this shows the differences in expression between gamma and delta
g_log1p_scatter <- ggplot(data = means_df, aes(x = log1p(delta_expr), y = log1p(gamma_expr))) +
  geom_point() +
  xlab("log(1 + Mean Expression in Delta Cells)") +
  ylab("log(1 + Mean Expression in Gamma Cells)") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 5.5,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.25,
    nudge_y = 0.5
  ) +
  cowplot::theme_cowplot() +
  ggtitle("Delta Vs. Gamma Expression - log1p Transformed") +
  theme(
    plot.title = element_text(size = 10),      # Title size
    axis.title.x = element_text(size = 8),     # X-axis label size
    axis.title.y = element_text(size = 8)      # Y-axis label size
  )

g_raw_scatter <- ggplot(data = means_df, aes(x = delta_expr, y = gamma_expr)) +
  geom_point() +
  xlab("Mean Expression in Delta Cells") +
  ylab("Mean Expression in Gamma Cells") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 5.5,                # adjust text size as desired
    box.padding = 1,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 5
  ) +
  cowplot::theme_cowplot() +
  ggtitle("Delta Vs. Gamma Expression - Raw") +
  theme(
    plot.title = element_text(size = 10),      # Title size
    axis.title.x = element_text(size = 8),     # X-axis label size
    axis.title.y = element_text(size = 8)      # Y-axis label size
  )

# overall, the two cell types are quite similar

# for the topic model, I should make an annotation heatmap of 2, 3, and 4
# for the log1p fit, I should make an annotation heatmap for 1, 4, 7, and 8, 
# and maybe 5, but I'm not sure it really makes sense.
colnames(fit_list[["Inf"]]$F) <- paste0("k", 1:K)
rownames(fit_list[["Inf"]]$F) <- colnames(counts)

colnames(fit_list[[as.character(1)]]$FF) <- paste0(
  "k", c(1, 3, 6, 5, 2, 7, 8, 9, 4)
)

colnames(fit_list[[as.character(1)]]$LL) <- paste0(
  "k", c(1, 3, 6, 5, 2, 7, 8, 9, 4)
)

colnames(fit_list[[as.character(Inf)]]$F) <- paste0(
  "k", 1:K
)

colnames(fit_list[[as.character(Inf)]]$L) <- paste0(
  "k", 1:K
)

fit_list[[as.character(1)]]$LL <- fit_list[[as.character(1)]]$LL[,paste0("k", 1:K)]
fit_list[[as.character(1)]]$FF <- fit_list[[as.character(1)]]$FF[,paste0("k", 1:K)]
fit_list[[as.character(Inf)]]$L <- fit_list[[as.character(Inf)]]$L[,paste0("k", 1:K)]
fit_list[[as.character(Inf)]]$F <- fit_list[[as.character(Inf)]]$F[,paste0("k", 1:K)]

normalize_bars <- function(LL, FF) {
  
  max_col <- apply(LL, 2, max)
  sweep(FF, 2, max_col, FUN = "*")
  
}


g_tm_high <- annotation_heatmap(
  effects_matrix = normalize_bars(fit_list[["Inf"]]$L, fit_list[["Inf"]]$F)[genes_to_label,paste0("k", c(5, 8))],
  feature_sign = "positive",
  n = 5
) + ggtitle("Marker Gene Factor Importance - Topic Model") +  theme(
  plot.title = element_text(size = 10)
)

g_log1p_high <- annotation_heatmap(
  effects_matrix = normalize_bars(fit_list[["1"]]$LL, fit_list[["1"]]$FF)[genes_to_label, paste0("k", c(1, 2, 5, 8))],
  feature_sign = "positive",
  n = 5
) + ggtitle("Marker Gene Factor Importance - log1p Link with c = 1") + theme(
  plot.title = element_text(size = 10)
)
rownames(fit_list[["0.1"]]$LL) <- rownames(counts)

# si_islet <- sample_info %>%
#   dplyr::filter(celltype %in% c("alpha", "beta", "delta", "gamma"))
# 
# topic_order <- rev(paste0(
#   "k",
#   c(1, 2, 4, 9, 5, 7, 3, 8, 6)
# ))
# 
# fit_list[[as.character(1)]]$LL <- fit_list[[as.character(1)]]$LL[si_islet$id, ]
# 
# log1p_sp <- normalized_structure_plot(
#   fit_list[[as.character(1)]],
#   grouping = paste(si_islet$celltype, substr(si_islet$id, 1, 3)),gap = 35,perplexity = 70,n = Inf, font.size = 12,
#   topics = topic_order
# ) 

g_log1p_sp <- log1p_sp$plot + 
  ggtitle("Loadings of log1p Link Poisson NMF With c = 1") +
  theme(
    plot.title = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )

g_tm_sp <- structure_plot(
  fit_list[[as.character(Inf)]], 
  grouping = celltype, topics = topic_order,
  gap = 35,perplexity = 70, font.size = 12, loadings_order = log1p_sp$loadings_order
)$plot + ggtitle("Loadings of Identity Link Poisson NMF / Topic Model") +
  theme(
    plot.title = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  )

g1 <- ggarrange(
  g_log1p_sp, g_tm_sp,
  nrow = 2, ncol = 1,
  labels = c("A", "B")
)

ggsave(
  "../pdfs/talk_pancreas_structure.png",
  g1,
  device = "png",
  width = 11,
  height = 6
)

g2 <- ggarrange(
  g_log1p_scatter, g_raw_scatter,
  nrow = 1, ncol = 2,
  labels = c("C", "D")
)

g3 <- ggarrange(
  g_log1p_high, g_tm_high,
  nrow = 1, ncol = 2,
  labels = c("E", "F")
)

g <- ggarrange(
  g1, g2, g3, nrow = 3, ncol = 1,
  heights = c(2, 1, 1)
)

g <- annotate_figure(g,
                     top = text_grob(glue::glue("Pancreas K = {K}"), size = 20, face = "bold"))

means_df <- means_df %>%
  dplyr::mutate(
    log1p_k5 = normalize_bars(fit_list[["1"]]$LL, fit_list[["1"]]$FF)[, "k5"],
    log1p_k8 = normalize_bars(fit_list[["1"]]$LL, fit_list[["1"]]$FF)[, "k8"],
    tm_k5 = normalize_bars(fit_list[["Inf"]]$L, fit_list[["Inf"]]$F)[, "k5"],
    tm_k8 = normalize_bars(fit_list[["Inf"]]$L, fit_list[["Inf"]]$F)[, "k8"]
  )


plot(log1p(means_df$gamma_expr), means_df$log1p_k8)
plot(log1p(means_df$delta_expr), means_df$log1p_k5)

plot(log1p(means_df$gamma_expr), log1p(means_df$tm_k8))
plot(log1p(means_df$delta_expr), log1p(means_df$tm_k5))

plot(means_df$log1p_k5, means_df$log1p_k8)
plot(log1p(means_df$tm_k5), log1p(means_df$tm_k8))

g_log1p_scatter <- ggplot(data = means_df, aes(x = log1p(delta_expr), y = log1p(gamma_expr))) +
  geom_point(size = 0.75, alpha = 0.33) +
  xlab("log(1 + Mean Expression in Delta Cells)") +
  ylab("log(1 + Mean Expression in Gamma Cells)") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 7,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.25
  ) +
  cowplot::theme_cowplot() +
  ggtitle("Delta vs. Gamma Expression - log1p Transformed") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),      # Title size
    axis.title.x = element_text(size = 14, face = "bold"),     # X-axis label size
    axis.title.y = element_text(size = 14, face = "bold")      # Y-axis label size
  )

g_raw_scatter <- ggplot(data = means_df, aes(x = delta_expr, y = gamma_expr)) +
  geom_point(size = 0.75, alpha = 0.33) +
  xlab("Mean Expression in Delta Cells") +
  ylab("Mean Expression in Gamma Cells") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 7,                # adjust text size as desired
    box.padding = 1,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 5
  ) +
  cowplot::theme_cowplot() +
  ggtitle("Delta vs. Gamma Expression - Raw") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),      # Title size
    axis.title.x = element_text(size = 14, face = "bold"),     # X-axis label size
    axis.title.y = element_text(size = 14, face = "bold")      # Y-axis label size
  )

g_log1p_gamma <- ggplot(data = means_df, aes(x = log1p(gamma_expr), y = log1p_k8)) +
  geom_point(size = 0.75, alpha = 0.33) +
  xlab("log(1 + Mean Expression in Gamma Cells)") +
  ylab("Gene Scores - Factor 8 (Pink)") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 7,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.5,
    nudge_y = 1.5
  ) +
  cowplot::theme_cowplot() +
  ggtitle("Gamma - Mean Expression Vs. log1p Model Score") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),      # Title size
    axis.title.x = element_text(size = 14, face = "bold"),     # X-axis label size
    axis.title.y = element_text(size = 14, face = "bold")      # Y-axis label size
  )

g_log1p_delta <- ggplot(data = means_df, aes(x = log1p(delta_expr), y = log1p_k5)) +
  geom_point(size = 0.75, alpha = 0.33) +
  xlab("log(1 + Mean Expression in Delta Cells)") +
  ylab("Gene Scores - Factor 5 (Orange)") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 7,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.75,
    nudge_y = 0.5
  ) +
  cowplot::theme_cowplot() +
  ggtitle("Delta - Mean Expression Vs. log1p Model Score") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),      # Title size
    axis.title.x = element_text(size = 14, face = "bold"),     # X-axis label size
    axis.title.y = element_text(size = 14, face = "bold")      # Y-axis label size
  )

g_tm_gamma <- ggplot(data = means_df, aes(x = log1p(gamma_expr), y = log1p(tm_k8))) +
  geom_point(size = 0.75, alpha = 0.33) +
  xlab("log(1 + Mean Expression in Gamma Cells)") +
  ylab("log(1 + Gene Scores) - Factor 8 (Pink)") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 7,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.5,
    nudge_y = 0.5
  ) +
  cowplot::theme_cowplot() +
  ggtitle("Gamma - Mean Expression Vs. Topic Model Score") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),      # Title size
    axis.title.x = element_text(size = 14, face = "bold"),     # X-axis label size
    axis.title.y = element_text(size = 14, face = "bold")      # Y-axis label size
  )

g_tm_delta <- ggplot(data = means_df, aes(x = log1p(delta_expr), y = log1p(tm_k5))) +
  geom_point(size = 0.75, alpha = 0.33) +
  xlab("log(1 + Mean Expression in Delta Cells)") +
  ylab("log(1 + Gene Scores) - Factor 5 (Orange)") +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 7,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.5,
    nudge_y = 0.5
  ) +
  cowplot::theme_cowplot() +
  ggtitle("Delta - Mean Expression Vs. Topic Model Score") +
  theme(
    plot.title = element_text(size = 12, face = "bold"),      # Title size
    axis.title.x = element_text(size = 14, face = "bold"),     # X-axis label size
    axis.title.y = element_text(size = 14, face = "bold")      # Y-axis label size
  )

i_dg <- which(sample_info$celltype %in% c("delta", "gamma"))
lo <- log1p_sp$loadings_order[log1p_sp$loadings_order %in% i_dg]
sample_info_dg <- sample_info[i_dg,]

log1p_sp_dg <- normalized_structure_plot(
  fit_list[[as.character(1)]],
  grouping = celltype,gap = 35,perplexity = 70,n = Inf, font.size = 12,
  topics = topic_order, loadings_order = lo
)

g_structure_dg_log1p <- log1p_sp_dg$plot + 
  ggtitle("log1p Model") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

g_structure_dg_tm <- structure_plot(
  fit_list[[as.character(Inf)]], 
  grouping = celltype, topics = topic_order,
  gap = 35,perplexity = 70, font.size = 12, loadings_order = lo,
)$plot + ggtitle("Topic Model") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

g_sp_dg <- ggarrange(
  g_structure_dg_log1p,
  g_structure_dg_tm,
  nrow = 2, ncol = 1,
  labels = c("A", "B")
  )

g_scatter_talk <- ggarrange(
  NULL, g_log1p_scatter, g_raw_scatter,
  g_structure_dg_log1p, g_log1p_gamma, g_log1p_delta,
  g_structure_dg_tm, g_tm_gamma, g_tm_delta,
  nrow = 3, ncol = 3,
  labels = c("", "A", "B", "C", "D", "E", "F", "G", "H")
)

ggsave(
  "../pdfs/scatter_talk.png",
  g_scatter_talk,
  device = "png",
  width = 14,
  height = 14
)

g1 <- ggarrange(
  g_log1p_sp, g_tm_sp,
  nrow = 2, ncol = 1,
  labels = c("A", "B")
)

g2 <- ggarrange(
  g_log1p_delta, g_log1p_gamma,
  nrow = 1, ncol = 2,
  labels = c("C", "D")
)

g3 <- ggarrange(
  g_tm_delta, g_tm_gamma,
  nrow = 1, ncol = 2,
  labels = c("E", "F")
)

g <- ggarrange(
  g1, g2, g3, nrow = 3, ncol = 1,
  heights = c(2, 1, 1)
)

ggsave(
  "../pdfs/pancreas_structure.pdf",
  g,
  device = "pdf",
  width = 8.75,
  height = 11.5
)

# here, I want to make correlation plots

gcor_a <- ggplot(data = means_df, aes(x = log1p_k5, y = log1p_k8)) +
  geom_point(size = 0.75, alpha = 0.33) +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 5.5,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.5,
    nudge_y = 0.5
  ) +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 10),      # Title size
    axis.title.x = element_text(size = 8),     # X-axis label size
    axis.title.y = element_text(size = 8)      # Y-axis label size
  ) +
  xlab("Log1p model Factor 5 (Orange)") +
  ylab("Log1p model Factor 8 (Pink)")

gcor_b <- ggplot(data = means_df, aes(x = log1p(tm_k5), y = log1p(tm_k8))) +
  geom_point(size = 0.75, alpha = 0.33) +
  geom_text_repel(
    data = subset(means_df, gene %in% genes_to_label),
    aes(label = gene),
    size = 5.5,                # adjust text size as desired
    box.padding = 0.5,       # adjust padding around labels
    max.overlaps = Inf,       # ensure all labels appear
    segment.color = 'grey40',  # color of connecting lines
    segment.size = 0.4,    
    min.segment.length = 0,
    nudge_x = 0.5,
    nudge_y = 0.5
  ) +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(size = 10),      # Title size
    axis.title.x = element_text(size = 8),     # X-axis label size
    axis.title.y = element_text(size = 8)      # Y-axis label size
  ) +
  xlab("log(1 + Topic Model Factor 5 (Orange))") +
  ylab("log(1 + Topic Model Factor 8 (Pink))")

ggarrange(gcor_a, gcor_b, nrow = 1, ncol = 2)

library(reshape2)
df <- melt(fit_list$`1`$cor_mat)

df$Var1 <- factor(df$Var1, levels = rev(rownames(fit_list$`1`$cor_mat)))
df$Var2 <- factor(df$Var2, levels = colnames(fit_list$`1`$cor_mat))

hma <- ggplot(df, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  # White for 0 and red for 1 (assuming all correlations are within [0, 1])
  scale_fill_gradient(limits = c(0, 1), low = "white", high = "red") +
  theme_minimal() +
  # Optional improvements to remove extra grid lines, etc.
  theme(panel.grid = element_blank()) +
  labs(x = "Column", y = "Row", fill = "Corr") +
  ggtitle("Separman Correlation log1p model")

rownames(fit_list$`Inf`$cor_mat) <- paste0("k", 1:9)
colnames(fit_list$`Inf`$cor_mat) <- paste0("k", 1:9)
df <- melt(fit_list$`Inf`$cor_mat)
df$Var1 <- factor(df$Var1, levels = rev(rownames(fit_list$`1`$cor_mat)))
df$Var2 <- factor(df$Var2, levels = colnames(fit_list$`1`$cor_mat))


hmb <- ggplot(df, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  # White for 0 and red for 1 (assuming all correlations are within [0, 1])
  scale_fill_gradient(limits = c(0, 1), low = "white", high = "red") +
  theme_minimal() +
  # Optional improvements to remove extra grid lines, etc.
  theme(panel.grid = element_blank()) +
  labs(x = "Column", y = "Row", fill = "Corr") +
  ggtitle("Separman Correlation topic model")

ggarrange(hma, hmb, nrow = 1, ncol = 2)
