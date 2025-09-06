library(ggplot2)
library(fastTopics)
library(log1pNMF)
library(ggpubr)

high_expressed_genes <- 500
low_expressed_genes <- 500

p <- high_expressed_genes + low_expressed_genes
n <- 1000

LL <- matrix(
  data = c(
    rep(1, 333), rep(0, 1000 - 333),
    rep(0, 333), rep(1, 333), rep(0, 334),
    rep(0, 1000 - 334), rep(1, 334)
  ),
  nrow = n,
  ncol = 3
)


FF <- matrix(
  data = c(
    rep(3, round(low_expressed_genes / 2)), 
    rep(0, round(low_expressed_genes / 2)),
    rep(1100, round(high_expressed_genes / 2)),
    rep(1000, round(high_expressed_genes / 2)),
    rep(0, round(low_expressed_genes / 2)), 
    rep(3, round(low_expressed_genes / 2)),
    rep(1000, round(high_expressed_genes / 2)),
    rep(1100, round(high_expressed_genes / 2)),
    rep(0, round(low_expressed_genes / 2)), 
    rep(0, round(low_expressed_genes / 2)),
    rep(1051.5, round(high_expressed_genes / 2)),
    rep(1051.5, round(high_expressed_genes / 2))
  ),
  nrow = p,
  ncol = 3
)

Lambda <- tcrossprod(LL, FF)
set.seed(1)
Y <- matrix(
  data = rpois(n = n * p, lambda = as.vector(Lambda)),
  nrow = n,
  ncol = p
)

library(fastTopics)
rownames(Y) <- paste0("cell", 1:n)
colnames(Y) <- paste0("gene", 1:p)

expr_df <- data.frame(
  group = c(
    rep("Group A", 1000), rep("Group B", 1000), rep("Group C", 1000)
  ),
  gene_id = c(
    1:1000, 1:1000, 1:1000
  ),
  expr = c(
    rep(3, round(low_expressed_genes / 2)), 
    rep(0, round(low_expressed_genes / 2)),
    rep(1100, round(high_expressed_genes / 2)),
    rep(1000, round(high_expressed_genes / 2)),
    rep(0, round(low_expressed_genes / 2)), 
    rep(3, round(low_expressed_genes / 2)),
    rep(1000, round(high_expressed_genes / 2)),
    rep(1100, round(high_expressed_genes / 2)),
    rep(0, round(low_expressed_genes / 2)), 
    rep(0, round(low_expressed_genes / 2)),
    rep(1051.5, round(high_expressed_genes / 2)),
    rep(1051.5, round(high_expressed_genes / 2))
  )
)

g_expr <- ggplot(expr_df, aes(x = gene_id, y = expr)) +
  geom_col(width = 1, colour = NA, linewidth = 0) +  # Use bars to represent lambda values
  facet_wrap(~ group) +          # Create a panel for each group
  labs(x = "Gene Index", y = "True Expression") +
  scale_y_continuous(trans = "log1p", breaks = c(0, 3, 100, 1000)) + 
  cowplot::theme_cowplot() + # Use a minimal theme for a clean look
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0.5, size = 13)
    )
  

ft_r1 <- fastTopics:::fit_pnmf_rank1(Y)
init_LL <- cbind(
  ft_r1$L,
  matrix(
    data = 1e-3,
    nrow = n,
    ncol = 2
  )
)
rownames(init_LL) <- rownames(Y)

init_FF <- cbind(
  ft_r1$F,
  matrix(
    data = 1e-3,
    nrow = p,
    ncol = 2
  )
)
rownames(init_FF) <- colnames(Y)

ft_init <- init_poisson_nmf(X = Y, F = init_FF, L = init_LL)

ft_mod <- fit_poisson_nmf(X = Y, fit0 = ft_init, verbose = "none")

group <- c(rep("A", 333), rep("B", 333), rep("C", 334))


g_tm_sp <- structure_plot(log1pNMF:::normalize_bars(ft_mod$L), grouping = group, gap = 10, loadings_order = 1:n, topics = rev(1:3)) + 
  theme(
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 13)
    ) + 
  ylab("Membership") +
  ggtitle("Topic Model Loadings") +
  guides(fill=guide_legend(title="Factor"))

set.seed(1)
log1p_mod <- fit_poisson_log1p_nmf(
  Y = Y, K = 3, loglik = "exact",
  control = list(maxiter = 250, verbose = FALSE)
)

colnames(log1p_mod$LL) <- paste0("k", 1:3)
colnames(log1p_mod$FF) <- paste0("k", 1:3)

LL <- log1p_mod$LL
FF <- log1p_mod$FF

log1p_mod$LL[,"k2"] <- LL[,"k3"]
log1p_mod$LL[,"k3"] <- LL[,"k2"]

log1p_mod$FF[,"k2"] <- FF[,"k3"]
log1p_mod$FF[,"k3"] <- FF[,"k2"]

g_log1p_sp <- normalized_structure_plot(
  log1p_mod, 
  grouping = group, 
  gap = 10,
  loadings_order = 1:n,
  topics = rev(1:3)
) + 
  theme(
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 13)
    ) + 
  ylab("Membership") +
  ggtitle("log1p Model Loadings (c = 1)") +
  guides(fill=guide_legend(title="Factor"))

g_sp <- ggarrange(
  g_log1p_sp, g_tm_sp,
  nrow = 1, ncol = 2,
  common.legend = TRUE,
  legend = "right",
  labels = c("B", "C")
)

max_col <- apply(log1p_mod$FF, 2, max)
log1p_FF <- sweep(log1p_mod$FF, 2, max_col, FUN = "*")

max_col <- apply(ft_mod$F, 2, max)
tm_FF <- sweep(ft_mod$F, 2, max_col, FUN = "*")
colnames(tm_FF) <- paste0("k", 1:3)

k2_df <- data.frame(
  model = c(
    rep("Topic Model Factor 2", 1000),
    rep("log1p Model Factor 2", 1000)
  ),
  factor_value = c(
    tm_FF[,"k2"],
    log1p_FF[,"k2"]
  ),
  gene_type = c(
    rep("Low On", 250),
    rep("Low Off", 250),
    rep("High Up", 250),
    rep("High Down", 250)
  )
)

k3_df <- data.frame(
  model = c(
    rep("Topic Model Factor 3", 1000),
    rep("log1p Model Factor 3", 1000)
  ),
  factor_value = c(
    tm_FF[,"k3"],
    log1p_FF[,"k3"]
  ),
  gene_type = c(
    rep("Low Off", 250),
    rep("Low On", 250),
    rep("High Down", 250),
    rep("High Up", 250)
  )
)

k_df <- rbind(k2_df, k3_df)

k_df$model <- factor(
  k_df$model, 
  levels = c(
    "log1p Model Factor 2", 
    "Topic Model Factor 2",
    "log1p Model Factor 3",
    "Topic Model Factor 3"
  )
)

g_k_a <- ggplot(
    k_df %>% dplyr::filter(model == "log1p Model Factor 2"), 
    aes(x = factor_value, fill = gene_type)
  ) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 60) +
  cowplot::theme_cowplot() +
  scale_x_continuous(trans = "log1p") +
  labs(x = "Factor Value",
       y = "Count",
       fill = "Gene Type",
       title = "log1p Model Factor 2"
       ) +
  theme(plot.title = element_text(hjust = 0.5, size = 13))

g_k_b <- ggplot(
  k_df %>% dplyr::filter(model == "Topic Model Factor 2"), 
  aes(x = factor_value, fill = gene_type)
) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 60) +
  cowplot::theme_cowplot() +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 10, 100)
  ) +
  labs(x = "Factor Value",
       y = "Count",
       fill = "Gene Type",
       title = "Topic Model Factor 2"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 13))

g_k_c <- ggplot(
  k_df %>% dplyr::filter(model == "log1p Model Factor 3"), 
  aes(x = factor_value, fill = gene_type)
) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 60) +
  cowplot::theme_cowplot() +
  scale_x_continuous(trans = "log1p") +
  labs(x = "Factor Value",
       y = "Count",
       fill = "Gene Type",
       title = "log1p Model Factor 3"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 13))

g_k_d <- ggplot(
  k_df %>% dplyr::filter(model == "Topic Model Factor 3"), 
  aes(x = factor_value, fill = gene_type)
) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 60) +
  cowplot::theme_cowplot() +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 10, 100)
  ) +
  labs(x = "Factor Value",
       y = "Count",
       fill = "Gene Type",
       title = "Topic Model Factor 3"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 13))


g_k <- ggarrange(
  g_k_a, g_k_b,
  g_k_c, g_k_d,
  nrow = 2, ncol = 2,
  common.legend = TRUE,
  legend = "right",
  labels = c("D", "E", "F", "G")
)

g <- ggarrange(
  g_expr,
  g_sp,
  g_k,
  nrow = 3,
  ncol = 1,
  labels = c("A", "", ""),
  heights = c(1, 1, 2)
)

ggsave(
  plot = g,
  device = "png",
  filename = "../images/factor_sim_fig.png",
  width = 7,
  height = 10
)
