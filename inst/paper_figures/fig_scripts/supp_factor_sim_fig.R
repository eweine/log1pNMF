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
  labs(x = "Feature Index", y = "True Rate") +
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


log1p_fit_list <- list()
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)

for (cc in cc_vec) {
  
  set.seed(1)
  log1p_fit_list[[as.character(cc)]] <- fit_poisson_log1p_nmf(
    Y = Y, K = 3, loglik = "exact", init_method = "rank1",
    control = list(maxiter = 250, verbose = FALSE), cc = cc
  )
  
}

colnames(log1p_fit_list$`1`$LL) <- paste0("k", 1:3)
colnames(log1p_fit_list$`1`$FF) <- paste0("k", 1:3)

LL <- log1p_fit_list$`1`$LL
FF <- log1p_fit_list$`1`$FF

log1p_fit_list$`1`$LL[,"k2"] <- LL[,"k3"]
log1p_fit_list$`1`$LL[,"k3"] <- LL[,"k2"]

log1p_fit_list$`1`$FF[,"k2"] <- FF[,"k3"]
log1p_fit_list$`1`$FF[,"k3"] <- FF[,"k2"]

plot_list <- list()

for (cc in cc_vec) {
  
  set.seed(1)
  plot_list[[as.character(cc)]] <- normalized_structure_plot(
    log1p_fit_list[[as.character(cc)]], 
    loadings_order = 1:1000, 
    grouping = group, 
    gap = 10, 
    topics = rev(1:3)
  ) +
    theme(
      axis.text.x = element_text(angle = 0,hjust = 0.5, size = 11),
      axis.title.y = element_text(size = 13),
      plot.title = element_text(size = 14)
    ) + 
    ylab("Membership") + 
    ggtitle(glue::glue("log1p Model Loadings (c = {cc})")) +
    guides(fill=guide_legend(title="Factor"))
  
}

plot_list[[as.character(Inf)]] <- structure_plot(
  log1pNMF:::normalize_bars(
    diag(1 / log1p_fit_list[[as.character(1)]]$s) %*% ft_mod$L
  ), 
  loadings_order = 1:1000, 
  grouping = group, 
  gap = 10
) +
  theme(
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 11),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(size = 14)
  ) +
  ylab("Membership") + 
  ggtitle("Topic Model Loadings") +
  guides(fill="none")

g <- ggarrange(
  plotlist = plot_list, 
  nrow = 4, 
  ncol = 2, 
  common.legend = TRUE, 
  legend = "right"
)

ggsave(
  "../images/supp_factor_sim_structure.png",
  g,
  device = "png",
  width = 10,
  height = 11
)

