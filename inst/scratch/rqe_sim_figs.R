library(ggplot2)
library(fastTopics)
library(log1pNMF)
library(ggpubr)
library(tidyr)
library(cowplot)
library(dplyr)

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
    rep(1000, round(high_expressed_genes / 2)),
    rep(1000, round(high_expressed_genes / 2))
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
    rep("Group C", 1000), rep("Group B", 1000), rep("Group A", 1000)
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
    rep(1000, round(high_expressed_genes / 2)),
    rep(1000, round(high_expressed_genes / 2))
  )
)

plot_df <- expr_df %>%
  mutate(
    # baseline heights
    baseline = case_when(
      gene_id <= 500 ~ pmin(expr, 3),
      gene_id > 500  ~ pmin(expr, 1000)
    ),
    excess = pmax(expr - baseline, 0),

    # log1p-transformed cumulative heights
    baseline_log = log1p(baseline),
    excess_log   = log1p(baseline + excess) - log1p(baseline)
  ) %>%
  select(group, gene_id, baseline_log, excess_log) %>%
  pivot_longer(
    c(baseline_log, excess_log),
    names_to = "segment",
    values_to = "height_log"
  ) %>%
  filter(height_log > 0) %>%
  mutate(
    fill = case_when(
      segment == "baseline_log" & gene_id <= 250 ~ "low_green",
      segment == "baseline_log" & gene_id <= 500 ~ "low_blue",
      segment == "baseline_log" & gene_id >  500 ~ "high_red",

      segment == "excess_log"   & gene_id <= 750 ~ "excess_green",
      segment == "excess_log"   & gene_id >  750 ~ "excess_blue"
    ),
    fill = factor(
      fill,
      levels = c("low_green", "low_blue", "excess_green",
                 "excess_blue", "high_red")
    )
  )

g_expr <- ggplot(plot_df, aes(x = gene_id, y = height_log, fill = fill)) +
  geom_col(
    width = 1,
    colour = NA,
    linewidth = 0,
    position = position_stack(reverse = FALSE)
  ) +
  facet_wrap(~ group) +
  labs(x = "Feature Index", y = "True Rate") +
  scale_y_continuous(
    breaks = log1p(c(0, 3, 100, 1000)),
    labels = c(0, 3, 100, 1000)
  ) +
  scale_fill_manual(
    values = c(
      low_green    = "#740001",
      low_blue     = "#005F73",
      high_red     = "#D3A625",
      excess_green = "#740001",
      excess_blue  = "#005F73"
    ),
    breaks = c("high_red", "excess_blue", "excess_green"),
    labels = c(
      expression(f^{(A)}),
      expression(f^{(Delta*B)}),
      expression(f^{(Delta*C)})
    ),
    guide = guide_legend(
      title = NULL,
      override.aes = list(colour = NA)
    )
  ) +
  theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13)
  )

ggsave(
  plot = g_expr,
  device = "png",
  filename = "~/Documents/log1pNMF/inst/paper_figures/images/factor_sim_fig_expr.png",
  width = 8,
  height = 6
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

group <- c(rep("C", 333), rep("B", 333), rep("A", 334))


set.seed(1)
log1p_mod <- fit_poisson_log1p_nmf(
  Y = Y, K = 3, loglik = "exact",
  control = list(maxiter = 250, verbose = FALSE)
)

g_tm_sp <- structure_plot(
  log1pNMF:::normalize_bars( diag(1 / log1p_mod$s) %*% ft_mod$L),
  grouping = group,
  gap = 10,
  loadings_order = 1:n,
  topics = rev(1:3),
  colors = c("#D3A625", "#740001", "#005F73")
) +
  theme(
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 12)
  ) +
  ylab("Membership") +
  ggtitle("Standard NMF Loadings") +
  guides(fill=guide_legend(title="Factor"))

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
  topics = rev(1:3),
  colors = c("#D3A625", "#740001", "#005F73")
) +
  theme(
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 12)
  ) +
  ylab("Membership") +
  ggtitle("log1p Model Loadings (c = 1)") +
  guides(fill=guide_legend(title="Factor"))


# --- 1. Keep one copy of g_expr with the desired legend ----------------------

g_expr_with_legend <- g_expr +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 14)
  )

expr_legend <- cowplot::get_legend(g_expr_with_legend)

# --- 2. Suppress legends inside the actual panels used in the layout ---------

g_expr_panel <- g_expr +
  theme(legend.position = "none")

g_tm_sp_panel <- g_tm_sp +
  theme(legend.position = "none")

g_log1p_sp_panel <- g_log1p_sp +
  theme(legend.position = "none")

# --- 3. Bottom row: two structure plots side by side -------------------------

g_sp_row <- ggarrange(
  g_log1p_sp_panel, g_tm_sp_panel,
  nrow = 1, ncol = 2,
  labels = c("", ""),
  align = "hv"
)

# --- 4. Stack top and bottom -------------------------------------------------

main_panels <- ggarrange(
  g_expr_panel,
  g_sp_row,
  ncol = 1, nrow = 2,
  heights = c(1.15, 1),   # adjust to taste
  align = "v"
)

# --- 5. Add the single legend on the right ----------------------------------

final_plot <- cowplot::plot_grid(
  main_panels,
  NULL,
  expr_legend,
  nrow = 1,
  rel_widths = c(1, 0.04, 0.12)
)

ggsave(
  plot = final_plot,
  device = "png",
  filename = "~/Documents/log1pNMF/inst/paper_figures/images/factor_sim_fig1_rqe.png",
  width = 8,
  height = 6
)

# factor plot ---------------------------------------------------------------

library(grid)

max_col <- apply(log1p_mod$FF, 2, max)
log1p_FF <- sweep(log1p_mod$FF, 2, max_col, FUN = "*")

max_col <- apply(ft_mod$F, 2, max)
tm_FF <- sweep(ft_mod$F, 2, max_col, FUN = "*")
colnames(tm_FF) <- paste0("k", 1:3)

k2_df <- data.frame(
  model = c(
    rep("Topic Model f^(ΔB) Factor", 1000),
    rep("log1p Model f^(ΔB) Factor", 1000)
  ),
  factor_value = c(
    tm_FF[, "k3"],
    log1p_FF[, "k3"]
  ),
  gene_type = c(
    rep("Off (1-250)", 250),
    rep("Low On (251-500)", 250),
    rep("Unchanged (501-750)", 250),
    rep("High Up (751-1000)", 250)
  )
)

g_k_a_with_legend <- ggplot(
  k2_df %>% dplyr::filter(model == "log1p Model f^(ΔB) Factor"),
  aes(x = factor_value, fill = gene_type)
) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 60) +
  cowplot::theme_cowplot() +
  scale_x_continuous(trans = "log1p") +
  labs(
    x = "Factor Value",
    y = "Count",
    fill = "Feature Type",
    title = expression(paste("log1p Model ", f^{(Delta*B)}, " Factor"))
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13),
    legend.position = "right",
    legend.title = element_text(size = 11.5),
    legend.text = element_text(size = 11.5),
    legend.key.size = unit(0.45, "cm")
  )

g_k_b_with_legend <- ggplot(
  k2_df %>% dplyr::filter(model == "Topic Model f^(ΔB) Factor"),
  aes(x = factor_value, fill = gene_type)
) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 60) +
  cowplot::theme_cowplot() +
  scale_x_continuous(
    trans = "log1p",
    breaks = c(0, 10, 100)
  ) +
  labs(
    x = "Factor Value",
    y = "Count",
    fill = "Feature Type",
    title = expression(paste("Standard NMF ", f^{(Delta*B)}, " Factor"))
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13),
    legend.position = "right",
    legend.title = element_text(size = 11.5),
    legend.text = element_text(size = 11.5),
    legend.key.size = unit(0.45, "cm")
  )

# --- 1. Extract legend for top panel ---------------------------------------

g_expr_with_legend <- g_expr +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 14)
  )

expr_legend <- cowplot::get_legend(g_expr_with_legend)

g_expr_panel <- g_expr +
  theme(legend.position = "none")

top_row <- cowplot::plot_grid(
  g_expr_panel,
  expr_legend,
  nrow = 1,
  rel_widths = c(1, 0.18),
  align = "h"
)

# --- 2. Extract shared legend for bottom row -------------------------------

factor_legend <- cowplot::get_legend(g_k_a_with_legend)

g_k_a <- g_k_a_with_legend +
  theme(legend.position = "none")

g_k_b <- g_k_b_with_legend +
  theme(legend.position = "none")

bottom_panels <- ggpubr::ggarrange(
  g_k_a, g_k_b,
  nrow = 1, ncol = 2,
  labels = c("", ""),
  align = "hv"
)

bottom_row <- cowplot::plot_grid(
  bottom_panels,
  factor_legend,
  nrow = 1,
  rel_widths = c(1, 0.30),
  align = "h"
)

# --- 3. Stack top and bottom -----------------------------------------------

final_plot_2 <- ggpubr::ggarrange(
  top_row,
  bottom_row,
  ncol = 1, nrow = 2,
  heights = c(1.15, 1),
  align = "v"
)

ggsave(
  plot = final_plot_2,
  device = "png",
  filename = "~/Documents/log1pNMF/inst/paper_figures/images/factor_sim_fig2_rqe.png",
  width = 8,
  height = 6
)