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
    rep(0, round(1000 - 334)), rep(1, 334)
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

plot_df_base <- expr_df %>%
  mutate(
    baseline = case_when(
      gene_id <= 500 ~ pmin(expr, 3),
      gene_id > 500  ~ pmin(expr, 1000)
    ),
    excess = pmax(expr - baseline, 0),
    baseline_log = log1p(baseline),
    excess_log   = log1p(baseline + excess) - log1p(baseline)
  ) %>%
  mutate(
    fill = case_when(
      gene_id <= 250 ~ "low_green",
      gene_id <= 500 ~ "low_blue",
      gene_id >  500 ~ "high_red"
    ),
    fill_excess = case_when(
      gene_id <= 750 ~ "excess_green",
      gene_id >  750 ~ "excess_blue"
    )
  )

plot_df_raw <- bind_rows(
  plot_df_base %>%
    transmute(group, gene_id, segment = "baseline", height = baseline, fill = fill),
  plot_df_base %>%
    transmute(group, gene_id, segment = "excess", height = excess, fill = fill_excess)
) %>%
  filter(height > 0) %>%
  mutate(
    fill = factor(
      fill,
      levels = c("low_green", "low_blue", "excess_green",
                 "excess_blue", "high_red")
    )
  )

plot_df_log <- bind_rows(
  plot_df_base %>%
    transmute(group, gene_id, segment = "baseline", height = baseline_log, fill = fill),
  plot_df_base %>%
    transmute(group, gene_id, segment = "excess", height = excess_log, fill = fill_excess)
) %>%
  filter(height > 0) %>%
  mutate(
    fill = factor(
      fill,
      levels = c("low_green", "low_blue", "excess_green",
                 "excess_blue", "high_red")
    )
  )

common_scale_fill <- scale_fill_manual(
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
)

g_expr_raw <- ggplot(plot_df_raw, aes(x = gene_id, y = height, fill = fill)) +
  geom_col(
    width = 1,
    colour = NA,
    linewidth = 0,
    position = position_stack(reverse = FALSE)
  ) +
  facet_wrap(~ group) +
  labs(x = "Feature Index", y = "True Rate", title = "Raw scale") +
  common_scale_fill +
  theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(hjust = 0.5)
  )

g_expr_log <- ggplot(plot_df_log, aes(x = gene_id, y = height, fill = fill)) +
  geom_col(
    width = 1,
    colour = NA,
    linewidth = 0,
    position = position_stack(reverse = FALSE)
  ) +
  facet_wrap(~ group) +
  labs(x = "Feature Index", y = "True Rate", title = "log1p scale") +
  scale_y_continuous(
    breaks = log1p(c(0, 3, 100, 1000)),
    labels = c(0, 3, 100, 1000)
  ) +
  common_scale_fill +
  theme_cowplot() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    plot.title = element_text(hjust = 0.5)
  )

library(grid)

legend_theme <- theme(
  legend.text = element_text(size = 18),
  legend.title = element_text(size = 12),   # can keep small if you want
  legend.key.size = unit(1.2, "cm"),
  legend.box.margin = margin(l = 20, unit = "pt")
)

g_expr_raw <- g_expr_raw + legend_theme
g_expr_log <- g_expr_log + legend_theme

g_expr_both <- ggpubr::ggarrange(
  g_expr_raw,
  g_expr_log,
  ncol = 1,
  nrow = 2,
  common.legend = TRUE,
  legend = "right"
)

ggsave(
  plot = g_expr_both,
  device = "png",
  filename = "~/Documents/log1pNMF/inst/paper_figures/images/factor_sim_fig_expr_both.png",
  width = 8,
  height = 7
)
