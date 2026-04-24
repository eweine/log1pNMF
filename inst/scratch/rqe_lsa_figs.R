library(Matrix)
library(log1pNMF)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)

load("~/Documents/log1pNMF/inst/paper_figures/data/raw_data/pancreas_cytokine_lsa.Rdata")

s <- Matrix::rowSums(counts)
s <- s / mean(s)
barcodes$s <- s

celltype_counts <- barcodes %>%
  count(celltype) %>%
  arrange(desc(n)) %>%
  mutate(celltype = factor(celltype, levels = celltype))

condition_counts <- barcodes %>%
  count(condition) %>%
  mutate(
    condition = if_else(condition == "IL-1B_IFNg", "IL-1B + IFNg", condition),
    condition = factor(condition, levels = c("Untreated", "IL-1B", "IFNg", "IL-1B + IFNg"))
  )

base_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

g_ct <- ggplot(celltype_counts, aes(x = celltype, y = n)) +
  geom_col(fill = "dodgerblue") +
  labs(y = "Number of Cells", title = "Cell Type") +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

g_cond <- ggplot(condition_counts, aes(x = condition, y = n)) +
  geom_col(fill = "olivedrab") +
  labs(y = "Number of Cells", title = "Condition") +
  base_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Convert to gtables
gt1 <- ggplotGrob(g_ct)
gt2 <- ggplotGrob(g_cond)

# Make widths identical across the two plots
max_widths <- unit.pmax(gt1$widths, gt2$widths)
gt1$widths <- max_widths
gt2$widths <- max_widths

# Make heights identical across the two plots
max_heights <- unit.pmax(gt1$heights, gt2$heights)
gt1$heights <- max_heights
gt2$heights <- max_heights

g_overview <- grid.arrange(gt1, gt2, ncol = 2)

ggsave(
  plot = g_overview,
  device = "png",
  filename = "~/Documents/log1pNMF/inst/paper_figures/images/lsa_overview.png",
  width = 6,
  height = 3.5
)

# expr_summary_all <- lapply(colnames(counts), function(gene) {
# 
#   barcodes$gene <- as.numeric(counts[, gene])
# 
#   expr_summary <- barcodes %>%
#     dplyr::group_by(condition, celltype) %>%
#     dplyr::summarise(
#       bulk_expr = sum(gene) / sum(s),
#       .groups = "drop"
#     )
# 
#   expr_summary$gene <- gene
#   expr_summary
# }) %>%
#   dplyr::bind_rows()
# 
# readr::write_csv(
#   expr_summary_all, "~/Documents/log1pNMF/inst/paper_figures/data/lsa_bulk.csv"
# )

expr_summary_all <- readr::read_csv(
  unzip(
    "~/Documents/log1pNMF/inst/paper_figures/data/lsa_bulk.csv.zip",
    "lsa_bulk.csv"
    )
  )

beta_df <- expr_summary_all %>%
  dplyr::filter(celltype == "Beta")

min_expr <- beta_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    min_expr = min(bulk_expr)
  )

high_enough_genes <- min_expr %>%
  dplyr::filter(min_expr >= 0.1) %>%
  dplyr::pull(gene)

beta_df <- beta_df %>%
  dplyr::filter(gene %in% high_enough_genes)


bad_genes <- c("Aldoa", "Gcat", "Pick1")

beta_df <- beta_df %>%
  dplyr::filter(!(gene %in% bad_genes))

treated_conditions <- c("IFNg", "IL-1B", "IL-1B_IFNg")

result_df <- beta_df %>%
  select(gene, condition, bulk_expr) %>%
  filter(condition %in% c("Untreated", treated_conditions)) %>%
  pivot_wider(
    names_from = condition,
    values_from = bulk_expr
  ) %>%
  pivot_longer(
    cols = all_of(treated_conditions),
    names_to = "condition",
    values_to = "treated_expr"
  ) %>%
  mutate(
    untreated_expr = Untreated,
    abs_diff = abs(treated_expr - untreated_expr),
    higher_condition = case_when(
      treated_expr > untreated_expr ~ condition,
      treated_expr < untreated_expr ~ "Untreated",
      TRUE ~ "Tie"
    ),
    higher_expr = pmax(treated_expr, untreated_expr),
    lower_expr  = pmin(treated_expr, untreated_expr),
    lop1pfc = log1p(higher_expr) - log1p(lower_expr),
    baseline_raw = lower_expr,
    baseline_log1p = log1p(lower_expr)
  ) %>%
  dplyr::select(
    gene,
    condition,
    lop1pfc,
    abs_diff,
    baseline_raw,
    baseline_log1p,
    higher_condition
  )

# now, it's worth plotting the differences on both a linear and a log scale

library(ggplot2)

g_log1p <- ggplot(data = (
  result_df %>% 
    dplyr::filter(condition == "IFNg")
)
) +
  geom_point(
    aes(x = baseline_log1p, y = lop1pfc), size = 1, alpha = 0.25
  ) +
  geom_smooth(
    aes(x = baseline_log1p, y = lop1pfc), se = FALSE, linewidth = 0.5
  ) +
  xlab("log1p(Lower Cond.)") +
  ylab("log1p(Higher Cond.) - log1p(Lower Cond.)") +
  ggtitle("log1p")  +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

g_baseline <- ggplot(data = (
  result_df %>% 
    dplyr::filter(condition == "IFNg")
)
) +
  geom_point(
    aes(x = baseline_raw, y = abs_diff), size = 1, alpha = 0.25
  ) +
  geom_smooth(
    aes(x = baseline_raw, y = abs_diff), se = FALSE, linewidth = 0.5
  ) +
  xlab("Lower Cond.") +
  ylab("|Higher Cond. - Lower Cond.|") +
  ggtitle("Raw")  +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

g_log1p_zoomed <- ggplot(data = (
  result_df %>% 
    dplyr::filter(condition == "IFNg" & baseline_raw < 500)
)
) +
  geom_point(
    aes(x = baseline_log1p, y = lop1pfc), size = 1, alpha = 0.25
  ) +
  geom_smooth(
    aes(x = baseline_log1p, y = lop1pfc), se = FALSE, linewidth = 0.5
  ) +
  xlab("log1p(Lower Cond.)") +
  ylab("log1p(Higher Cond.) - log1p(Lower Cond.)") +
  ggtitle("log1p (Zoomed in)")  +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

g_baseline_zoomed <- ggplot(data = (
  result_df %>% 
    dplyr::filter(condition == "IFNg" & baseline_raw < 500)
)
) +
  geom_point(
    aes(x = baseline_raw, y = abs_diff), size = 1, alpha = 0.25
  ) +
  geom_smooth(
    aes(x = baseline_raw, y = abs_diff), se = FALSE, linewidth = 0.5
  ) +
  xlab("Lower Cond") +
  ylab("|Higher Cond. - Lower Cond.|") +
  ggtitle("Raw (Zoomed in)")  +
  cowplot::theme_cowplot() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

g_diff <- ggarrange(
  g_baseline, g_log1p,
  g_baseline_zoomed, g_log1p_zoomed,
  nrow = 2, ncol = 2
)

ggsave(
  plot = g_diff,
  device = "png",
  filename = "~/Documents/log1pNMF/inst/paper_figures/images/lsa_diff.png",
  width = 9,
  height = 8
)

# another potential thing to plot:
# A + B vs. C
# log(A) + log(B) vs. log(C)

combo_df <- beta_df %>%
  dplyr::select(gene, condition, bulk_expr) %>%
  dplyr::filter(condition %in% c("Untreated", "IFNg", "IL-1B", "IL-1B_IFNg")) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = bulk_expr
  ) %>%
  dplyr::mutate(
    # raw-scale single-treatment differences
    diff_ifng_raw   = IFNg - Untreated,
    diff_il1b_raw   = `IL-1B` - Untreated,
    diff_combo_raw  = `IL-1B_IFNg` - Untreated,
    
    # raw-scale additive expectation from singles
    sum_single_raw = diff_ifng_raw + diff_il1b_raw,
    
    # log1p-scale condition values
    untreated_log1p = log1p(Untreated),
    ifng_log1p      = log1p(IFNg),
    il1b_log1p      = log1p(`IL-1B`),
    combo_log1p     = log1p(`IL-1B_IFNg`),
    
    # log1p-scale single-treatment differences
    diff_ifng_log1p  = ifng_log1p - untreated_log1p,
    diff_il1b_log1p  = il1b_log1p - untreated_log1p,
    diff_combo_log1p = combo_log1p - untreated_log1p,
    
    # log1p-scale additive expectation from singles
    sum_single_log1p = diff_ifng_log1p + diff_il1b_log1p,
    
    # optional: deviation from additivity
    interaction_raw   = diff_combo_raw - sum_single_raw,
    interaction_log1p = diff_combo_log1p - sum_single_log1p
  )


ggplot(data = combo_df) +
  geom_point(aes(x = sum_single_raw, y = diff_combo_raw)) +
  geom_abline(
    slope = 1, intercept = 0, color = "red", linetype = "dashed"
  )


ggplot(data = combo_df) +
  geom_point(aes(x = sum_single_log1p, y = diff_combo_log1p)) +
  geom_abline(
    slope = 1, intercept = 0, color = "red", linetype = "dashed"
  )
