counts <- Matrix::t(readRDS("scdata.rds"))
md <- read.csv("whole_cyto_annot.csv")

library(dplyr)

### specifiy the 45 conditions measured on the same batch
trts <- c("Ctrl_2", "CCL20", "CXCL1", "CCL22", "CXCL5", "CCL11", "CCL4", "CCL17", "CCL5", "CXCL13", "CXCL10", "CXCL9",
          "CXCL12", "GCSF", "MCSF", "GMCSF", "IFNg", "IL10", "IL12p70", "IL17a", "IL13", "IL15", "IL17f", "IL22",
          "IL18", "IL1a", "IL2", "IL3", "IL1b", "IL23", "IL21", "IL33", "IL25", "IL34", "IL36a",
          "IL4", "IL6", "IL5", "IL7", "IL9", "IL11", "TGFb", "CCL2", "CCL3", "TSLP")

md_sum <- md %>%
  dplyr::group_by(cell_type) %>%
  summarize(tot = n())

# for simplicity, I think that it would be useful to just start
# with the neutrophils. Later, I can think about multiple
# cell types with multiple conditions

md$s <- Matrix::rowSums(counts)

md <- md %>%
  dplyr::filter(cell_type == "Neutrophils" & (sample == "IL1a" | sample == "Ctrl_2"))

counts <- counts[md$X0, ]

counts <- counts[, Matrix::colSums(counts) > 0]
counts <- counts[Matrix::rowSums(counts) > 0, ]

est_vec <- c()
se_vec <- c()

for (gene in colnames(counts)) {

  md$cts <- counts[, gene]
  mod <- glm(
    cts ~ offset(log(s)) + sample,
    data = md,
    family = poisson()
  )

  est_vec <- c(est_vec, coef(summary(mod))["sampleIL1a", "Estimate"])
  se_vec <- c(se_vec, coef(summary(mod))["sampleIL1a", "Std. Error"])

}

res_df <- data.frame(
  gene = colnames(counts),
  est = est_vec,
  se = se_vec
)

a <- ashr::ash(res_df$est, res_df$se)

res_df$pm <- ashr::get_pm(a)
res_df$lfsr <- ashr::get_lfsr(a)

sig_df <- res_df %>% dplyr::filter(pm > 0 & lfsr < .05)

ctrl_expr_vec <- c()
Il_expr_vec <- c()

for (gene in sig_df$gene) {

  md$cts <- counts[, gene]
  md_ctrl <- md %>% dplyr::filter(sample == "Ctrl_2")
  md_Il <- md %>% dplyr::filter(sample != "Ctrl_2")

  ctrl_expr <- sum(md_ctrl$cts) / sum(md_ctrl$s)
  Il_expr <- sum(md_Il$cts) / sum(md_Il$s)

  ctrl_expr_vec <- c(ctrl_expr_vec, ctrl_expr)
  Il_expr_vec <- c(Il_expr_vec, Il_expr)

}

sig_df$ctrl_expr <- ctrl_expr_vec
sig_df$Il_expr <- Il_expr_vec

sig_df$abs_diff <- sig_df$Il_expr - sig_df$ctrl_expr

ggplot(data = sig_df) +
  geom_point(aes(x = 10000 * ctrl_expr, y = 10000 * abs_diff)) +
  xlab("Control Count / 10K") +
  ylab("Difference Between Stim and Control Count / 10K") +
  cowplot::theme_cowplot()

s2 <- sig_df %>% dplyr::filter(ctrl_expr < .01 & abs_diff < .025)

# Now, I want to run the passPCA algorithm on this
# it shouldn't take too long
set.seed(1)
log1p_fit <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = 15,
  maxiter = 100,
  approx_range = c(0, 1.25),
  s = as.vector(Matrix::rowSums(counts) / mean(Matrix::rowSums(counts)))
)

library(fastTopics)


nmf_fit <- fastTopics::fit_poisson_nmf(
  X = counts,
  k = 15
)

U_normalized <- apply(log1p_fit$U, 2, function(col) col / max(col))

md <- md %>%
  dplyr::mutate(ct_t = glue::glue("{cell_type} {sample}"))

structure_plot(nmf_fit, grouping = md$ct_t)
structure_plot(U_normalized, grouping = md$ct_t)
#structure_plot(U_normalized, grouping = md$sample)
