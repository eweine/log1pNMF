load("../data/raw_data/pancreas_cytokine_lsa.Rdata")

library(glmGamPoi)

# there are probably multiple different ways to do DE here, but I think 
# the simplest is likely to just go celltype by celltype instead of 
# estimating a model with interaction effects. However, a relevant
# question is how much of an impact the interaction effects have vs.
# the marginal effects

barcodes <- barcodes %>%
  dplyr::mutate(
    condition = case_when(
      condition == "IL-1B" ~ "IL_1B",
      condition == "IL-1B_IFNg" ~ "IL_1B_IFNg",
      TRUE ~ condition
    )
  )

cts <- unique(barcodes$celltype)
cts <- cts[cts != "Acinar"]

gp_out_list <- list()

for (ct in rev(cts)) {
  
  print(ct)
  bc_dat <- barcodes %>%
    dplyr::filter(
      celltype == ct
    )
  
  celltype_counts <- counts[bc_dat$cell_bc, ]
  # genes must be expressed in at least 5 cells
  genes_to_use <- which(Matrix::colSums(celltype_counts>0)>9)
  
  celltype_counts <- as.matrix(Matrix::t(celltype_counts[,genes_to_use]))
  gp_mod <- glm_gp(
    data = celltype_counts,
    design = ~ condition,
    col_data = bc_dat,
    reference_level = "Untreated",
    verbose = TRUE
  )
  
  gp_out_list[[ct]] <- list()
  
  for (test_cond in c("conditionIFNg", "conditionIL_1B", "conditionIL_1B_IFNg")) {
    
    gp_out_list[[ct]][[test_cond]] <- test_de(
      fit = gp_mod,
      contrast = test_cond,
      verbose = TRUE
    )
    
  }
  
}

load("../data/experiment_results.Rdata")

nmf_fit <- res_list$pancreas$`Inf`

F <- nmf_fit$F
colnames(F) <- paste0(
  "k", 
  c(11,12, 7, 6, 5, 1, 9, 4, 3, 2, 10, 13, 8)
)
F <- F[,paste0("k", 1:13)]


k10_df <- data.frame(
  f = unname(F[,"k10"]),
  gene = names(F[,"k10"])
)

lfc_df <- gp_out_list$Alpha$conditionIL_1B %>%
  dplyr::select(
    name, lfc
  ) %>%
  dplyr::rename(gene = name)

alpha_df <- barcodes %>%
  dplyr::filter(celltype == "Alpha")

alpha_counts <- counts[alpha_df$cell_bc, ]
gene_means <- Matrix::colSums(alpha_counts) / sum(alpha_counts)
tpm <- log2(1 + gene_means * 1e6)
tpm_df <- data.frame(
  gene = names(tpm),
  avg_expr = unname(tpm)
)

tpm_df <- tpm_df %>%
  dplyr::inner_join(lfc_df) %>%
  dplyr::inner_join(k10_df)

plot(tpm_df$lfc, log1p(tpm_df$f))
plot(tpm_df$avg_expr, log1p(tpm_df$f))

log1p_mod <- res_list$pancreas$`1`

log1p_k10_df <- data.frame(
  f_log1p = unname(log1p_mod$FF[,"k_10"]),
  gene = names(log1p_mod$FF[,"k_10"])
)

tpm_df <- tpm_df %>%
  dplyr::inner_join(log1p_k10_df)

plot(tpm_df$lfc, tpm_df$f_log1p)

barcodes <- barcodes %>%
  dplyr::filter(celltype != "Acinar")
counts <- counts[barcodes$cell_bc, ]

# each gene must be expressed in at least 25 cells
genes_to_use <- which(Matrix::colSums(counts>0)>24)
counts <- counts[, genes_to_use]

gp_full_mod <- glm_gp(
  data = as.matrix(Matrix::t(counts)),
  design = ~ condition * celltype,
  col_data = barcodes,
  reference_level = "Untreated",
  verbose = TRUE,
  subsample = TRUE
)


