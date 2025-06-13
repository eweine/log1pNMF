# I think that the correct way to cluster the other datasets is to start with
# (a) the same vector of factors and
# (b) A non-negative least squares estimate (minimizing the KL) of loadings

# hopefully this will give a reasonably good initialization which should make
# it fairly easy to figure out what the cell types are in the other datasets

tm_fit <- readr::read_rds("~/Downloads/pancreas_cytokine_tm_k12.rds")

load("~/Downloads/pancreas_cytokine.RData")

i       <- which(samples$mouse == "S4")
samples <- samples[i,]
counts  <- counts[i,]

outliers <- rownames(counts)[Matrix::rowSums(counts) > 80000]
i        <- which(!is.element(samples$barcode,outliers))
samples  <- samples[i,]
counts   <- counts[i,]

#j      <- which(colSums(counts > 0) > 2)
#genes  <- genes[j,]
#counts <- counts[,j]

#set.seed(1)
#tm <- fit_poisson_nmf(X = counts, k = 14, control = list(nc = 7))

# 
counts <- counts[, colnames(counts) %in% rownames(tm_fit$F)]
j      <- which(colSums(counts > 0) > 0)
genes  <- genes[j,]
counts <- counts[,j]

tm_fit$F <- tm_fit$F[colnames(counts),]

library(fastTopics)
set.seed(1)
nm_pred <- predict(tm_fit, counts)

# 
fit0 <- init_poisson_nmf(counts, tm_fit$F, nm_pred)
# 
structure_plot(fit0, n = Inf)

# factor associations
# k1: Beta Cells
# k2: Macrophage Cells
# k3: Beta Cells
# k4: Beta Cells
# k5: Acinar Cells
# k6: Gamma (PP) Cells
# k7: Beta Cells
# k8: Alpha Cells
# k9: Delta Cells
# k10: Beta Cells
# k11: Ductal Cells
# k12: Endothelial / Mesnchymal Cells
library(dplyr)

L <- poisson2multinom(fit0)$L

L_df <- as.data.frame(L)

L_df <- L_df %>%
  dplyr::mutate(
    cluster = case_when(
      k2 > (1/3) ~ "Macrophage",
      k5 > (1/3) ~ "Acinar",
      k6 > (1/3) ~ "Gamma",
      k8 > (1/3) ~ "Alpha",
      k9 > (1/3) ~ "Delta",
      k11 > (1/3) ~ "Ductal",
      k12 > (1/3) ~ "Endothelial/Mesnchymal",
      TRUE ~ "Beta"
    )
  )

structure_plot(fit0, n = Inf, grouping = L_df$cluster, gap = 20)

celltype_df <- data.frame(
  barcode = rownames(L_df),
  celltype = L_df$cluster
)

readr::write_csv(celltype_df, "~/Downloads/pancreas_cytokine_S4_celltypes.csv")

