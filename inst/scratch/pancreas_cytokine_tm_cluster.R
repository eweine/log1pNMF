library(Matrix)
load("~/Downloads/pancreas_cytokine.RData")

i       <- which(samples$mouse == "S1")
samples <- samples[i,]
counts  <- counts[i,]

outliers <- c("TTTGTTGTCGTTAGTG-1","TTTGTTGGTAGAGCTG-1")
i        <- which(!is.element(samples$barcode,outliers))
samples  <- samples[i,]
counts   <- counts[i,]

j      <- which(colSums(counts > 0) > 2)
genes  <- genes[j,]
counts <- counts[,j]

library(fastTopics)
set.seed(1)
tm <- fit_poisson_nmf(counts,k = 12, control = list(numiter = 100,nc = 7))

structure_plot(tm)

scale_rows <- function (A)
  A / apply(A,1,max)
# using marker genes, I would like to connect each topic to a celltype
marker_genes <- c("Ins1","Ins2","Gcg","Sst","Ghrl",
                  "Ppy","Chga","Iapp","Krt19","Ccr5","Pecam1","Esam",
                  "Col1a1","Ghrl", "Cd34", "Krt17", "Cpa1")
j <- match(marker_genes,genes$symbol)
F <- poisson2multinom(tm)$F
F <- F[j,]
F <- scale_rows(F)
rownames(F) <- marker_genes
p <- annotation_heatmap(F,select_features = "all",verbose = FALSE)
print(p)
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

L <- poisson2multinom(tm)$L

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

structure_plot(tm, grouping = L_df$cluster, gap = 20)
