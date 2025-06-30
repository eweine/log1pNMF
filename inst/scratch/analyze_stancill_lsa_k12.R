# here, I want to analyze the different factor models and see how they
# represent both the relevant celltypes and treatments

load("~/Downloads/panc_cyto_lsa.Rdata")

log1p_fit12 <- readr::read_rds(
  "~/Downloads/panc_cyto_lsa_res/lsa_pancreas_cytokine_log1p_c1_rank1_init_K12.rds"
  )

library(log1pNMF)

normalized_structure_plot(
  log1p_fit12, 
  grouping = barcodes$condition, 
  gap = 20, 
  n = Inf
  )

# now, want to identify cell type specific and treatment specific factors
# celltypes:
# k1 (beta)
# k2 (alpha)
# k3 (delta, somewhat mixed)
# k5 (macrophage)
# k6 (ductal)
# k7 (EM)
# k8 (acinar)
# k9 (delta)

# treatments:
# k4 (IFNg)
# k10 (not IL-1B)
# k11 (IL-1B)
# k12 (?)

normalized_structure_plot(
  log1p_fit12,
  grouping = barcodes$celltype,
  topics = c(
    paste0("k", "_", c(1, 2, 3, 5, 6, 7, 8, 9))
  ),
  gap = 20
)

# this looks good, though some pockets of purple and yellow could be misclassification

normalized_structure_plot(
  log1p_fit12,
  grouping = barcodes$celltype,
  topics = c(
    paste0("k", "_", c(1, 2, 3, 5, 6, 7, 8, 9))
  ),
  gap = 20
)



normalized_structure_plot(
  log1p_fit12,
  grouping = barcodes$condition,
  topics = c(
    paste0("k", "_", c(4, 10, 11, 12))
  ),
  gap = 20
)

# now, I'd like to checkout the topic model fits

tm_fit <- readr::read_rds("~/Downloads/panc_cyto_lsa_res/tm_k12.rds")

library(fastTopics)

structure_plot(
  tm_fit,
  n = Inf,
  grouping = barcodes$celltype,
  gap = 20
)

# celltype factors
# k1 (?)
# k2 (acinar)
# k3 (EM)
# k4 (?)
# k5 (Ductal)
# k6 (? is beta specific but could be treatment related)
# k7 (?)
# k8 (Delta)
# k9 (Alpha)
# k10 (Gamma)
# k11 (Macrophage)
# k12 (?)

structure_plot(
  tm_fit,
  n = Inf,
  grouping = barcodes$celltype,
  gap = 20,
  topics = paste0("k", c(2, 3, 5, 6, 7, 8, 9, 10, 11, 12))
)


structure_plot(
  tm_fit,
  n = Inf,
  grouping = barcodes$condition,
  gap = 40,
  topics = paste0("k", c(1, 4))
)

# treatments:
# k4 (IL-1B)
# 

