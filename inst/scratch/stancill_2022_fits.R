library(Matrix)
library(log1pNMF)

load("~/Downloads/pancreas_cytokine.RData")

set.seed(1)
outliers <- rownames(counts)[Matrix::rowSums(counts) > 80000]
i        <- which(!is.element(samples$barcode,outliers))
samples  <- samples[i,]
counts   <- counts[i,]

j      <- which(colSums(counts > 0) > 2)
genes  <- genes[j,]
counts <- counts[,j]

celltypes_S1 <- readr::read_csv(
  "~/Documents/single-cell-jamboree/data/pancreas_cytokine_S1_celltypes.csv"
  )

celltypes_S1$mouse <- "S1"

celltypes_S2 <- readr::read_csv(
  "~/Documents/single-cell-jamboree/data/pancreas_cytokine_S2_celltypes.csv"
)

celltypes_S2$mouse <- "S2"

celltypes_S3 <- readr::read_csv(
  "~/Documents/single-cell-jamboree/data/pancreas_cytokine_S3_celltypes.csv"
)

celltypes_S3$mouse <- "S3"

celltypes_S4 <- readr::read_csv(
  "~/Documents/single-cell-jamboree/data/pancreas_cytokine_S4_celltypes.csv"
)

celltypes_S4$mouse <- "S4"

celltypes <- rbind(
  celltypes_S1,
  celltypes_S2,
  celltypes_S3,
  celltypes_S4
)

celltypes$barcode <- gsub("\\.", "-", celltypes$barcode)

samples <- samples %>%
  dplyr::inner_join(celltypes)

tm_fit15 <- readr::read_rds("~/Downloads/pancreas_cytokine_full_res/pancreas_cytokine_tm_K15.rds")

library(fastTopics)

structure_plot(
  tm_fit15, 
  grouping = paste0(samples$celltype, samples$mouse),
  n = Inf,
  gap = 20
  )

log1p_fit15 <- readr::read_rds("~/Downloads/pancreas_cytokine_full_res/pancreas_cytokine_log1p_c1_rank1_init_K15.rds")

log1pNMF::normalized_structure_plot(
  log1p_fit15, 
  grouping = paste0(samples$celltype, samples$mouse),
  n = Inf,
  gap = 20
)




