counts1 <- Matrix::t(
  Matrix::readMM("~/Documents/single-cell-jamboree/data/GSE156175_RAW/GSM5842388_Rep1_S1_matrix.mtx")
  )

genes1 <- readr::read_tsv(
  "~/Documents/single-cell-jamboree/data/GSE156175_RAW/GSM5842388_Rep1_S1_features.tsv",
  col_names = c("ensembl", "symbol", "type")
  )

barcodes1 <- readr::read_tsv(
  "~/Documents/single-cell-jamboree/data/GSE156175_RAW/GSM5842388_Rep1_S1_barcodes.tsv",
  col_names = c("barcode")
  )

barcodes1$barcode <- substr(barcodes1$barcode, 1, nchar(barcodes1$barcode) - 2)

barcodes1$condition <- "Untreated"

rownames(counts1) <- paste0(barcodes1$barcode, "_", barcodes1$condition)
colnames(counts1) <- genes1$symbol


counts2 <- Matrix::t(
  Matrix::readMM("~/Documents/single-cell-jamboree/data/GSE156175_RAW/GSM4726017_S2_matrix.mtx")
)

genes2 <- readr::read_tsv(
  "~/Documents/single-cell-jamboree/data/GSE156175_RAW/GSM4726017_S2_features.tsv",
  col_names = c("ensembl", "symbol", "type")
)

barcodes2 <- readr::read_tsv(
  "~/Documents/single-cell-jamboree/data/GSE156175_RAW/GSM4726017_S2_barcodes.tsv",
  col_names = c("barcode")
)

barcodes2$barcode <- substr(barcodes2$barcode, 1, nchar(barcodes2$barcode) - 2)

barcodes2$condition <- "IL-1B"

rownames(counts2) <- paste0(barcodes2$barcode, "_", barcodes2$condition)
colnames(counts2) <- genes2$symbol

counts3 <- Matrix::t(
  Matrix::readMM("~/Documents/single-cell-jamboree/data/GSE156175_RAW/GSM4726018_S3_matrix.mtx")
)

genes3 <- readr::read_tsv(
  "~/Documents/single-cell-jamboree/data/GSE156175_RAW/GSM4726018_S3_features.tsv",
  col_names = c("ensembl", "symbol", "type")
)

barcodes3 <- readr::read_tsv(
  "~/Documents/single-cell-jamboree/data/GSE156175_RAW/GSM4726018_S3_barcodes.tsv",
  col_names = c("barcode")
)

barcodes3$barcode <- substr(barcodes3$barcode, 1, nchar(barcodes3$barcode) - 2)

barcodes3$condition <- "IFNg"

rownames(counts3) <- paste0(barcodes3$barcode, "_", barcodes3$condition)
colnames(counts3) <- genes3$symbol

counts4 <- Matrix::t(
  Matrix::readMM("~/Documents/single-cell-jamboree/data/GSE156175_RAW/GSM4726019_S4_matrix.mtx")
)

genes4 <- readr::read_tsv(
  "~/Documents/single-cell-jamboree/data/GSE156175_RAW/GSM4726019_S4_features.tsv",
  col_names = c("ensembl", "symbol", "type")
)

barcodes4 <- readr::read_tsv(
  "~/Documents/single-cell-jamboree/data/GSE156175_RAW/GSM4726019_S4_barcodes.tsv",
  col_names = c("barcode")
)

barcodes4$barcode <- substr(barcodes4$barcode, 1, nchar(barcodes4$barcode) - 2)

barcodes4$condition <- "IL-1B_IFNg"

rownames(counts4) <- paste0(barcodes4$barcode, "_", barcodes4$condition)
colnames(counts4) <- genes4$symbol

counts <- rbind(
  as.matrix(counts1),
  as.matrix(counts2),
  as.matrix(counts3),
  as.matrix(counts4)
)

counts <- as(counts, "CsparseMatrix")

barcodes <- rbind(
  barcodes1, barcodes2, barcodes3, barcodes4
)

barcodes$cell_bc <- paste0(barcodes$barcode, "_", barcodes$condition)

rm(counts1, counts2, counts3, counts4, barcodes1, barcodes2, barcodes3, barcodes4,
   genes1, genes2, genes3, genes4
   )
gc()

counts <- counts[, Matrix::colSums(counts) > 0]

x <- rowSums(counts > 0)
i <- which(x > 2000)
counts  <- counts[i,]

outliers <- rownames(counts)[Matrix::rowSums(counts) > 60000]
i        <- which(!is.element(rownames(counts),outliers))
counts   <- counts[i,]

library(fastTopics)

mito_genes <- which(substr(colnames(counts),1,2) == "mt")
s          <- rowSums(counts)
s_mito     <- counts[,mito_genes]
prop_mito  <- rowSums(s_mito)/s
i          <- which(prop_mito < 0.1)
counts  <- counts[i,]

j <- which(!(grepl("^mt-",colnames(counts)) |
               grepl("^Rp[sl]",colnames(counts))))
counts <- counts[,j]

j      <- which(colnames(counts) != "Malat1")
counts <- counts[,j]

j      <- which(colSums(counts > 0) > 2)
counts <- counts[,j]

rm(s_mito, i, j, prop_mito, s, x)
gc()

barcodes <- barcodes %>%
  dplyr::filter(cell_bc %in% rownames(counts))

set.seed(1)
# tm_fit <- fit_poisson_nmf(
#   X = as(counts, "CsparseMatrix"),
#   k = 12,
#   control = list(nc = 7)
# )

tm_fit <- readr::read_rds("~/Documents/single-cell-jamboree/data/panc_cyto_lsa_tm_k12.rds")

structure_plot(tm_fit, n = Inf)


scale_rows <- function (A)
  A / apply(A,1,max)
# using marker genes, I would like to connect each topic to a celltype
marker_genes <- c("Ins1","Ins2","Gcg","Sst",
                  "Ppy","Krt17","Ccr5",
                  "Col1a1", "Cd34", "Cpa1")
F <- poisson2multinom(tm_fit)$F
F <- F[marker_genes,]
F <- scale_rows(F)
rownames(F) <- marker_genes
p <- annotation_heatmap(F,select_features = "all",verbose = FALSE)
print(p)
# factor associations
# k1: ?
# k2: Acinar
# k3: Endothelial/Mesnchymal
# k4: ? (presumably Beta)
# k5: Ductal
# k6: Beta
# k7: ? (presumably Beta)
# k8: Delta
# k9: Alpha
# k10: Gamma (PP) 
# k11: Macrophage
# k12: ? (presumably Beta)

L <- poisson2multinom(tm_fit)$L

L_df <- as.data.frame(L)

L_df <- L_df %>%
  dplyr::mutate(
    cluster = case_when(
      k2 > (1/3) ~ "Acinar",
      k3 > (1/3) ~ "Endothelial/Mesnchymal",
      k5 > (1/3) ~ "Ductal",
      k8 > (1/3) ~ "Delta",
      k9 > (1/3) ~ "Alpha",
      k10 > (1/3) ~ "Gamma",
      k11 > (1/3) ~ "Macrophage",
      TRUE ~ "Beta"
    )
  )

structure_plot(tm_fit, grouping = L_df$cluster, gap = 20, n = Inf)

barcodes$celltype <- L_df$cluster

library(Seurat)

counts <- counts[,!duplicated(colnames(counts))]

so <- CreateSeuratObject(counts = Matrix::t(counts), meta.data = barcodes)

so <- subset(x = so, subset = celltype %in% c("Beta", "Alpha", "Delta"))

so <- SCTransform(so, verbose = TRUE)

so <- RunPCA(so)

so <- RunUMAP(so, dims = 1:7)

DimPlot(so, group.by = "celltype", reduction = "pca")

save(counts, barcodes, file="~/Downloads/panc_cyto_lsa.Rdata")

structure_plot(tm_fit, grouping = paste0(barcodes$celltype, "_", barcodes$condition), gap = 20, n = Inf)

readr::write_rds(tm_fit, "~/Downloads/panc_cyto_lsa_res/tm_k12.rds")

set.seed(1)
tm_fit15 <- fit_poisson_nmf(
  X = as(counts, "CsparseMatrix"),
  k = 15,
  control = list(nc = 7)
)

readr::write_rds(tm_fit15, "~/Downloads/panc_cyto_lsa_res/tm_k15.rds")

structure_plot(tm_fit15, grouping = paste0(barcodes$celltype, "_", barcodes$condition), gap = 20, n = Inf)

