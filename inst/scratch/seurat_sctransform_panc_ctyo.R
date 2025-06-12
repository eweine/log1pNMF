library(data.table)
library(Matrix)
library(tools)
library(rsvd)
library(uwot)
library(ggplot2)
library(cowplot)

set.seed(1)

read_geo_data <- function (i, r, prefix = "GSM000000", dir = ".") {
  infile   <- sprintf("%s_Rep%d_S%d_barcodes.tsv.gz",prefix,r,i)
  barcodes <- fread(file.path(dir,infile),quote = FALSE,header = FALSE,
                    stringsAsFactors = FALSE)
  class(barcodes) <- "data.frame"
  barcodes <- barcodes[,1]
  infile <- sprintf("%s_Rep%d_S%d_features.tsv.gz",prefix,r,i)
  genes  <- fread(file.path(dir,infile),sep = "\t",quote = FALSE,
                  header = FALSE,stringsAsFactors = FALSE)
  class(genes) <- "data.frame"
  names(genes) <- c("ensembl","symbol","type")
  genes        <- transform(genes,type = factor(type))
  infile <- sprintf("%s_Rep%d_S%d_matrix.mtx.gz",prefix,r,i)
  counts <- fread(file.path(dir,infile),quote = FALSE,header = FALSE,
                  skip = 2)
  class(counts) <- "data.frame"
  names(counts) <- c("row","col","value")
  n <- max(counts$row)
  m <- max(counts$col)
  counts <- sparseMatrix(i = counts$row,j = counts$col,x = counts$value,
                         dims = c(n,m))
  rownames(counts) <- genes$ensembl
  colnames(counts) <- barcodes
  return(list(genes  = genes,
              counts = counts))
}
dat     <- vector("list",8)
dataset <- 0
samples <- NULL
for (r in 1:2) {
  for (i in 1:4) {
    dataset <- dataset + 1
    dat[[dataset]] <-
      read_geo_data(i,r,prefix = paste0("GSM55486",23 + dataset),
                    dir = "~/Downloads/GSE183010_RAW/")
    samples <- rbind(samples,
                     data.frame(barcode   = colnames(dat[[dataset]]$counts),
                                mouse     = paste0("S",i),
                                replicate = r,
                                stringsAsFactors = FALSE))
  }
}
samples <- transform(samples,
                     mouse     = factor(mouse),
                     replicate = factor(replicate))
features <- Reduce(intersect,lapply(dat,function (x) x$genes$ensembl))
genes    <- subset(dat[[1]]$genes,is.element(ensembl,features))
counts   <- do.call("cbind",lapply(dat,function (x) x$counts[features,]))
counts   <- t(counts)

par(mar = c(4,4,1,1))
x <- rowSums(counts > 0)
i <- which(x > 2000)
samples <- samples[i,]
counts  <- counts[i,]
hist(x,n = 64,main = "",xlab = "number of genes",ylab = "number of cells")

par(mar = c(4,4,1,1))
mito_genes <- which(substr(genes$symbol,1,2) == "mt")
s          <- rowSums(counts)
s_mito     <- counts[,mito_genes]
prop_mito  <- rowSums(s_mito)/s
i          <- which(prop_mito < 0.1)
samples <- samples[i,]
counts  <- counts[i,]
hist(prop_mito,n = 64,main = "",xlab = "proportion mitochondrial",
     ylab = "number of cells")


j      <- which(genes$symbol != "Malat1")
genes  <- genes[j,]
counts <- counts[,j]

x      <- colSums(counts > 0)
j      <- which(x > 0)
genes  <- genes[j,]
counts <- counts[,j]

library(Seurat)

j <- which(!duplicated(genes$symbol))
genes  <- genes[j,]
counts <- counts[,j]

colnames(counts) <- genes$symbol

rownames(counts) <- paste0(samples$barcode, "-", samples$mouse)
i          <- which(!duplicated(rownames(counts)))
samples <- samples[i,]
counts  <- counts[i,]

so <- CreateSeuratObject(counts = Matrix::t(counts), meta.data = samples)

so <- subset(x = so, subset = mouse == "S1")

rm(counts, samples)
gc()

so <- PercentageFeatureSet(so, pattern = "^mt-", col.name = "percent.mt")

so <- SCTransform(so, vars.to.regress = "percent.mt", verbose = TRUE)

so <- RunPCA(so, verbose = TRUE)
ElbowPlot(so)

so <- RunUMAP(so, dims = 1:10, verbose = TRUE)

so <- FindNeighbors(so, dims = 1:10, verbose = TRUE)
so <- FindClusters(so, verbose = TRUE)

DimPlot(so, reduction = "umap", label = TRUE)

# VlnPlot(so, features = c(
#   "Ppy", "Sst","Rbp4", "Gcg"
#   ), idents = c(7, 9),
#         pt.size = 0)
# 
# # 5, 6, 7, and 9 appear to be some mix of 
# # alpha, PP (gamma), and delta cells
# 
# adg_cells <- subset(so, subset = seurat_clusters %in% c(7, 9))
# adg_md <- adg_cells@meta.data
# adg_md$ppy <- adg_cells@assays$RNA$counts["Ppy", rownames(adg_md)]
# adg_md$sst <- adg_cells@assays$RNA$counts["Sst", rownames(adg_md)]
# adg_md$rbp4 <- adg_cells@assays$RNA$counts["Rbp4", rownames(adg_md)]
# 
# 
# plot(log1p(adg_md$ppy), log1p(adg_md$rbp4 + adg_md$sst))
# 
# # clusters 5 and 6 are alpha, if Gcg > 4 
# 
# 
# 
# 
# VlnPlot(so, features = c(
# "Pecam1"
# ),
# pt.size = 0)
# 
# em_cells <- subset(so, subset = seurat_clusters %in% c(11, 13))
# md_em <- em_cells@meta.data
# 
# md_em$cd34 <- em_cells@assays$RNA$counts["Cd34", rownames(md_em)]
# md_em$col1a <- em_cells@assays$RNA$counts["Col1a1", rownames(md_em)] + em_cells@assays$RNA$counts["Col1a2", rownames(md_em)]
# 
# plot(log1p(md_em$cd34), log1p(md_em$col1a))

md <- so@meta.data
md$gcg <- so@assays$RNA$counts["Gcg", ]

# here, I want to build a map from cluster to cell type label
celltype_df <- data.frame(
  cluster = as.factor(c(
    0, 
    1, 
    2,
    3,
    4,
    5,
    6,
    7,
    8, 
    9,
    10,
    11,
    12,
    13,
    14
  )),
  celltype_init = c(
    "Beta", # 0
    "Beta", # 1
    "Beta", # 2
    "Beta", # 3
    "Beta", # 4
    "Alpha", # 5 IFF log1p(Gcg) > 4, Otherwise D/G
    "Alpha", # 6 IFF log1p(Gcg) > 4, Otherwise D/G
    "D/G", # 7
    "Beta", # 8
    "D/G", # 9
    "Macrophage", # 10
    "E/M", # 11
    "E/M", # 12
    "E/M", # 13
    "Ductal" # 14
  )
)

library(dplyr)
md <- md %>%
  dplyr::inner_join(
    celltype_df,
    by = c("seurat_clusters" = "cluster")
  )

md <- md %>%
  dplyr::mutate(
    celltype_final = case_when(
      seurat_clusters %in% c(5, 6) & log1p(gcg) > 4 ~ "Alpha",
      seurat_clusters %in% c(5, 6) ~ "D/G",
      TRUE ~ celltype_init
    )
  )

md$umap1 <- so@reductions$umap@cell.embeddings[,"umap_1"]
md$umap2 <- so@reductions$umap@cell.embeddings[,"umap_2"]

library(ggplot2)

ggplot(data = md) +
  geom_point(aes(x = umap1, y = umap2, colour = celltype_final))

md <- md %>%
  dplyr::select(barcode, mouse, celltype_final)

readr::write_rds(md, "~/Downloads/panc_ctyo_S1_celltypes.rds")

