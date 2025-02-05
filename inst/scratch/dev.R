md <- readr::read_csv("~/Downloads/GSE189161_metadata.csv")
md <- md %>% dplyr::select(-`...1`)

cells <- readr::read_tsv("~/Downloads/GSE189161_cells.txt", col_names = FALSE)
colnames(cells) <- c("bc")

genes <- readr::read_tsv("~/Downloads/GSE189161_features.txt", col_names = FALSE)
colnames(genes) <- c("gene")

counts <- Matrix::readMM("~/Downloads/GSE189161_matrix.mtx")
counts <- Matrix::t(counts)
rownames(counts) <- cells$bc
colnames(counts) <- genes$gene
counts <- as(counts, "CsparseMatrix")


# it could be reasonable here to subset to only bone marrow, or
# only fetal liver...
genes_to_use <- which(Matrix::colSums(counts>0)>9)
counts <- counts[, genes_to_use]

readr::write_rds(counts, "~/Downloads/hspcs.rds")

