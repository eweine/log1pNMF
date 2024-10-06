load("~/Documents/data/fastglmpca/raw_data/pbmc_68k.RData")

samp <- names(sort(log1p_mod$V[,1], decreasing = TRUE))[1:10]

immune_genes <- c(
  "TNF", "IFNG", "IL6", "IL10",
  "IL1B", "CD4", "CD8A", "CD28",
  "CD40", "FOXP3"
)

library(clusterProfiler)
library(fgsea)
library(AnnotationDbi)
library(org.Mm.eg.db)

entrez_ids <- mapIds(
  org.Mm.eg.db,
  keys = samp,
  column = "ENTREZID",
  keytype = "SYMBOL"
)

go_result <- enrichGO(gene = entrez_ids,
                      OrgDb = org.Mm.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP", # Biological Process. Can also use "MF" for Molecular Function or "CC" for Cellular Component
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.1)

barplot(go_result, showCategory = 10)

