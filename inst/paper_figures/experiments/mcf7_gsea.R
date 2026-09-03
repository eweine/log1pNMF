# GO over-representation of the top-10 genes of each MCF-7 factor, for
# the c = 1 and c = Inf (topic model) fits shown in the paper -- a trial
# of the factor_gsea.R pipeline on the first biological example.
#
#   Rscript mcf7_gsea.R        (run from this directory)
#
# The factors are the same ones as in the paper's Figure/table: the topic
# model's columns are matched to the c = 1 reference fit with the
# Hungarian/Spearman procedure of factor_matching.R, so "k2" means the
# RA-associated factor in both models. Gene identifiers in the stored
# fits are NCBI GeneIDs; display symbols come from the same GEO
# annotation file used to build the count matrix.
#
# Writes mcf7_gsea.csv (one row per factor x enriched GO BP term).

source("factor_gsea.R")
source("../fig_scripts/factor_matching.R")

load("../data/experiment_results.Rdata")
log1p_fit <- res_list$mcf7$`1`
tm_fit    <- res_list$mcf7$`Inf`

F1 <- as.matrix(log1p_fit$FF)
colnames(F1) <- paste0("k", 1:3)
Ft <- as.matrix(tm_fit$F)
colnames(Ft) <- paste0("k", match_display_labels(F1, Ft))
Ft <- Ft[, paste0("k", 1:3)]

## GeneID -> symbol, from the GEO annotation used to build the matrix
ann <- read.delim(gzfile("../data/raw_data/Human.GRCh38.p13.annot.tsv.gz"),
                  sep = "\t", quote = "", check.names = FALSE)
sym <- setNames(as.character(ann$Symbol), as.character(ann$GeneID))

run_one <- function(FF, model) {
  message("=== ", model, " ===")
  d <- factor_ora(FF, n_top = 10, id_type = "entrez", symbols = sym)
  if (nrow(d) > 0) d <- cbind(model = model, d)
  d
}

res <- rbind(run_one(F1, "c = 1"), run_one(Ft, "c = Inf"))
write.csv(res, "mcf7_gsea.csv", row.names = FALSE)
message("Wrote mcf7_gsea.csv (", nrow(res), " rows)")

## ranked GSEA on the full gene-score vectors, for comparison
run_ranked <- function(FF, model) {
  message("=== ranked, ", model, " ===")
  set.seed(1)
  d <- factor_gsea_ranked(FF, id_type = "entrez", symbols = sym)
  if (nrow(d) > 0) d <- cbind(model = model, d)
  d
}
resr <- rbind(run_ranked(F1, "c = 1"), run_ranked(Ft, "c = Inf"))
write.csv(resr, "mcf7_gsea_ranked.csv", row.names = FALSE)
message("Wrote mcf7_gsea_ranked.csv (", nrow(resr), " rows)")

## console summary: top terms per (model, factor), both methods
show <- function(res, gene_col) {
  for (m in unique(res$model)) for (k in unique(res$factor[res$model == m])) {
    d <- res[res$model == m & res$factor == k, ]
    cat(sprintf("\n%s  %s  (%d terms at FDR < 0.05)\n", m, k, nrow(d)))
    for (i in seq_len(min(4, nrow(d))))
      cat(sprintf("  %-55s p_adj=%.2g  [%s]\n",
                  substr(d$term[i], 1, 55), d$p_adj[i], d[[gene_col]][i]))
  }
}
cat("\n================ TOP-10 ORA ================\n"); show(res, "genes")
cat("\n================ RANKED GSEA ===============\n"); show(resr, "leading_edge")
