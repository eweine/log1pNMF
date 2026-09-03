# GO over-representation of the top-50 genes of each pancreas factor
# (K = 13), for the c = 1 and c = Inf fits shown in the paper. Murine
# data, so annotations come from org.Mm.eg.db; gene identifiers in the
# fits are mouse symbols.
#
#   Rscript pancreas_gsea.R      (run from this directory)
#
# The topic model's columns are matched to the c = 1 reference fit with
# the Hungarian/Spearman procedure of factor_matching.R, so factor labels
# agree with the paper's figures (c = 1: k10 = IL-1beta, k4 = IFN-gamma,
# k3/k13 = delta/gamma cells, ...).
#
# Writes pancreas_gsea.csv (one row per model x factor x enriched term).

source("factor_gsea.R")
source("../fig_scripts/factor_matching.R")
suppressMessages(library(org.Mm.eg.db))

load("../data/experiment_results.Rdata")
log1p_fit <- res_list$pancreas$`1`
tm_fit    <- res_list$pancreas$`Inf`

K  <- ncol(log1p_fit$FF)
F1 <- as.matrix(log1p_fit$FF)
colnames(F1) <- paste0("k", seq_len(K))
Ft <- as.matrix(tm_fit$F)
colnames(Ft) <- paste0("k", match_display_labels(F1, Ft))
Ft <- Ft[, paste0("k", seq_len(K))]

N_TOP <- 50

run_one <- function(FF, model) {
  message("=== ", model, " ===")
  d <- factor_ora(FF, n_top = N_TOP, id_type = "symbol",
                  orgdb = org.Mm.eg.db)
  if (nrow(d) > 0) d <- cbind(model = model, d)
  d
}

res <- rbind(run_one(F1, "c = 1"), run_one(Ft, "c = Inf"))
write.csv(res, "pancreas_gsea.csv", row.names = FALSE)
message("Wrote pancreas_gsea.csv (", nrow(res), " rows)")

for (m in unique(res$model)) for (k in paste0("k", seq_len(K))) {
  d <- res[res$model == m & res$factor == k, ]
  if (nrow(d) == 0) next
  cat(sprintf("\n%s  %s  (%d terms at FDR < 0.05)\n", m, k, nrow(d)))
  for (i in seq_len(min(3, nrow(d))))
    cat(sprintf("  %-52s p_adj=%.2g  [%s%s]\n",
                substr(d$term[i], 1, 52), d$p_adj[i],
                paste(strsplit(d$genes[i], ", ")[[1]][1:min(6, d$k[i])],
                      collapse = ", "),
                ifelse(d$k[i] > 6, ", ...", "")))
}
