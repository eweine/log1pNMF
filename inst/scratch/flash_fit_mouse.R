load("~/Downloads/mouse_brain_stim.Rdata")

library(dplyr)

cells <- cells %>% dplyr::filter(!is.na(maintype))


counts <- counts[rownames(counts) %in% cells$...1, ]
#counts <- counts[, Matrix::colSums(counts) > 0]
counts <- counts[Matrix::rowSums(counts) > 0, ]

counts <- counts[sample(x = rownames(counts), size = 7500, replace = FALSE), ]
counts <- counts[, Matrix::colSums(counts) > 0]
counts <- as(counts, "CsparseMatrix")

ff <- passPCA::run_flash_log1p_with_greedy_init(
  Y = counts,
  var_type = 2,
  greedy_Kmax = 20
)

out <- list(
  LL = ff$L_pm,
  FF = ff$F_pm
)

readr::write_rds(out, "mouse_brain_flash_fit_out.rds")
