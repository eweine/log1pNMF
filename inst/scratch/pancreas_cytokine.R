library(Matrix)
library(log1pNMF)

load("pancreas_cytokine.RData")

set.seed(1)

i       <- which(samples$mouse == "S1")
samples <- samples[i,]
counts  <- counts[i,]

outliers <- c("TTTGTTGTCGTTAGTG-1","TTTGTTGGTAGAGCTG-1")
i        <- which(!is.element(samples$barcode,outliers))
samples  <- samples[i,]
counts   <- counts[i,]

j      <- which(colSums(counts > 0) > 4)
genes  <- genes[j,]
counts <- counts[,j]

# cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)
# 
# K_vec <- c(7, 10)
# 
# init_method_vec <- c("rank1", "random")
# 
# for (init_method in init_method_vec) {
#   
#   for (K in K_vec) {
#     
#     for (cc in cc_vec) {
#       
#       set.seed(1)
#       fit <- fit_poisson_log1p_nmf(
#         Y = counts, 
#         K = K,
#         cc = cc,
#         init_method = init_method,
#         loglik = "exact",
#         control = list(maxiter = 250)
#       )
#       
#       readr::write_rds(
#         fit,
#         glue::glue("pancreas_cytokine_log1p_c{cc}_{init_method}_init_K{K}.rds")
#       )
#       
#     }
#     
#   }
#   
# }

md <- readr::read_rds("~/Downloads/panc_ctyo_S1_celltypes.rds")

fit1 <- readr::read_rds("~/Downloads/pancreas_cytokine/pancreas_cytokine_log1p_c1_rank1_init_K7.rds")
md <- md %>%
  dplyr::filter(barcode %in% rownames(fit1$LL))
fit1$LL <- fit1$LL[md$barcode, ]

normalized_structure_plot(fit1, grouping = md$celltype_final, gap = 10)


fit2 <- readr::read_rds("~/Downloads/pancreas_cytokine/pancreas_cytokine_log1p_c1_random_init_K7.rds")
fit2$LL <- fit2$LL[md$barcode, ]
normalized_structure_plot(fit2, grouping = md$celltype_final, gap = 10)


