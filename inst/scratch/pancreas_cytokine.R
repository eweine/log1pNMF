library(Matrix)
load("~/Downloads/pancreas_cytokine.RData")

i       <- which(samples$mouse == "S1")
samples <- samples[i,]
counts  <- counts[i,]

outliers <- c("TTTGTTGTCGTTAGTG-1","TTTGTTGGTAGAGCTG-1")
i        <- which(!is.element(samples$barcode,outliers))
samples  <- samples[i,]
counts   <- counts[i,]

j      <- which(colSums(counts > 0) > 2)
genes  <- genes[j,]
counts <- counts[,j]

celltypes <- readr::read_csv("~/Downloads/pancreas_cytokine_S1_celltypes.csv")
celltypes$barcode <- gsub("\\.", "-", celltypes$barcode)

tm_fit <- readr::read_rds("~/Downloads/pancreas_cytokine_tm_k12.rds")

log1p_random <- readr::read_rds("~/Downloads/pancreas_cytokine/pancreas_cytokine_log1p_c1_random_init_K12.rds")

library(log1pNMF)

normalized_structure_plot(
  log1p_random, grouping = celltypes$celltype, n = Inf, gap = 20
  )

log1p_rank1 <- readr::read_rds("~/Downloads/pancreas_cytokine/pancreas_cytokine_log1p_c1_rank1_init_K12.rds")

normalized_structure_plot(
  log1p_rank1, grouping = celltypes$celltype, n = Inf, gap = 20
)

rownames(log1p_rank1$FF) <- genes$symbol

library(fastTopics)
set.seed(1)
#tm <- fit_poisson_nmf(counts,k = 12, control = list(numiter = 100,nc = 7))

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

celltypes <- readr::read_csv("~/Downloads/pancreas_cytokine_S1_celltypes.csv")
celltypes$barcode <- gsub("\\.", "-", celltypes$barcode)

plot_list <- list()

cc_vec <- c(1)

for (cc in cc_vec) {
  
  print(cc)
  
  fit <- readr::read_rds(
    glue::glue(
      "~/Downloads/pancreas_cytokine/pancreas_cytokine_log1p_c{cc}_rank1_init_K10.rds"
    )
  )
  
  plot_list[[as.character(cc)]] <- normalized_structure_plot(
    fit, grouping = celltypes$celltype, n = Inf, gap = 20
  )$plot
  
}

plot_list$`1000`
