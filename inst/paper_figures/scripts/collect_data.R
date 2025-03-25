sim_quality_res_list <- list()

load("~/Documents/data/fit_list_sim_c1e-3.Rdata")
sim_quality_res_list[["c = 1e-3"]]$fit_list_approx <- fit_list_approx
sim_quality_res_list[["c = 1e-3"]]$fit_list_exact <- fit_list_exact

load("~/Documents/data/fit_list_sim_c1.Rdata")
sim_quality_res_list[["c = 1"]]$fit_list_approx <- fit_list_approx
sim_quality_res_list[["c = 1"]]$fit_list_exact <- fit_list_exact

load("~/Documents/data/fit_list_sim_tm.Rdata")
sim_quality_res_list[["tm"]]$fit_list_approx <- fit_list_approx
sim_quality_res_list[["tm"]]$fit_list_exact <- fit_list_exact

res_list <- list()
res_list$sim_quality <- sim_quality_res_list

bbc_fit_list <- list()

cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)
K <- 10

for (cc in cc_vec) {
  
  bbc_fit_list[[as.character(cc)]] <- readr::read_rds(
    glue::glue(
      "~/Documents/data/bbc_log1p_c{cc}_k{K}_exact_1000_iter.rds"
    )
  )
  
}

bbc_fit_list[["Inf"]] <- readr::read_rds(
  glue::glue("~/Documents/data/bbc_nmf_k{K}_exact_1000_iter.rds")
)

res_list$bbc <- bbc_fit_list

pancreas_fit_list <- list()
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)
K <- 9

for (cc in cc_vec) {
  
  pancreas_fit_list[[as.character(cc)]] <- readr::read_rds(
    glue::glue(
      "~/Documents/data/pancreas_cellseq2_log1p_c{cc}_k{K}_exact_250_iter.rds"
    )
  )

}

pancreas_fit_list[["Inf"]] <- readr::read_rds(
  glue::glue("~/Documents/data/passPCA/pancreas/pancreas_pois_nmf_k{K}_exact_250_iter.rds")
)

res_list$pancreas <- pancreas_fit_list

save(res_list, file = "../data/experiment_results.Rdata")


