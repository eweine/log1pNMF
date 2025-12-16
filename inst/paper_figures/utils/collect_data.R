sim_quality_res_list <- list()

load("~/Documents/data/fit_list_sim_c1e-3.Rdata")
sim_quality_res_list[["c = 1e-3"]]$fit_list_approx <- fit_list_approx
sim_quality_res_list[["c = 1e-3"]]$fit_list_approx_cheb <- fit_list_approx_cheb
sim_quality_res_list[["c = 1e-3"]]$fit_list_exact <- fit_list_exact

load("~/Documents/data/fit_list_sim_c1.Rdata")
sim_quality_res_list[["c = 1"]]$fit_list_approx <- fit_list_approx
sim_quality_res_list[["c = 1"]]$fit_list_approx_cheb <- fit_list_approx_cheb
sim_quality_res_list[["c = 1"]]$fit_list_exact <- fit_list_exact

load("~/Documents/data/fit_list_sim_tm.Rdata")
sim_quality_res_list[["tm"]]$fit_list_approx <- fit_list_approx
sim_quality_res_list[["tm"]]$fit_list_approx_cheb <- fit_list_approx_cheb
sim_quality_res_list[["tm"]]$fit_list_exact <- fit_list_exact

res_list <- list()
res_list$sim_quality <- sim_quality_res_list

res_list$bbc <- readr::read_rds("~/Documents/data/bbc_fit_list.rds")

pancreas_fit_list <- list()
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)
K <- 13

for (cc in cc_vec) {
  
  pancreas_fit_list[[as.character(cc)]] <- readr::read_rds(
    glue::glue(
      "~/Documents/data/passPCA/log1p_k13_c_{cc}_lsa.rds"
    )
  )

}

pancreas_fit_list[["Inf"]] <- readr::read_rds(
  glue::glue("~/Documents/single-cell-jamboree/output/panc_cyto_lsa_res/stancill_lsa_k13_r1_init_250_iter.rds")
)

pancreas_fit_list[["c = 1, cheby"]] <- readr::read_rds(
  "~/Documents/data/passPCA/lsa_k13_cheby_approx.rds"
)

pancreas_fit_list[["c = 1, frob"]] <- readr::read_rds(
  "~/Documents/data/passPCA/lsa_k13_frob_fit.rds"
)

pancreas_fit_list[["c = 1, frob, k = 12"]] <- readr::read_rds(
  "~/Documents/data/passPCA/lsa_k12_frob_fit.rds"
)

res_list$pancreas <- pancreas_fit_list

res_list$mcf7 <- readr::read_rds("~/Documents/data/mcf7_fit_list.rds")

save(res_list, file = "../data/experiment_results.Rdata")


