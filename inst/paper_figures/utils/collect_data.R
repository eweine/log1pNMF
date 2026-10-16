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

pancreas_fit_list[["c = 1, frob, k = 11"]] <- readr::read_rds(
  "~/Documents/data/passPCA/lsa_k11_frob_fit.rds"
)


res_list$pancreas <- pancreas_fit_list

res_list$mcf7 <- readr::read_rds("~/Documents/data/mcf7_fit_list.rds")



# Cluster-experiment summary tables (small, read directly from the
# experiment outputs in ../experiments; see the README there):
#   - approx_sim_metrics: per-fit metrics of the approx-sim pipeline
#     (analyze_approx_sim.R) -- Figs 3, 5, K.1, K.3
#   - exact_vs_approx_traces: rank-1 trace tables of the exact-vs-approx
#     timing experiment (run_experiment.R) -- Fig F.1
#   - bbc_init_stability: summary rows of the BBC initialization
#     experiment (run_bbc_init.R; extracted from the per-task outputs
#     into bbc_init_metrics.rds) -- Fig K.2

res_list$approx_sim_metrics <- readRDS("../experiments/approx_sim_metrics.rds")

tags <- c(moca = "exact_vs_approx_moca_K25_c1",
          pancreas = "exact_vs_approx_pancreas_K13_c1",
          bbc = "exact_vs_approx_bbc_K10_c1")
res_list$exact_vs_approx_traces <- lapply(tags, function(tag)
  sapply(c("exact", "approx"), function(m)
    readRDS(file.path("../experiments",
                      sprintf("%s_rank1_%s_traces.rds", tag, m))),
    simplify = FALSE))

res_list$bbc_init_stability <- readRDS("../experiments/bbc_init_metrics.rds")

save(res_list, file = "../data/experiment_results.Rdata")
