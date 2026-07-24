# Generic runner for the exact-vs-approx timing experiment.
#
#   Rscript run_experiment.R <dataset> <scheme> <method>
#     dataset : pancreas | moca | bbc
#     scheme  : random | rank1
#     method  : approx | exact
#
# Sources the dataset's _common.R (configuration + data load + shared helpers),
# then fits the single (scheme, method) combination and writes its fitted object
# and tidy trace table. One combination per invocation so each is checkpointed
# independently: a job hitting the wall-clock limit loses only its own fit.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L)
  stop("usage: run_experiment.R <pancreas|moca|bbc> <random|rank1> <approx|exact>")
dataset <- args[1]; scheme <- args[2]; method <- args[3]

common <- switch(dataset,
  pancreas = "exact_vs_approx_common.R",
  moca     = "exact_vs_approx_moca_common.R",
  bbc      = "exact_vs_approx_bbc_common.R",
  stop("unknown dataset: ", dataset))

message(sprintf("=== %s / %s / %s ===", dataset, scheme, method))
source(common)
## make sure the output directory exists so the fit's final saveRDS cannot fail
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
run_job(scheme, method)
