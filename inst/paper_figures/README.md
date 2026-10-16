# Manuscript figures

This directory contains the code used to run the experiments and produce the
figures in "A New Family of Poisson Non-negative Matrix Factorization Methods
Using the Shifted Log Link". Figure numbers below refer to the current
revision of the manuscript.

Layout:

- `experiments/` — code to run the experiments (model fitting, simulations,
  and post-processing), plus the small result files the figure scripts read.
  Scripts with `.sbatch`/`.sh` extensions run the heavier experiments on a
  Slurm cluster.
- `fig_scripts/` — code to make the plots. Each script writes its figure to
  `images/`. Run from inside `fig_scripts/` (paths are relative).
- `images/` — the rendered figures, as included in the manuscript.
- `data/` — fitted-model results (`experiment_results.Rdata`, assembled by
  `utils/collect_data.R`) and, in `data/raw_data/`, the raw datasets.
- `utils/` — `collect_data.R`, which gathers the real-data fits into
  `data/experiment_results.Rdata`.

## Figures

For each figure: (a) the code to run the underlying experiment, and (b) the
code to make the plot. When a single self-contained script does both, one
file is listed.

### Main text

| Figure | (a) Experiment code | (b) Plot code |
|---|---|---|
| Fig 1 (link function) | — | `fig_scripts/link_fig.R` |
| Fig 2 (computational scaling) | — | `fig_scripts/computational_complexity_fig.R` |
| Fig 3 (approximation accuracy) | approx-sim pipeline (below) | `fig_scripts/approximation_quality_fig.R` |
| Fig 4 (toy example) | `fig_scripts/factor_sim_identifiable_reduc_fig.R` (fits the models in-script) | same file |
| Fig 5 (sparsity / log-likelihood vs. c) | approx-sim pipeline (below) | `fig_scripts/sim2_main_fig.R` |
| Fig 6 (MCF-7) | `experiments/fit_mcf7_models.R`; GO terms: `experiments/mcf7_gsea.R` | `fig_scripts/mcf7_fig.R` |
| Figs 7, 9 (pancreas cell-type / treatment factors) | `experiments/fit_lsa_models.R`; GO terms: `experiments/pancreas_gsea.R` | `fig_scripts/lsa_fig.R` |
| Fig 8 (delta vs. gamma cells) | `experiments/fit_lsa_models.R` | `fig_scripts/lsa_delta_gamma_fig.R` |
| Fig 10 (BBC) | `experiments/fit_bbc_models.R` | `fig_scripts/bbc_fig.R` |
| Fig 11 (fit summaries across c) | `experiments/fit_mcf7_models.R`, `fit_lsa_models.R`, `fit_bbc_models.R` | `fig_scripts/fit_summary_fig.R` |

### Appendices

| Figure | (a) Experiment code | (b) Plot code |
|---|---|---|
| Fig F.1 (exact vs. approximate on real data) | `experiments/run_experiment.R` with `exact_vs_approx_common.R`, `exact_vs_approx_helpers.R`, `exact_vs_approx_bbc_common.R`, `exact_vs_approx_moca_common.R`; submitted via `run_pancreas.sbatch`, `run_bbc.sbatch`, `run_moca.sbatch` | `fig_scripts/exact_vs_approx_all_datasets_fig.R` |
| Figs G.1, G.2 (exact vs. approximations, pancreas) | `experiments/fit_lsa_approximate_supp.R` | `fig_scripts/approx_comp_supp_fig.R` |
| Fig K.1 (iterations to convergence) | approx-sim pipeline (below) | `fig_scripts/approx_sim_iterations_fig.R` |
| Fig K.2 (BBC initialization stability) | `experiments/run_bbc_init.R` with `bbc_init_common.R`, submitted via `run_bbc_init.sbatch` (QC: `check_bbc_init.R`) | `fig_scripts/bbc_init_stability_fig.R` |
| Fig K.3 (RMSE / factor correlation vs. c) | approx-sim pipeline (below) | `fig_scripts/sim2_rel_rmse_fig.R` |
| Fig K.4 (MCF-7 structure plots across c) | `experiments/fit_mcf7_models.R` | `fig_scripts/supp_mcf7_fig.R` |
| Fig K.5 (pancreas structure plots across c) | `experiments/fit_lsa_models.R` | `fig_scripts/supp_lsa_fig.R` |
| Fig K.6 (BBC structure plots across c) | `experiments/fit_bbc_models.R` | `fig_scripts/supp_bbc_fig.R` |

### The approx-sim pipeline (Figs 3, 5, K.1, K.3)

A single simulation experiment underlies the approximation-accuracy and
misspecification figures. In order:

1. `experiments/calibrate_approx_sim.R` — calibrates the simulation
   parameters to the pancreas data (reads
   `pancreas_count_histogram.rds`, writes `approx_sim_calibration.rds`).
2. `experiments/run_approx_sim.R` (shared settings in
   `approx_sim_common.R`) — fits every (seed, c_true, c_fit, method)
   combination; submitted on the cluster via `submit_approx_sim.sh` /
   `run_approx_sim.sbatch`.
3. `experiments/analyze_approx_sim.R` — collects the fits into the flat
   table `approx_sim_metrics.rds`, which the four figure scripts read.
   (QC: `check_approx_sim.R`.)

### Shared helpers

- `fig_scripts/factor_matching.R` — matches factors across fits (Hungarian
  algorithm); sourced by most real-data figure scripts.
- `experiments/factor_gsea.R` — GO-enrichment helpers sourced by
  `mcf7_gsea.R` and `pancreas_gsea.R` (GO terms shown in the tables of
  Figs 6, 7, and 9).
