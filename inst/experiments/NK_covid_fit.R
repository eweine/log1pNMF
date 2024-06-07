library(Seurat)
library(SeuratData)
library(SeuratDisk)

Convert("nk_covid.h5ad", dest = "h5seurat", overwrite = FALSE)


covid_dat <- LoadH5Seurat("nk_covid.h5seurat", meta.data = FALSE, misc = FALSE)

covid_md <- readr::read_csv(
  "~/Downloads/nk_md.csv"
)

counts <- covid_dat@assays$RNA@counts

#covid_dat <- FindVariableFeatures(covid_dat, nfeatures = 6000)

#counts_samp <- counts[VariableFeatures(covid_dat), ]
library(dplyr)
samples_df <- covid_md %>%
  dplyr::rename(
    severity = `CoVID-19 severity`,
    type = `Sample type`,
    time_sampled = `Sample time`
  ) %>%
  select(severity, type, sampleID, time_sampled)

library(passPCA)

counts <- Matrix::t(counts)

counts <- counts[, Matrix::colSums(counts) > 10]

fit <- fit_factor_model_log1p_quad_approx_sparse(
  Y = counts,
  K = 10,
  maxiter = 750,
  approx_range = c(0, 1.25),
  s = Matrix::rowSums(counts) / mean(Matrix::rowSums(counts))
)

readr::write_rds(fit, "~/Documents/passPCA/inst/experiments/results/fit_log1p_covid.rds")

library(fastTopics)

nmf_fit <- fit_poisson_nmf(
  X = counts, 
  k = 10,
  numiter = 750
)

readr::write_rds(nmf_fit, "~/Documents/passPCA/inst/experiments/results/fit_nmf_pois_covid.rds")


samples_df <- samples_df %>%
  dplyr::mutate(ts = glue::glue("{time_sampled} {severity}"))

samples_df <- samples_df %>%
  dplyr::mutate(
    ts = case_when(
      ts == "control control" ~ "control",
      TRUE ~ ts
    )
  )

fit$U <- t(t(fit$U) / apply(fit$U,2,max))

library(fastTopics)

structure_plot(
  fit$U, grouping = samples_df$ts
)

