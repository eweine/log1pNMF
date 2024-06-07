ids <- readr::read_csv("~/Downloads/GSE102827_cell_type_assignments.csv")
library(dplyr)
ids <- ids %>%
  dplyr::filter(!is.na(maintype))

counts <- readr::read_csv("~/Downloads/GSE102827_merged_all_raw.csv")
genes <- counts$...1
counts <- counts %>% dplyr::select(-`...1`)
counts <- counts %>% dplyr::select(ids$`...1`)
rm(ids)
gc()
counts <- as(counts, "matrix")
