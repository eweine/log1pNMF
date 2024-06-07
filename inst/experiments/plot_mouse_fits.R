nmf12 <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/pois_nmf_k12_mouse_light.rds")

#nmf15 <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/pois_nmf_k15_mouse_light.rds")

counts <- readr::read_rds("~/Downloads/counts.rds")

library(fastTopics)

ids <- readr::read_csv("~/Downloads/GSE102827_cell_type_assignments.csv")

library(dplyr)

ids <- ids %>%
  dplyr::rename(cell_bc = `...1`)

counts_ids <- data.frame(
  row_idx = 1:nrow(counts),
  cell_bc = rownames(counts)
)

ids <- ids %>%
  dplyr::inner_join(counts_ids, by = "cell_bc")

ids <- ids %>% arrange(row_idx)

ids <- ids %>%
  dplyr::mutate(
    type_stim = glue::glue("{maintype}_{stim}"))

g1 <- structure_plot(
  nmf12,
  grouping = ids$type_stim,
  n = 15000,
  gap = 300
) + ggplot2::ggtitle("Poisson NMF")

log1p_mod <- readr::read_rds(
  "../passPCA/inst/experiments/results/log1p_mouse_light_K12.rds"
)

rownames(log1p_mod$V) <- colnames(counts)

U_normalized <- apply(log1p_mod$U, 2, function(col) col / max(col))

g2 <- structure_plot(
  U_normalized,
  grouping = ids$type_stim,
  n = 15000,
  gap = 300
) + ggplot2::ggtitle("log1p")

# Now, it would be interesting to just look at the excitatory cells

ids_excit <- ids %>%
  dplyr::filter(maintype == "Excitatory")

nmf12_e <- nmf12
nmf12_e$L <- nmf12$L[ids_excit$cell_bc, ]
nmf12_e$Ln <- nmf12$Ln[ids_excit$cell_bc, ]
nmf12_e$Ly <- nmf12$Ly[ids_excit$cell_bc, ]

ids_excit <- ids_excit %>% dplyr::mutate(stim_sub = glue::glue("{celltype}_{stim}"))

g3 <- structure_plot(
  nmf12_e,
  grouping = ids_excit$stim,
  n = nrow(ids_excit),
  gap = 300
) + ggplot2::ggtitle("Poisson NMF")

rownames(U_normalized) <- ids$cell_bc

g4 <- structure_plot(
  U_normalized[ids_excit$cell_bc, ],
  grouping = ids_excit$stim,
  n = nrow(ids_excit),
  gap = 300
) + ggplot2::ggtitle("log1p")

ids_excit_0h <- ids_excit %>% dplyr::filter(stim == "0h")
ids_excit_1h <- ids_excit %>% dplyr::filter(stim == "1h")
ids_excit_4h <- ids_excit %>% dplyr::filter(stim == "4h")

U0 <- U_normalized[ids_excit_0h$cell_bc, ]
U1 <- U_normalized[ids_excit_1h$cell_bc, ]
U4 <- U_normalized[ids_excit_4h$cell_bc, ]

cm_U0 <- colMeans(U0)
cm_U1 <- colMeans(U1)
cm_U4 <- colMeans(U4)
