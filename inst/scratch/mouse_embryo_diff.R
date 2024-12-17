library(dplyr)
dat <- readr::read_rds("~/Documents/data/lung_and_airway.rds")

cells <- dat@meta.data
genes <- dat@assays$RNA@meta.features

counts <- Matrix::t(dat@assays$RNA@counts)

genes <- genes %>% dplyr::filter(feature_type == "protein_coding")

counts <- counts[, colnames(counts) %in% rownames(genes)]

counts <- counts[,Matrix::colSums(counts) > 0]
counts <- counts[Matrix::rowSums(counts) > 0, ]

cells_sum <- cells %>%
  dplyr::group_by(author_day, author_cell_type) %>%
  dplyr::summarise(
    tot = n()
  )

counts <- as(counts, "CsparseMatrix")

cells_av2 <- cells %>%
  dplyr::filter(author_cell_type == "Alveolar Type 2 cells")

df_list <- list()

for (day in unique(cells_av2$author_day)) {

  print(day)
  cells_av2_day <- cells_av2 %>%
    dplyr::filter(author_day == day)

  counts_av2_day <- counts[rownames(cells_av2_day),,drop=FALSE ]
  gene_rates <- Matrix::colSums(counts_av2_day) / sum(counts_av2_day)

  df_list[[day]] <- data.frame(
    gene = names(gene_rates),
    expr_rate = unname(gene_rates),
    day = day,
    n_cells = nrow(counts_av2_day)
  )

}

df <- dplyr::bind_rows(df_list)

df$expr_rate <- df$expr_rate * mean(Matrix::rowSums(counts))

# now, the question is how should I select which set of genes to look at?
df <- df %>%
  dplyr::filter(day %in% c("E1875", "P0000"))

df <- df %>%
  dplyr::mutate(sqrt_expr_rate = 2 * sqrt((3/8) + expr_rate))


start_df <- df %>%
  dplyr::filter(day == "E1875") %>%
  dplyr::select(-c(day, n_cells)) %>%
  dplyr::rename(expr_rate_start = expr_rate, sqrt_expr_rate_start = sqrt_expr_rate)

end_df <- df %>%
  dplyr::filter(day == "P0000") %>%
  dplyr::select(-c(day, n_cells)) %>%
  dplyr::rename(expr_rate_end = expr_rate, sqrt_expr_rate_end = sqrt_expr_rate)


bookend_df <- start_df %>%
  dplyr::inner_join(end_df, by = "gene")

bookend_df <- bookend_df %>%
  dplyr::mutate(
    abs_diff = abs(sqrt_expr_rate_end - sqrt_expr_rate_start),
    diff = sqrt_expr_rate_end - sqrt_expr_rate_start
  )

cc <- 0.1

bookend_df <- bookend_df %>%
  dplyr::mutate(
    base_rate = case_when(
      expr_rate_end > expr_rate_start ~ expr_rate_start,
      TRUE ~ expr_rate_end
    ),
    abs_diff_rate = abs(expr_rate_start - expr_rate_end),
    log1p_diff_rate = abs(log1p(expr_rate_start/cc) - log1p(expr_rate_end/cc))
  )


library(ggplot2)

ggplot(data = bookend_df) +
  geom_point(aes(x = log1p(base_rate/cc), y = log1p_diff_rate)) +
  geom_smooth(aes(x = log1p(base_rate/cc), y = log1p_diff_rate))

# it may also be useful to look at positive and negative signals here


