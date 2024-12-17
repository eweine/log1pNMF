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
  dplyr::filter(day != "E1625")

df <- df %>% 
  dplyr::mutate(sqrt_expr_rate = 2 * sqrt((3/8) + expr_rate))

# I think first I would like to look at big changes from start to finish

start_df <- df %>%
  dplyr::filter(day == "E1675") %>%
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

mover_df <- bookend_df %>% dplyr::filter(abs_diff > 0.5)

mover_genes <- mover_df$gene

# I think it might be useful to plot expression over time for these genes

day_str <- substr(df$day, 2, 5)
day_str <- ifelse(day_str == "0000", "1900", day_str)
day <- as.numeric(day_str)
df$day <- day

df_g <- df %>%
  dplyr::filter(gene == "ENSMUSG00000048583")

# a lot of genes may look like this...
plot(df_g$day, df_g$expr_rate)

# for each gene, I'm interested in the rank correlation between
# the difference and the point in time

cor_vec <- c()
gene_vec <- c()

for (g in unique(mover_df$gene)) {
  
  df_g <- df %>%
    dplyr::filter(gene == g) %>%
    dplyr::arrange(day)
  
  cor_g <- cor(x = df_g$day[-1], y = diff(df_g$expr_rate), method = "spearman")
  
  cor_vec <- c(cor_vec, cor_g)
  gene_vec <- c(gene_vec, g)
  
}

cor_df <- data.frame(
  gene = gene_vec,
  cor = cor_vec
)



df_g <- df %>%
  dplyr::filter(gene == "ENSMUSG00000039286")

# a lot of genes may look like this...
plot(df_g$day, df_g$expr_rate)


# it would also be useful to look at differences between before and after




