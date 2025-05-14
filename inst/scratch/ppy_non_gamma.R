load("../data/raw_data/pancreas.RData")

i           <- which(sample_info$tech == "celseq2")
sample_info <- sample_info[i,]
counts      <- counts[i,]

s <- Matrix::rowSums(counts)

ppy_expr <- counts[,"PPY"]
xist_df <- data.frame(d = substr(sample_info$id, 1, 3), xist = counts[,"XIST"])

xist_df <- xist_df %>%
  dplyr::group_by(d) %>%
  dplyr::summarise(xist = mean(log1p(xist)))

r_df <- data.frame(d = substr(sample_info$id, 1, 3), r = counts[,"RPS4Y1"])

r_df <- r_df %>%
  dplyr::group_by(d) %>%
  dplyr::summarise(r = mean(log1p(r)))

ppy_df <- data.frame(
  ppy = ppy_expr,
  celltype = sample_info$celltype
)

n_cts <- ppy_df %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(tot = dplyr::n())

n_cts <- n_cts %>%
  dplyr::filter(tot >= 25)

# ppy_df <- ppy_df %>%
#   dplyr::filter(celltype %in% n_cts$celltype)

counts <- counts[, Matrix::colSums(counts) > 0]
genes_to_use <- which(Matrix::colSums(counts>0)>9)
counts <- counts[,genes_to_use]

gene_vec <- character((ncol(counts) - 1) * (length(unique(n_cts$celltype)) - 1))
ct_vec <- character(length(gene_vec))
rank_cor_vec <- numeric(length(gene_vec))
p_vec <- numeric(length(gene_vec))
i <- 1
j <- 1

for (gene in colnames(counts)) {

  if (gene != "PPY") {
    j <- j + 1
    print(j)

    ppy_df$other_gene <- counts[, gene]

    for (ct in unique(n_cts$celltype)) {

      if (ct != "gamma")

      ct_df <- ppy_df %>%
        dplyr::filter(celltype == ct)

      out <- cor.test(ct_df$ppy / s, ct_df$other_gene / s, method = "spearman")

      rank_cor <- out$estimate

      gene_vec[i] <- gene
      ct_vec[i] <- ct
      rank_cor_vec[i] <- rank_cor
      p_vec[i] <- out$p.value
      i <- i + 1

    }

  }
}

cor_df <- data.frame(
  gene = gene_vec,
  celltype = ct_vec,
  cor = rank_cor_vec,
  p = p_vec
)

cor_df$q <- p.adjust(cor_df$p, method = "BH")

cor_df <- cor_df %>%
  dplyr::filter(q < 0.05)

ggplot(ppy_df, aes(x = log1p(ppy))) +
  geom_histogram(colour = "black", fill = "steelblue") +
  facet_wrap(~ celltype, scales = "free_y") +
  theme_minimal() +
  labs(
    x = "",
    y = "log1p(Count)",
    title = "Distribution of PPY by Cell Type"
  ) +
  cowplot::theme_cowplot()

ggplot(ppy_df, aes(x = ppy)) +
  geom_histogram(colour = "black", fill = "steelblue", bins = 100) +
  facet_wrap(~ celltype, scales = "free") +
  theme_minimal() +
  labs(
    x = "",
    y = "Count",
    title = "Distribution of PPY by Cell Type"
  ) +
  cowplot::theme_cowplot()


i_alpha <- which(sample_info$celltype == "acinar")
counts_alpha <- counts[i_alpha, ]
plot(log1p(counts_alpha[,"RPS4Y1"]), log1p(counts_alpha[,"PPY"]))
