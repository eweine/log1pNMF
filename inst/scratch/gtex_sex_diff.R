df <- read.delim(file="~/Downloads/gene_reads_2017-06-05_v8_artery_coronary.gct", skip = 2)

counts <- df %>%
  dplyr::select(-c(id, Name, Description)) %>%
  as.matrix()

rownames(counts) <- df$Description

# filter out very low expression genes
counts <- counts[Matrix::rowSums(counts) >= 100, ]

md <- read.delim("~/Downloads/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

col_obs_ids <- sub("\\..*", "", sub("\\.", "-", colnames(counts)))

md <- md %>% dplyr::filter(SUBJID %in% col_obs_ids)
colnames(counts) <- col_obs_ids

md$gender <- as.factor(ifelse(md$SEX == 1, "M", "F"))

counts <- counts[,md$SUBJID]
tot <- Matrix::colSums(counts)
md$tot <- tot
b0_vec <- c()
b1_vec <- c()
b1_se_vec <- c()
gene_vec <- c()

it <- 1

for (gene in rownames(counts)) {

  print(round(it / nrow(counts), 2))

  md$expr <- counts[gene, ]

  # Wrap the model fitting and extraction in tryCatch
  tryCatch({
    mod <- MASS::glm.nb(
      expr ~ gender + offset(log(tot)),
      data = md
    )

    b0_vec <- c(b0_vec, coef(summary(mod))["(Intercept)", "Estimate"])
    b1_vec <- c(b1_vec, coef(summary(mod))["genderM", "Estimate"])
    b1_se_vec <- c(b1_se_vec, coef(summary(mod))["genderM", "Std. Error"])
    gene_vec <- c(gene_vec, gene)

  }, error = function(e) {
    message(paste("Error with gene:", gene, "- Skipping this gene"))
  })

  it <- it + 1
}

out_df <- data.frame(
  gene = gene_vec,
  b1 = b1_vec,
  b0 = b0_vec,
  b1_se = b1_se_vec
)

library(ashr)
a_out <- ash(betahat = out_df$b1, sebetahat = out_df$b1_se)

out_df$b1_pm <- ashr::get_pm(a_out)

off_m <- mean(tot)

out_df <- out_df %>%
  dplyr::mutate(
    pred_diff = off_m * exp(b0) - off_m * exp(b0 + b1_pm)
  )

out_df <- out_df %>%
  dplyr::mutate(
    baseline = if_else(
      pred_diff > 0, off_m * exp(b0), off_m * exp(b0 + b1)
    )
  )

out_df$lfsr <- ashr::get_lfsr(a_out)

readr::write_rds(out_df, "~/Documents/data/passPCA/experiment_results/gtex_artery_sex_diff.rds")



# I think the two plots below are useful
plot(out_df$baseline, abs(out_df$pred_diff))
plot(log(out_df$baseline), abs(out_df$b1_pm))
