reticulate::use_condaenv("nmf")

ad <- reticulate::import("anndata")
dat <- ad$read_h5ad("/Users/eweine/Downloads/healthy_human_4chamber_map_unnormalized_V3.h5ad")

tot <- dat$obs %>%
  dplyr::group_by(biological.individual) %>%
  dplyr::summarise(tot = dplyr::n())

d1_obs <- dat$obs %>%
  dplyr::filter(biological.individual == 1221)

counts <- dat$X
rownames(counts) <- rownames(dat$obs)
d1_counts <- counts[rownames(d1_obs), ]

colnames(d1_counts) <- rownames(dat$var)
d1_counts <- as(d1_counts, "CsparseMatrix")
d1_counts <- d1_counts[,Matrix::colSums(d1_counts) > 0]
readr::write_rds(d1_counts, "~/Documents/data/heart.rds")


#genes_to_use <- which(Matrix::colSums(d1_counts>0)>4)
#d1_counts <- d1_counts[,genes_to_use]
