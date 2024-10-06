# First, I'm interested in comparing the fits from the different model
# in terms of lambda

load("~/Documents/data/fastglmpca/raw_data/droplet.RData")

fit_c1 <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/log1p_quad_approx_c1_droplets_7_factors_2500_iter.rds")
fit_c2 <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/log1p_quad_approx_c2_droplets_7_factors_2500_iter.rds")
fit_c4 <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/log1p_quad_approx_c4_droplets_7_factors_2500_iter.rds")
fit_c8 <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/log1p_quad_approx_c8_droplets_7_factors_2500_iter.rds")
fit_c16 <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/log1p_quad_approx_c16_droplets_7_factors_2500_iter.rds")


s <- Matrix::rowSums(counts) / mean(Matrix::rowSums(counts))

Lambda_c1 <- Matrix::Diagonal(x=s) %*% (exp(tcrossprod(fit_c1$U, fit_c1$V)) - 1)
Lambda_c2 <- Matrix::Diagonal(x=2 * s) %*% (exp(tcrossprod(fit_c2$U, fit_c2$V)) - 1)

lc1 <- as.vector(as.matrix(Lambda_c1))
lc2 <- as.vector(as.matrix(Lambda_c2))

lambda_df <- data.frame(
  lc1 = lc1,
  lc2 = lc2
)

l_samp <- lambda_df %>%
  sample_frac(.005)

library(ggplot2)

# ggplot(data = l_samp) +
#   geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
#   geom_point(aes(x = lc1, y = lc2), alpha = 0.25) +
#   cowplot::theme_cowplot() +
#   xlab("UV' c = 1") +
#   ylab("UV' c = 2")

# Now, it would be nice to make structure plots for each of these fits
# this will give us a sense of different values of c
library(fastTopics)

# Now, I think it would be interesting to examine different
# structure plots

fit_c1$U <- as.matrix(Matrix::colScale(fit_c1$U, 1/apply(fit_c1$U,2,max)))
fit_c2$U <- as.matrix(Matrix::colScale(fit_c2$U, 1/apply(fit_c2$U,2,max)))

fit_c4$U <- as.matrix(Matrix::colScale(fit_c4$U, 1/apply(fit_c4$U,2,max)))
fit_c8$U <- as.matrix(Matrix::colScale(fit_c8$U, 1/apply(fit_c8$U,2,max)))

fit_c16$U <- as.matrix(Matrix::colScale(fit_c16$U, 1/apply(fit_c16$U,2,max)))



p1 <- structure_plot(fit_c1$U, grouping = samples$tissue, gap = 75) + ggtitle("log1p c = 1")
p2 <- structure_plot(fit_c2$U, grouping = samples$tissue, gap = 75) + ggtitle("log1p c = 2")

p3 <- structure_plot(fit_c4$U, grouping = samples$tissue, gap = 75) + ggtitle("log1p c = 4")
p4 <- structure_plot(fit_c8$U, grouping = samples$tissue, gap = 75) + ggtitle("log1p c = 8")

p5 <- structure_plot(fit_c16$U, grouping = samples$tissue, gap = 75) + ggtitle("log1p c = 16")

# now, want to look at plots for other models
# first, the actual NMF fit

nmf_pois_fit <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/pois_nmf_droplets_7_factors_2500_iter.rds")
p6 <- structure_plot(nmf_pois_fit, grouping = samples$tissue, gap = 75) + ggtitle("Poisson NMF")

nmf_frob_fit <- readr::read_rds("~/Documents/passPCA/inst/experiments/results/log1p_frobenius_nmf_droplets.rds")

nmf_frob_fit$U <- nmf_frob_fit$w %*% diag(nmf_frob_fit$d)
nmf_frob_fit$U <- as.matrix(Matrix::colScale(nmf_frob_fit$U, 1/apply(nmf_frob_fit$U,2,max)))
p7 <- structure_plot(nmf_frob_fit$U, grouping = samples$tissue, gap = 75) + ggtitle("log1p Trans + Frobenius")

library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, p6, p7, nrow = 4, ncol = 2)
