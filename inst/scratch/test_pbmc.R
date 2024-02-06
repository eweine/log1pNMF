library(fastglmpca)
library(Seurat)

set.seed(1)
seurat_obj <- CreateSeuratObject(counts = pbmc_facs$counts)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

cts <- pbmc_facs$counts[VariableFeatures(seurat_obj), ]
cts <- Matrix::t(cts)

tic()
mod1 <- passPCA::fit_factor_model_log1p_exact(
  cts,
  K = 5,
  maxiter = 100,
  s = Matrix::rowSums(cts) / mean(Matrix::rowSums(cts))
)
toc()

s_d <- Matrix::Diagonal(x = Matrix::rowSums(cts) / mean(Matrix::rowSums(cts)))

H <- tcrossprod(mod1$U, mod1$V)
Lambda <- as.matrix(s_d %*% (exp(H) - 1))

tic()
mod2 <- passPCA::fit_factor_model_log1p_quad_approx_sparse(
  cts,
  K = 5,
  approx_range = c(0, 1.25),
  maxiter = 100,
  s = Matrix::rowSums(cts) / mean(Matrix::rowSums(cts))
)
toc()

H2 <- tcrossprod(mod2$U, mod2$V)
Lambda2 <- as.matrix(s_d %*% (exp(H2) - 1))

mod3 <- passPCA::fit_factor_model_log1p_lin_approx_sparse(
  cts,
  K = 5,
  a = 1,
  maxiter = 100,
  s = Matrix::rowSums(cts) / mean(Matrix::rowSums(cts))
)

H3 <- tcrossprod(mod3$U, mod3$V)
Lambda3 <- as.matrix(s_d %*% (exp(H2) - 1))

h_df <- data.frame(
  h = as.vector(H),
  h2 = as.vector(H2),
  h3 = as.vector(H3)
)

library(ggplot2)

g1 <- ggplot(data = h_df) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_point(aes(x = h, y = h2), alpha = 0.25) +
  cowplot::theme_cowplot() +
  xlab("Exact MLE UV'") +
  ylab("Quadratic Approximation UV'")

g2 <- ggplot(data = h_df) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_point(aes(x = h, y = h3), alpha = 0.25) +
  cowplot::theme_cowplot() +
  xlab("Exact MLE UV'") +
  ylab("Linear Approximation UV'")
