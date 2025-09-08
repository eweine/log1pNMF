library(ggplot2)
library(fastTopics)
library(log1pNMF)
library(ggh4x)
library(ggpubr)
library(dplyr)
library(fastglmpca)

n_cells <- 1000

grouping <- c(
  rep("A1", 100), rep("A2", 100),
  rep("B1", 100), rep("B2", 100),
  rep("C1", 100), rep("C2", 100),
  rep("D1", 100), rep("D2", 100),
  rep("E1", 100), rep("E2", 100)
)

set.seed(1)
lambda_group1_c <- c(rexp(200), rep(0, 800))
lambda_group2_c <- c(rep(0, 200), rexp(200), rep(0, 600))
lambda_group3_c <- c(rep(0, 400), rexp(200), rep(0, 400))
lambda_group4_c <- c(rep(0, 600), rexp(200), rep(0, 200))
lambda_group5_c <- c(rep(0, 800), rexp(200))

mult_vec <- rep(2.5, 1000)
mult_vec[sample(1:1000, size = 800)] <- 0

lambda_group1_t <- if_else(
  mult_vec == 0,
  lambda_group1_c,
  exp(mult_vec + lambda_group1_c)
)

lambda_group2_t <- if_else(
  mult_vec == 0,
  lambda_group2_c,
  exp(mult_vec + lambda_group2_c)
)

lambda_group3_t <- if_else(
  mult_vec == 0,
  lambda_group3_c,
  exp(mult_vec + lambda_group3_c)
)

lambda_group4_t <- if_else(
  mult_vec == 0,
  lambda_group4_c,
  exp(mult_vec + lambda_group4_c)
)

lambda_group5_t <- if_else(
  mult_vec == 0,
  lambda_group5_c,
  exp(mult_vec + lambda_group5_c)
)

lambda_list <- list()
Y_list <- list()

lambda_list[[1]] <- lambda_group1_c
lambda_list[[2]] <- lambda_group1_t
lambda_list[[3]] <- lambda_group2_c
lambda_list[[4]] <- lambda_group2_t
lambda_list[[5]] <- lambda_group3_c
lambda_list[[6]] <- lambda_group3_t
lambda_list[[7]] <- lambda_group4_c
lambda_list[[8]] <- lambda_group4_t
lambda_list[[9]] <- lambda_group5_c
lambda_list[[10]] <- lambda_group5_t

set.seed(1)
n_cells_per_group <- 100
n_genes <- 1000

for (group in 1:10) {
  
  Lambda <- matrix(
    data = rep(lambda_list[[group]], n_cells_per_group),
    nrow = n_cells_per_group,
    ncol = n_genes,
    byrow = TRUE
  )
  
  Y <- matrix(
    data = rpois(n = n_cells_per_group * n_genes, lambda = as.vector(Lambda)),
    nrow = n_cells_per_group,
    ncol = n_genes
  )
  
  Y_list[[group]] <- Y
  
}

Y <- do.call(rbind, Y_list)

Y <- as(Y, "CsparseMatrix")
Y <- Y[,Matrix::colSums(Y) > 0]

library(fastTopics)

ft_r1 <- fastTopics:::fit_pnmf_rank1(Y)

init_LL <- cbind(
  ft_r1$L,
  matrix(data = 1e-5, nrow = nrow(Y), ncol = 5)
)

init_FF <- cbind(
  ft_r1$F,
  matrix(data = 1e-5, nrow = ncol(Y), ncol = 5)
)

ft_init <- init_poisson_nmf(X = Y, F = init_FF, L = init_LL)

ft_r1_init <- fit_poisson_nmf(
  X = Y, 
  fit0 = ft_init, 
  control = list(nc = 7),
  verbose = "none"
)

log1p_fit_list <- list()
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)

for (cc in cc_vec) {
  
  set.seed(1)
  log1p_fit_list[[as.character(cc)]] <- fit_poisson_log1p_nmf(
    Y = Y, K = 6, loglik = "exact", init_method = "rank1",
    control = list(maxiter = 250, verbose = FALSE), cc = cc
  )
  
}

colnames(ft_r1_init$L) <- paste0("k", c(1, 2, 4, 3, 5, 6))
colnames(ft_r1_init$F) <- paste0("k", c(1, 2, 4, 3, 5, 6))

ft_r1_init$L <- ft_r1_init$L[,paste0("k", 1:6)]
ft_r1_init$F <- ft_r1_init$F[,paste0("k", 1:6)]

sp1 <- structure_plot(
  log1pNMF:::normalize_bars(
    diag(1 / log1p_fit_list[[as.character(1)]]$s) %*% ft_r1_init$L
  ), 
  loadings_order = 1:n_cells, 
  grouping = grouping, 
  gap = 10
) +
  theme(
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 11),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(size = 14)
  ) +
  ylab("Membership") + 
  ggtitle("Topic Model Loadings") +
  guides(fill="none")

sp2 <- normalized_structure_plot(
    log1p_fit_list[[as.character(10)]], 
    loadings_order = 1:n_cells, 
    grouping = grouping, 
    gap = 10, 
    topics = rev(1:6)
  ) +
  theme(
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 11),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(size = 14)
    ) + 
  ylab("Membership") + 
  ggtitle("log1p Model Loadings (c = 10)") +
  guides(fill=guide_legend(title="Factor"))

sp3 <- normalized_structure_plot(
  log1p_fit_list[[as.character(0.001)]], 
  loadings_order = 1:n_cells, 
  grouping = grouping, 
  gap = 10, 
  topics = rev(1:6)
) +
  theme(
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 11),
    axis.title.y = element_text(size = 13),
    plot.title = element_text(size = 14)
  ) + 
  ylab("Membership") + 
  ggtitle("log1p Model Loadings (c = 1e-3)") +
  guides(
    fill=guide_legend(title="Factor")
    )

g_sp <- ggarrange(
  sp1, sp2, sp3, 
  nrow = 3, ncol = 1, 
  common.legend = TRUE, legend = "right",
  labels = c("B", "C", "D")
  )

hoyer_sparsity <- function(x) {
  
  n <- length(x)
  (1 / (sqrt(n) - 1)) * (sqrt(n) - (sum(x) / (sqrt(sum(x ^ 2)))))
  
}

for (cc in cc_vec) {
  
  log1p_fit_list[[as.character(cc)]]$l_sparsity <- apply(
    log1p_fit_list[[as.character(cc)]]$LL, 2, hoyer_sparsity
  )
  
  log1p_fit_list[[as.character(cc)]]$f_sparsity <- apply(
    log1p_fit_list[[as.character(cc)]]$FF, 2, hoyer_sparsity
  )
  
}

l_sparsity_vec <- unlist(lapply(log1p_fit_list, function(x) {mean(x$l_sparsity)}))

f_sparsity_vec <- unlist(lapply(log1p_fit_list, function(x) {mean(x$f_sparsity)}))

s <- Matrix::rowSums(Y)
s <- s / mean(s)

l_tm_sparsity <- mean(apply(diag(1/s) %*% ft_r1_init$L, 2, hoyer_sparsity))
f_tm_sparsity <- mean(apply(ft_r1_init$F, 2, hoyer_sparsity))

df_sparsity_l <- data.frame(
  cc = cc_vec,
  sparsity = l_sparsity_vec
)

g_l_sparsity <- ggplot(data = df_sparsity_l, aes(x = cc, y = sparsity)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-3, 1, 1e3), transform = "log10") +
  xlab("c") +
  ylab("Mean Loading Sparsity") +
  geom_hline(yintercept = l_tm_sparsity, color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.006, y=l_tm_sparsity + 0.05, label="Topic Model", color="red"
  )

df_sparsity_f <- data.frame(
  cc = cc_vec,
  sparsity = f_sparsity_vec
)

g_f_sparsity <- ggplot(data = df_sparsity_f, aes(x = cc, y = sparsity)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(1e-3, 1, 1e3), transform = "log10") +
  xlab("c") +
  ylab("Mean Factor Sparsity") +
  geom_hline(yintercept = f_tm_sparsity, color = "red", linetype = "dashed") +
  ggplot2::annotate(
    geom="text", x=0.006, y=f_tm_sparsity + 0.05, label="Topic Model", color="red"
  )

g_sparsity <- ggarrange(
  g_l_sparsity, g_f_sparsity, nrow = 1, ncol = 2,
  labels = c("E", "F")
  )

gene_order_df <- data.frame(
  gene_id = 1:1000,
  group = c(
    rep("A", 200),
    rep("B", 200),
    rep("C", 200),
    rep("D", 200),
    rep("E", 200)
  ),
  mult = mult_vec > 0
  ) %>%
  dplyr::arrange(group, desc(mult_vec))

gene_order_df$gene_order <- 1:1000
gene_order_df <- gene_order_df %>%
  dplyr::select(gene_id, gene_order)

expr_df <- data.frame(
  group = c(
    rep("A1", 1000), rep("A2", 1000),
    rep("B1", 1000), rep("B2", 1000),
    rep("C1", 1000), rep("C2", 1000),
    rep("D1", 1000), rep("D2", 1000),
    rep("E1", 1000), rep("E2", 1000)
  ),
  expr = unlist(lambda_list),
  gene_id = rep(1:1000, 10)
) %>%
  dplyr::inner_join(gene_order_df)

g_expr <- ggplot(
  expr_df %>% dplyr::filter(
    group %in% c("A1", "A2", "B1", "B2")
    ), 
    aes(x = gene_order, y = expr)
  ) +
  geom_col(width = 1, colour = NA, linewidth = 0) +  # Use bars to represent lambda values
  facet_manual(
    ~group, #scales = "free", 
    design = c("ABCD"),
    labeller = labeller(group = function(x) paste("Group", x))) +
  labs(x = "Gene Index", y = "True Expression") +
  scale_y_continuous(trans = "log1p", breaks = c(0, 10, 100, 2000)) + 
  cowplot::theme_cowplot() + # Use a minimal theme for a clean look
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0.5, size = 13)
  ) 

set.seed(1)
fgpca_fit <- fit_glmpca_pois(Y = Matrix::t(Y), K = 2, control = list(maxiter = 350))
gpca_df <- as.data.frame(fgpca_fit$V)
gpca_df$Group <- grouping
gpca <- ggplot(data = gpca_df, aes(x = k_1, y = k_2)) +
  geom_point(aes(color = Group)) +
  cowplot::theme_cowplot() +
  scale_color_manual(values = fastTopics:::kelly()[2:11]) +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("GLM-PCA Loadings") +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

g_gpca <- ggarrange(
  NULL, gpca, NULL,
  nrow = 1, ncol = 3,
  widths = c(0.33, 1, 0.33),
  labels = c("", "G", "")
)

g <- ggarrange(
  g_expr,
  g_sp,
  g_sparsity,
  g_gpca,
  nrow = 4,
  ncol = 1,
  heights = c(0.5, 1, 0.7, 0.7),
  labels = c("A", "", "", "")
)

ggsave(
  plot = g,
  device = "png",
  filename = "../images/loading_sim_fig.png",
  width = 8,
  height = 14
)

