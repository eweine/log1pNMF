library(dplyr)
library(fastTopics)
library(log1pNMF)
library(ggplot2)
library(cowplot)
library(Matrix)
library(readr)

set.seed(1)
K <- 13

load("../data/raw_data/pancreas_cytokine_lsa.Rdata")
load("../data/experiment_results.Rdata")
barcodes   <- as.data.frame(barcodes)


barcodes <- barcodes %>%
  dplyr::mutate(
    condition = if_else(
      condition == "IL-1B_IFNg",
      "IL-1B + IFNg",
      condition
    )
  )

clusters   <- factor(barcodes$celltype,
                     c("Acinar","Ductal","Endothelial/Mesenchymal","Macrophage",
                       "Alpha","Beta","Delta","Gamma"))

conditions <- factor(barcodes$condition,
                     c("Untreated","IL-1B","IFNg","IL-1B + IFNg"))

i <- c(sample(which(clusters == "Beta"),900),
       which(clusters != "Beta"))

fit_list <- list()
cc_vec <- c(1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3)

for (cc in cc_vec) {
  
  fit_list[[as.character(cc)]] <- res_list$pancreas[[as.character(cc)]]
  
  colnames(fit_list[[as.character(cc)]]$LL) <- paste0("k", 1:K)
  
}

fit_list[["Inf"]] <- res_list$pancreas$`Inf`

log1p_loadings_order <- c(11, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 12, 13)
tm_loadings_order <- c(11,12, 7, 6, 5, 1, 9, 4, 3, 2, 10, 13, 8)


colnames(fit_list[[as.character(1)]]$LL) <- paste0(
  "k", log1p_loadings_order
)

fit_list[[as.character(1)]]$LL <- fit_list[[as.character(1)]]$LL[,paste0("k", 1:13)]

colnames(fit_list[[as.character(1)]]$FF) <- paste0(
  "k", log1p_loadings_order
)

fit_list[[as.character(1)]]$FF <- fit_list[[as.character(1)]]$FF[,paste0("k", 1:13)]



colnames(fit_list[[as.character(Inf)]]$L) <- paste0(
  "k", tm_loadings_order
)

fit_list[[as.character(Inf)]]$L <- fit_list[[as.character(Inf)]]$L[,paste0("k", 1:13)]

colnames(fit_list[[as.character(Inf)]]$F) <- paste0(
  "k", tm_loadings_order
)

fit_list[[as.character(Inf)]]$F <- fit_list[[as.character(Inf)]]$F[,paste0("k", 1:13)]

plot_list <- list()

for (cc in cc_vec) {
  
  set.seed(1)
  plot_list[[as.character(cc)]] <- structure_plot(
    log1pNMF:::normalize_bars(fit_list[[as.character(cc)]]$LL[i,]), 
    grouping = clusters[i],gap = 25,perplexity = 70,n = Inf, font.size = 14
  )  + ggtitle(glue::glue("log1p Model Loadings (c = {cc})")) +
    theme(
      plot.title = element_text(size = 12),
      axis.title.y = element_text(size = 11),
      axis.text.x = element_text(size = 8)
    ) + ylab("Membership") +
    guides(fill=guide_legend(title="Factor", ncol=1)) +
    guides(colour = "none")
  
}

L <- log1pNMF:::normalize_bars(diag(1 / fit_list[[as.character(1)]]$s) %*% fit_list[[as.character(Inf)]]$L)

set.seed(1)
plot_list[[as.character(Inf)]] <- structure_plot(
  L[i,], 
  grouping = clusters[i],gap = 25,perplexity = 70,n = Inf, font.size = 12
)  + ggtitle("log1p Model Loadings (c = \u221E)") +
  theme(
    plot.title = element_text(size = 12),
    axis.title.y = element_text(size = 11),
    axis.text.x = element_text(size = 8)
  ) + ylab("Membership") +
  guides(fill=guide_legend(title="Factor", ncol=1)) +
  guides(colour = "none")

# 2. create a single blank plot
blank <- ggplot() + theme_void()

# 3. interleave: real, blank, real, blank, … but don’t add a blank after the last
combined <- vector("list", length(plot_list)*2 - 1)
combined[seq(1, length(combined), by = 2)] <- plot_list
combined[seq(2, length(combined)-1, by = 2)] <- list(blank)

# 4. set relative heights: 1 unit for each real plot, 0.1 for each blank
heights <- rep(c(1, -0.085), length(plot_list))
#    → this gives c(1,0.1, 1,0.1, …) of total length 15

# 5. arrange
g <- ggarrange(
  plotlist     = combined,
  ncol         = 1,
  nrow         = length(combined),
  heights      = heights,
  common.legend = TRUE,
  legend       = "right",
  align        = "v"
)


ggsave(
  "../images/supp_lsa_structure.png",
  g,
  device = "png",
  width = 10,
  height = 17
)
