library(ggpubr)
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

## the raw Rdata spells this cell type "Endothelial/Mesnchymal";
## fix it so the factor level below matches (otherwise those cells
## silently become NA and vanish from the plots)
barcodes$celltype[barcodes$celltype == "Endothelial/Mesnchymal"] <-
  "Endothelial/Mesenchymal"
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

## Factor matching for display (see factor_matching.R): the c = 1 fit is the
## reference, shown in the manually chosen semantic order below (identical to
## the main-text figure). Every other panel -- each c and the topic model --
## is matched to it by Hungarian assignment on Spearman correlations between
## the paired columns of F, so a factor keeps one color everywhere.
source("factor_matching.R")
ref_labels_pancreas <- c(11, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 12, 13)
ref_FF <- res_list$pancreas[["1"]]$FF

for (cc in cc_vec) {
  lab <- if (cc == 1) ref_labels_pancreas else
    match_display_labels(ref_FF, fit_list[[as.character(cc)]]$FF,
                         ref_labels_pancreas)
  fit_list[[as.character(cc)]]$LL <-
    relabel_and_sort(fit_list[[as.character(cc)]]$LL, lab)
  fit_list[[as.character(cc)]]$FF <-
    relabel_and_sort(fit_list[[as.character(cc)]]$FF, lab)
}

tm_labels <- match_display_labels(ref_FF, fit_list[[as.character(Inf)]]$F,
                                  ref_labels_pancreas)
fit_list[[as.character(Inf)]]$L <-
  relabel_and_sort(fit_list[[as.character(Inf)]]$L, tm_labels)
fit_list[[as.character(Inf)]]$F <-
  relabel_and_sort(fit_list[[as.character(Inf)]]$F, tm_labels)

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
