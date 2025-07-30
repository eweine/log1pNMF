library(rsvd)
library(fastglmpca)
library(fastTopics)
library(flashier)
library(data.table)
library(GEOquery)
library(ggplot2)
library(cowplot)
library(log1pNMF)
library(ggpubr)

load("../data/experiment_results.Rdata")
topic_colors <- c("tomato","darkblue","dodgerblue")
set.seed(10)
counts <- fread("../data/raw_data/GSE152749_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                sep = "\t",header = TRUE,stringsAsFactors = FALSE)
class(counts) <- "data.frame"
rownames(counts) <- counts$GeneID
counts <- counts[,-1]
counts <- as.matrix(counts)
storage.mode(counts) <- "double"
counts <- t(counts)
ids <- rownames(counts)

genes <- fread("../data/raw_data/Human.GRCh38.p13.annot.tsv.gz",sep = "\t",
               header = TRUE,stringsAsFactors = FALSE)
class(genes) <- "data.frame"
genes <- genes[1:10]
genes <- transform(genes,
                   GeneType = factor(GeneType),
                   Status   = factor(Status))

geo <- getGEO(filename = "../data/raw_data/GSE152749_family.soft.gz")
samples <- data.frame(id = names(GSMList(geo)),
                      treatment = sapply(GSMList(geo),
                                         function (x) Meta(x)$title))
samples <- samples[ids,]
rownames(samples) <- NULL
samples <- transform(samples,
                     EtOH = grepl("EtOH",treatment,fixed = TRUE),
                     RA   = grepl("RA",treatment,fixed = TRUE),
                     TGFb = grepl("TGFb",treatment,fixed = TRUE))
samples$label                 <- "EtOH"
samples[samples$RA,"label"]   <- "RA"
samples[samples$TGFb,"label"] <- "TGFb"
samples[with(samples,RA & TGFb),"label"] <- "RA+TGFb"
samples <- transform(samples,
                     label = factor(label,c("EtOH","RA","TGFb","RA+TGFb")))

K <- 3
x <- colSums(counts > 0)
i <- which(x > 3 &
             genes$GeneType == "protein-coding" &
             genes$Status == "active")
genes  <- genes[i,]
counts <- counts[,i]
n <- nrow(counts)

fit_list <- list()
cc_vec <- c(1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4)

for (cc in cc_vec) {
  
  fit_list[[as.character(cc)]] <- res_list$mcf7[[as.character(cc)]]
  
  colnames(fit_list[[as.character(cc)]]$LL) <- paste0("k", 1:K)
  
}

fit_list[["Inf"]] <- res_list$mcf7$`Inf`

log1p_loadings_order <- 1:3
tm_loadings_order <- c(1,3,2)

colnames(fit_list[[as.character(1)]]$LL) <- paste0(
  "k", log1p_loadings_order
)

fit_list[[as.character(1)]]$LL <- fit_list[[as.character(1)]]$LL[,paste0("k", 1:3)]

colnames(fit_list[[as.character(1e-3)]]$FF) <- paste0(
  "k", log1p_loadings_order
)

fit_list[[as.character(1e-3)]]$FF <- fit_list[[as.character(1e-3)]]$FF[,paste0("k", 1:3)]

colnames(fit_list[[as.character(Inf)]]$L) <- paste0(
  "k", tm_loadings_order
)

fit_list[[as.character(Inf)]]$L <- fit_list[[as.character(Inf)]]$L[,paste0("k", 1:3)]

colnames(fit_list[[as.character(Inf)]]$F) <- paste0(
  "k", tm_loadings_order
)

fit_list[[as.character(Inf)]]$F <- fit_list[[as.character(Inf)]]$F[,paste0("k", 1:3)]

plot_list <- list()

for (cc in cc_vec) {
  
  set.seed(1)
  plot_list[[as.character(cc)]] <- normalized_structure_plot(
    fit_list[[as.character(cc)]], 
    colors = topic_colors,
    grouping = samples$label,
    loadings_order = 1:n,
    topics = rev(1:3)
  )  + ggtitle(glue::glue("log1p Model Loadings (c = {cc})")) +
    theme(
      plot.title = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text.x = element_text(angle = 0,hjust = 0.5, size = 11)
    ) + ylab("Membership") +
    guides(fill=guide_legend(title="Factor")) +
    guides(colour = "none")
  
}

set.seed(1)
plot_list[[as.character(Inf)]] <- structure_plot(
  fit_list[[as.character(Inf)]], 
  colors = topic_colors,
  grouping = samples$label,
  loadings_order = 1:n,
  topics = rev(1:3)
)  + ggtitle("Topic Model Loadings") +
  theme(
    plot.title = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(angle = 0,hjust = 0.5, size = 11)
  ) + ylab("Membership") +
  guides(fill=guide_legend(title="Factor")) +
  guides(colour = "none")

g <- ggarrange(
  plotlist = plot_list, 
  nrow = 5, 
  ncol = 2, 
  common.legend = TRUE, 
  legend = "right"
)

ggsave(
  "../images/supp_mcf7_structure.png",
  g,
  device = "png",
  width = 10,
  height = 12
)

