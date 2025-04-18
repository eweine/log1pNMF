g1 <- readr::read_rds("../data/bbc_sparsity_ggplot.rds")
g2 <- readr::read_rds("../data/pancreas_sparsity_ggplot.rds")

library(ggpubr)

g <- ggarrange(g1, g2, nrow = 2, ncol = 1)

ggplot2::ggsave(
  "../pdfs/sparsity.pdf",
  g,
  device = "pdf",
  width = 11.5,
  height = 7.5
)