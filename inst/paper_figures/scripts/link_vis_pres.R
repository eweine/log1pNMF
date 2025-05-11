out_list <- list()

for (cc in c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1e3)) {
  
  lambda <- seq(0, 100, 0.01)
  b <- max(1, cc) * log1p(lambda / cc)
  out_list[[as.character(cc)]] <- data.frame(
    lambda = lambda, 
    b = b,
    cc = cc
  )
  
}

out_df <- do.call(rbind, out_list)

library(ggplot2)
library(ggh4x)

g <- ggplot(data = out_df, aes(x = lambda, y = b)) +
  geom_line() +
  facet_manual(
    ~cc, scales = "free", 
    design = c("
               ABC
               DEF
               #G#
               "),
    labeller = labeller(cc = function(x) paste("c =", x))) +
  cowplot::theme_cowplot() +
  xlab(bquote("   " ~ lambda)) +
  ylab(expression(alpha ~ "×" ~ log(1 + lambda / c)))

ggsave(
  "/Users/eweine/Documents/passPCA/inst/paper_figures/pdfs/link_vis_pres.pdf",
  g,
  device = "pdf",
  width = 5,
  height = 5
)