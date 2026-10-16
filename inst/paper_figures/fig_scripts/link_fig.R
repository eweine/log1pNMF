# here, for a fixed grid of values of lambda, I would like to plot
# lambda vs. log1p(lambda / c)

# I think that I could start with c = 1e-2, 1e-1, 1, 10, 100,
# and vary lambda between 0 and 100

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

library(ggplot2)
library(patchwork)

# build one square, strip-labelled panel per value of c; because every panel is
# forced square and patchwork aligns them, all 7 come out the same size
make_panel <- function(df) {
  ggplot(data = df, aes(x = lambda, y = b)) +
    geom_line() +
    facet_wrap(
      ~cc,
      labeller = labeller(cc = function(x) paste("c =", x))) +
    cowplot::theme_cowplot() +
    theme(
      aspect.ratio = 1,
      axis.title = element_blank()
    )
}

plots <- lapply(out_list, make_panel)

# put the y-axis title natively on the top-left panel only, so ggplot places it
# correctly beside that panel's axis (always visible, no manual positioning)
plots[[1]] <- plots[[1]] +
  ylab(bquote(g(lambda * ";" ~ c))) +
  theme(axis.title.y = element_text(size = 14))

# top row: 4 panels; bottom row: 3 panels centered with half-width spacers
row1 <- wrap_plots(plots[1:4], nrow = 1)
# left/right spacers are slightly asymmetric to counteract the rightward shift
# that the top-left panel's y-axis title adds to the top row, keeping the
# bottom row visually centred under the top row
row2 <- wrap_plots(
  c(list(plot_spacer()), plots[5:7], list(plot_spacer())),
  nrow = 1,
  widths = c(0.40, 1, 1, 1, 0.60)
)

# shared x-axis title
xlab <- wrap_elements(grid::textGrob(
  bquote(lambda), gp = grid::gpar(fontsize = 14)))

# stack both panel rows plus the x-axis title in a single column
g <- row1 / row2 / xlab + plot_layout(heights = c(1, 1, 0.12))

ggsave(
  "../images/link.png",
  g,
  device = "png",
  width = 8,
  height = 4.4
)
