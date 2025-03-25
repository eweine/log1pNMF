# here, I want to plot the differences in runtime based on the computational
# approximation. I'm not completely clear on what scale to plot things,
# but I will try a few things out and see what is most helpful to visualize

get_glmpca_o <- function(n, m, k) {

  n * m * k

}

get_tm_o <- function(n, m, k, pct_0) {

  size_nz <- n * m * (1 - pct_0)
  ((n + m) * k) + size_nz * k

}

get_log1p_approx_o <- function(n, m, k, pct_0) {

  size_nz <- n * m * (1 - pct_0)
  ((n + m) * (k + k^2)) + size_nz * k

}

m <- 10000
n_vec <- seq(26, 1000000, 10)
k <- 25
pct_0 <- 0.95

o_glmpca <- numeric(length(n_vec))
o_tm <- numeric(length(n_vec))
o_log1p_approx <- numeric(length(n_vec))

for (i in 1:length(n_vec)) {

  o_glmpca[i] <- get_glmpca_o(n_vec[i], m, k)
  o_tm[i] <- get_tm_o(n_vec[i], m, k, pct_0)
  o_log1p_approx[i] <- get_log1p_approx_o(n_vec[i], m, k, pct_0)

}

o_df <- data.frame(
  o = c(
    o_glmpca / o_tm,
    o_log1p_approx / o_tm
  ),
  Calculation = c(
    rep("log1p Exact", length(o_glmpca)),
    rep("log1p Approximation", length(o_tm))
  ),
  n = rep(n_vec, 2)
)
o_df$Calculation <- factor(x = o_df$Calculation, levels = c("log1p Exact", "log1p Approximation"))

library(ggplot2)

ga <- ggplot(data = o_df, aes(x = n, y = o, color = Calculation)) +
  geom_line(linewidth=1) +
  scale_x_continuous(
    labels = function(x) format(x, scientific = TRUE),
    trans = "log10"
  ) +
  scale_y_continuous(
    breaks = c(1, 5, 10, 20),
    trans = "log10"  # Ensure 1 is included
  ) +
  ylab("Complexity Relative to Topic Model (log10 scale)") +
  xlab("n (log10 scale)") +
  cowplot::theme_cowplot() +
  ggtitle("m = 10,000, K = 25, Sparsity = 95%") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values=c("blue", "red")) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 10)
  )


m <- 10000
n <- 10000
k_vec <- seq(1, 250, 1)
pct_0 <- 0.95

o_glmpca <- numeric(length(k_vec))
o_tm <- numeric(length(k_vec))
o_log1p_approx <- numeric(length(k_vec))

for (i in 1:length(k_vec)) {

  o_glmpca[i] <- get_glmpca_o(n, m, k_vec[i])
  o_tm[i] <- get_tm_o(n, m, k_vec[i], pct_0)
  o_log1p_approx[i] <- get_log1p_approx_o(n, m, k_vec[i], pct_0)

}

o_df <- data.frame(
  o = c(
    o_glmpca / o_tm,
    o_log1p_approx / o_tm
  ),
  Calculation = c(
    rep("log1p Exact", length(o_glmpca)),
    rep("log1p Approximation", length(o_tm))
  ),
  K = rep(k_vec, 2)
)
o_df$Calculation <- factor(x = o_df$Calculation, levels = c("log1p Exact", "log1p Approximation"))


library(ggplot2)

gb <- ggplot(data = o_df, aes(x = K, y = o, color = Calculation)) +
  geom_line(linewidth=1) +
  scale_y_continuous(
    transform = "log10", # Ensure 1 is included
    breaks = c(1, 5, 10, 20)
  ) +
  ylab("Complexity Relative to Topic Model (log10 scale)") +
  cowplot::theme_cowplot() +
  ggtitle("n = 10,000, m = 10,000, Sparsity = 95%") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values=c("blue", "red")) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 10)
  )

# finally, I think it would be useful to look at this
# as a function of how sparse the data are

m <- 10000
n <- 10000
k <- 25
pct_0_vec <- seq(0.01, 0.99, 0.001)

o_glmpca <- numeric(length(pct_0_vec))
o_tm <- numeric(length(pct_0_vec))
o_log1p_approx <- numeric(length(pct_0_vec))

for (i in 1:length(pct_0_vec)) {

  o_glmpca[i] <- get_glmpca_o(n, m, k)
  o_tm[i] <- get_tm_o(n, m, k, pct_0_vec[i])
  o_log1p_approx[i] <- get_log1p_approx_o(n, m, k, pct_0_vec[i])

}

o_df <- data.frame(
  o = c(
    o_glmpca / o_tm,
    o_log1p_approx / o_tm
  ),
  Calculation = c(
    rep("log1p Exact", length(o_glmpca)),
    rep("log1p Approximation", length(o_tm))
  ),
  pct_0 = rep(pct_0_vec, 2)
)
o_df$Calculation <- factor(x = o_df$Calculation, levels = c("log1p Exact", "log1p Approximation"))


gc <- ggplot(data = o_df, aes(x = pct_0, y = o, color = Calculation)) +
  geom_line(linewidth=1) +
  scale_y_continuous(
    transform = "log10" # Ensure 1 is included
  ) +
  scale_x_continuous(
    labels = scales::percent,      # Format x-axis labels as percentages
    breaks = scales::pretty_breaks(n = 5) # Optional: Adjust the number of breaks
  ) +
  ylab("Complexity Relative to Topic Model (log10 scale)") +
  xlab("Sparsity") +
  cowplot::theme_cowplot() +
  ggtitle("n = 10,000, m = 10,000, K = 25") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(values=c("blue", "red")) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 10),
    axis.text = element_text(size = 10)
  )

library(ggpubr)

blank <- ggplot() + theme_void()
f1 <- ggarrange(
  ggarrange(ga, gb, blank, nrow = 1,
            widths = c(1, 1, 0), legend = "none", labels = c("A", "B")),
  ggarrange(blank, gc, blank, nrow = 1,
            widths = c(0.25, 1.5, 0.25), legend = "right", labels = c("", "C", "")),
  nrow = 2
)

ggsave(
  "/Users/eweine/Documents/passPCA/inst/paper_figures/pdfs/computational_scaling.pdf",
  f1,
  device = "pdf",
  width = 8.25,
  height = 8
)

