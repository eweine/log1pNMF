load("../data/experiment_results.Rdata")

ll_approx_c_small <- unlist(
  lapply(
    res_list$sim_quality$`c = 1e-3`$fit_list_approx, function(x){x$ll}
    )
  )
ll_exact_c_small <- unlist(
  lapply(
    res_list$sim_quality$`c = 1e-3`$fit_list_exact, function(x){x$ll}
    )
  )

ll_approx_c1 <- unlist(
  lapply(
    res_list$sim_quality$`c = 1`$fit_list_approx, function(x){x$ll}
    )
  )
ll_exact_c1 <- unlist(
  lapply(
    res_list$sim_quality$`c = 1`$fit_list_exact, function(x){x$ll}
    )
  )

ll_approx_tm <- unlist(
  lapply(
    res_list$sim_quality$tm$fit_list_approx, function(x){x$ll}
  )
)
ll_exact_tm <- unlist(
  lapply(
    res_list$sim_quality$tm$fit_list_exact, function(x){x$ll}
  )
)

library(ggplot2)

ll_diff_df <- data.frame(
  cc = as.numeric(names(ll_approx_c_small)),
  ll_ratio = exp(ll_approx_c_small - ll_exact_c_small)
)

g1 <- ggplot(data = ll_diff_df, aes(x = cc, y = ll_ratio)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  xlab("c (log10 scale)") +
  ylab("Likelihood Ratio (approx / exact)") +
  ggtitle("log1p Link With c = 1e-3") +
  theme(
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12)
    ) +
  ylim(0, 1.05) +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10")

ll_diff_df <- data.frame(
  cc = as.numeric(names(ll_approx_c1)),
  ll_ratio = exp(ll_approx_c1 - ll_exact_c1)
)

g2 <- ggplot(data = ll_diff_df, aes(x = cc, y = ll_ratio)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  xlab("c (log10 scale)") +
  ylab("Likelihood Ratio (approx / exact)") +
  ggtitle("log1p Link With c = 1") +
  theme(
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12)
    ) +
  ylim(0, 1.05) +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10")

ll_diff_df <- data.frame(
  cc = as.numeric(names(ll_approx_tm)),
  ll_ratio = exp(ll_approx_tm - ll_exact_tm)
)

g3 <- ggplot(data = ll_diff_df, aes(x = cc, y = ll_ratio)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  xlab("c (log10 scale)") +
  ylab("Likelihood Ratio (approx / exact)") +
  ggtitle("Identity Link") +
  theme(
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12)
    ) +
  ylim(0, 1.05) +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10")

library(ggpubr)
g <- ggarrange(g1, g2, g3, nrow = 1, labels = "AUTO")

ggsave(
  "../pdfs/approx_quality.pdf",
  g,
  device = "pdf",
  width = 11.2,
  height = 3.5
)
