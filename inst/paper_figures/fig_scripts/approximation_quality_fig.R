load("../data/experiment_results.Rdata")

ll_approx_c_small <- unlist(
  lapply(
    res_list$sim_quality$`c = 1e-3`$fit_list_approx, function(x){x$ll}
    )
  )
ll_approx_cheb_c_small <- unlist(
  lapply(
    res_list$sim_quality$`c = 1e-3`$fit_list_approx_cheb, function(x){x$ll}
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
ll_approx_cheb_c1 <- unlist(
  lapply(
    res_list$sim_quality$`c = 1`$fit_list_approx_cheb, function(x){x$ll}
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
ll_approx_cheb_tm <- unlist(
  lapply(
    res_list$sim_quality$tm$fit_list_approx_cheb, function(x){x$ll}
  )
)
ll_exact_tm <- unlist(
  lapply(
    res_list$sim_quality$tm$fit_list_exact, function(x){x$ll}
  )
)

library(ggplot2)

ll_diff_df <- data.frame(
  cc = rep(as.numeric(names(ll_approx_c_small)), 2),
  ll_ratio = c(
    exp(ll_approx_c_small - ll_exact_c_small), 
    exp(ll_approx_cheb_c_small - ll_exact_c_small)
  ),
  approx_type = c(
    rep("Taylor", length(as.numeric(names(ll_approx_c_small)))),
    rep("Chebyshev", length(as.numeric(names(ll_approx_c_small))))
  )
)

g1 <- ggplot(data = ll_diff_df, aes(x = cc, y = ll_ratio, color = approx_type)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  xlab("c") +
  ylab("Likelihood Ratio (approx / exact)") +
  ggtitle("Data from log1p Model with c = 1e-3") +
  theme(
    axis.title.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    plot.title = element_text(size = 12.33)
    ) +
  ylim(0, 1.05) +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10") +
  labs(color = "Approximation") +
  geom_hline(yintercept = 1, linetype = "dashed")

ll_diff_df <- data.frame(
  cc = rep(as.numeric(names(ll_approx_c1)), 2),
  ll_ratio = c(
    exp(ll_approx_c1 - ll_exact_c1), 
    exp(ll_approx_cheb_c1 - ll_exact_c1)
  ),
  approx_type = c(
    rep("Taylor", length(as.numeric(names(ll_approx_c1)))),
    rep("Chebyshev", length(as.numeric(names(ll_approx_c1))))
  )
)

g2 <- ggplot(data = ll_diff_df, aes(x = cc, y = ll_ratio, color = approx_type)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  xlab("c") +
  ylab("Likelihood Ratio (approx / exact)") +
  ggtitle("Data from log1p Model with c = 1") +
  theme(
    axis.title.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    plot.title = element_text(size = 12.33)
  ) +
  ylim(0, 1.05) +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10") +
  labs(color = "Approximation") +
  geom_hline(yintercept = 1, linetype = "dashed")

ll_diff_df <- data.frame(
  cc = rep(as.numeric(names(ll_approx_tm)), 2),
  ll_ratio = c(
    exp(ll_approx_tm - ll_exact_tm), 
    exp(ll_approx_cheb_tm - ll_exact_tm)
  ),
  approx_type = c(
    rep("Taylor", length(as.numeric(names(ll_approx_tm)))),
    rep("Chebyshev", length(as.numeric(names(ll_approx_tm))))
  )
)

g3 <- ggplot(data = ll_diff_df, aes(x = cc, y = ll_ratio, color = approx_type)) +
  geom_point() +
  geom_line() +
  cowplot::theme_cowplot() +
  xlab("c") +
  ylab("Likelihood Ratio (approx / exact)") +
  ggtitle("Data from Topic Model") +
  theme(
    axis.title.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    plot.title = element_text(size = 12.33)
  ) +
  ylim(0, 1.05) +
  scale_x_continuous(breaks = c(1e-4, 1e-2, 1, 1e2, 1e4), transform = "log10") +
  labs(color = "Approximation") +
  geom_hline(yintercept = 1, linetype = "dashed")

library(ggpubr)
g <- ggarrange(g1, g2, g3, nrow = 1, common.legend = TRUE, legend = "right", labels = "AUTO")

ggsave(
  "../images/approx_quality.png",
  g,
  device = "png",
  width = 12.25,
  height = 3.5
)
