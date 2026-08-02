# Shared loading and panel construction for the "one panel per simulating c"
# family of c-grid figures. Sourced by:
#
#   c_grid_loglik_by_ctrue_fig.R    log-likelihood shortfall
#   c_grid_rankcor_F_fig.R          mean rank correlation between factors
#   c_grid_hoyer_L_fig.R            mean Hoyer sparsity of the loadings
#   c_grid_hoyer_F_fig.R            mean Hoyer sparsity of the factors
#
# All four read c_grid_sim_metrics.rds (experiments/analyze_c_grid_sim.R) and
# share the same layout: a 2 x 4 grid of panels, one per c_true, with all 10
# seeds drawn individually as thin lines plus their mean, and the panel's
# matching c_fit circled. Keeping the seeds visible matters -- the aggregate
# hides whether a flat region is genuinely flat for every seed or only flat on
# average.
#
# c is a discrete axis throughout because c_fit = Inf is a real level of the
# grid (the c -> Inf limit, plain Poisson NMF) and has no place on a log axis.

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot(font_size = 11))

exp_dir <- "../experiments"

## ---- data ------------------------------------------------------------------
## Each cell was fit from a rank-1 and a random start; take the better of the
## two by log-likelihood. In this simulation that is immaterial -- the two agree
## to within 0.008 log-likelihood on all 640 cells -- but it is the right
## summary for a nonconvex objective.
c_grid_load <- function() {
  d   <- readRDS(file.path(exp_dir, "c_grid_sim_metrics.rds"))
  key <- paste(d$seed, d$c_true, d$c_fit)
  b   <- d[order(key, -d$loglik), ]
  b   <- b[!duplicated(paste(b$seed, b$c_true, b$c_fit)), ]

  b$subopt <- ave(b$loglik, paste(b$seed, b$c_true), FUN = max) - b$loglik

  clev <- sort(unique(b$c_fit))
  clab <- ifelse(is.infinite(clev), "∞",
                 format(clev, scientific = FALSE, trim = TRUE, drop0trailing = TRUE))
  fc   <- function(x) factor(x, levels = clev, labels = clab)

  b$cf     <- fc(b$c_fit)
  b$cf_lab <- fc(b$c_true)                       # matching c_fit, as a level
  b$ct     <- factor(paste("simulated at c =", fc(b$c_true)),
                     levels = paste("simulated at c =", clab))
  b
}

## ---- panel figure ----------------------------------------------------------
#' @param b     table from c_grid_load()
#' @param yvar  column to plot on y
#' @param ref   optional column holding the TRUE value of that quantity; drawn
#'              as a dashed horizontal line per panel, averaged over seeds so it
#'              is comparable to the mean curve
#' @param ...   passed to scale_y_continuous (trans, breaks, labels)
c_grid_panels <- function(b, yvar, ylab, title, subtitle = NULL, ref = NULL, ...) {

  b$.y <- b[[yvar]]
  avg  <- aggregate(.y ~ ct + cf, data = b, FUN = mean)

  ## the point where the fitted c equals the simulating c
  dg <- unique(b[as.character(b$cf) == as.character(b$cf_lab), c("ct", "cf")])
  dg <- merge(dg, avg, by = c("ct", "cf"))

  p <- ggplot(b, aes(cf, .y)) +
    geom_line(aes(group = seed), colour = "#9BB7BD", linewidth = 0.35, alpha = 0.8)

  if (!is.null(ref)) {
    ## mean over seeds, matching the thick mean curve. Using the first seed's
    ## value instead would put the line somewhere no plotted average lives.
    rf <- aggregate(stats::as.formula(paste(ref, "~ ct")), data = b, FUN = mean)
    names(rf)[2] <- ".ref"
    p <- p + geom_hline(data = rf, aes(yintercept = .ref),
                        linetype = "dashed", colour = "#C2410C", linewidth = 0.5)
  }

  p +
    geom_line(data = avg, aes(group = 1), colour = "#0E5C6B", linewidth = 1) +
    geom_point(data = avg, colour = "#0E5C6B", size = 1.6) +
    geom_point(data = dg, shape = 21, size = 3.6, stroke = 1.1,
               fill = NA, colour = "#C2410C") +
    facet_wrap(~ ct, nrow = 2) +
    scale_y_continuous(...) +
    labs(x = "c used for fitting", y = ylab, title = title, subtitle = subtitle) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 10),
          plot.title = element_text(size = 13, face = "plain"),
          plot.subtitle = element_text(size = 10, colour = "grey35"),
          panel.border = element_rect(colour = "grey80", fill = NA),
          axis.text.x = element_text(size = 8))
}

SUB_SEEDS <- "thin lines: 10 individual seeds;  thick line: mean;  circled: c_fit = c_true"
SUB_REF   <- paste(SUB_SEEDS, "\ndashed: the true value for that simulation")
