# Identifiable "relative vs. absolute" simulation, individual-factor version
# (the reduced 2x2 design from inst/scratch/factor_sim_identifiable_reduc). The
# figure shows, on a single recovered factor, that log1p NMF prioritizes RELATIVE
# (fold) changes while the topic model (log1p as c -> Inf) prioritizes ABSOLUTE
# changes -- with no identifiability ambiguity, because a "+A and +B" group forbids
# a hard-cluster representation at rank 3.
#
# The data are small, so the models are fit here in the script.
#
# Composite:
#   A  ground-truth rates by group (baseline + additive program caps)
#   B  structure plots of sample scores: truth, topic (c = Inf), log1p (c = 1)
#   C  factor-value histograms of the recovered A/B factors, split by feature role
#
# Palettes (validated with the dataviz colorblind-safe checker):
#   FACCOL  factor colors -- the paper-wide Baseline/A/B (grey / red / teal).
#   ROLECOL 4-way feature-role palette -- ALL checks pass (worst adjacent CVD
#           dE 17.4; lighter fills get their contrast relief from the legend).

library(ggplot2)
library(fastTopics)
library(log1pNMF)
library(cowplot)
theme_set(theme_cowplot())
set.seed(1)

out_png <- "../images/factor_sim_identifiable_reduc_fig.png"

## ---- design: two programs (A, B), each acting on a low- and a high-baseline panel
nP <- 100                          # genes per program panel
p  <- 4 * nP

off_rel <- 2;    on_rel <- 20      # relative-change genes (10x fold, +18)
off_abs <- 1000; on_abs <- 1100    # absolute-change genes (1.1x fold, +100)

panel <- factor(rep(c("A-rel", "A-abs", "B-rel", "B-abs"), c(nP, nP, nP, nP)),
                levels = c("A-rel", "A-abs", "B-rel", "B-abs"))

baseP <- c(rep(off_rel, nP), rep(off_abs, nP), rep(off_rel, nP), rep(off_abs, nP))
devA  <- c(rep(on_rel - off_rel, nP), rep(on_abs - off_abs, nP), rep(0, nP), rep(0, nP))
devB  <- c(rep(0, nP), rep(0, nP), rep(on_rel - off_rel, nP), rep(on_abs - off_abs, nP))

## ---- simulate the four groups (Baseline, +A, +B, +A and B)
npop      <- 250
La        <- c(0, 1, 0, 1)         # program-A indicator per group
Lb        <- c(0, 0, 1, 1)         # program-B indicator per group
groupname <- c("Group 1 (Baseline)", "Group 2 (Baseline + A)",
               "Group 3 (Baseline + B)", "Group 4 (Baseline + A + B)")

prof   <- sapply(1:4, function(i) baseP + La[i] * devA + Lb[i] * devB)  # p x 4
Lambda <- t(prof[, rep(1:4, each = npop)])                             # n x p
n      <- nrow(Lambda)

Y <- matrix(rpois(n * p, as.vector(t(Lambda))), nrow = n, ncol = p, byrow = TRUE)
rownames(Y) <- paste0("cell", 1:n)
colnames(Y) <- paste0("gene", 1:p)

Ltrue <- cbind(Baseline = 1, A = La[rep(1:4, each = npop)], B = Lb[rep(1:4, each = npop)])

## ---- fit both models
# topic model (= log1p as c -> Inf): build the rank-1 init manually
r1 <- fastTopics:::fit_pnmf_rank1(Y)
init_LL <- cbind(r1$L, matrix(1e-3, n, 2)); rownames(init_LL) <- rownames(Y)
init_FF <- cbind(r1$F, matrix(1e-3, p, 2)); rownames(init_FF) <- colnames(Y)
tm <- fit_poisson_nmf(Y, fit0 = init_poisson_nmf(Y, F = init_FF, L = init_LL),
                      verbose = "none")

# log1p model at c = 1 (rank-1 is the default init_method)
set.seed(1)
lp <- fit_poisson_log1p_nmf(Y, K = 3, loglik = "exact", init_method = "rank1",
                            control = list(maxiter = 200, verbose = FALSE))

# put fitted factors in (Baseline, A, B) order so every panel agrees
ord_of <- function(L) {
  cc <- function(t) suppressWarnings(apply(L, 2, function(x) cor(x, t)))
  jA <- which.max(cc(Ltrue[, "A"])); jB <- which.max(cc(Ltrue[, "B"]))
  c(setdiff(seq_len(ncol(L)), c(jA, jB)), jA, jB)
}
relabel <- function(L) { M <- L[, ord_of(L), drop = FALSE]; colnames(M) <- c("Baseline", "A", "B"); M }

## ---- palettes (see header for validation)
FACCOL  <- c("Baseline" = "grey60", "A" = "#AE2012", "B" = "#005F73")
ROLECOL <- c("High-Baseline, Up"   = "#B36A00", "High-Baseline, Flat" = "#E0A81E",
             "Low-Baseline, Up"    = "#574494", "Low-Baseline, Flat"  = "#9B7FD4")

## ---- Panel A: ground-truth rates by group (grey baseline + additive colored caps)
seg1 <- function(k) {
  out <- data.frame(gene = 1:p, ymin = 0, ymax = baseP, factor = "Baseline")
  if (La[k] == 1) out <- rbind(out, data.frame(gene = 1:p, ymin = baseP,
      ymax = baseP + devA, factor = "A")[devA > 0, ])
  if (Lb[k] == 1) out <- rbind(out, data.frame(gene = 1:p, ymin = baseP,
      ymax = baseP + devB, factor = "B")[devB > 0, ])
  out$factor <- factor(out$factor, levels = c("Baseline", "A", "B")); out
}
segA <- do.call(rbind, lapply(1:4, function(k)
  transform(seg1(k), group = factor(groupname[k], levels = groupname))))

leg_factor <- get_legend(
  ggplot(data.frame(x = 1, f = factor(c("Baseline", "A", "B"),
                                      levels = c("Baseline", "A", "B"))),
         aes(x, 1, fill = f)) + geom_col() +
    scale_fill_manual(values = FACCOL, name = "Factor") +
    theme(legend.position = "right"))

panelA <- ggplot(segA, aes(xmin = gene - 0.5, xmax = gene + 0.5,
                           ymin = ymin, ymax = ymax, fill = factor)) +
  geom_rect() +
  facet_wrap(~ group, nrow = 1) +
  scale_fill_manual(values = FACCOL, name = "Factor", drop = FALSE) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 2, 20, 1000), limits = c(0, 1100)) +
  labs(x = "Feature", y = "True Rate") +
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_text(size = 13, face = "bold"))

## ---- Panel B: structure plots of sample scores (baseline on the bottom)
grp <- factor(rep(1:4, each = npop))
SP <- function(L, title, ylab = FALSE) {
  g <- structure_plot(log1pNMF:::normalize_bars(pmax(L, 0)),
                      loadings_order = seq_len(nrow(L)), grouping = grp, gap = 24,
                      topics = c("B", "A", "Baseline"), colors = FACCOL) +
    ggtitle(title) + labs(x = "Group", y = if (ylab) "Membership" else NULL) +
    theme(plot.title = element_text(size = 13, hjust = 0.5), legend.position = "none",
          axis.title.y = element_text(size = 14),
          axis.title.x = element_text(size = 15),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 13))
  if (!ylab) g <- g + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  g
}
panelB <- plot_grid(SP(Ltrue, "True Sample Scores", ylab = TRUE),
                    SP(relabel(tm$L),  "c = ∞ Sample Scores"),
                    SP(relabel(lp$LL), "c = 1 Sample Scores"), nrow = 1)

## ---- Panel C: factor-value histograms by feature role
role_for <- function(prog) {
  other <- setdiff(c("A", "B"), prog)
  factor(ifelse(panel == paste0(prog,  "-abs"), "High-Baseline, Up",
         ifelse(panel == paste0(other, "-abs"), "High-Baseline, Flat",
         ifelse(panel == paste0(prog,  "-rel"), "Low-Baseline, Up", "Low-Baseline, Flat"))),
         levels = c("High-Baseline, Up", "High-Baseline, Flat", "Low-Baseline, Up", "Low-Baseline, Flat"))
}
fac_df <- function(F, L, model) {
  Fo <- F[, ord_of(L)]                               # (Baseline, A, B)
  rbind(data.frame(value = Fo[, 2] / max(Fo[, 2]), role = role_for("A"),
                   factor = "Factor A", model = model),
        data.frame(value = Fo[, 3] / max(Fo[, 3]), role = role_for("B"),
                   factor = "Factor B", model = model))
}
facC <- rbind(fac_df(tm$F, tm$L, "c = ∞"), fac_df(lp$FF, lp$LL, "c = 1"))
facC$model <- factor(facC$model, levels = c("c = ∞", "c = 1"))

panelC <- ggplot(facC, aes(value, fill = role)) +
  geom_histogram(position = "identity", alpha = 0.85, bins = 60) +
  facet_grid(model ~ factor) +
  scale_fill_manual(values = ROLECOL, name = "Feature Type") +
  labs(x = "Factor Value (Normalized to Max 1)", y = "Count") +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7))

## ---- assemble: A over B share one Factor legend (right); C centered below
AB   <- plot_grid(panelA, panelB, ncol = 1, labels = c("A", "B"), rel_heights = c(1, 1))
top  <- plot_grid(AB, leg_factor, ncol = 2, rel_widths = c(1, 0.10))
Crow <- plot_grid(NULL, panelC, NULL, ncol = 3, rel_widths = c(0.12, 1, 0.12),
                  labels = c("", "C", ""))
final <- plot_grid(top, Crow, ncol = 1, rel_heights = c(2, 1.05))

ggsave(out_png, final, width = 12, height = 10.5, dpi = 150, bg = "white")
message("Wrote ", normalizePath(out_png, mustWork = FALSE))
