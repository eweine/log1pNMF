# Factor matching for display: when two fits of the same dataset are shown
# together, permute one fit's factors to line up with the other's.
#
# PROCEDURE (described in the manuscript): the factors of every non-reference
# fit are put in one-to-one correspondence with the reference fit's factors
# by maximizing the sum of Spearman correlations between the paired columns
# of F, computed as an optimal assignment (Hungarian algorithm,
# clue::solve_LSAP). Rank correlation is used because the scale of F changes
# non-linearly across values of c. The matching affects visualization only.
#
# Each dataset has ONE reference fit -- the log1p fit featured in the main
# text (mcf7: c = 1; bbc: c = 1e-3; pancreas: c = 1) -- displayed in a fixed,
# manually chosen semantic order. Every other fit shown anywhere (the topic
# model, and each panel of the across-c supplementary figures) is matched
# DIRECTLY to that reference, so a factor keeps one color everywhere it
# appears, in the main and supplementary figures alike.
#
# Validation against the previous manual alignments (see the session notes):
# identical for MCF-7 (3/3) and BBC (10/10); for the pancreas, 11/13
# identical, and the remaining two pairs -- the topic model's two
# Endothelial/Mesenchymal factors, whose leading candidates' rank
# correlations differ by < 0.01 -- are swapped relative to the old manual
# choice.

library(clue)

#' Display labels for the columns of a fit, matched to a reference fit.
#'
#' @param ref_F p x K factor matrix of the reference fit.
#' @param F     p x K factor matrix of the fit to be aligned.
#' @param ref_labels integer display labels of the reference columns
#'   (default 1:K). Column j of F receives the label of the reference column
#'   it is paired with.
#' @return integer vector `labs` of length K: column j of F should be
#'   displayed as factor `labs[j]`.
match_display_labels <- function(ref_F, F, ref_labels = seq_len(ncol(ref_F))) {
  stopifnot(ncol(ref_F) == ncol(F), nrow(ref_F) == nrow(F))
  S <- suppressWarnings(stats::cor(as.matrix(ref_F), as.matrix(F),
                                   method = "spearman"))
  S[!is.finite(S)] <- 0
  pm <- as.integer(clue::solve_LSAP(S - min(S) + 1e-9, maximum = TRUE))
  labs <- integer(ncol(F))
  labs[pm] <- ref_labels           # reference column l is paired with F column pm[l]
  labs
}

#' Rename a fit's loading/factor matrices with matched labels and sort the
#' columns into k1..kK display order. `lab` comes from match_display_labels().
relabel_and_sort <- function(M, lab) {
  colnames(M) <- paste0("k", lab)
  M[, paste0("k", seq_len(ncol(M)))]
}
