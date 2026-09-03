# Gene-set over-representation analysis (ORA) for fitted factors, to help
# interpret the biological examples (reviewer request). For each factor,
# the TOP N genes by gene score (the same genes shown in the paper's
# top-genes tables; ranking within a column of F is invariant to the
# display normalizations) are tested for over-representation in GO
# Biological Process terms against the background of all genes in the
# dataset, with a hypergeometric test and BH correction across terms.
#
# With only ~10 genes per list this is deliberately over-representation
# analysis, not ranked GSEA: a term is reported when the top genes hit it
# far more often than a random draw of the same size from the background
# would. Everything runs locally from Bioconductor annotation packages
# (org.Hs.eg.db + GO.db); no web APIs.
#
# Main entry point:
#   factor_ora(FF, n_top = 10, universe = rownames(FF), ...)
# returns one tidy data.frame (factor x enriched term), see below.

suppressMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
})

## symbols -> Entrez ids: exact symbol match first, alias fallback for the
## rest; unmapped genes are dropped (reported via attr "n_unmapped").
symbols_to_entrez <- function(symbols, orgdb = org.Hs.eg.db) {
  eg <- suppressMessages(mapIds(orgdb, keys = symbols,
                                column = "ENTREZID", keytype = "SYMBOL",
                                multiVals = "first"))
  miss <- names(eg)[is.na(eg)]
  if (length(miss) > 0) {
    ali <- suppressMessages(mapIds(orgdb, keys = miss,
                                   column = "ENTREZID", keytype = "ALIAS",
                                   multiVals = "first"))
    eg[names(ali)] <- ali
  }
  eg
}

## GO BP term -> Entrez gene sets (GO2ALLEGS includes descendant terms),
## restricted to a universe of Entrez ids and to sets of a sensible size.
go_bp_sets <- function(universe_eg, min_size = 5, max_size = 500,
                       orgdb = org.Hs.eg.db) {
  m <- suppressMessages(AnnotationDbi::select(
    orgdb, keys = unique(unname(universe_eg)),
    columns = c("GOALL", "ONTOLOGYALL"), keytype = "ENTREZID"))
  m <- unique(m[!is.na(m$GOALL) & m$ONTOLOGYALL == "BP",
                c("ENTREZID", "GOALL")])
  sets <- split(m$ENTREZID, m$GOALL)
  sets[lengths(sets) >= min_size & lengths(sets) <= max_size]
}

#' ORA of each factor's top genes against GO Biological Process.
#'
#' @param FF p x K matrix of gene scores. Rownames are gene identifiers:
#'   symbols by default, or NCBI/Entrez GeneIDs with id_type = "entrez"
#'   (the MCF-7 fits store GeneIDs).
#' @param n_top number of top genes per factor (default 10, matching the
#'   paper's top-genes tables).
#' @param universe character vector of background gene ids (default:
#'   all rownames of FF, i.e. every gene the model saw).
#' @param id_type "symbol" (mapped to Entrez via org.Hs.eg.db) or
#'   "entrez" (used as-is).
#' @param symbols optional named vector (names = ids in `universe`) of
#'   display symbols for the output's `genes` column; defaults to the ids
#'   themselves when id_type = "entrez".
#' @param fdr keep terms with BH-adjusted p below this (default 0.05).
#' @param max_terms keep at most this many terms per factor (default 10).
#' @return data.frame: factor, go_id, term, k (top genes in term),
#'   K (top genes mapped), m (background genes in term), N (background
#'   size), p, p_adj, genes (the overlapping genes, as display symbols).
factor_ora <- function(FF, n_top = 10, universe = rownames(FF),
                       id_type = c("symbol", "entrez"), symbols = NULL,
                       fdr = 0.05, max_terms = 10,
                       min_size = 5, max_size = 500,
                       orgdb = org.Hs.eg.db) {
  stopifnot(!is.null(rownames(FF)), !is.null(colnames(FF)))
  id_type <- match.arg(id_type)
  universe <- unique(universe)
  if (id_type == "symbol") {
    uni_eg <- symbols_to_entrez(universe, orgdb) # values Entrez, names = symbol
  } else {
    uni_eg <- structure(universe, names = if (is.null(symbols)) universe
                        else ifelse(is.na(symbols[universe]) |
                                    symbols[universe] == "",
                                    universe, symbols[universe]))
    uni_eg[!uni_eg %in% AnnotationDbi::keys(orgdb)] <- NA
  }
  n_unmapped <- sum(is.na(uni_eg))
  uni_eg <- uni_eg[!is.na(uni_eg)]
  message(sprintf("universe: %d ids, %d usable Entrez (%d unmapped)",
                  length(universe), length(uni_eg), n_unmapped))
  sets <- go_bp_sets(uni_eg, min_size, max_size, orgdb)
  message(sprintf("GO BP sets of size %d-%d in universe: %d",
                  min_size, max_size, length(sets)))
  N <- length(uni_eg)

  out <- list()
  for (kname in colnames(FF)) {
    ## uni_eg: values are Entrez ids, names are display symbols. FF's
    ## rownames are symbols (id_type = "symbol") or the Entrez ids
    ## themselves (id_type = "entrez"), so match on whichever side holds
    ## the FF ids.
    top_ids <- names(sort(FF[, kname], decreasing = TRUE))[seq_len(n_top)]
    top_eg  <- uni_eg[if (id_type == "symbol") names(uni_eg) %in% top_ids
                      else uni_eg %in% top_ids]
    K <- length(top_eg)
    hit <- vapply(sets, function(g) sum(top_eg %in% g), integer(1))
    m   <- lengths(sets)
    p   <- stats::phyper(hit - 1, m, N - m, K, lower.tail = FALSE)
    keep <- hit > 0
    d <- data.frame(factor = kname, go_id = names(sets)[keep],
                    k = hit[keep], K = K, m = m[keep], N = N,
                    p = p[keep], row.names = NULL)
    d$p_adj <- stats::p.adjust(d$p, method = "BH")
    d <- d[d$p_adj < fdr, ]
    d <- d[order(d$p), ]
    if (nrow(d) > max_terms) d <- d[seq_len(max_terms), ]
    if (nrow(d) > 0) {
      d$term  <- suppressMessages(mapIds(GO.db, keys = d$go_id,
                                         column = "TERM", keytype = "GOID"))
      d$genes <- vapply(d$go_id, function(g) paste(
        names(top_eg)[top_eg %in% sets[[g]]], collapse = ", "), "")
      out[[kname]] <- d
    } else {
      message("  ", kname, ": no GO BP terms at FDR < ", fdr)
    }
  }
  res <- do.call(rbind, out)
  if (is.null(res)) return(data.frame())
  rownames(res) <- NULL
  res[, c("factor", "go_id", "term", "k", "K", "m", "N",
          "p", "p_adj", "genes")]
}

#' Ranked GSEA of each factor's FULL gene-score vector against GO BP.
#'
#' Instead of testing a top-N cutoff, ranks every gene by its factor
#' score and asks (via fgsea) whether each GO term's genes concentrate
#' near the top of the ranking. Factor scores are non-negative, so the
#' test is one-sided (scoreType = "pos"); the many zero/near-zero scores
#' of sparse factors form ties at the bottom, which fgsea breaks
#' arbitrarily (harmless for the top of the ranking, where the signal
#' lives). Same id conventions as factor_ora().
#'
#' @return data.frame: factor, go_id, term, size, ES, NES, p, p_adj,
#'   leading_edge (the leading-edge genes, as display symbols).
factor_gsea_ranked <- function(FF, universe = rownames(FF),
                               id_type = c("symbol", "entrez"),
                               symbols = NULL, fdr = 0.05, max_terms = 10,
                               min_size = 5, max_size = 500,
                               orgdb = org.Hs.eg.db) {
  stopifnot(!is.null(rownames(FF)), !is.null(colnames(FF)))
  id_type <- match.arg(id_type)
  universe <- unique(universe)
  if (id_type == "symbol") {
    uni_eg <- symbols_to_entrez(universe, orgdb)
  } else {
    uni_eg <- structure(universe, names = if (is.null(symbols)) universe
                        else ifelse(is.na(symbols[universe]) |
                                    symbols[universe] == "",
                                    universe, symbols[universe]))
    uni_eg[!uni_eg %in% AnnotationDbi::keys(orgdb)] <- NA
  }
  uni_eg <- uni_eg[!is.na(uni_eg)]
  sets <- go_bp_sets(uni_eg, min_size, max_size, orgdb)
  message(sprintf("ranked GSEA: %d genes, %d GO BP sets",
                  length(uni_eg), length(sets)))
  sym_of <- setNames(names(uni_eg), uni_eg)    # entrez -> display symbol

  out <- list()
  for (kname in colnames(FF)) {
    v <- FF[, kname]
    ids <- if (id_type == "symbol") uni_eg[names(uni_eg) %in% names(v)] else
             uni_eg[uni_eg %in% names(v)]
    stats <- v[if (id_type == "symbol") names(ids) else ids]
    names(stats) <- ids
    fr <- suppressWarnings(fgsea::fgsea(
      pathways = sets, stats = stats, scoreType = "pos",
      minSize = min_size, maxSize = max_size, eps = 0))
    fr <- fr[fr$padj < fdr, ]
    fr <- fr[order(fr$pval), ]
    if (nrow(fr) > max_terms) fr <- fr[seq_len(max_terms), ]
    if (nrow(fr) > 0)
      out[[kname]] <- data.frame(
        factor = kname, go_id = fr$pathway,
        term = suppressMessages(mapIds(GO.db, keys = fr$pathway,
                                       column = "TERM", keytype = "GOID")),
        size = fr$size, ES = fr$ES, NES = fr$NES,
        p = fr$pval, p_adj = fr$padj,
        leading_edge = vapply(fr$leadingEdge, function(g)
          paste(sym_of[utils::head(g, 10)], collapse = ", "), ""),
        row.names = NULL)
    else message("  ", kname, ": no GO BP terms at FDR < ", fdr)
  }
  res <- do.call(rbind, out)
  if (is.null(res)) return(data.frame())
  rownames(res) <- NULL
  res
}
