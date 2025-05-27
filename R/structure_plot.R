# normalizes matrix A to have maximum value 1 in each column
normalize_bars <- function(A) {
  
  max_col <- apply(A, 2, max)
  sweep(A, 2, max_col, FUN = "/")
  
}

#' Structure Plot
#' 
#' Plot a "structure plot" of the loadings of a Poisson NMF fit. Each column of
#' the loadings is normalized to have a maximum value of one before plotting.
#' See function \code{structure_plot} in package \code{fastTopics} for more 
#' details.
#'
#' @param fit object of type \code{log1p_nmf_fit}.
#' @param ... additional argument to be passed to 
#' \code{fastTopics::structure_plot}.
#'
#' @return a \code{ggplot} object.
#' @export
#'
#' @examples
#' 
#' set.seed(1)
#' dat <- generate_log1p_pois_data(n = 500, p = 250, K = 4)
#'
#' fit_out <- fit_poisson_log1p_nmf(
#'  Y = dat$Y,
#'  K = 4,
#'  init_method = "random",
#'  control = list(verbose = FALSE)
#'  )
#'  
#' normalized_structure_plot(fit_out)
#' 
normalized_structure_plot <- function(fit, ...) {
  
  fastTopics::structure_plot(
    normalize_bars(pmax(fit$LL, 0)), ...
  )
  
}
