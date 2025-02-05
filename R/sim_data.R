#' @title Generate Data from a log1p Poisson NMF Model
#'
#' @description Generate data from log1p Poisson NMF Model with a specified
#'   rank.
#'
#' @details Generates data from the model
#' \deqn{y_{ij} \sim \textrm{Poisson}(\lambda_{ij})}
#' \deqn{\log(1 + \lambda_{ij} / c) = \sum_{k = 1}^{K} \ell_{ik}f_{jk},} where
#' \eqn{Y} is an \eqn{n \times p} matrix of non-negative counts, \eqn{L} is a
#' non-negative \eqn{n \times K} matrix of "loadings," \eqn{F} is a
#' non-negative \eqn{p \times K} matrix of "factors," and \eqn{c} is a
#' tuning parameter controlling the shape of the link function.
#'
#' @param n Number of rows.
#'
#' @param p Number of columns.
#'
#' @param K Rank of the underlying mean structure.
#'
#' @param cc Value of \eqn{c} in link function.
#'
#' @return list with the following components
#' \describe{
#'   \item{LL}{loadings of underlying mean structure. An \eqn{n \times K}
#'   matrix}
#'   \item{FF}{factors of underlying mean structure. A \eqn{p \times K} matrix.}
#'   \item{Y}{\eqn{n \times p} \code{dgCMatrix} matrix of generated data.}
#' }
#'
#' @examples
#' set.seed(1)
#' dat <- generate_log1p_pois_data(n = 1000, p = 500, K = 2)
#'
#' @importFrom stats rpois
#' @importFrom distr UnivarMixingDistribution
#' @importFrom distr Unif
#' @importFrom distr Dirac
#' @importFrom distr r
#' @importFrom methods as
#'
#' @export
#'
generate_log1p_pois_data <- function(n, p, K, cc = 1) {
  if (!is.scalar(K) || K < 1)
    stop("\"K\" must be an integer greater than or equal to 1")
  if (!is.scalar(n) || n < 1)
    stop("\"n\" must be an integer greater than or equal to 1")
  if (!is.scalar(p) || p < 1)
    stop("\"p\" must be an integer greater than or equal to 1")
  if (!is.scalar(cc) || cc <= 0)
    stop("\"cc\" must be a scalar greater than 0")

  LL <- matrix(nrow = n,ncol = K)
  FF <- matrix(nrow = p,ncol = K)

  l_dist <- UnivarMixingDistribution(
    Unif(0,0.1),
    Dirac(0.5),
    Dirac(1),
    mixCoeff = c(0.6, 0.2, 0.2)
  )

  f_dist <- UnivarMixingDistribution(
    Unif(0,0.1),
    Dirac(0.5),
    Dirac(1),
    mixCoeff = c(0.6, 0.2, 0.2)
  )

  l_sampler <- distr::r(l_dist)
  f_sampler <- distr::r(f_dist)

  for (k in 1:K) {
    LL[,k] <- l_sampler(n)
    FF[,k] <- f_sampler(p)
  }

  Lambda <- cc * (exp(tcrossprod(LL, FF)) - 1)

  Y <- matrix(
    data = rpois(n * p, as.vector(Lambda)),
    nrow = n, ncol = p
  )

  rownames(Y) <- paste0("c_",1:n)
  colnames(Y) <- paste0("g_",1:p)
  return(list(Y = as(Y, "CsparseMatrix"), LL = LL, FF = FF))
}
