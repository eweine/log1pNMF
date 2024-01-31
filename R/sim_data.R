#' @title Generate Data from a log1p Factor Model
#'
#' @param n Number of rows.
#'
#' @param p Number of columns.
#'
#' @param size_factor boolean indicating if size factor should be generated
#'
#' @param K Rank of the underlying mean structure
#' @importFrom stats rpois
#' @importFrom distr UnivarMixingDistribution
#' @importFrom distr Unif
#' @importFrom distr r
#' @importFrom distr Dirac
#'
#' @return list with underlying mean structure and observations
#'
#' @export
generate_data_simple <- function(n, p, K, size_factor = FALSE) {
  if (!is.scalar(K) || K < 1)
    stop("\"K\" must be an integer greater than or equal to 1")
  if (!is.scalar(n) || n < 1)
    stop("\"n\" must be an integer greater than or equal to 1")
  if (!is.scalar(p) || p < 1)
    stop("\"p\" must be an integer greater than or equal to 1")

  U <- matrix(nrow = n, ncol = K)
  V <- matrix(nrow = p, ncol = K)

  u_dist <- UnivarMixingDistribution(
    Unif(0,0.0001),
    Unif(0, 0.15),
    Dirac(.5),
    mixCoeff = c(.75, .2, .05)
  )

  v_dist <- UnivarMixingDistribution(
    Unif(0,0.75),
    Unif(0.75,7),
    mixCoeff = c(.98,.02)
  )

  u_sampler <- distr::r(u_dist)
  v_sampler <- distr::r(v_dist)

  for (k in 1:K) {
    U[, k] <- u_sampler(n)
    V[, k] <- v_sampler(p)
  }

  Lambda <- exp(tcrossprod(U, V)) - 1

  if (size_factor) {

    s_dist <- UnivarMixingDistribution(
      Unif(0.75,1.25),
      Dirac(0.5),
      Dirac(0.25),
      Dirac(2),
      Dirac(5),
      mixCoeff = c(0.5, rep(.15, 3), .05)
    )

    s_sampler <- distr::r(s_dist)

    s <- s_sampler(n)

    Lambda <- diag(s) %*% Lambda

  } else {

    s <- NULL

  }

  Y_dat <- rpois(n = n * p,lambda = as.vector(Lambda))
  Y <- matrix(data = Y_dat, nrow = n, ncol = p)

  return(
    list(
      Y = as(Y, "sparseMatrix"),
      U = U,
      V = V,
      s = s
    )
  )

}
