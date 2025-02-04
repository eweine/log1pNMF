#' Fit a Poisson NMF Model with log1p Link Function
#'
#' @description Fit a Poisson non-negative matrix factorization (NMF) model with
#' a log1p link using (approximate) maximum likelihood.
#'
#' @details Here, we fit the model
#' \deqn{y_{ij} \sim \textrm{Poisson}(s_{i}\lambda_{ij})}
#' \deqn{\log(1 + \lambda_{ij} / c) = \sum_{k = 1}^{K} \ell_{ik}f_{jk},} where
#' \eqn{Y} is an \eqn{n \times p} matrix of non-negative counts, \eqn{L} is a
#' non-negative \eqn{n \times K} matrix of "loadings," \eqn{F} is a
#' non-negative \eqn{p \times K} matrix of "factors," \eqn{s_{i}} is a
#' "size-factor" controlling for the total size of the \eqn{i^{\textrm{th}}}
#' row of \eqn{Y}, and \eqn{c} is a tuning parameter controlling the shape of
#' the link function.
#'
#' Exact optimization of the log-likelihood of this function can be slow, owing
#' mostly to the fact that \eqn{\exp(\sum_{k = 1}^{K} \ell_{ik}f_{jk})} must be
#' calculated for all \eqn{np} entries in \eqn{Y}. Especially for very sparse
#' data \eqn{Y}, the value of \eqn{\sum_{k = 1}^{K} \hat{\ell}_{ik}\hat{f}_{jk}}
#' evaluated at the maximum likelihood estimates \eqn{\hat{L}} and \eqn{\hat{F}}
#' is likely to be small for indices \eqn{(i, j)} where \eqn{y_{ij} = 0}. Thus,
#' we provide an approximate log-likelihood technique which approximates only
#' terms in the log-likelihood corresponding to indices \eqn{(i, j)} where
#' \eqn{y_{ij} = 0}. This approximation involves assumes
#' \eqn{\exp(x) \approx a_{0} + a_{1}x + a_{2}x^{2}} for some constants
#' \eqn{a_{0}}, \eqn{a_{1}}, and \eqn{a_{2}}. Selection of these constants can
#' be done either by a second order Taylor expansion about \eqn{x = 0} or by
#' optimizing a Chebyshev polynomial in a specified interval. Use of the
#' approximate technique can substantially increase the speed of computation,
#' but due to loss of accuracy we do not encourage it's use for values of
#' \eqn{c} smaller than \eqn{1}.
#'
#'
#' @param Y an \eqn{n \times p} data matrix. For single cell data, this should
#' be a cells by genes matrix. Must be convertible to an object of
#' type \code{dgCMatrix}.
#' @param K rank of the factorization. Must be specified unless
#' \code{init_LL} and \code{init_FF} are provided.
#' @param s size factor for each row of \code{Y} used to scale the link function
#' (see Details). The default value of \code{NULL} will
#' use the row sums of \eqn{Y} scaled to have an average value of \eqn{1}.
#' Custom size factors can be provided with a numeric vector of length \eqn{n}.
#' If set to \code{FALSE}, no size factors will be used.
#' @param cc value of \code{c} in the link function. See Details.
#' @param init_LL initial value of \eqn{L}.
#' @param init_FF initial value of \eqn{F}.
#' @param loglik character string specifying the log-likelihood to optimize.
#' Setting this to \code{"approx"} will approximate the log-likelihood
#' corresponding to the \eqn{0} values of \eqn{Y} with a quadratic.
#' \code{"exact"} optimizes the true log-likelihood. Generally, this
#' approximation works poorly for values of \eqn{c} small relative to \eqn{1}
#' and well for values of \eqn{c} large relative to \eqn{1}. Thus, the
#' \code{"default"} method will use the approximate log likelihood for values
#' of \eqn{c \geq 1}.
#' @param init_method Method for initialization (ignored if \code{init_LL}
#' and \code{init_FF} are specified). \code{"random"} will initialize \eqn{L}
#' and \eqn{F} to small positive values generated from an exponential
#' distribution. \code{"rank1"} will initialize the first factor by fitting
#' a model with \eqn{K = 1} to the data, with all other factors initialized to
#' small positive values from an exponential distribution.
#' @param approx_technique technique used for quadratic approximation of terms
#' corresponding to \eqn{0} values of \eqn{Y} in the log likelihood. Method
#' \code{"taylor"} uses a second order Taylor approximation of \eqn{\exp(x)}
#' about \eqn{x = 0}. Method \code{"chebyshev"} uses the Chebyshev coefficients
#' to find the best quadratic approximation of \eqn{\exp(x)} in the
#' interval \code{chebyshev_interval}.
#' @param chebyshev_interval a sorted numeric vector of length 2
#' (e.g. \code{c(0, 1)}) specifying the interval for Chebyshev approximation of
#' the terms of the log likelihood correspond to \eqn{0} values of \eqn{Y}. Only
#' used when \code{loglik} is set accordingly and when \code{approx_technique}
#' is set to \code{"chebyshev"}. By default, the interval is set to
#' \eqn{[0, \log(1 + 1/c)]}, which best approximates values of \eqn{\lambda}
#' between \eqn{0} and \eqn{1}. See Details for more information.
#' @param control a list of control parameters.
#'
#' @return a list containing the fitted values of \eqn{L} and \eqn{F}.
#' @export
#'
fit_poisson_log1p_nmf <- function(
  Y,
  K = NULL,
  s = NULL,
  cc = 1,
  init_LL = NULL,
  init_FF = NULL,
  loglik = c("default", "approx", "exact"),
  init_method = c("rank1", "random"),
  approx_technique = c("taylor", "chebyshev"),
  chebyshev_interval = NULL,
  control = list()
) {

  # TODO:
  # (1) Finish documentation above.
  # (2) Add control parameters.
  # (3) Test code.
  # (4) Add examples.

  # eventually, these should be turned into control parameters
  init_maxiter <- 5
  fit_maxiter <- 100

  verify.count.matrix(Y)
  loglik <- match.arg(loglik)
  init_method <- match.arg(init_method)
  approx_technique <- match.arg(approx_technique)

  # convert to dgCMatrix is Y isn't already
  if (!inherits(Y, "dgCMatrix")) {

    Y <- as(Y, "CsparseMatrix")

  }

  # default to rowSums on the order of 1 for size factor
  if (is.null(s)) {

    s <- Matrix::rowSums(Y)
    s <- s / mean(s)

  } else if (is.logical(s)) {

    if (s) {

      stop(
        "s cannot be set to TRUE. Either provide size factors or set s to FALSE"
        )

    } else {

      s <- rep(1, nrow(Y))

    }

  } else {

    if (!is.positive.numeric(s)) {

      stop("s must be a positive numeric vector.")

    }

    if (length(s) != nrow(Y)) {

      stop("s must have the same number of elements as there are rows of Y.")

    }

  }

  fit$s <- s

  if (!is.positive.scalar(cc)) {

    stop("cc must be a positive value")

  }

  if (loglik == "default") {

    # use the approximate optimization method if cc is at least 1
    if (cc >= 1) {

      loglik <- "approx"

    } else {

      loglik <- "exact"

    }

  }

  fit$loglik <- loglik

  if (loglik == "approx") {

    fit_fn <- fit_factor_model_log1p_quad_approx_sparse

    # find quadratic approximation coefficients
    if (approx_technique == "taylor") {

      a1 <- 1
      a2 <- 0.5

    }

    else if (approx_technique == "chebyshev") {

      # by default approximate exp(x) well for x in [0, 1]
      if (is.null(chebyshev_interval)) {

        chebyshev_interval <- c(0, log1p(1/cc))

      } else {

        if (
          !is.positive.numeric(chebyshev_interval) ||
          is.unsorted(chebyshev_interval) ||
          length(chebyshev_interval) != 2
        ) {

          stop(
            "chebyshev_interval must be a length 2 non-negative sorted numeric vector"
            )

        }

        poly_approx <- pracma::polyApprox(
          exp,
          chebyshev_interval[1],
          chebyshev_interval[2],
          2
        )

        a1 <- poly_approx$p[2]
        a2 <- poly_approx$p[1]

      }

    }

    fit$approx_technique <- approx_technique
    fit$a1 <- a1
    fit$a2 <- a2

  } else {

    fit_fn <- fit_factor_model_log1p_exact

  }

  sc <- Matrix::summary(Y)
  sc_t <- Matrix::summary(Matrix::t(Y))

  if (!is.null(init_LL) || !(is.null(init_FF))) {

    if (is.null(init_LL) || is.null(init_FF)) {

      stop(
        "If init_LL or init_FF is provided, BOTH init_LL and init_FF must be provided"
        )

    }

    fit <- list()
    if (ncol(init_LL) != ncol(init_FF)) {

      stop("Number of columns of init_LL and init_FF must match.")

    }

    if (nrow(init_LL) != nrow(Y)) {

      stop("The number of rows of init_LL must match the number of rows of Y.")

    }

    if (!is.nonneg.numeric(init_LL)) {

      stop("init_LL must be a non-negative numeric matrix.")

    }

    # add small positive constant to prevent initial NA log-likelihood
    fit$LL <- pmax(init_LL, 1e-16)

    if (nrow(init_FF) != ncol(Y)) {

      stop(
        "The number of rows of init_FF must match the number of columns of Y."
        )

    }

    if (!is.nonneg.numeric(init_FF)) {

      stop("init_FF must be a non-negative numeric matrix.")

    }

    # add small positive constant to prevent initial NA log-likelihood
    fit$FF <- pmax(init_FF, 1e-16)

  } else {

    if(is.null(K)) {

      stop("K must be specified if init_LL and init_FF are NULL.")

    }
    else if (!is.above1.scalar(K)) {

      stop("K must be a positive integer that is at least 1")

    } else {

      K <- as.integer(K)

    }

    if (init_method == "random") {

      exp_rate <- sqrt((1 / mean(Y@x)) * K * cc)

      fit$LL <- matrix(
        data = rexp(
          n = nrow(Y) * K,
          rate = exp_rate
        ),
        nrow = nrow(Y),
        ncol = K
      )

      fit$FF <- matrix(
        data = rexp(
          n = ncol(Y) * K,
          rate = exp_rate
        ),
        nrow = ncol(Y),
        ncol = K
      )

    } else if (init_method == "rank1") {

      exp_rate <- sqrt((1 / mean(Y@x)) * cc)

      fit$LL <- matrix(
        data = rexp(
          n = nrow(Y),
          rate = exp_rate
        ),
        nrow = nrow(Y),
        ncol = 1
      )

      fit$FF <- matrix(
        data = rexp(
          n = ncol(Y),
          rate = exp_rate
        ),
        nrow = ncol(Y),
        ncol = 1
      )

      fit <- fit_fn(
        sc = sc,
        sc_t = sc_t,
        n = nrow(Y),
        p = ncol(Y),
        fit = fit,
        maxiter = init_maxiter
      )

      remaining_exp_rate <- sqrt(cc * 15)

      fit$LL <- cbind(
        fit$LL,
          matrix(
            data = rexp(
              n = nrow(Y) * (K - 1),
              rate = remaining_exp_rate
            ),
            nrow = nrow(Y),
            ncol = K - 1
          )
        )

      fit$FF <- cbind(
        fit$FF,
        matrix(
          data = rexp(
            n = ncol(Y) * (K - 1),
            rate = remaining_exp_rate
          ),
          nrow = ncol(Y),
          ncol = K - 1
        )
      )

    }

  }

  # main fit
  fit <- fit_fn(
    sc = sc,
    sc_t = sc_t,
    n = nrow(Y),
    p = ncol(Y),
    fit = fit,
    maxiter = fit_maxiter
  )

  class(fit) <- "log1p_nmf"
  return(fit)

}
