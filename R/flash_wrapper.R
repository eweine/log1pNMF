make_log1p_init_fn <- function(Y, s) {

  log1p_init_fn <- function(f) {

    if (is.null(f$EF)) {

      L_init <- matrix(
        data = rexp(n = nrow(Y), rate = 15),
        nrow = nrow(Y)
      )
      F_init <- matrix(
        data = rexp(n = ncol(Y), rate = 15),
        nrow = ncol(Y)
      )
      new_factor_num <- 1

    } else {

      L_init <- cbind(f$EF[[1]], rexp(n = nrow(Y), rate = 15))
      F_init <- cbind(f$EF[[2]], rexp(n = ncol(Y), rate = 15))
      new_factor_num <- 1 + ncol(f$EF[[1]])

    }

    new_fit <- fit_factor_model_log1p_quad_approx_sparse(
      Y = Y,
      K = new_factor_num,
      maxiter = 15,
      init_U = L_init,
      init_V = F_init,
      update_idx = c(new_factor_num),
      approx_range = c(0, 1.25),
      s = s
    )

    return(
      list(
        new_fit$U[,new_factor_num], new_fit$V[,new_factor_num]
      )
    )

  }

  return(log1p_init_fn)

}

#' Wrapper for flash with greedy initialization
#'
#' @param Y Original data Y
#' @param greedy_Kmax max number of greedy factors to add
#' @param var_type variance type for flash init
#' @param s size factor, on the order of 1.
#'
#' @return flashier fit
#' @export
#'
run_flash_log1p_with_greedy_init <- function(
    Y, greedy_Kmax = 50, var_type = 2, s = NULL
  ) {

  if (is.null(s)) {

    Y_tilde <- MatrixExtra::mapSparse(Y, log1p)
    S_tilde <- MatrixExtra::mapSparse(Y, function(x){sqrt(x / ((1 + x) ^ 2))})
    S_tilde_vals <- as.vector(as.matrix(S_tilde))
    S_tilde <- matrix(
      data = ifelse(
        S_tilde_vals == 0, 1, S_tilde_vals
      ),
      nrow = nrow(S_tilde),
      ncol = ncol(S_tilde)
    )

    rm(S_tilde_vals)
    gc()

    Y_tilde_vals <- as.vector(as.matrix(Y_tilde))

    Y_tilde <- matrix(
      data = ifelse(
        Y_tilde_vals == 0, -1, Y_tilde_vals
      ),
      nrow = nrow(S_tilde),
      ncol = ncol(S_tilde)
    )
    rm(Y_tilde_vals)
    gc()

  } else {

    Y_scaled <- Matrix::Diagonal(x = 1/s) %*% Y
    Y_tilde <- MatrixExtra::mapSparse(Y_scaled, log1p)

    Y_tilde_vals <- as.vector(as.matrix(Y_tilde))

    Y_tilde <- matrix(
      data = ifelse(
        Y_tilde_vals == 0, -1, Y_tilde_vals
      ),
      nrow = nrow(Y_tilde),
      ncol = ncol(Y_tilde)
    )
    rm(Y_tilde_vals)
    gc()

    S_nz_denom <- MatrixExtra::mapSparse(Y_scaled, function(x){1 / ((1 + x)^2)})
    S_nz_num <- Matrix::Diagonal(x = 1/(s^2)) %*% Y
    S_tilde <- S_nz_num * S_nz_denom




    S_tilde <- MatrixExtra::mapSparse(Y, function(x){sqrt(x / ((1 + x) ^ 2))})
    S_tilde_vals <- as.vector(as.matrix(S_tilde))
    S_tilde <- matrix(
      data = ifelse(
        S_tilde_vals == 0, 1, S_tilde_vals
      ),
      nrow = nrow(S_tilde),
      ncol = ncol(S_tilde)
    )

    rm(S_tilde_vals)
    gc()

    Y_tilde_vals <- as.vector(as.matrix(Y_tilde))

    Y_tilde <- matrix(
      data = ifelse(
        Y_tilde_vals == 0, -1, Y_tilde_vals
      ),
      nrow = nrow(S_tilde),
      ncol = ncol(S_tilde)
    )
    rm(Y_tilde_vals)
    gc()

  }


  f_log1p_init <- flashier::flash_init(
    data = Y_tilde,
    S = S_tilde,
    var_type = var_type
  )

  f_log1p_init <- flashier::flash_greedy(
    flash = f_log1p_init,
    Kmax = greedy_Kmax,
    ebnm_fn = ebnm::ebnm_point_exponential,
    init_fn = make_log1p_init_fn(Y, s)
  )

  f_log1p_init <- flashier::flash_backfit(f_log1p_init)
  return(f_log1p_init)

}

#' Wrapper for flash with MLE initialization
#'
#' @param Y data
#' @param K Number of factors in the MLE initialization
#' @param var_type variance type for flash backfit
#' @param mle_iter number of iterations for the MLE
#'
#' @return flash_object
#' @export
#'
run_flash_log1p_with_MLE_init <- function(Y, mle_iter = 100, K = 50, var_type = 2) {

  mle_init <- fit_factor_model_log1p_quad_approx_sparse(
    Y = Y,
    K = K,
    maxiter = mle_iter,
    approx_range = c(0, 1.25)
  )

  Y_tilde <- MatrixExtra::mapSparse(Y, log1p)
  S_tilde <- MatrixExtra::mapSparse(Y, function(x){sqrt(x / ((1 + x) ^ 2))})
  S_tilde_vals <- as.vector(as.matrix(S_tilde))
  S_tilde <- matrix(
    data = ifelse(
      S_tilde_vals == 0, 1, S_tilde_vals
    ),
    nrow = nrow(S_tilde),
    ncol = ncol(S_tilde)
  )

  rm(S_tilde_vals)
  gc()

  Y_tilde_vals <- as.vector(as.matrix(Y_tilde))

  Y_tilde <- matrix(
    data = ifelse(
      Y_tilde_vals == 0, -1, Y_tilde_vals
    ),
    nrow = nrow(S_tilde),
    ncol = ncol(S_tilde)
  )
  rm(Y_tilde_vals)
  gc()

  f_log1p_init <- flashier::flash_init(
    data = Y_tilde,
    S = S_tilde,
    var_type = var_type
  )

  f_log1p_init <- flashier::flash_factors_init(
    flash = f_log1p_init,
    init = list(
      mle_init$U, mle_init$V
    ),
    ebnm_fn = ebnm::ebnm_point_exponential
  )

  f_log1p_init <- flashier::flash_backfit(f_log1p_init)
  return(f_log1p_init)

}

