get_loglik_quad_approx_full <- function(
  U,
  V,
  y_nz_vals,
  y_nz_rows_idx,
  y_nz_cols_idx,
  a1,
  a2
) {

  sparse_term <- get_sparse_term_loglik(
    t(U),
    t(V),
    y_nz_vals,
    y_nz_rows_idx,
    y_nz_cols_idx,
    length(y_nz_vals)
  )

  lin_term <- a1 * sum(colSums(U) * colSums(V))

  U_T_U <- crossprod(U)
  U_T_U_V_T <- tcrossprod(U_T_U, V)
  quad_term <- a2 * sum(t(V) * U_T_U_V_T)

  loglik <- sparse_term - lin_term - quad_term

  return(loglik)

}
