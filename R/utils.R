prep_data <- function(Y) {

  proc_data <- list()

  proc_data$sc <- Matrix::summary(Y)

  col_num_repeats <- get_num_repeats(
    proc_data$sc$j,
    ncol(Y),
    length(proc_data$sc$j)
  )

  proc_data$y_cols_data <- create_vals_list(
    ncol(Y),
    col_num_repeats,
    proc_data$sc$x
  )

  proc_data$y_cols_idx <- create_vals_list(
    ncol(Y),
    col_num_repeats,
    proc_data$sc$i - 1 # account for 0 indexing in C++
  )

  sc_t <- Matrix::summary(Matrix::t(Y))

  row_num_repeats <- get_num_repeats(
    sc_t$j,
    nrow(Y),
    length(sc_t$j)
  )

  proc_data$y_rows_data <- create_vals_list(
    nrow(Y),
    row_num_repeats,
    sc_t$x
  )

  proc_data$y_rows_idx <- create_vals_list(
    nrow(Y),
    row_num_repeats,
    sc_t$i - 1 # account for 0 indexing in C++
  )

  return(proc_data)

}
