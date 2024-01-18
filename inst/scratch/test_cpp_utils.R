# Here, I want to explore different ways of dealing with
# sparse matrices in R so that I can easily feed them to C++
library(Matrix)
load("~/Documents/data/fastglmpca/raw_data/pbmc_68k.RData")

sc <- summary(counts)

col_num_repeats <- passPCA:::get_num_repeats(
  sc$j,
  ncol(counts),
  length(sc$j)
)

y_cols_list <- passPCA:::create_vals_list(
  ncol(counts),
  col_num_repeats,
  sc$x
)

y_cols_idx <- passPCA:::create_vals_list(
  ncol(counts),
  col_num_repeats,
  sc$i
)


sc_t <- summary(Matrix::t(counts))

row_num_repeats <- passPCA:::get_num_repeats(
  sc_t$j,
  nrow(counts),
  length(sc_t$j)
)

y_rows_list <- passPCA:::create_vals_list(
  nrow(counts),
  row_num_repeats,
  sc_t$x
)

y_rows_idx <- passPCA:::create_vals_list(
  nrow(counts),
  row_num_repeats,
  sc_t$i
)
