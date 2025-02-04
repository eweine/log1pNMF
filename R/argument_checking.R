# Return TRUE if x is a finite scalar with no missing entries.
is.scalar <- function (x)
  is.numeric(x) &
  length(x) == 1 &
  all(!is.na(x)) &
  all(is.finite(x))

# Verify that x is a valid count matrix.
#
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#'
verify.count.matrix <- function (x, arg.name = deparse(substitute(x))) {
  arg.name <- sprintf("\"%s\"",arg.name)
  msg <- paste("Input argument",arg.name,"should be a non-negative,",
               "numeric matrix with at least 2 rows and 2 columns,",
               "and all entries should be finite and non-missing")
  verify.matrix(x,arg.name)
  if (!(nrow(x) > 1 & ncol(x) > 1))
    stop(msg)
  else if (any(x < 0))
    stop(msg)

  if (any(Matrix::rowSums(x) == 0))
    stop(
      paste(
        "Input argument",arg.name,"cannot have any rows with all 0 entries"
      )
    )
  if (any(Matrix::colSums(x) == 0))
    stop(
      paste(
        "Input argument",arg.name,"cannot have any columns with all 0 entries"
      )
    )

  return(TRUE)
}

is.positive.numeric <- function(x) {
  is.numeric(x) && all(x > 0)
}

is.nonneg.numeric <- function(x) {
  is.numeric(x) && all(x >= 0)
}

is.positive.scalar <- function(x) {

  is.numeric(x) && (length(x) == 1)

}

is.above1.scalar <- function(x) {

  is.numeric(x) && (length(x) == 1) && (x >= 1)

}
