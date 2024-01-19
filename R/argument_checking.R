# Return TRUE if x is a finite scalar with no missing entries.
is.scalar <- function (x)
  is.numeric(x) &
  length(x) == 1 &
  all(!is.na(x)) &
  all(is.finite(x))
