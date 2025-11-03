# Used to check that x is a vector in which x[i+1] >= x[i] for all i.
expect_nondecreasing <- function (x, tol=1e-7)
  expect_equal(diff(x) >= -tol,rep(TRUE,length(x) - 1))

expect_all_equal <- function(x, tolerance = 1e-8) {
  expect_true(all(abs(x - x[1]) < tolerance),
              info = sprintf("Not all elements of x are equal within tolerance = %g", tolerance))
}
