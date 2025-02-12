.onAttach <- function(libname, pkgname) {
  if (!openmp_available()) {
    packageStartupMessage(
      "Note: This package was compiled without OpenMP support.
      For faster performance, please download OpenMP and reinstall the package."
      )
  }
}
