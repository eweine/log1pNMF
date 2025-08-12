.onAttach <- function(libname, pkgname) {
  if (!openmp_available()) {
    packageStartupMessage(
      "Note: This package was compiled without OpenMP support.
      For faster performance, please download OpenMP and reinstall the package.
      For installation instructions, see the README of our github repo: https://github.com/eweine/log1pNMF"
      )
  }
}
