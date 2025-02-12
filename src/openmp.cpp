#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
bool openmp_available() {
#ifdef _OPENMP
  return true;
#else
  return false;
#endif
}

