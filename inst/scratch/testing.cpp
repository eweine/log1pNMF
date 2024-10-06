#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void timesTwo(double &x) {
  x = 4;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
x <- 3
timesTwo(x)
print(x)
*/
