#ifndef UPDATE_ADD_CONST_H
#define UPDATE_ADD_CONST_H

#include <RcppArmadillo.h>

double update_a(
    double a,
    const arma::mat U_T,
    const arma::mat V_T,
    const std::vector<int> nonzero_y,
    const std::vector<int> nonzero_y_i_idx,
    const std::vector<int> nonzero_y_j_idx,
    const int num_nonzero_y,
    const arma::vec s,
    const int n,
    const int p,
    double alpha,
    double beta
);

#endif
