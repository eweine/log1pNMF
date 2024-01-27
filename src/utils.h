#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>

std::vector<int> get_num_repeats_cpp(
    const std::vector<int> idx,
    const int p,
    const int total_idx
);

std::vector<arma::vec> create_vals_vector_arma(
    const int num_vectors,
    const std::vector<int> vector_sizes,
    const std::vector<int> values
);

std::vector<arma::uvec> create_vals_vector_uvec(
    const int num_vectors,
    const std::vector<int> vector_sizes,
    const std::vector<int> values
);

#endif
