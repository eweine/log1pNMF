#include <RcppArmadillo.h>
#include "utils.h"
using namespace Rcpp;
using namespace arma;

std::vector<int> get_num_repeats_cpp(
    const std::vector<int> idx,
    const int p,
    const int total_idx
) {

  std::vector<int> idx_cts(p, 1);

  int prev_val = idx[0];
  int ctr = 0;

  for (int i = 1; i < total_idx; i++) {

    if (idx[i] != prev_val) {

      ctr += 1;
      prev_val = idx[i];

    } else {

      idx_cts[ctr] += 1;

    }

  }

  return(idx_cts);

}

std::vector<arma::vec> create_vals_vector_arma(
    const int num_vectors,
    const std::vector<int> vector_sizes,
    const std::vector<int> values
) {
  // Preallocate memory for the vector
  std::vector<arma::vec> vectorOfVectors;
  vectorOfVectors.reserve(num_vectors);

  int values_ctr = 0;

  // Create vectors in a for loop and add them to the vector
  for (int i = 0; i < num_vectors; i++) {
    // Create an Armadillo vector
    arma::vec vec(vector_sizes[i]);

    // Populate the Armadillo vector with values (e.g., incrementing values)
    for (int j = 0; j < vector_sizes[i]; j++) {
      vec(j) = values[values_ctr];
      values_ctr++;
    }

    // Add the Armadillo vector to the vector
    vectorOfVectors.push_back(vec);
  }

  return vectorOfVectors;
}

std::vector<arma::uvec> create_vals_vector_uvec(
    const int num_vectors,
    const std::vector<int> vector_sizes,
    const std::vector<int> values
) {
  // Preallocate memory for the vector
  std::vector<arma::uvec> vectorOfUvecs;
  vectorOfUvecs.reserve(num_vectors);

  int values_ctr = 0;

  // Create vectors in a for loop and add them to the vector
  for (int i = 0; i < num_vectors; i++) {
    // Create an Armadillo unsigned integer vector
    arma::uvec uvec(vector_sizes[i]);

    // Populate the Armadillo unsigned integer vector with values (e.g., incrementing values)
    for (int j = 0; j < vector_sizes[i]; j++) {
      uvec(j) = static_cast<unsigned int>(values[values_ctr]);
      values_ctr++;
    }

    // Add the Armadillo unsigned integer vector to the vector
    vectorOfUvecs.push_back(uvec);
  }

  return vectorOfUvecs;
}

