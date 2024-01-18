#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<int> get_num_repeats(
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

// Function to create a list of vectors using a for loop
// [[Rcpp::export]]
List create_vals_list(
    const int num_vectors,
    const std::vector<int> vector_sizes,
    const std::vector<int> values
  ) {
  // Create an empty list
  List myList(num_vectors);

  int values_ctr = 0;

  // Create vectors in a for loop and add them to the list
  for (int i = 0; i < num_vectors; i++) {
    // Create a numeric vector
    IntegerVector vec(vector_sizes[i]);

    // Populate the vector with values (e.g., incrementing values)
    for (int j = 0; j < vector_sizes[i]; j++) {
      vec[j] = values[values_ctr];
      values_ctr++;
    }

    // Add the vector to the list
    myList[i] = vec;
  }

  return myList;
}

