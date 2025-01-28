#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppParallel.h>
#include "ll.h"

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;


struct SparseTermLogLikExactWorker : public Worker
{
  // Pointers to the underlying data (read-only)
  const double* U_T_data;
  const double* V_T_data;
  const double* s_data;

  // Read-only references to vectors of indices/values
  const std::vector<int> &nonzero_y;
  const std::vector<int> &nonzero_y_i_idx;
  const std::vector<int> &nonzero_y_j_idx;

  // Dimensions (number of rows) for U_T and V_T
  // (they must share the same n_rows if you are doing dot products)
  const int n_rows_U;
  const int n_rows_V;

  // Accumulator for partial sums
  double sum;

  // Constructor
  SparseTermLogLikExactWorker(const arma::mat &U_T,
                              const arma::mat &V_T,
                              const arma::vec &s,
                              const std::vector<int> &nonzero_y_,
                              const std::vector<int> &nonzero_y_i_idx_,
                              const std::vector<int> &nonzero_y_j_idx_)
    : U_T_data(U_T.memptr()),
      V_T_data(V_T.memptr()),
      s_data(s.memptr()),
      nonzero_y(nonzero_y_),
      nonzero_y_i_idx(nonzero_y_i_idx_),
      nonzero_y_j_idx(nonzero_y_j_idx_),
      n_rows_U(U_T.n_rows),
      n_rows_V(V_T.n_rows),
      sum(0.0)
  {
    // (Optionally) you could add checks here to ensure dimensions match:
    // e.g.,  Rcpp::Rcout << "U_T: " << U_T.n_rows << "x" << U_T.n_cols << std::endl;
    //        Rcpp::Rcout << "V_T: " << V_T.n_rows << "x" << V_T.n_cols << std::endl;
  }

  // Split constructor for parallelReduce
  SparseTermLogLikExactWorker(const SparseTermLogLikExactWorker &w, Split)
    : U_T_data(w.U_T_data),
      V_T_data(w.V_T_data),
      s_data(w.s_data),
      nonzero_y(w.nonzero_y),
      nonzero_y_i_idx(w.nonzero_y_i_idx),
      nonzero_y_j_idx(w.nonzero_y_j_idx),
      n_rows_U(w.n_rows_U),
      n_rows_V(w.n_rows_V),
      sum(0.0)
  {}

  // The main loop for each thread-chunk
  void operator()(std::size_t begin, std::size_t end)
  {
    double localSum = 0.0;
    for (std::size_t r = begin; r < end; r++)
    {
      const int i = nonzero_y_i_idx[r];   // column index in U_T
      const int j = nonzero_y_j_idx[r];   // column index in V_T
      const int yval = nonzero_y[r];

      // Compute the dot product of U_T.col(i) and V_T.col(j).
      // Recall Armadillo is column-major, so element (row, col) is at
      //   pointer[row + col * n_rows].
      double dot_val = 0.0;
      for (int k = 0; k < n_rows_U; k++)
      {
        dot_val += U_T_data[k + i * n_rows_U] *
          V_T_data[k + j * n_rows_V];
      }

      // Evaluate the term yval * log(exp(dot) - s[i])
      double tmp = std::exp(dot_val) - s_data[i];
      // IMPORTANT: you might want to add a small check to ensure tmp>0
      // if there's any possibility it goes non-positive.
      localSum += yval * std::log(tmp);
    }
    sum += localSum;
  }

  // How partial sums are combined
  void join(const SparseTermLogLikExactWorker &rhs)
  {
    sum += rhs.sum;
  }
};

double get_sparse_term_loglik_exact_parallel(
    const arma::mat &U_T,
    const arma::mat &V_T,
    const std::vector<int> &nonzero_y,
    const std::vector<int> &nonzero_y_i_idx,
    const std::vector<int> &nonzero_y_j_idx,
    const arma::vec &s,
    const int num_nonzero_y
)
{
  // Instantiate the worker
  SparseTermLogLikExactWorker worker(U_T, V_T, s,
                                     nonzero_y,
                                     nonzero_y_i_idx,
                                     nonzero_y_j_idx);

  // Run the parallel reduction
  parallelReduce(0, num_nonzero_y, worker);

  // Return the accumulated sum
  return worker.sum;
}

struct DenseTermLogLikExactWorker : public Worker
{
  // Pointers to the underlying data
  const double* U_T_data;
  const double* V_T_data;

  // Dimensions (number of rows) for U_T and V_T
  // We assume these match so that dot products make sense
  const int n_rows_U;
  const int n_rows_V;

  // Size of outer loops
  const int n;
  const int p;

  // Accumulator for partial sums
  double sum;

  // Constructor
  DenseTermLogLikExactWorker(const arma::mat &U_T,
                             const arma::mat &V_T,
                             const int n_,
                             const int p_)
    : U_T_data(U_T.memptr()),
      V_T_data(V_T.memptr()),
      n_rows_U(U_T.n_rows),
      n_rows_V(V_T.n_rows),
      n(n_),
      p(p_),
      sum(0.0)
  {
  }

  // Split constructor for parallelReduce
  DenseTermLogLikExactWorker(const DenseTermLogLikExactWorker &w, Split)
    : U_T_data(w.U_T_data),
      V_T_data(w.V_T_data),
      n_rows_U(w.n_rows_U),
      n_rows_V(w.n_rows_V),
      n(w.n),
      p(w.p),
      sum(0.0)
  {
  }

  // The main loop for each chunk (1D indexing of i,j pairs)
  void operator()(std::size_t begin, std::size_t end)
  {
    double localSum = 0.0;
    for (std::size_t idx = begin; idx < end; idx++)
    {
      // Decode i, j from single index
      int i = idx / p;  // row-block index
      int j = idx % p;  // column-block index

      // Dot product of U_T.col(i) and V_T.col(j)
      double dot_val = 0.0;
      for (int k = 0; k < n_rows_U; k++)
      {
        // element (k,i) in U_T is at [k + i*n_rows_U]
        // element (k,j) in V_T is at [k + j*n_rows_V]
        dot_val += U_T_data[k + i * n_rows_U] *
          V_T_data[k + j * n_rows_V];
      }

      // Subtract exp(dot_val)
      localSum -= std::exp(dot_val);
    }
    // Accumulate into this functor's sum
    sum += localSum;
  }

  // Join partial sums
  void join(const DenseTermLogLikExactWorker &rhs)
  {
    sum += rhs.sum;
  }
};


double get_dense_term_loglik_exact_parallel(const arma::mat &U_T,
                                            const arma::mat &V_T,
                                            const int n,
                                            const int p)
{
  // Instantiate the worker
  DenseTermLogLikExactWorker worker(U_T, V_T, n, p);

  // We process [0, n*p) as a 1D range
  const std::size_t total = static_cast<std::size_t>(n) * p;

  // Execute the parallel reduce over that range
  parallelReduce(0, total, worker);

  // Return the accumulated result
  return worker.sum;
}

struct SparseTermLogLikQuadWorker : public Worker
{
  // Pointers to raw matrix data
  const double* U_T_data;
  const double* V_T_data;
  // Pointer to s
  const double* s_data;

  // Read-only references to vectors
  const std::vector<int> &nonzero_y;
  const std::vector<int> &nonzero_y_i_idx;
  const std::vector<int> &nonzero_y_j_idx;

  // Dimensions for U_T, V_T (number of rows)
  const int n_rows_U;
  const int n_rows_V;

  // Additional constants
  const double a1;
  const double a2;

  // Partial sums
  double sp_term;
  double lin_correction;
  double quad_correction;

  // Constructor
  SparseTermLogLikQuadWorker(
    const arma::mat &U_T,
    const arma::mat &V_T,
    const arma::vec &s,
    const std::vector<int> &nonzero_y_,
    const std::vector<int> &nonzero_y_i_idx_,
    const std::vector<int> &nonzero_y_j_idx_,
    double a1_,
    double a2_
  )
    : U_T_data(U_T.memptr()),
      V_T_data(V_T.memptr()),
      s_data(s.memptr()),
      nonzero_y(nonzero_y_),
      nonzero_y_i_idx(nonzero_y_i_idx_),
      nonzero_y_j_idx(nonzero_y_j_idx_),
      n_rows_U(U_T.n_rows),
      n_rows_V(V_T.n_rows),
      a1(a1_),
      a2(a2_),
      sp_term(0.0),
      lin_correction(0.0),
      quad_correction(0.0)
  {
  }

  // Split constructor for parallelReduce
  SparseTermLogLikQuadWorker(const SparseTermLogLikQuadWorker &w, Split)
    : U_T_data(w.U_T_data),
      V_T_data(w.V_T_data),
      s_data(w.s_data),
      nonzero_y(w.nonzero_y),
      nonzero_y_i_idx(w.nonzero_y_i_idx),
      nonzero_y_j_idx(w.nonzero_y_j_idx),
      n_rows_U(w.n_rows_U),
      n_rows_V(w.n_rows_V),
      a1(w.a1),
      a2(w.a2),
      sp_term(0.0),
      lin_correction(0.0),
      quad_correction(0.0)
  {
  }

  // The parallel loop
  void operator()(std::size_t begin, std::size_t end)
  {
    double local_sp_term = 0.0;
    double local_lin_correction = 0.0;
    double local_quad_correction = 0.0;

    for (std::size_t r = begin; r < end; r++)
    {
      const int i = nonzero_y_i_idx[r];
      const int j = nonzero_y_j_idx[r];
      const int yval = nonzero_y[r];

      // Dot product of U_T.col(i) and V_T.col(j)
      double cp = 0.0;
      for (int k = 0; k < n_rows_U; k++)
      {
        cp += U_T_data[k + i * n_rows_U] * V_T_data[k + j * n_rows_V];
      }

      // sp_term piece
      // sp_term += y[r]*log(exp(cp) - 1) - s[i]*exp(cp)
      double exp_cp = std::exp(cp);

      // NOTE: If cp < 0, it's possible exp_cp - 1 <= 0 => log() invalid.
      // We assume the original code ensures cp>0 or handles the domain properly.
      local_sp_term += yval * std::log(exp_cp - 1.0)
        - s_data[i] * exp_cp;

      // lin_correction piece
      local_lin_correction += s_data[i] * cp;

      // quad_correction piece
      local_quad_correction += s_data[i] * (cp * cp);
    }

    // Accumulate partial sums into the functor
    sp_term         += local_sp_term;
    lin_correction  += local_lin_correction;
    quad_correction += local_quad_correction;
  }

  // Join partial sums from different threads
  void join(const SparseTermLogLikQuadWorker &rhs)
  {
    sp_term         += rhs.sp_term;
    lin_correction  += rhs.lin_correction;
    quad_correction += rhs.quad_correction;
  }
};

// [[Rcpp::export]]
double get_sparse_term_loglik_quad_sparse_approx_parallel(
    const arma::mat &U_T,
    const arma::mat &V_T,
    const std::vector<int> &nonzero_y,
    const std::vector<int> &nonzero_y_i_idx,
    const std::vector<int> &nonzero_y_j_idx,
    const int num_nonzero_y,
    const arma::vec &s,
    const double a1,
    const double a2
)
{
  // Instantiate worker
  SparseTermLogLikQuadWorker worker(
      U_T, V_T, s, nonzero_y, nonzero_y_i_idx, nonzero_y_j_idx, a1, a2
  );

  // Parallel reduce over the range [0, num_nonzero_y)
  parallelReduce(0, num_nonzero_y, worker);

  // Combine results into the final
  double ll = worker.sp_term
  + a1 * worker.lin_correction
  + a2 * worker.quad_correction;

  return ll;
}


double get_loglik_quad_approx_sparse(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const std::vector<int>& y_nz_vals,
    const std::vector<int>& y_nz_rows_idx,
    const std::vector<int>& y_nz_cols_idx,
    const arma::vec& s,
    const double a1,
    const double a2
) {

  double loglik = get_sparse_term_loglik_quad_sparse_approx_parallel(
    U_T,
    V_T,
    y_nz_vals,
    y_nz_rows_idx,
    y_nz_cols_idx,
    y_nz_vals.size(),
    s,
    a1,
    a2
  );

  double lin_term = a1 * arma::dot(U_T * s, arma::sum(V_T, 1));

  arma::mat U_T_U_V_T = (U_T.each_row() % s.t()) * U_T.t() * V_T;
  double quad_term = a2 * arma::accu(V_T % U_T_U_V_T);

  loglik = loglik - lin_term - quad_term;
  return(loglik);

}


double get_loglik_exact(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const std::vector<int>& y_nz_vals,
    const std::vector<int>& y_nz_rows_idx,
    const std::vector<int>& y_nz_cols_idx,
    const arma::vec& s,
    const int n,
    const int p
) {

  double loglik_sparse_term = get_sparse_term_loglik_exact_parallel(
    U_T,
    V_T,
    y_nz_vals,
    y_nz_rows_idx,
    y_nz_cols_idx,
    s,
    y_nz_vals.size()
  );

  double loglik_dense_term = get_dense_term_loglik_exact_parallel(
    U_T,
    V_T,
    n,
    p
  );

  double loglik = loglik_sparse_term + loglik_dense_term;

  return(loglik);

}
