#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <omp.h>
#include "ll.h"

using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

double get_sparse_term_loglik_exact(
    const arma::mat U_T,
    const arma::mat V_T,
    const std::vector<int> nonzero_y,
    const std::vector<int> nonzero_y_i_idx,
    const std::vector<int> nonzero_y_j_idx,
    const arma::vec s,
    const int num_nonzero_y
) {

  double sum = 0.0;

  //#pragma omp parallel for reduction(+:sum)
  for (int r = 0; r < num_nonzero_y; r++) {

    sum += nonzero_y[r] * log(
      exp(
        dot(
          U_T.col(nonzero_y_i_idx[r]), V_T.col(nonzero_y_j_idx[r])
        )
      ) - s[nonzero_y_i_idx[r]]
    );
  }

  return sum;
}

// A worker struct that processes subsets of the data
struct SparseTermLogLikWorker : public Worker
{
  // References to input data
  const arma::mat &U_T;
  const arma::mat &V_T;
  const std::vector<int> &nonzero_y;
  const std::vector<int> &nonzero_y_i_idx;
  const std::vector<int> &nonzero_y_j_idx;
  const arma::vec &s;

  // Shared accumulator for partial sums
  double sum;

  // Constructor
  SparseTermLogLikWorker(const arma::mat &U_T,
                         const arma::mat &V_T,
                         const std::vector<int> &nonzero_y,
                         const std::vector<int> &nonzero_y_i_idx,
                         const std::vector<int> &nonzero_y_j_idx,
                         const arma::vec &s)
    : U_T(U_T), V_T(V_T),
      nonzero_y(nonzero_y),
      nonzero_y_i_idx(nonzero_y_i_idx),
      nonzero_y_j_idx(nonzero_y_j_idx),
      s(s),
      sum(0.0)
  {}

  // Splitting constructor (used internally by parallelFor)
  SparseTermLogLikWorker(const SparseTermLogLikWorker &w, Split)
    : U_T(w.U_T), V_T(w.V_T),
      nonzero_y(w.nonzero_y),
      nonzero_y_i_idx(w.nonzero_y_i_idx),
      nonzero_y_j_idx(w.nonzero_y_j_idx),
      s(w.s),
      sum(0.0)
  {}

  // Process the range [begin, end)
  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t r = begin; r < end; ++r) {
      const int i_idx = nonzero_y_i_idx[r];
      const int j_idx = nonzero_y_j_idx[r];

      // Compute the dot product of U_T.col(i_idx) and V_T.col(j_idx)
      double dot_val = arma::dot(U_T.col(i_idx), V_T.col(j_idx));

      // Accumulate sum: nonzero_y[r] * log(exp(dot_val) - s(i_idx))
      sum += nonzero_y[r] * std::log(std::exp(dot_val) - s[i_idx]);
    }
  }

  // Reduce partial sums
  void join(const SparseTermLogLikWorker &rhs) {
    sum += rhs.sum;
  }
};

// [[Rcpp::export]]
double get_sparse_term_loglik_exact_parallel(
    const arma::mat &U_T,
    const arma::mat &V_T,
    const std::vector<int> &nonzero_y,
    const std::vector<int> &nonzero_y_i_idx,
    const std::vector<int> &nonzero_y_j_idx,
    const arma::vec &s
) {
  // Create worker
  SparseTermLogLikWorker worker(U_T, V_T, nonzero_y,
                                nonzero_y_i_idx,
                                nonzero_y_j_idx,
                                s);

  // Dispatch parallelFor over the range [0, nonzero_y.size())
  parallelFor(0, nonzero_y.size(), worker);

  // The accumulated sum is in worker.sum
  return worker.sum;
}

double get_dense_term_loglik_exact(
    const arma::mat U_T,
    const arma::mat V_T,
    const int n,
    const int p
) {

  double sum = 0.0;

  //#pragma omp parallel for reduction(+:sum)
  for (int i = 0; i < n; i++) {

    for (int j = 0; j < p; j++) {

      sum -= exp(dot(U_T.col(i), V_T.col(j)));

    }
  }

  return sum;
}

struct DenseTermLogLikWorker : public Worker
{
  // References to input
  const arma::mat &U_T;
  const arma::mat &V_T;
  const int n;
  const int p;

  // Accumulator
  double sum;

  // Constructor
  DenseTermLogLikWorker(const arma::mat &U_T_,
                        const arma::mat &V_T_,
                        int n_,
                        int p_)
    : U_T(U_T_), V_T(V_T_),
      n(n_), p(p_),
      sum(0.0)
  {}

  // Splitting constructor
  DenseTermLogLikWorker(const DenseTermLogLikWorker &w, Split)
    : U_T(w.U_T), V_T(w.V_T),
      n(w.n), p(w.p),
      sum(0.0)
  {}

  // operator() processes i in [begin, end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      for (int j = 0; j < p; j++) {
        double dot_val = arma::dot(U_T.col(i), V_T.col(j));
        sum -= std::exp(dot_val);
      }
    }
  }

  // join() combines partial sums
  void join(const DenseTermLogLikWorker &rhs) {
    sum += rhs.sum;
  }
};

// [[Rcpp::export]]
double get_dense_term_loglik_exact_parallel(
    const arma::mat &U_T,
    const arma::mat &V_T,
    const int n,
    const int p
) {
  // Create worker
  DenseTermLogLikWorker worker(U_T, V_T, n, p);

  // Parallelize only over i in [0, n)
  parallelFor(0, n, worker);

  return worker.sum;
}

double get_sparse_term_loglik_quad_sparse_approx(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const std::vector<int>& nonzero_y,
    const std::vector<int>& nonzero_y_i_idx,
    const std::vector<int>& nonzero_y_j_idx,
    const int num_nonzero_y,
    const arma::vec& s,
    const double a1,
    const double a2
) {

  double sp_term = 0.0;
  double lin_correction = 0.0;
  double quad_correction = 0;
  double cp;

  //#pragma omp parallel for reduction(+:sp_term, lin_correction, quad_correction)
  for (int r = 0; r < num_nonzero_y; r++) {

    cp = dot(U_T.col(nonzero_y_i_idx[r]), V_T.col(nonzero_y_j_idx[r]));

    sp_term += nonzero_y[r] * log(exp(cp) - 1) -
      s[nonzero_y_i_idx[r]] * exp(cp);
    lin_correction += s[nonzero_y_i_idx[r]] * cp;
    quad_correction += s[nonzero_y_i_idx[r]] * cp * cp;

  }

  double ll = sp_term + a1 * lin_correction + a2 * quad_correction;

  return(ll);

}

struct SparseTermQuadApproxWorker : public Worker
{
  // Inputs
  const arma::mat &U_T;
  const arma::mat &V_T;
  const std::vector<int> &nonzero_y;
  const std::vector<int> &nonzero_y_i_idx;
  const std::vector<int> &nonzero_y_j_idx;
  const arma::vec &s;
  const double a1;
  const double a2;

  // Single partial sum for each worker
  double ll_sum;

  // Constructor
  SparseTermQuadApproxWorker(const arma::mat &U_T_,
                             const arma::mat &V_T_,
                             const std::vector<int> &nonzero_y_,
                             const std::vector<int> &nonzero_y_i_idx_,
                             const std::vector<int> &nonzero_y_j_idx_,
                             const arma::vec &s_,
                             double a1_,
                             double a2_)
    : U_T(U_T_), V_T(V_T_),
      nonzero_y(nonzero_y_),
      nonzero_y_i_idx(nonzero_y_i_idx_),
      nonzero_y_j_idx(nonzero_y_j_idx_),
      s(s_),
      a1(a1_), a2(a2_),
      ll_sum(0.0)
  {}

  // Splitting constructor (used internally by TBB)
  SparseTermQuadApproxWorker(const SparseTermQuadApproxWorker &w, Split)
    : U_T(w.U_T), V_T(w.V_T),
      nonzero_y(w.nonzero_y),
      nonzero_y_i_idx(w.nonzero_y_i_idx),
      nonzero_y_j_idx(w.nonzero_y_j_idx),
      s(w.s),
      a1(w.a1), a2(w.a2),
      ll_sum(0.0) // each split copy starts at 0
  {}

  // Process indices [begin, end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t r = begin; r < end; ++r) {
      double cp = arma::dot(
        U_T.col(nonzero_y_i_idx[r]),
        V_T.col(nonzero_y_j_idx[r])
      );

      // Add everything in one shot:
      //  = nonzero_y[r] * log(exp(cp) - 1)
      //    - s[i_idx]*exp(cp)
      //    + a1 * s[i_idx] * cp
      //    + a2 * s[i_idx] * cp^2
      double si = s[ nonzero_y_i_idx[r] ];

      ll_sum +=
        nonzero_y[r] * std::log(std::exp(cp) - 1.0)
        - si * std::exp(cp)
        + a1 * si * cp
        + a2 * si * cp * cp;
    }
  }

  // Combine partial results from another worker
  void join(const SparseTermQuadApproxWorker &rhs) {
    ll_sum += rhs.ll_sum;
  }
};

// [[Rcpp::export]]
double get_sparse_term_loglik_quad_sparse_approx_parallel(
    const arma::mat &U_T,
    const arma::mat &V_T,
    const std::vector<int> &nonzero_y,
    const std::vector<int> &nonzero_y_i_idx,
    const std::vector<int> &nonzero_y_j_idx,
    const arma::vec &s,
    const double a1,
    const double a2
) {
  // Number of elements in nonzero_y
  int num_nonzero = nonzero_y.size();

  // Instantiate the worker
  SparseTermQuadApproxWorker worker(
      U_T, V_T,
      nonzero_y,
      nonzero_y_i_idx,
      nonzero_y_j_idx,
      s,
      a1,
      a2
  );

  // Parallel loop over [0, num_nonzero)
  parallelFor(0, num_nonzero, worker);

  // The final sum is in worker.ll_sum
  return worker.ll_sum;
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

  double loglik = get_sparse_term_loglik_quad_sparse_approx(
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

  double loglik_sparse_term = get_sparse_term_loglik_exact(
    U_T,
    V_T,
    y_nz_vals,
    y_nz_rows_idx,
    y_nz_cols_idx,
    s,
    y_nz_vals.size()
  );

  Rprintf("loglik_sparse_term = %f\n", loglik_sparse_term);

  double loglik_dense_term = get_dense_term_loglik_exact(
    U_T,
    V_T,
    n,
    p
  );

  Rprintf("loglik_dense_term = %f\n", loglik_dense_term);


  double loglik = loglik_sparse_term + loglik_dense_term;

  return(loglik);

}

double get_loglik_quad_approx_sparse_parallel(
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

double get_loglik_exact_parallel(
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
    s
  );
  Rprintf("loglik_sparse_term = %f\n", loglik_sparse_term);

  double loglik_dense_term = get_dense_term_loglik_exact_parallel(
    U_T,
    V_T,
    n,
    p
  );

  Rprintf("loglik_dense_term = %f\n", loglik_dense_term);

  double loglik = loglik_sparse_term + loglik_dense_term;

  return(loglik);

}
