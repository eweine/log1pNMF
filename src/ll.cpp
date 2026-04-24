#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <RcppParallel.h>
#include "ll.h"

using namespace Rcpp;
using namespace arma;

struct SparseTermLoglikExact : public RcppParallel::Worker {
  const arma::mat& U_T;
  const arma::mat& V_T;
  const std::vector<int>& nonzero_y;
  const std::vector<int>& nonzero_y_i_idx;
  const std::vector<int>& nonzero_y_j_idx;
  const arma::vec& s;
  double sum;

  SparseTermLoglikExact(
    const arma::mat& U_T, const arma::mat& V_T,
    const std::vector<int>& nonzero_y,
    const std::vector<int>& nonzero_y_i_idx,
    const std::vector<int>& nonzero_y_j_idx,
    const arma::vec& s
  ) : U_T(U_T), V_T(V_T), nonzero_y(nonzero_y),
      nonzero_y_i_idx(nonzero_y_i_idx), nonzero_y_j_idx(nonzero_y_j_idx),
      s(s), sum(0.0) {}

  SparseTermLoglikExact(const SparseTermLoglikExact& other, RcppParallel::Split)
    : U_T(other.U_T), V_T(other.V_T), nonzero_y(other.nonzero_y),
      nonzero_y_i_idx(other.nonzero_y_i_idx), nonzero_y_j_idx(other.nonzero_y_j_idx),
      s(other.s), sum(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t r = begin; r < end; r++) {
      sum += nonzero_y[r] * std::log(
        std::exp(arma::dot(U_T.col(nonzero_y_i_idx[r]), V_T.col(nonzero_y_j_idx[r])))
        - s[nonzero_y_i_idx[r]]
      );
    }
  }

  void join(const SparseTermLoglikExact& rhs) { sum += rhs.sum; }
};

double get_sparse_term_loglik_exact(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const std::vector<int>& nonzero_y,
    const std::vector<int>& nonzero_y_i_idx,
    const std::vector<int>& nonzero_y_j_idx,
    const arma::vec& s,
    const int num_nonzero_y
) {
  SparseTermLoglikExact worker(
    U_T, V_T, nonzero_y, nonzero_y_i_idx, nonzero_y_j_idx, s
  );
  RcppParallel::parallelReduce(0, num_nonzero_y, worker);
  return worker.sum;
}

struct DenseTermLoglikExact : public RcppParallel::Worker {
  const arma::mat& U_T;
  const arma::mat& V_T;
  const int p;
  double sum;

  DenseTermLoglikExact(const arma::mat& U_T, const arma::mat& V_T, int p)
    : U_T(U_T), V_T(V_T), p(p), sum(0.0) {}

  DenseTermLoglikExact(const DenseTermLoglikExact& other, RcppParallel::Split)
    : U_T(other.U_T), V_T(other.V_T), p(other.p), sum(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      arma::vec dots = V_T.t() * U_T.col(i);
      sum -= arma::sum(arma::exp(dots));
    }
  }

  void join(const DenseTermLoglikExact& rhs) { sum += rhs.sum; }
};

double get_dense_term_loglik_exact(
    const arma::mat& U_T,
    const arma::mat& V_T,
    const int n,
    const int p
) {
  DenseTermLoglikExact worker(U_T, V_T, p);
  RcppParallel::parallelReduce(0, n, worker);
  return worker.sum;
}

struct SparseTermLoglikQuadApprox : public RcppParallel::Worker {
  const arma::mat& U_T;
  const arma::mat& V_T;
  const std::vector<int>& nonzero_y;
  const std::vector<int>& nonzero_y_i_idx;
  const std::vector<int>& nonzero_y_j_idx;
  const arma::vec& s;
  const double a1;
  const double a2;
  double ll;

  SparseTermLoglikQuadApprox(
    const arma::mat& U_T, const arma::mat& V_T,
    const std::vector<int>& nonzero_y,
    const std::vector<int>& nonzero_y_i_idx,
    const std::vector<int>& nonzero_y_j_idx,
    const arma::vec& s, double a1, double a2
  ) : U_T(U_T), V_T(V_T), nonzero_y(nonzero_y),
      nonzero_y_i_idx(nonzero_y_i_idx), nonzero_y_j_idx(nonzero_y_j_idx),
      s(s), a1(a1), a2(a2), ll(0.0) {}

  SparseTermLoglikQuadApprox(const SparseTermLoglikQuadApprox& other, RcppParallel::Split)
    : U_T(other.U_T), V_T(other.V_T), nonzero_y(other.nonzero_y),
      nonzero_y_i_idx(other.nonzero_y_i_idx), nonzero_y_j_idx(other.nonzero_y_j_idx),
      s(other.s), a1(other.a1), a2(other.a2), ll(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t r = begin; r < end; r++) {
      double cp = arma::dot(
        U_T.col(nonzero_y_i_idx[r]), V_T.col(nonzero_y_j_idx[r])
      );
      ll += nonzero_y[r] * std::log(std::expm1(cp)) -
        s[nonzero_y_i_idx[r]] * std::exp(cp) +
        a1 * s[nonzero_y_i_idx[r]] * cp +
        a2 * s[nonzero_y_i_idx[r]] * cp * cp;
    }
  }

  void join(const SparseTermLoglikQuadApprox& rhs) { ll += rhs.ll; }
};

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
  SparseTermLoglikQuadApprox worker(
    U_T, V_T, nonzero_y, nonzero_y_i_idx, nonzero_y_j_idx, s, a1, a2
  );
  RcppParallel::parallelReduce(0, num_nonzero_y, worker);
  return worker.ll;
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


// [[Rcpp::export]]
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

  double loglik_dense_term = get_dense_term_loglik_exact(
    U_T,
    V_T,
    n,
    p
  );

  double loglik = loglik_sparse_term + loglik_dense_term;

  return(loglik);

}
