#include <vector>
#include <algorithm>
#include <cmath>

#include "threads.h"

#include <R.h>
#include <Rinternals.h>

using namespace std;

/*
 Diagonal of the inverse of Q = L L^T, computed with the Takahashi recursion
 on the sparsity pattern of the Cholesky factor.

 Input is a triangular factor in CSC format (0-based). If lower is TRUE it is
 the lower triangular L, otherwise the upper triangular R = L^T. The upper
 case is transposed to lower CSC once, which is the same as reading R in CSR.

 The selected inverse Z is stored flat, sharing L's pattern exactly, so that
 Z(r, c) for r >= c is z[q] at the offset q that holds L(r, c). The current
 column is scattered into dense workspaces (pos, lj), which turns what would
 be a search per access into an array index.

 Column j only needs the columns of the rows below its diagonal, which are
 its ancestors in the elimination tree, so the columns can be computed in
 any order where each column comes after its parent, the first row below its
 diagonal. The columns are therefore computed in parallel by walking the tree
 from the roots: after a column is computed, each child with a large subtree
 becomes a new task, and the smaller subtrees are computed directly by the
 same thread. Each column is computed by one thread in the same way, so the
 results do not depend on the number of threads. The longest path from a
 root limits the speed-up, to about four for two-dimensional fields.
*/

namespace {

// Workspace for computing a column, see qinv_diag.
struct QinvWork {
  std::vector<int> pos;
  std::vector<double> lj, ljc, acc;
  explicit QinvWork(int n) : pos(n, -1), lj(n, 0.0) {}
};

struct QinvFactor {
  const int *Lp, *Li;
  const double *Lx;
  double *z;
};

// Scatter column j into the workspace, and return the number of rows below
// the diagonal.
inline int qinv_scatter(const QinvFactor &F, int j, QinvWork &w) {
  const int start = F.Lp[j], end = F.Lp[j + 1];
  const double Ljj = F.Lx[start];
  const int m = end - start - 1;
  w.ljc.resize(m);
  for (int t = 0; t < m; ++t) {
    const int r = F.Li[start + 1 + t];
    w.pos[r] = t;
    w.lj[r] = w.ljc[t] = F.Lx[start + 1 + t] / Ljj;
  }
  return m;
}

// acc[a] -= sum over the rows b in [b0, b1) of column j. Walking column
// c = r_b yields Z(r, c) for r >= c, and that value serves both the
// (i = r, k = c) and the (i = c, k = r) term, so one pass over the already
// computed columns covers every pair without a search. If the rows of column
// c below its diagonal are the rows of column j after r_b, which is the
// common case for a supernodal pattern, the rows are matched by position, and
// the pass is a contiguous update and inner product.
inline void qinv_accumulate(const QinvFactor &F, int j, int b0, int b1,
                            const QinvWork &w, double *acc) {
  const int start = F.Lp[j];
  const int m = F.Lp[j + 1] - start - 1;
  const int *rows = F.Li + start + 1;
  const double *ljc = w.ljc.data();
  for (int b = b0; b < b1; ++b) {
    const int c = rows[b];
    const double lc = ljc[b];
    const int q0 = F.Lp[c], q1 = F.Lp[c + 1];
    const int rem = m - b - 1;
    if (q1 - q0 - 1 == rem && std::equal(F.Li + q0 + 1, F.Li + q1, rows + b + 1)) {
      const double *zc = F.z + q0 + 1;
      double *acc_b = acc + b + 1;
      const double *lj_b = ljc + b + 1;
      double s = 0.0;
      for (int t = 0; t < rem; ++t) {
        const double v = zc[t];
        acc_b[t] -= lc * v;
        s += lj_b[t] * v;
      }
      acc[b] -= lc * F.z[q0] + s;
      continue;
    }
    for (int q = q0; q < q1; ++q) {
      const int a = w.pos[F.Li[q]];
      if (a < 0)
        continue;
      const double v = F.z[q];
      acc[a] -= lc * v;
      if (F.Li[q] != c)
        acc[b] -= w.lj[F.Li[q]] * v;
    }
  }
}

// Store column j from acc, and clear the workspace.
inline void qinv_finish(const QinvFactor &F, int j, int m, QinvWork &w,
                        const double *acc, double *variances) {
  const int start = F.Lp[j];
  const double Ljj = F.Lx[start];
  double accd = 0.0;
  for (int t = 0; t < m; ++t) {
    F.z[start + 1 + t] = acc[t];
    accd += w.lj[F.Li[start + 1 + t]] * acc[t];
  }
  F.z[start] = 1.0 / (Ljj * Ljj) - accd;
  variances[j] = F.z[start];
  for (int t = 0; t < m; ++t) {
    const int r = F.Li[start + 1 + t];
    w.pos[r] = -1;
    w.lj[r] = 0.0;
  }
}

// Column j, sequentially.
inline void qinv_column(const QinvFactor &F, int j, QinvWork &w,
                        double *variances) {
  const int m = qinv_scatter(F, j, w);
  w.acc.assign(m, 0.0);
  qinv_accumulate(F, j, 0, m, w, w.acc.data());
  qinv_finish(F, j, m, w, w.acc.data(), variances);
}

// The columns of the subtree of column j, each after its parent, where the
// children with a subtree of at least grain work are new tasks.
void qinv_subtree(const QinvFactor &F, int j, const std::vector<int> &cp,
                  const std::vector<int> &ci, const std::vector<double> &sub,
                  double grain, std::vector<QinvWork *> &work_t,
                  double *variances) {
  int t = 0;
#ifdef _OPENMP
  t = omp_get_thread_num();
#endif
  QinvWork &w = *work_t[t];
  std::vector<int> stack(1, j);
  while (!stack.empty()) {
    const int k = stack.back();
    stack.pop_back();
    qinv_column(F, k, w, variances);
    for (int q = cp[k]; q < cp[k + 1]; ++q) {
      const int c = ci[q];
      if (sub[c] >= grain) {
        #pragma omp task firstprivate(c) shared(F, cp, ci, sub, work_t)
        qinv_subtree(F, c, cp, ci, sub, grain, work_t, variances);
      } else {
        stack.push_back(c);
      }
    }
  }
}

} // namespace

// Returns 0 on success, or the 1-based column with a missing diagonal. If
// selected is not NULL, the selected inverse is also copied to it, in the
// pattern of L (only used for lower triangular input).
static int qinv_diag(int n, const int *Ap, const int *Ai, const double *Ax,
                     bool lower, double *variances, double *selected = NULL,
                     int n_threads = 1) {
  const int nnz = Ap[n];

  const int *Lp, *Li;
  const double *Lx;
  vector<int> tp, ti;
  vector<double> tx;
  if (lower) {
    Lp = Ap;
    Li = Ai;
    Lx = Ax;
  } else {
    // Transpose the upper triangular R to lower triangular L = R^T. A counting
    // sort by row keeps the row indices in each column of L sorted, so the
    // diagonal comes first in each column.
    tp.assign(n + 1, 0);
    ti.resize(nnz);
    tx.resize(nnz);
    for (int q = 0; q < nnz; ++q)
      ++tp[Ai[q] + 1];
    for (int j = 0; j < n; ++j)
      tp[j + 1] += tp[j];
    vector<int> next(tp.begin(), tp.end() - 1);
    for (int c = 0; c < n; ++c) {
      for (int q = Ap[c]; q < Ap[c + 1]; ++q) {
        const int r = Ai[q];
        ti[next[r]] = c;
        tx[next[r]++] = Ax[q];
      }
    }
    Lp = tp.data();
    Li = ti.data();
    Lx = tx.data();
  }

  for (int j = 0; j < n; ++j) {
    if (Lp[j] >= Lp[j + 1] || Li[Lp[j]] != j)
      return j + 1;
  }

  vector<double> z(nnz, 0.0);
  QinvFactor F = {Lp, Li, Lx, z.data()};

  // The elimination tree, the work of each column, and the work of each
  // subtree. The parent of a column has a larger index.
  vector<int> parent(n, -1);
  vector<double> work(n), sub(n);
  double total = 0.0;
  for (int j = 0; j < n; ++j) {
    const int start = Lp[j], end = Lp[j + 1];
    if (end - start > 1)
      parent[j] = Li[start + 1];
    double wj = end - start;
    for (int q = start + 1; q < end; ++q)
      wj += Lp[Li[q] + 1] - Lp[Li[q]];
    work[j] = wj;
    total += wj;
  }
  for (int j = 0; j < n; ++j) {
    sub[j] += work[j];
    if (parent[j] >= 0)
      sub[parent[j]] += sub[j];
  }
  // The children of each column
  vector<int> cp(n + 1, 0), ci(n);
  for (int j = 0; j < n; ++j)
    if (parent[j] >= 0)
      ++cp[parent[j] + 1];
  for (int j = 0; j < n; ++j)
    cp[j + 1] += cp[j];
  {
    vector<int> next(cp.begin(), cp.end() - 1);
    for (int j = 0; j < n; ++j)
      if (parent[j] >= 0)
        ci[next[parent[j]]++] = j;
  }

  // The work of the longest path from a root bounds the speed-up, and more
  // threads than twice the bound only add overhead, and can put the columns
  // of the longest path on slower cores.
  double longest = 0.0;
  {
    vector<double> path(n);
    for (int j = n - 1; j >= 0; --j) {
      path[j] = work[j] + (parent[j] >= 0 ? path[parent[j]] : 0.0);
      longest = std::max(longest, path[j]);
    }
  }
  const int max_useful = (int) std::ceil(2.0 * total / std::max(longest, 1.0));
  const int nP = std::max(1, std::min(excursions_threads(n_threads), max_useful));
  // Subtrees with less work than this are computed by the thread that
  // computed their parent.
  const double grain = std::max(total / 4096.0, 20000.0);
  vector<QinvWork *> work_t(nP, NULL);
  for (int t = 0; t < nP; ++t)
    work_t[t] = new QinvWork(n);
  #pragma omp parallel num_threads(nP)
  {
    #pragma omp single
    {
      for (int j = n - 1; j >= 0; --j) {
        if (parent[j] < 0) {
          #pragma omp task firstprivate(j)
          qinv_subtree(F, j, cp, ci, sub, grain, work_t, variances);
        }
      }
    }
  }
  for (int t = 0; t < nP; ++t)
    delete work_t[t];

  if (selected)
    std::copy(z.begin(), z.end(), selected);
  return 0;
}

extern "C" SEXP Qinv(SEXP Rp, SEXP Ri, SEXP Rx, SEXP Rlower, SEXP Rthreads) {
  const int n = Rf_length(Rp) - 1;
  if (n < 0 || Rf_length(Ri) < INTEGER(Rp)[n] || Rf_length(Rx) < INTEGER(Rp)[n])
    Rf_error("Qinv: invalid sparse matrix.");
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  const int bad = n > 0 ? qinv_diag(n, INTEGER(Rp), INTEGER(Ri), REAL(Rx),
                                    Rf_asLogical(Rlower) == TRUE, REAL(out),
                                    NULL, Rf_asInteger(Rthreads))
                        : 0;
  UNPROTECT(1);
  if (bad)
    Rf_error("Qinv: the Cholesky factor has no diagonal element in column %d.", bad);
  return out;
}

// The selected inverse of Q = L L^T on the pattern of the lower triangular
// factor L, in the same order as L@x.
extern "C" SEXP Qinv_selected(SEXP Rp, SEXP Ri, SEXP Rx, SEXP Rthreads) {
  const int n = Rf_length(Rp) - 1;
  if (n < 0 || Rf_length(Ri) < INTEGER(Rp)[n] || Rf_length(Rx) < INTEGER(Rp)[n])
    Rf_error("Qinv_selected: invalid sparse matrix.");
  SEXP var = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n > 0 ? INTEGER(Rp)[n] : 0));
  const int bad = n > 0 ? qinv_diag(n, INTEGER(Rp), INTEGER(Ri), REAL(Rx),
                                    true, REAL(var), REAL(out),
                                    Rf_asInteger(Rthreads))
                        : 0;
  UNPROTECT(2);
  if (bad)
    Rf_error("Qinv_selected: the Cholesky factor has no diagonal element in column %d.", bad);
  return out;
}
