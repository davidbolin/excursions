#include <vector>
#include <algorithm>

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
 be a search per access into an array index. This replaces the earlier
 vector-of-vectors version, which searched sorted rows for every product.
*/
// Returns 0 on success, or the 1-based column with a missing diagonal.
static int qinv_diag(int n, const int *Ap, const int *Ai, const double *Ax,
                     bool lower, double *variances) {
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

  vector<double> z(nnz, 0.0);
  vector<int> pos(n, -1);
  vector<double> lj(n, 0.0);
  vector<double> acc;

  for (int j = n - 1; j >= 0; --j) {
    const int start = Lp[j], end = Lp[j + 1];
    if (start >= end || Li[start] != j)
      return j + 1;
    const double Ljj = Lx[start];
    const int m = end - start - 1;

    for (int t = 0; t < m; ++t) {
      const int r = Li[start + 1 + t];
      pos[r] = t;
      lj[r] = Lx[start + 1 + t] / Ljj;
    }
    acc.assign(m, 0.0);

    // acc[a] accumulates -sum_k Ltil(k, j) Z(r_a, k) over the rows k of this
    // column. Walking column c yields Z(r, c) for r >= c, and that value
    // serves both the (i = r, k = c) and the (i = c, k = r) term, so one pass
    // over the already computed columns covers every pair without a search.
    for (int b = 0; b < m; ++b) {
      const int c = Li[start + 1 + b];
      const double lc = lj[c];
      for (int q = Lp[c]; q < Lp[c + 1]; ++q) {
        const int a = pos[Li[q]];
        if (a < 0)
          continue;
        const double v = z[q];
        acc[a] -= lc * v;
        if (Li[q] != c)
          acc[b] -= lj[Li[q]] * v;
      }
    }

    double accd = 0.0;
    for (int t = 0; t < m; ++t) {
      z[start + 1 + t] = acc[t];
      accd += lj[Li[start + 1 + t]] * acc[t];
    }
    z[start] = 1.0 / (Ljj * Ljj) - accd;
    variances[j] = z[start];

    for (int t = 0; t < m; ++t) {
      const int r = Li[start + 1 + t];
      pos[r] = -1;
      lj[r] = 0.0;
    }
  }
  return 0;
}

extern "C" SEXP Qinv(SEXP Rp, SEXP Ri, SEXP Rx, SEXP Rlower) {
  const int n = Rf_length(Rp) - 1;
  if (n < 0 || Rf_length(Ri) < INTEGER(Rp)[n] || Rf_length(Rx) < INTEGER(Rp)[n])
    Rf_error("Qinv: invalid sparse matrix.");
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  const int bad = n > 0 ? qinv_diag(n, INTEGER(Rp), INTEGER(Ri), REAL(Rx),
                                    Rf_asLogical(Rlower) == TRUE, REAL(out))
                        : 0;
  UNPROTECT(1);
  if (bad)
    Rf_error("Qinv: the Cholesky factor has no diagonal element in column %d.", bad);
  return out;
}
