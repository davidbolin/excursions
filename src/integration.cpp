#include <fcntl.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <memory>
#include <algorithm>
#include <limits>
#include <time.h>
#include "gsl_fix.h"
#include "threads.h"

/* Needed on Linux: */
#include <unistd.h>

#ifdef _OPENMP
#include<omp.h>
#endif

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

extern "C"{
  #include "RngStream.h"
}

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#define used_with_openmp(X) (void)X
using namespace std;

// See excursions_threads in threads.h.
static int setup_threads(int n_threads) {
  return excursions_threads(n_threads);
}

// Summary of the OpenMP support, for checking the installation: whether the
// package was compiled with OpenMP, and the default and maximal number of
// threads.
extern "C" SEXP excursions_openmp_info() {
  SEXP out = PROTECT(Rf_allocVector(INTSXP, 4));
  int *v = INTEGER(out);
  #ifdef _OPENMP
    v[0] = 1;
    v[1] = omp_get_max_threads();
    v[2] = omp_get_thread_limit();
    v[3] = omp_get_num_procs();
  #else
    v[0] = 0;
    v[1] = v[2] = v[3] = 1;
  #endif
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(names, 0, Rf_mkChar("openmp"));
  SET_STRING_ELT(names, 1, Rf_mkChar("max.threads"));
  SET_STRING_ELT(names, 2, Rf_mkChar("thread.limit"));
  SET_STRING_ELT(names, 3, Rf_mkChar("num.procs"));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

// Seed the RngStream package, from seed_in if provided and from R otherwise.
static void setup_seed(int seed_provided, int * seed_in) {
  unsigned long m_1 = 4294967087U;
  unsigned long m_2 = 4294944443U;
  unsigned long seed[6];

  if(seed_provided == 1){
    for(int i=0;i<6;i++){
      seed[i] = (unsigned long) seed_in[i];
    }
  } else {
    ssize_t seed_read = 0;
    #if defined (__APPLE__) && defined (__linux__)
      int randomSrc = open("/dev/urandom", O_RDONLY);
      if (randomSrc > 0) {
      seed_read = read(randomSrc, seed, sizeof(seed));
      close(randomSrc);
      }
    #endif
    if (seed_read != (ssize_t) sizeof(seed)) {
      GetRNGstate();
      for(int i=0;i<6;i++){
        seed[i] = round(RAND_MAX*unif_rand());
      }
      PutRNGstate();
    }
  }

  seed[0] = seed[0] % m_1;
  seed[1] = seed[1] % m_1;
  seed[2] = seed[2] % m_1;
  seed[3] = seed[3] % m_2;
  seed[4] = seed[4] % m_2;
  seed[5] = seed[5] % m_2;

  RngStream_SetPackageSeed(seed);
}

// Number of samples in each chunk of a batch with Kb samples, see shape_int.
// Small chunks give more chunks to distribute over the threads, and larger
// chunks less overhead for each row. It only depends on Kb, so that the
// results do not depend on the number of threads.
static int chunk_size(int Kb) {
  return (Kb < 2048) ? 32 : ((Kb < 4096) ? 64 : 128);
}

/*
 The rows of R in a group of rows g, ..., ge of the sequential integration,
 restricted to the columns k > g, which are the rows of x that are computed
 before the group. If these are dense enough, they are stored as a dense
 matrix W with the union U of their columns, so that their contributions to
 the conditional means of all rows of the group are computed together, see
 group_product. split[r] is the position in row i = g - r of R where the
 columns k > g start, and the columns i < k <= g before it are added row by
 row, since they are computed in the group.
*/
struct RowGroup {
  bool dense = false;
  int nU = 0;
  int Bp = 0;                // number of rows rounded up to a multiple of 4
  vector<int> split;
  vector<int> U;
  vector<size_t> off;        // offsets of the rows of U in the samples of a chunk
  vector<double> W;          // W[(r / 4) * nU * 4 + u * 4 + r % 4]
  vector<int> upos;          // position of a column in U, -1 if not in U

  void build(int g, int ge, const vector<int> & rp, const vector<int> & ri,
             const vector<double> & rv, int lo, int C) {
    const int B = g - ge + 1;
    split.resize(B);
    U.clear();
    size_t nnz = 0;
    for (int r = 0; r < B; r++) {
      const int i = g - r;
      int q = rp[i];
      while (q < rp[i + 1] && ri[q] <= g) {
        q++;
      }
      split[r] = q;
      for (; q < rp[i + 1]; q++) {
        if (upos[ri[q]] < 0) {
          upos[ri[q]] = 0;
          U.push_back(ri[q]);
        }
      }
      nnz += (size_t) (rp[i + 1] - split[r]);
    }
    nU = (int) U.size();
    // The dense product does B * nU multiplications for nnz entries, and is
    // several times faster for each multiplication than the sparse loop.
    dense = B >= 4 && nU > 0 && (double) nnz >= 0.3 * (double) B * (double) nU;
    if (dense) {
      sort(U.begin(), U.end());
      for (int u = 0; u < nU; u++) {
        upos[U[u]] = u;
      }
      Bp = (B + 3) / 4 * 4;
      W.assign((size_t) Bp * nU, 0.0);
      off.resize(nU);
      for (int u = 0; u < nU; u++) {
        off[u] = (size_t) (U[u] - lo) * C;
      }
      for (int r = 0; r < B; r++) {
        const int i = g - r;
        for (int q = split[r]; q < rp[i + 1]; q++) {
          W[(size_t) (r / 4) * nU * 4 + (size_t) upos[ri[q]] * 4 + r % 4] = rv[q];
        }
      }
    }
    for (int u = 0; u < nU; u++) {
      upos[U[u]] = -1;
    }
  }
};

// Vectors of two doubles, with the vector extensions of GCC and Clang, which
// are the SIMD registers of NEON and SSE2. The loads and stores do not
// assume alignment.
typedef double v2d __attribute__((vector_size(16)));
static inline v2d load2(const double * p) {
  v2d v;
  memcpy(&v, p, sizeof v);
  return v;
}
static inline void store2(double * p, v2d v) {
  memcpy(p, &v, sizeof v);
}

/*
 S[r * C + j] = sum_u W(r, u) x_{U[u]}[j] for the rows r < G.Bp of a group
 and the samples j < m of the chunk with samples xc, see RowGroup. Each block
 of four rows and eight samples is accumulated in sixteen vector registers,
 so that each sample that is loaded is used for four rows. The registers are
 written out explicitly, since compilers otherwise do not reliably keep the
 sums in registers. The columns u are taken in blocks of kc, so that the
 samples of a block are reused from the L1 cache for all rows, and the sums
 are kept in S between the blocks.
*/
static void group_product(const RowGroup & G, const double * xc, int C, int m,
                          double * __restrict S) {
  const int nU = G.nU;
  const int kc = 32;
  const size_t * off = G.off.data();
  for (int r = 0; r < G.Bp; r++) {
    for (int j = 0; j < m; j++) {
      S[(size_t) r * C + j] = 0.0;
    }
  }
  const int m8 = m / 8 * 8;
  for (int u0 = 0; u0 < nU; u0 += kc) {
    const int u1 = min(nU, u0 + kc);
    for (int jb = 0; jb < m8; jb += 8) {
      for (int rb = 0; rb < G.Bp; rb += 4) {
        const double * w = &G.W[(size_t) (rb / 4) * nU * 4];
        double * S0 = &S[(size_t) rb * C + jb];
        double * S1 = S0 + C;
        double * S2 = S1 + C;
        double * S3 = S2 + C;
        v2d a00 = load2(S0), a01 = load2(S0 + 2), a02 = load2(S0 + 4), a03 = load2(S0 + 6);
        v2d a10 = load2(S1), a11 = load2(S1 + 2), a12 = load2(S1 + 4), a13 = load2(S1 + 6);
        v2d a20 = load2(S2), a21 = load2(S2 + 2), a22 = load2(S2 + 4), a23 = load2(S2 + 6);
        v2d a30 = load2(S3), a31 = load2(S3 + 2), a32 = load2(S3 + 4), a33 = load2(S3 + 6);
        for (int u = u0; u < u1; u++) {
          const double * xu = xc + off[u] + jb;
          const v2d x0 = load2(xu), x1 = load2(xu + 2), x2 = load2(xu + 4), x3 = load2(xu + 6);
          const double w0 = w[4 * u], w1 = w[4 * u + 1], w2 = w[4 * u + 2], w3 = w[4 * u + 3];
          a00 += w0 * x0; a01 += w0 * x1; a02 += w0 * x2; a03 += w0 * x3;
          a10 += w1 * x0; a11 += w1 * x1; a12 += w1 * x2; a13 += w1 * x3;
          a20 += w2 * x0; a21 += w2 * x1; a22 += w2 * x2; a23 += w2 * x3;
          a30 += w3 * x0; a31 += w3 * x1; a32 += w3 * x2; a33 += w3 * x3;
        }
        store2(S0, a00); store2(S0 + 2, a01); store2(S0 + 4, a02); store2(S0 + 6, a03);
        store2(S1, a10); store2(S1 + 2, a11); store2(S1 + 4, a12); store2(S1 + 6, a13);
        store2(S2, a20); store2(S2 + 2, a21); store2(S2 + 4, a22); store2(S2 + 6, a23);
        store2(S3, a30); store2(S3 + 2, a31); store2(S3 + 4, a32); store2(S3 + 6, a33);
      }
    }
    // The samples after the last block of eight
    if (m8 < m) {
      for (int rb = 0; rb < G.Bp; rb += 4) {
        const double * w = &G.W[(size_t) (rb / 4) * nU * 4];
        for (int u = u0; u < u1; u++) {
          const double * xu = xc + off[u];
          for (int t = 0; t < 4; t++) {
            const double wt = w[4 * u + t];
            for (int j = m8; j < m; j++) {
              S[(size_t) (rb + t) * C + j] += wt * xu[j];
            }
          }
        }
      }
    }
  }
}

// Merge the sum Sb and the sum of squared deviations from the mean M2b of nb
// values into those of na values.
static inline void merge_moments(double & na, double & Sa, double & M2a,
                                 double nb, double Sb, double M2b) {
  if (nb <= 0) {
    return;
  }
  if (na <= 0) {
    na = nb;
    Sa = Sb;
    M2a = M2b;
    return;
  }
  const double delta = Sb / nb - Sa / na;
  M2a += M2b + delta * delta * na * nb / (na + nb);
  Sa += Sb;
  na += nb;
}

/*
 Sequential importance sampling for P(a < X < b), where X has precision
 Q = R^T R and R is upper triangular given in CSC format (Mp, Mi, Mv). The
 rows are integrated from the last to the first, and Pv[i] and Ev[i] are the
 estimate and its standard error for the rows i, ..., n-1.

 Rows with probe[i] != 0 are probes: they should have infinite limits a and
 b, so that they are sampled without constraint, and Pp[i] is then the
 estimate of the probability for rows i+1, ..., n-1 together with the probe
 limits pa[i] < X_i < pb[i], with standard error Pe[i]. Since the probe rows
 are sampled from their conditional distribution, the estimate for each probe
 only involves the constraints of the other rows. probe can be NULL.

 The samples are drawn in batches of independent samples, and the estimates
 are the averages over all samples of the batches that reached the row. With
 tol <= 0 there is a single batch of K samples. Otherwise the first batch has
 min(K, K0) samples, and further batches are added until the standard error
 at the target is at most tol, or K samples have been used. The target is
 the first row, from the end, where the estimate is below level, and the row
 before it, or the last computed row if there is no such row or level <= 0.
 Each batch is sized to reach tol from the current standard error, which
 decreases as one over the square root of the number of samples. A batch
 stops at the row where the estimate over all batches goes below lim. The
 number of samples that are used is returned in K_used.

 With size_tol > 0 instead of tol, and level > 0, the target is the
 standard error of the number k of rows, from the end, before the estimate
 goes below level, relative to k, which for excursions() is the size of the
 excursion set. It is the standard error of the estimate at the k-th row
 divided by the slope of the estimates there, per row, which is estimated
 from the max(3, k / 50) rows before it, since the rows after it may not be
 computed. The target is at least half a row, since the row where the
 estimate passes level is uncertain for any number of samples. If k = 0
 there is no error, since the estimate for the first row is exact, and if
 the error cannot be estimated, all K samples are used. The adaptive batches
 have at most max_batch samples, which bounds the memory for the samples.

 The samples of a batch are split into chunks of consecutive samples, where
 each chunk draws from its own random stream, and the chunks are distributed
 dynamically over the threads. The samples of a chunk are stored contiguously,
 x[c][i][j], so that the conditional means s[j] = sum_k R(i, k) x[c][k][j]
 are computed with unit stride and the rows that a chunk uses stay in the
 cache. The chunks only depend on the batch size, and the sums are added in
 the order of the chunks, so for a given seed the results do not depend on
 the number of threads.

 Since the samples of a chunk do not depend on the other chunks, the threads
 only synchronise after groups of rows, where the sums of the rows are added
 and the estimates are compared to lim. Rows after the row where the
 estimate goes below lim are computed but not used. The groups grow from one
 row to max_group rows, so that at most as many rows are wasted as have been
 used. For the rows of a group, the contributions of the rows before the
 group to the conditional means are computed together, see RowGroup, which
 changes the order of the sums. The groups do not depend on the number of
 threads, so the results do not either, and each batch uses new random
 streams.
*/
static void shape_int(int * Mp, int * Mi, double * Mv, double * a,double * b, int * opts, double * lim_in, double * Pv, double * Ev,int * seed_in,
                      const int * probe, const double * pa, const double * pb, double * Pp, double * Pe,
                      double tol, double level, int K0, int * K_used,
                      double size_tol = 0.0){

  const int n = opts[0];
  const int K = opts[1];
  const int max_size = opts[2];
  const int n_threads = opts[3];
  const int seed_provided = opts[4];
  const double lim = lim_in[0];

  if (K_used != NULL) {
    *K_used = 0;
  }
  if (n <= 0 || K <= 0) {
    return;
  }

  // Rows of R without the diagonal, in CSR format, and the diagonal.
  vector<int> rp(n + 1, 0), ri;
  vector<double> rv, Li(n, 0.0);
  for (int c = 0; c < n; c++) {
    for (int q = Mp[c]; q < Mp[c + 1]; q++) {
      if (Mi[q] == c) {
        Li[c] = Mv[q];
      } else {
        rp[Mi[q] + 1]++;
      }
    }
  }
  for (int i = 0; i < n; i++) {
    rp[i + 1] += rp[i];
  }
  ri.resize(rp[n]);
  rv.resize(rp[n]);
  {
    vector<int> next(rp.begin(), rp.end() - 1);
    for (int c = 0; c < n; c++) {
      for (int q = Mp[c]; q < Mp[c + 1]; q++) {
        const int r = Mi[q];
        if (r != c) {
          ri[next[r]] = c;
          rv[next[r]++] = Mv[q];
        }
      }
    }
  }

  // Only the rows lo, ..., n-1 are needed, see max_size below.
  const int lo = max(n - max_size, 0);
  const size_t nrows = (size_t) (n - lo);

  vector<double> al(n), bl(n);
  for (int i = 0; i < n; i++) {
    al[i] = Li[i]*a[i];
    bl[i] = Li[i]*b[i];
  }

  const int nP = setup_threads(n_threads);
  setup_seed(seed_provided, seed_in);

  // Sums over the samples of all batches, the sums of squared deviations
  // from the mean (S2 and PS2), and the number of samples, for each row.
  vector<double> S1(n, 0.0), S2(n, 0.0), PS1(n, 0.0), PS2(n, 0.0), N(n, 0.0);
  // The samples x are not initialised: row i only reads the rows k > i, which
  // are written first, and the memory of the rows that are never reached is
  // then not touched.
  std::unique_ptr<double[]> x;
  size_t x_size = 0;
  vector<double> f;
  vector<RngStream> rng;
  // Sums and sums of squared deviations for each chunk and row of the group.
  const int max_group = 32;
  vector<double> cs1, cs2, cp1, cp2;
  RowGroup grp;
  grp.upos.assign(n, -1);
  bool nan_found = false;
  int nan_row = 0;
  long long K_tot = 0;
  const bool adaptive = tol > 0 || (size_tol > 0 && level > 0);
  // The largest adaptive batch, see above
  const long long max_batch = 10000;
  int Kb = adaptive ? min(K, max(K0, 1)) : K;

  while (Kb > 0) {
    const int C = chunk_size(Kb);
    const int nC = (Kb + C - 1) / C;
    const size_t cstride = nrows * (size_t) C;
    if (cstride * (size_t) nC > x_size) {
      x_size = cstride * (size_t) nC;
      x.reset(new double[x_size]);
    }
    f.assign(Kb, 1.0);
    rng.resize(nC);
    for (int c = 0; c < nC; c++) {
      rng[c] = RngStream_CreateStream("chunk");
    }
    cs1.assign((size_t) nC * max_group, 0.0);
    cs2.assign((size_t) nC * max_group, 0.0);
    cp1.assign((size_t) nC * max_group, 0.0);
    cp2.assign((size_t) nC * max_group, 0.0);

    // The group is the rows g, ..., g - B + 1. g, B and done are only changed
    // in the single region, between barriers.
    int g = n - 1;
    int B = 1;
    bool done = false;
    grp.build(g, max(g - B + 1, lo), rp, ri, rv, lo, C);

    #pragma omp parallel num_threads(nP)
    {
      vector<double> s(C), pf(C);
      vector<double> Sg((size_t) max_group * C);
      while (!done && g >= lo) {
        const int ge = max(g - B + 1, lo);

        #pragma omp for schedule(dynamic, 1)
        for (int c = 0; c < nC; c++) {
          const int m = min(C, Kb - c * C);
          double * xc = &x[(size_t) c * cstride];
          double * fc = &f[(size_t) c * C];
          if (grp.dense) {
            group_product(grp, xc, C, m, Sg.data());
          }
          for (int i = g; i >= ge; i--) {
            double * xi = xc + (size_t) (i - lo) * C;
            const double ali = al[i];
            const double bli = bl[i];
            const double Lii = Li[i];
            const bool is_probe = probe != NULL && probe[i];
            const double pali = is_probe ? Lii * pa[i] : 0.0;
            const double pbli = is_probe ? Lii * pb[i] : 0.0;

            // The conditional means. For a dense group, the columns k > g
            // are in Sg, and the other columns are added row by row.
            double * __restrict sp = s.data();
            int q = rp[i];
            int q_end = rp[i + 1];
            if (grp.dense) {
              const double * Sr = &Sg[(size_t) (g - i) * C];
              for (int j = 0; j < m; j++) {
                sp[j] = Sr[j];
              }
              q_end = grp.split[g - i];
            } else {
              for (int j = 0; j < m; j++) {
                sp[j] = 0.0;
              }
            }
            // Four rows of x at a time, which reads and writes s a quarter
            // as often.
            for (; q + 3 < q_end; q += 4) {
              const double v0 = rv[q], v1 = rv[q + 1], v2 = rv[q + 2], v3 = rv[q + 3];
              const double * __restrict x0 = xc + (size_t) (ri[q] - lo) * C;
              const double * __restrict x1 = xc + (size_t) (ri[q + 1] - lo) * C;
              const double * __restrict x2 = xc + (size_t) (ri[q + 2] - lo) * C;
              const double * __restrict x3 = xc + (size_t) (ri[q + 3] - lo) * C;
              for (int j = 0; j < m; j++) {
                sp[j] += v0*x0[j] + v1*x1[j] + v2*x2[j] + v3*x3[j];
              }
            }
            for (; q < q_end; q++) {
              const double v = rv[q];
              const double * __restrict xk = xc + (size_t) (ri[q] - lo) * C;
              for (int j = 0; j < m; j++) {
                sp[j] += v*xk[j];
              }
            }

            double fsum = 0.0, psum = 0.0;
            for (int j = 0; j < m; j++) {
              double ai, bi, c_, d, rtmp = 0;

              if (is_probe) {
                const double pc = (pali == -numeric_limits<double>::infinity()) ? 0.0 :
                  gsl_cdf_ugaussian_P(pali + s[j]);
                const double pd = (pbli == numeric_limits<double>::infinity()) ? 1.0 :
                  gsl_cdf_ugaussian_P(pbli + s[j]);
                const double fp = fc[j] * max(pd - pc, 0.0);
                pf[j] = fp;
                psum += fp;
              }

              if (ali == -numeric_limits<double>::infinity()){
                ai = -numeric_limits<double>::infinity();
              } else {
                ai = ali + s[j];
              }

              if (bli == numeric_limits<double>::infinity()){
                bi = numeric_limits<double>::infinity();
              } else {
                bi = bli + s[j];
              }

              if (ai<-9) {
                c_ = 0;
              }else if(ai>9){
                c_ = 1;
              }else {
                c_ = gsl_cdf_ugaussian_P(ai);
              }
              if (bi<-9) {
                d = 0;
              }else if(bi>9){
                d = 1;
              }else {
                d = gsl_cdf_ugaussian_P(bi);
              }

              fc[j] = fc[j]*(d-c_);
              fsum += fc[j];

              if (d-c_<1e-12) { //no weight is given to this sample
                xi[j] = 0; //just set x to zero
              } else {
                rtmp = c_+(d-c_)* RngStream_RandU01(rng[c]);
                xi[j] = (gsl_cdf_ugaussian_Pinv(rtmp)-s[j])/Lii;
              }

              if (xi[j] == numeric_limits<double>::infinity()){
                xi[j] = 0;
              }
            }
            // The squared deviations from the mean of the chunk
            double fm2 = 0.0, pm2 = 0.0;
            const double fmean = fsum / m;
            for (int j = 0; j < m; j++) {
              fm2 += (fc[j] - fmean) * (fc[j] - fmean);
            }
            if (is_probe) {
              const double pmean = psum / m;
              for (int j = 0; j < m; j++) {
                pm2 += (pf[j] - pmean) * (pf[j] - pmean);
              }
            }
            const size_t k = (size_t) c * max_group + (g - i);
            cs1[k] = fsum;
            cs2[k] = fm2;
            cp1[k] = psum;
            cp2[k] = pm2;
          }
        }

        #pragma omp single
        {
          for (int i = g; i >= ge; i--) {
            double fn = 0.0, fsum = 0.0, fm2 = 0.0;
            double pn = 0.0, psum = 0.0, pm2 = 0.0;
            for (int c = 0; c < nC; c++) {
              const size_t k = (size_t) c * max_group + (g - i);
              const double mc = min(C, Kb - c * C);
              merge_moments(fn, fsum, fm2, mc, cs1[k], cs2[k]);
              merge_moments(pn, psum, pm2, mc, cp1[k], cp2[k]);
            }
            double Ni = N[i];
            merge_moments(Ni, S1[i], S2[i], fn, fsum, fm2);
            Ni = N[i];
            merge_moments(Ni, PS1[i], PS2[i], pn, psum, pm2);
            N[i] = Ni;

            const double Pi = S1[i]/N[i];
            if (Pi!=Pi) {
              nan_found = true;
              nan_row = i;
              done = true;
              break;
            }
            if (Pi<lim) {
              done = true;
              break;
            }
          }
          g = ge - 1;
          B = min(2 * B, max_group);
          if (!done && g >= lo) {
            grp.build(g, max(g - B + 1, lo), rp, ri, rv, lo, C);
          }
        }
      }
    }

    for (int c = 0; c < nC; c++) {
      RngStream_DeleteStream(&rng[c]);
    }
    K_tot += Kb;

    if (!adaptive || nan_found || K_tot >= K) {
      break;
    }

    // The error relative to the target, where more samples are needed if it
    // is larger than one.
    double ratio;
    if (tol > 0) {
      // Standard error at the target, from the end to the last computed row.
      double se = 0.0, se_prev = 0.0;
      for (int i = n-1; i >= lo && N[i] > 0; i--) {
        const double Pi = S1[i]/N[i];
        se_prev = se;
        se = sqrt(max(S2[i], 0.0)) / N[i];
        if (level > 0 && Pi < level) {
          if (i < n-1) {
            se = max(se, se_prev);
          }
          break;
        }
      }
      ratio = se / tol;
    } else {
      // Standard error of the number of rows before the estimate goes below
      // level, relative to size_tol times the number. Row n - p is the p-th
      // row from the end.
      int k = 0;
      for (int i = n-1; i >= lo && N[i] > 0 && S1[i] / N[i] >= level; i--) {
        k++;
      }
      if (k == 0) {
        ratio = 0.0;
      } else {
        const int w = max(3, (int) round(0.02 * k));
        const int top = max(1, k - w);
        const int ik = n - k, it = n - top;
        const double slope = (top < k) ? (S1[it] / N[it] - S1[ik] / N[ik]) / (k - top) : 0.0;
        const double se_k = sqrt(max(S2[ik], 0.0)) / N[ik];
        ratio = (slope > 0) ? se_k / slope / max(size_tol * k, 0.5)
                            : numeric_limits<double>::infinity();
      }
    }
    if (ratio <= 1.0) {
      break;
    }
    const double K_need = (ratio < 1e6) ? 1.1 * (double) K_tot * ratio * ratio : (double) K;
    Kb = (int) min((double) min(K - K_tot, max_batch),
                   max((double) K0, ceil(K_need - (double) K_tot)));
  }

  if (K_used != NULL) {
    *K_used = (int) K_tot;
  }

  for (int i = n-1; i >= lo && N[i] > 0; i--) {
    if (probe != NULL && probe[i]) {
      Pp[i] = PS1[i] / N[i];
      Pe[i] = sqrt(max(PS2[i], 0.0)) / N[i];
    }
    if (nan_found && i == nan_row) {
      Rprintf("%d Estimated probability is nan, stopping estimation\n",i);
      break;
    }
    const double Pi = S1[i]/N[i];
    if (Pi<lim) {
      break;
    }
    Pv[i] = Pi;
    Ev[i] = sqrt(max(S2[i], 0.0)) / N[i];
  }

}

extern "C" void shapeInt(int * Mp, int * Mi, double * Mv, double * a,double * b, int * opts, double * lim_in, double * Pv, double * Ev,int * seed_in){
  shape_int(Mp, Mi, Mv, a, b, opts, lim_in, Pv, Ev, seed_in, NULL, NULL, NULL, NULL, NULL,
            0.0, 0.0, 0, NULL);
}

/*
 shapeInt with adaptive number of samples, adapt = (tol, level, K0) or
 adapt = (tol, level, K0, size_tol), see shape_int. It is called with .Call,
 which does not copy the arguments, so that the Cholesky factor is not
 copied. Returns list(Pv, Ev, K_used), where K_used is the number of samples
 that are used.
*/
extern "C" SEXP shapeIntCall(SEXP Mp, SEXP Mi, SEXP Mv, SEXP a, SEXP b, SEXP opts, SEXP lim, SEXP seed_in, SEXP adapt){
  const int n = INTEGER(opts)[0];
  if (Rf_length(a) != n || Rf_length(b) != n || Rf_length(Mp) != n + 1 ||
      Rf_length(Mi) != Rf_length(Mv) || Rf_length(opts) < 5 ||
      Rf_length(seed_in) < 6 || Rf_length(adapt) < 3 || Rf_length(lim) < 1) {
    Rf_error("shapeIntCall: arguments of wrong lengths");
  }
  SEXP Pv = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP Ev = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP K_used = PROTECT(Rf_allocVector(INTSXP, 1));
  for (int i = 0; i < n; i++) {
    REAL(Pv)[i] = 0.0;
    REAL(Ev)[i] = 0.0;
  }
  // shape_int only reads the inputs.
  shape_int(INTEGER(Mp), INTEGER(Mi), REAL(Mv), REAL(a), REAL(b), INTEGER(opts),
            REAL(lim), REAL(Pv), REAL(Ev), INTEGER(seed_in), NULL, NULL, NULL, NULL, NULL,
            REAL(adapt)[0], REAL(adapt)[1], (int) REAL(adapt)[2], INTEGER(K_used),
            Rf_length(adapt) >= 4 ? REAL(adapt)[3] : 0.0);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, Pv);
  SET_VECTOR_ELT(out, 1, Ev);
  SET_VECTOR_ELT(out, 2, K_used);
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(names, 0, Rf_mkChar("Pv"));
  SET_STRING_ELT(names, 1, Rf_mkChar("Ev"));
  SET_STRING_ELT(names, 2, Rf_mkChar("K_used"));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(5);
  return out;
}

// shapeInt with probes, see shape_int, and with an adaptive number of
// samples as in shapeIntCall, with adapt = (tol, level, K0, size_tol), where
// K_used is the number of samples that are used.
extern "C" void shapeIntProbe(int * Mp, int * Mi, double * Mv, double * a,double * b, int * opts, double * lim_in, double * Pv, double * Ev,int * seed_in,
                              int * probe, double * pa, double * pb, double * Pp, double * Pe,
                              double * adapt, int * K_used){
  shape_int(Mp, Mi, Mv, a, b, opts, lim_in, Pv, Ev, seed_in, probe, pa, pb, Pp, Pe,
            adapt[0], adapt[1], (int) adapt[2], K_used, adapt[3]);
}


extern "C" void testRand( int * opts, double * x, int * seed_in){

  const int n = opts[0];
  const int nP = setup_threads(opts[1]);
  setup_seed(opts[2], seed_in);

  vector<RngStream> RngArray(nP);
  for (int t = 0; t < nP; t++) {
    RngArray[t] = RngStream_CreateStream("namehere");
  }

  #pragma omp parallel num_threads(nP)
  {
    int myrank = 0;
  #ifdef _OPENMP
    myrank = omp_get_thread_num();
  #endif

    #pragma omp for
    for (int i=0; i<n; i++) {
      x[i] = RngStream_RandU01(RngArray[myrank]);
    }
  }

  for (int t = 0; t < nP; t++) {
    RngStream_DeleteStream(&RngArray[t]);
  }
}
