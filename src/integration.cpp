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

// Number of threads to request for the parallel regions. With n_threads = 0
// this is the OpenMP default, which respects OMP_NUM_THREADS, and otherwise
// n_threads, at most the number of processors. In both cases it is at most
// OMP_THREAD_LIMIT. The runtime may still give fewer threads, so the parallel
// regions must use the size of the team they get. The number of threads is
// requested with num_threads() rather than omp_set_num_threads(), which would
// change the default for later calls and for other packages.
static int setup_threads(int n_threads) {
  used_with_openmp(n_threads);
  #ifdef _OPENMP
    int nP = (n_threads <= 0) ? omp_get_max_threads() : min(n_threads, omp_get_num_procs());
    nP = min(nP, omp_get_thread_limit());
    return max(nP, 1);
  #else
    return 1;
  #endif
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
 estimate goes below lim are computed but not used. The groups have one row
 for one thread, and grow from one row to max_group rows otherwise. The
 group sizes do not change the results, since the rows that are used are the
 same, and each batch uses new random streams.
*/
static void shape_int(int * Mp, int * Mi, double * Mv, double * a,double * b, int * opts, double * lim_in, double * Pv, double * Ev,int * seed_in,
                      const int * probe, const double * pa, const double * pb, double * Pp, double * Pe,
                      double tol, double level, int K0, int * K_used){

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

  // Sums over the samples of all batches, and the number of samples, for
  // each row.
  vector<double> S1(n, 0.0), S2(n, 0.0), PS1(n, 0.0), PS2(n, 0.0), N(n, 0.0);
  // The samples x are not initialised: row i only reads the rows k > i, which
  // are written first, and the memory of the rows that are never reached is
  // then not touched.
  std::unique_ptr<double[]> x;
  size_t x_size = 0;
  vector<double> f;
  vector<RngStream> rng;
  // Sums for each chunk and row of the group.
  const int max_group = 32;
  vector<double> cs1, cs2, cp1, cp2;
  bool nan_found = false;
  int nan_row = 0;
  long long K_tot = 0;
  const bool adaptive = tol > 0;
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

    #pragma omp parallel num_threads(nP)
    {
      vector<double> s(C);
      while (!done && g >= lo) {
        const int ge = max(g - B + 1, lo);

        #pragma omp for schedule(dynamic, 1)
        for (int c = 0; c < nC; c++) {
          const int m = min(C, Kb - c * C);
          double * xc = &x[(size_t) c * cstride];
          double * fc = &f[(size_t) c * C];
          for (int i = g; i >= ge; i--) {
            double * xi = xc + (size_t) (i - lo) * C;
            const double ali = al[i];
            const double bli = bl[i];
            const double Lii = Li[i];
            const bool is_probe = probe != NULL && probe[i];
            const double pali = is_probe ? Lii * pa[i] : 0.0;
            const double pbli = is_probe ? Lii * pb[i] : 0.0;

            // The conditional means, four rows of x at a time, which reads
            // and writes s a quarter as often.
            double * __restrict sp = s.data();
            for (int j = 0; j < m; j++) {
              sp[j] = 0.0;
            }
            int q = rp[i];
            const int q_end = rp[i + 1];
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

            double fsum = 0.0, fsum2 = 0.0, psum = 0.0, psum2 = 0.0;
            for (int j = 0; j < m; j++) {
              double ai, bi, c_, d, rtmp = 0;

              if (is_probe) {
                const double pc = (pali == -numeric_limits<double>::infinity()) ? 0.0 :
                  gsl_cdf_ugaussian_P(pali + s[j]);
                const double pd = (pbli == numeric_limits<double>::infinity()) ? 1.0 :
                  gsl_cdf_ugaussian_P(pbli + s[j]);
                const double fp = fc[j] * max(pd - pc, 0.0);
                psum += fp;
                psum2 += fp * fp;
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
              fsum2 += fc[j]*fc[j];

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
            const size_t k = (size_t) c * max_group + (g - i);
            cs1[k] = fsum;
            cs2[k] = fsum2;
            cp1[k] = psum;
            cp2[k] = psum2;
          }
        }

        #pragma omp single
        {
          for (int i = g; i >= ge; i--) {
            double fsum = 0.0, fsum2 = 0.0, psum = 0.0, psum2 = 0.0;
            for (int c = 0; c < nC; c++) {
              const size_t k = (size_t) c * max_group + (g - i);
              fsum += cs1[k];
              fsum2 += cs2[k];
              psum += cp1[k];
              psum2 += cp2[k];
            }
            S1[i] += fsum;
            S2[i] += fsum2;
            PS1[i] += psum;
            PS2[i] += psum2;
            N[i] += Kb;

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
          if (nP > 1) {
            B = min(2 * B, max_group);
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

    // Standard error at the target, from the end to the last computed row.
    double se = 0.0, se_prev = 0.0;
    for (int i = n-1; i >= lo && N[i] > 0; i--) {
      const double Pi = S1[i]/N[i];
      se_prev = se;
      se = sqrt(max((S2[i]-S1[i]*S1[i]/N[i])/N[i]/N[i],0));
      if (level > 0 && Pi < level) {
        if (i < n-1) {
          se = max(se, se_prev);
        }
        break;
      }
    }
    if (se <= tol) {
      break;
    }
    const double K_need = 1.1 * (double) K_tot * (se / tol) * (se / tol);
    Kb = (int) min((double) (K - K_tot), max((double) K0, ceil(K_need - (double) K_tot)));
  }

  if (K_used != NULL) {
    *K_used = (int) K_tot;
  }

  for (int i = n-1; i >= lo && N[i] > 0; i--) {
    if (probe != NULL && probe[i]) {
      Pp[i] = PS1[i] / N[i];
      Pe[i] = sqrt(max((PS2[i] - PS1[i] * PS1[i] / N[i]) / N[i] / N[i], 0.0));
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
    Ev[i] = sqrt(max((S2[i]-S1[i]*S1[i]/N[i])/N[i]/N[i],0));
  }

}

extern "C" void shapeInt(int * Mp, int * Mi, double * Mv, double * a,double * b, int * opts, double * lim_in, double * Pv, double * Ev,int * seed_in){
  shape_int(Mp, Mi, Mv, a, b, opts, lim_in, Pv, Ev, seed_in, NULL, NULL, NULL, NULL, NULL,
            0.0, 0.0, 0, NULL);
}

/*
 shapeInt with adaptive number of samples, adapt = (tol, level, K0), see
 shape_int. It is called with .Call, which does not copy the arguments, so that
 the Cholesky factor is not copied. Returns list(Pv, Ev, K_used), where
 K_used is the number of samples that are used.
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
            REAL(adapt)[0], REAL(adapt)[1], (int) REAL(adapt)[2], INTEGER(K_used));
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

extern "C" void shapeIntProbe(int * Mp, int * Mi, double * Mv, double * a,double * b, int * opts, double * lim_in, double * Pv, double * Ev,int * seed_in,
                              int * probe, double * pa, double * pb, double * Pp, double * Pe){
  shape_int(Mp, Mi, Mv, a, b, opts, lim_in, Pv, Ev, seed_in, probe, pa, pb, Pp, Pe,
            0.0, 0.0, 0, NULL);
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
