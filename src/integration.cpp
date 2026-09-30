#include <fcntl.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
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

/*
 Sequential importance sampling estimate of P(a < x < b) for x ~ N(0, Q^-1),
 with Q = R^T R and R upper triangular given in CSC format (Mp, Mi, Mv).
 Components are integrated from the last to the first, and Pv[i] is the
 probability for the components i, ..., n-1.

 The K samples are split into one contiguous block per thread, and each
 thread draws from its own random stream. The blocks are determined by the
 number of threads in the team, so for a given seed the results only depend
 on the number of threads that are used. The samples are stored with the
 sample index running fastest, x[i*K + j], so that the conditional means
 s[j] = sum_k R(i, k) x[k, j] are computed with unit stride over the block.
*/
extern "C" void shapeInt(int * Mp, int * Mi, double * Mv, double * a,double * b, int * opts, double * lim_in, double * Pv, double * Ev,int * seed_in){

  const int n = opts[0];
  const int K = opts[1];
  const int max_size = opts[2];
  const int n_threads = opts[3];
  const int seed_provided = opts[4];
  const double lim = lim_in[0];

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

  vector<double> al(n), bl(n), f(K, 1.0), s(K);
  vector<double> x(nrows * (size_t) K, 0.0);
  for (int i = 0; i < n; i++) {
    al[i] = Li[i]*a[i];
    bl[i] = Li[i]*b[i];
  }

  const int nP = setup_threads(n_threads);
  setup_seed(seed_provided, seed_in);

  vector<RngStream> RngArray(nP);
  for (int t = 0; t < nP; t++) {
    RngArray[t] = RngStream_CreateStream("namehere");
  }
  vector<double> fsum_t(nP), fsum2_t(nP);
  int team = 1;

  for (int i = n-1; i >= lo; i--) {
    double * xi = &x[(size_t) (i - lo) * K];
    const double ali = al[i];
    const double bli = bl[i];
    const double Lii = Li[i];

    #pragma omp parallel num_threads(nP)
    {
      int myrank = 0, nteam = 1;
      #ifdef _OPENMP
        myrank = omp_get_thread_num();
        nteam = omp_get_num_threads();
      #endif
      if (myrank == 0) {
        team = nteam;
      }
      const int j0 = (int) (((long long) K * myrank) / nteam);
      const int j1 = (int) (((long long) K * (myrank + 1)) / nteam);

      for (int j = j0; j < j1; j++) {
        s[j] = 0.0;
      }
      for (int q = rp[i]; q < rp[i + 1]; q++) {
        const double v = rv[q];
        const double * xk = &x[(size_t) (ri[q] - lo) * K];
        for (int j = j0; j < j1; j++) {
          s[j] += v*xk[j];
        }
      }

      double fsum = 0.0, fsum2 = 0.0;
      for (int j = j0; j < j1; j++) {
        double ai, bi, c, d, rtmp = 0;

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
          c = 0;
        }else if(ai>9){
          c = 1;
        }else {
          c = gsl_cdf_ugaussian_P(ai);
        }
        if (bi<-9) {
          d = 0;
        }else if(bi>9){
          d = 1;
        }else {
          d = gsl_cdf_ugaussian_P(bi);
        }

        f[j] = f[j]*(d-c);
        fsum += f[j];
        fsum2 += f[j]*f[j];

        if (d-c<1e-12) { //no weight is given to this sample
          xi[j] = 0; //just set x to zero
        } else {
          rtmp = c+(d-c)* RngStream_RandU01(RngArray[myrank]);
          xi[j] = (gsl_cdf_ugaussian_Pinv(rtmp)-s[j])/Lii;
        }

        if (xi[j] == numeric_limits<double>::infinity()){
          xi[j] = 0;
        }
      }
      fsum_t[myrank] = fsum;
      fsum2_t[myrank] = fsum2;
    }

    double fsum = 0.0, fsum2 = 0.0;
    for (int t = 0; t < team; t++) {
      fsum += fsum_t[t];
      fsum2 += fsum2_t[t];
    }

    const double Pi = fsum/K;
    if (Pi!=Pi) {
      Rprintf("%d Estimated probability is nan, stopping estimation\n",i);
      break;
    }
    const double Ei = sqrt(max((fsum2-fsum*fsum/K)/K/K,0));

    if (Pi<lim) {
      break;
    }

    Pv[i] = Pi;
    Ev[i] = Ei;
  }

  for (int t = 0; t < nP; t++) {
    RngStream_DeleteStream(&RngArray[t]);
  }
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
