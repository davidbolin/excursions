#include <fcntl.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <map>
#include <algorithm>
#include <limits>
#include <time.h>
#include "gsl_fix.h"

/* Needed on Linux: */
#include <unistd.h>

#ifdef _OPENMP
#include<omp.h>
#endif

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>

extern "C"{
  #include "RngStream.h"
}

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

using namespace std;

extern "C" void shapeInt(int * Mp, int * Mi, double * Mv, double * a,double * b, int * opts, double * lim_in, double * Pv, double * Ev,int * seed_in){

  int n = opts[0];
  int K = opts[1];
  int max_size = opts[2];
  int n_threads = opts[3];
  int seed_provided = opts[4];

  double lim = lim_in[0];

  vector< map<int,double> > L;

  int i,row,col;


  for (i=0; i<n; i++) {
    map<int,double> m;
    L.push_back(m);
  }

  col = 0;
  for(i=0;i<Mp[n];i++){
    row = Mi[i];
    if (i>=Mp[col+1]) {
      col++;
    }
    L[row][col] = Mv[i];
  }



  double *al, *bl,*f,*s,*Li;
  double ** x;
  double fsum,fsum2,Pi,Ei,ai,bi,c,d,rtmp;
  int j;

  al = new double[n];
  bl = new double[n];
  f = new double[K];
  s = new double[K];
  Li = new double[n];

  x = new double*[n];
  for (i=0; i<n; i++) {
    Li[i] = L[i][i];
    al[i] = Li[i]*a[i];
    bl[i] = Li[i]*b[i];
    x[i] = new double[K];
    for(j=0;j<K;j++){
      x[i][j] = 0.0;
    }
  }

  for (i=0; i<K; i++) {
    f[i] = 1.0;
  }


  #ifdef _OPENMP
    const int max_nP = omp_get_num_procs();
    int nPtmp;
    if(n_threads == 0){
      nPtmp = max_nP;
    } else {
      nPtmp = min(max_nP, max(n_threads,1));
    }
    const int nP = nPtmp;
    omp_set_num_threads(nP);
  #else
    const int nP = 1;
  #endif


  unsigned long m_1 = 4294967087U;
  unsigned long m_2 = 4294944443U;
  unsigned long seed[6];

  if(seed_provided == 1){
    for(i=0;i<6;i++){
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
      //srand(time(0));
      GetRNGstate();
      for(i=0;i<6;i++){
        //seed[i] = rand();
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
  RngStream * RngArray = new RngStream[nP];
  int myrank = 0;

  for (i=0; i<nP; i++) {
    RngArray[i] = RngStream_CreateStream("namehere");
  }

  for (i=n-1; i>=0; i--) {
    for (j=0; j<K; j++) {
      s[j] = 0;
    }
    fsum = 0;
    fsum2 = 0;

    for(map<int,double>::iterator iter = (L[i]).begin(); iter != (L[i]).end(); iter++ ) {
      for (j=0; j<K; j++) {
        s[j] += (iter->second)*x[iter->first][j];
      }
    }

    #pragma omp parallel private(myrank,ai,bi,c,d,j,rtmp)
    {
      #ifdef _OPENMP
        myrank = omp_get_thread_num();
      #endif

      #pragma omp for reduction(+:fsum,fsum2)
      for (j=0; j<K; j++) {
        ai = al[i] + s[j];
        bi = bl[i] + s[j];

        if (al[i] == -numeric_limits<double>::infinity()){
          ai = -numeric_limits<double>::infinity();
        } else {
          ai = al[i] + s[j];
        }

        if (bl[i] == numeric_limits<double>::infinity()){
          bi = numeric_limits<double>::infinity();
        } else {
          bi = bl[i] + s[j];
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
          rtmp = 0;
          x[i][j] = 0; //just set x to zero
        } else {
          rtmp = c+(d-c)* RngStream_RandU01(RngArray[myrank]);
          x[i][j] = (gsl_cdf_ugaussian_Pinv(rtmp)-s[j])/Li[i];
        }

        if (x[i][j] == numeric_limits<double>::infinity()){
          Rprintf("simulated infinite value, changing to zero\n",i);
          Rprintf("c= %f, d= %f, ,d-c= %f\n",c,d,d-c);
          x[i][j] = 0;
        }

        if (x[i][j]!=x[i][j]) {
          Rprintf("%d x is nan: rtmp= %f, c= %f, d= %f, ai=%f",rtmp,c,d,ai);
          Rprintf(", bi= %f, s[j] = %f, Li = %f",bi,s[j],Li[i]);
        }
      }
    }

    Pi = fsum/K;
    if (Pi!=Pi) {
      Rprintf("%d Estimated probability is nan, stopping estimation\n",i);
      break;
    }
    Ei = sqrt(max((fsum2-fsum*fsum/K)/K/K,0));

    if (Pi<lim) {
      break;
    }

    if(i<n-max_size){
      break;
    }

    Pv[i] = Pi;
    Ev[i] = Ei;
  }

  delete[] al;
  delete[] bl;
  delete[] f;
  delete[] s;
  delete[] Li;
  for (i=0; i<n; i++) {
    delete[] x[i];
  }
  delete[] x;
  delete[] RngArray;

}

extern "C" void testRand( int * opts, double * x, int * seed_in){

  int n = opts[0];
  int n_threads = opts[1];
  int seed_provided = opts[2];

  #ifdef _OPENMP
    const int max_nP = omp_get_num_procs();
    int nPtmp;
    if(n_threads == 0){
      nPtmp = max_nP;
    } else {
      nPtmp = min(max_nP, max(n_threads,1));
    }
    const int nP = nPtmp;
    omp_set_num_threads(nP);
  #else
    const int nP = 1;
  #endif


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
      //srand(time(0));
      GetRNGstate();
      for(int i=0;i<6;i++){
        //seed[i] = rand();
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
  RngStream * RngArray = new RngStream[nP];
  int myrank = 0;

  for (int i=0; i<nP; i++) {
    RngArray[i] = RngStream_CreateStream("namehere");
  }

  #pragma omp parallel private(myrank)
  {
  #ifdef _OPENMP
    myrank = omp_get_thread_num();
  #endif

    #pragma omp for
    for (int i=0; i<n; i++) {
      x[i] = RngStream_RandU01(RngArray[myrank]);
    }
  }
  delete[] RngArray;
}
