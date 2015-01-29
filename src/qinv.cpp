#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <map>
#include <algorithm>


#ifdef _OPENMP
	#include<omp.h>
#endif

using namespace std;

extern "C" void Qinv(int * Rir, int * Rjc, double * Rpr, double * variances, int * nin, int n_threads){
  int n = nin[0];

  typedef pair<size_t,double> Qpairtype;
  typedef vector< Qpairtype > Qvectype;
  typedef vector< Qvectype > Qtype;

  int i,j;

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

  /*
  Copy cholesky factor to more convenient format and extract diagonal elements
  */

  Qtype R(n);
  vector<double> D(n);
  //Extract the elements and store the sparse R-matrix
  //in a more convinient format.
  if(Rjc[n]-Rjc[n-1] == 1){
    //only one element in the last column, assume lower triangular matrix
    for(int c=0;c<n;++c){
      D[c] = Rpr[Rjc[c]];
      R[c].resize(Rjc[c+1]-Rjc[c]);
      for(j=Rjc[c],i=0;j<Rjc[c+1];++j,++i)
        R[c][i] = Qpairtype(Rir[j], Rpr[j]);
    }
  }else{
    //assume upper triangular matrix - first find number of element in each row
    vector<size_t> nRow(n), iRow(n);
    for(int c=0;c<n;++c){
      for(j=Rjc[c];j<Rjc[c+1];++j)
        ++nRow[Rir[j]];
      D[c] = Rpr[Rjc[c+1]-1];
    }
    for(int c=0;c<n;++c)
      R[c].resize( nRow[c] );
    for(int c=0;c<n;++c){
      for(j=Rjc[c];j<Rjc[c+1];++j)
        R[Rir[j]][iRow[Rir[j]]++] = Qpairtype(c, Rpr[j]);
    }
  }

  /* Calculate inverse */
  Qvectype::iterator pos;
  size_t Nmax=0;
  //divide all elemnts in R by the diagonal-elements
  for(i=0; i<n; ++i){
    //find the maximal number of non-zero elements in any row of R
    if(Nmax < R[i].size())
      Nmax = R[i].size();
    //compute R[i,j]/D[i]
    for(pos=R[i].begin(); pos!=R[i].end(); ++pos)
      (pos->second) /=D[i];
    //and compute 1/d^2
    D[i] = 1/(D[i]*D[i]);
  }

  //count number of elements that is going to end up in iQ
  vector<size_t> nnz(n,1);
  for(i=0; i<n; ++i){
    //first find the indices of the non-zero elements
    for(pos=R[i].begin(), ++pos; pos!=R[i].end(); ++pos){
      nnz[i]++;
      nnz[pos->first]++;
    }
  }

  //vectors containing the location and values within one column
  vector<size_t> ii(Nmax);
  vector<double> s(Nmax);
  vector< Qvectype::iterator > iQpos(Nmax);
  vector< Qvectype::iterator > iQstart(n);

  //create a structure holding the inverse matrix
  Qtype iQ(n);
  for(i=0; i<n; ++i){
    iQ[i].resize(nnz[i]);
    iQstart[i] = iQ[i].end();
  }

  //loop over the columns of the matrix
  i = n;
  while(i>0){
    --i;
    //first find the indices of the non-zero elements
    for(pos=R[i].begin(), ++pos, j=0; pos!=R[i].end(); ++pos, j++){
      ii[j] = pos->first; //index of elements
      s[j] = 0; //set values to zero
      iQpos[j] = iQstart[ii[j]]; //start of each iQ row
    }

    //multiply the row of R with the rows of iQ
    #pragma omp parallel for private(pos)
    for(int j2=0; j2<(R[i].size()-1); ++j2){
      Qvectype::iterator iQpos_tmp = iQpos[j2];
      Qvectype::iterator iQend = iQ[ii[j2]].end();
      for(pos=R[i].begin(), ++pos; pos!=R[i].end(); ++pos){
        for(;iQpos_tmp != iQend && iQpos_tmp->first < pos->first; ++iQpos_tmp){}
        if(iQpos_tmp != iQend && iQpos_tmp->first == pos->first)
          s[j2] += (iQpos_tmp->second) * (pos->second);
      }
    }

    //the diagonal elements
    double diag = D[i];
    for(pos=R[i].begin(), ++pos, j=0; pos!=R[i].end(); ++pos, ++j)
      diag += s[j] * (pos->second);

    //add the elements to iQ
    j = R[i].size()-1;
    while(j>0){
      --j;
      *(--iQstart[i]) = Qpairtype(ii[j], -s[j]);
      *(--iQstart[ ii[j] ]) = Qpairtype(i, -s[j]);
    }
    *(--iQstart[i]) = Qpairtype(i, diag);
  }
  for (i=0; i<n; i++) {
    for(pos=iQ[i].begin(); pos!=iQ[i].end(); ++pos){
      if (i==(int)pos->first) {
        variances[i] = pos->second;
      }
    }
  }
}