/* utils.cpp
 *
 *   Copyright (C) 2012, 2013 David Bolin, Finn Lindgren
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "utils.h"



void readSparseMatrix(vector< map<int,double> > * matrix, int * m_size, char* ifile, char* vfile){
    int i,j;
    FILE * pFile;
    long lSize;
    int * indexV;
    double * valV;
    size_t result;
	
    /*read index vector*/
	
    pFile = fopen(ifile, "rb");
    if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	
    // obtain file size:
    fseek (pFile , 0 , SEEK_END);
    lSize = ftell (pFile);
    rewind (pFile);
    long nbrelements = lSize/sizeof(int);
    
    // allocate memory to contain the whole file:
    indexV = (int*) malloc (lSize);
    if (indexV == NULL) {fputs ("Memory error",stderr); exit (2);}
	
    // copy the file into the buffer:
	
    result = fread (indexV,1,lSize,pFile);
    if ((long)result != lSize) {fputs ("Reading error",stderr); exit (3);}
    fclose (pFile);
	
    /*read values*/
	
    pFile = fopen(vfile, "rb");
    if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	
    // obtain file size:
    fseek (pFile , 0 , SEEK_END);
    lSize = ftell (pFile);
    rewind (pFile);
    // allocate memory to contain the whole file:
    valV = (double*) malloc (lSize);
    if (valV == NULL) {fputs ("Memory error",stderr); exit (2);}
	
    // copy the file into the buffer:
    result = fread (valV,1,lSize,pFile);
    if ((long)result != lSize) {fputs ("Reading error",stderr); exit (3);}
    fclose (pFile);
	
    /* values are now located in valV and corresponding indices in indexV */
    j=0;
    int row,col;
    double val;
    int current = -1;
	
    for(i=0;i<nbrelements;i+=2){
		
		row = indexV[i];   //node number
		col = indexV[i+1]; //neighbor number
		val = valV[j];     //value
		while(current < row){
			current += 1;
			map<int,double> m;
			(*matrix).push_back(m);
		}
		(*matrix)[current][col] = val;
		j++;
    }
	m_size[1] = nbrelements/2;
	m_size[0] = (*matrix).size();
    free (indexV);
    free(valV);
} 

void convert_to_ccs(vector< map<int,double> > * matrix, cholmod_sparse ** A,int n,int nz,cholmod_common * cm){
	// convert to ccs format and only save lower triangular part

	*A = cholmod_allocate_sparse(n,n,nz,1,1,-1,CHOLMOD_REAL,cm);
	(**A).nzmax = nz;
	int * Ap = (int*) (*A)->p; 
	int * Ai = (int*) (*A)->i;
	double * Ax = (double*) (*A)->x;
	int pos = 0;
	for (int i = 0; i < n; ++i) {
		Ap[i] = pos;
		for (map<int,double>::const_iterator iter = ((*matrix)[i]).begin(); iter != ((*matrix)[i]).end(); ++iter) {
			Ai[pos] = (*iter).first;
			Ax[pos] = (*iter).second;
			//cout << i << " " << (*iter).first << " " << (*iter).second << endl;
			++pos;
		}
	}
	Ap[n] = nz; //is this needed?
 }



void convert_from_ccs(vector< map<int,double> > * matrix, cholmod_sparse ** A,cholmod_common * cm){

	//convert and transpose!
	
    int i,row,col;

	// build vector
    for (i=0; i<(int)(**A).nrow; i++) {
		map<int,double> m;
		(*matrix).push_back(m);
	}
	
	//add elements
	col = 0;
	for(i=0;i<((int*) (**A).p)[(int)(**A).ncol];i++){
		
		row = ((int*)(**A).i)[i]; 
		
		if (i>=((int*) (**A).p)[col+1]) {
			col++;
		} 
		//(*matrix)[row][col] = (**A).values.d[i]; //no transpose
		(*matrix)[col][row] = ((double*)(**A).x)[i];
		
    }
}


void sort_vector(int n, double v[], double vs[], int ind[]){

	vec_ind * vi = (vec_ind *) malloc(sizeof(vec_ind) * n);
	int i;
	for(i = 0; i< n;i++){
		vi[i].val = v[i];
		vi[i].ind = i;
	}
	
	qsort(vi, n, sizeof(vi[0]), compare_structs);
	
	for(i = 0; i< n;i++){
		vs[i] = vi[i].val;
		ind[i] = vi[i].ind;
	}
	
	free(vi);
}

int compare_structs(const void *a, const void *b){
	
    vec_ind *struct_a = (vec_ind *) a;
    vec_ind *struct_b = (vec_ind *) b;
	
    if (struct_a->val < struct_b->val) return -1;
    else if (struct_a->val == struct_b->val) return 0;
    else return 1;
}

void Qinv(int * Rir, int * Rjc, double * Rpr, double * variances, int n){	
	typedef pair<size_t,double> Qpairtype;
	typedef vector< Qpairtype > Qvectype;
	typedef vector< Qvectype > Qtype;
				
	int i,j;
	
	/*
	 * Copy cholesky factor to more convenient format and extract diagonal elements
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
	
	/*
	 * Calculate inverse
	 */
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

void output_matrix(int m, int n, vector< map<int,double> > * matrix){

	for (int i=0; i<n; i++) {
		for (int j=0; j<m; j++) {
			cout << (*matrix)[i][j] << " ";
		}
		cout << endl;
	}

}

void vec_permute(int n, double * a,double * a_sort,int * reo){
	for(int i=0;i<n;i++){
		a_sort[i] = a[reo[i]];
	}
}
