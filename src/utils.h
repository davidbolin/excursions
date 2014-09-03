/* utils.h
 *
 *   Copyright (C) 2013 David Bolin, Finn Lindgren
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

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <map>
#include <algorithm>
#include "camd.h"

#ifdef _OPENMP
	#include<omp.h>
#endif

extern "C" {
#include "cholmod.h"
}



using namespace std;

typedef struct _amplitude_index{ 
	double val;
	int ind;
} vec_ind;

/*
 Internal comparisson function used for sorting vectors.
	@param a first structure to compare
	@param b second structure to compare
	@return 1 if a>b, 0 if a=b, and -1 if a<b
 */
int compare_structs (const void *a, const void *b);

/*
 Function used for reading sparce matrices into a vector< map<int,double> > format.
	@param matrix pointer to vector< map<int,double> > object to store the matrix in
	@param m_size vector where the size of the matrix and the number of non-zero elements are stored.
	@param ifile path to binary file containing indices to non-zero elements
	@param vfile path to binary file containing the corresponding non-zero values
 */
void readSparseMatrix(vector< map<int,double> > * matrix, int * m_size, char* ifile, char* vfile);

/*
 Function used to convert a matrix stored in vector< map<int,double> > format to CCS format.
	@param matrix pointer to the matrix to convert
	@param A double pointer to a taucs_ccs_matrix where the matrix is to be stored
	@param n size of matrix (the matrix is n x n)
	@param nz the number of non-zero elements in the matrix
 */
void convert_to_ccs(vector< map<int,double> > * matrix, cholmod_sparse ** A,int n,int nz, cholmod_common * cm);

/*
 Function used to convert a matrix stored in taucs_ccs_matrix format to vector< map<int,double> > format.
 Note that this function also transposes the matrix before converting.
	@param matrix pointer to vector< map<int,double> > object to store the matrix in
	@param A double pointer to a taucs_ccs_matrix to convert
 */
void convert_from_ccs(vector< map<int,double> > * matrix,cholmod_sparse ** A, cholmod_common * cm);

/*
 Function used for sorting a vector and extracting the index vector of the resulting permutation
	@param n the number of elements in the vector
	@param v the vector to be sorted
	@param vs vector to save the sorted reult in
	@param ind vector to save the permutation in
 */
void sort_vector(int n, double v[], double vs[], int ind[]);

/*
 Function used for calculating the diagonal elements of the inverse of a sparse matrix stored in CCS format
	@param Rir vector with row indices
	@param Rjc vector with column pointers
	@param Rpr vector with corresponding values
	@param variances vector to store the diagonal elements in
	@param n the size of the matrix (n x n)
*/
void Qinv(int * Rir, int * Rjc, double * Rpr, double * variances, int n);

void output_matrix(int m, int n, vector< map<int,double> > * matrix);

void vec_permute(int n, double * a,double * a_sort,int * reo);

