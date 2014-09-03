/* integration.h
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

extern "C"{
	#include "RngStream.h"
}

#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

using namespace std;

/*
 Function used for calculating Gaussian integrals P = |Q|/(2*pi)^n/2 int_a^b exp(-x'*Q*x/(2pi))dx.
 @param L pointer to vector< map<int,double> > object containing cholesky factor of Q 
 @param a vector containing lower limits of integral
 @param b vector containing upper limits of integral
 @param K number of particles in particle filter
 @param n dimension of integral
 @param lim limit value for integral, the integration is stopped if this value is reached
 @param Pv vector containing the resulting probabilities.
 @param Ev error estimate for the probabilities in Pv.
 @param max_size max size of the excursion set.
 @param n_threads the max number of threads the program can use.
 */
void shapeInt(vector< map<int,double> > * L, double * a,double * b,int K, const int n,double lim,double * Pv, double * Ev, int max_size, int n_threads);
