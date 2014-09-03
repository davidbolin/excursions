/* gaussint.cpp
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
#include "integration.h"

using namespace std;

int main (int argc, char * const argv[]) {
	
	string path;
	
	if (argc > 1) {
		path = argv[1];
    } else {
		path =  "/tmp/";
    }
	
	FILE * pFile;
	int i,n_itr,chol,ind_p,nz,max_size,verbose,n_threads;
	int m_size[2];

	double alpha,u, sigma;
	vector< map<int,double> > Q;
	vector< map<int,double> > R;
	cholmod_sparse * Qt;
		
	/* Read parameters */
	pFile = fopen((path + "initdata.bin").c_str(), "rb");
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	if(1!=fread(&alpha,sizeof(double),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1!=fread(&n_itr,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1!=fread(&chol,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1!=fread(&ind_p,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1!=fread(&max_size,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1!=fread(&verbose,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1!=fread(&n_threads,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	fclose (pFile);
	
	/* Read precision matrix */
	readSparseMatrix(&Q,m_size,(char*) (path + "precI.bin").c_str(), (char*) (path + "precV.bin").c_str());
	
	int n = m_size[0];
	nz = m_size[1];
	
	double * mu = new double[n]; 
	double * a = new double[n];
	double * b = new double[n];
	double * Pv = new double[n];
	double * Ev = new double[n];
	double * result = new double[2];
	int * reo = new int[n];
	int * ireo = new int[n];
	int * ind = new int[n];

	for(int i=0;i<n;i++){
		Pv[i] = 0.0;
		Ev[i] = 0.0;
	}
	
	cholmod_common cm ;
	cholmod_start (&cm);
		
	(&cm)->print = 5;
	(&cm)->final_ll = 1;
	
	convert_to_ccs(&Q,&Qt,n,nz,&cm);
	
	pFile = fopen((path+"mu.bin").c_str(), "rb");
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	if(n!=fread(mu,sizeof(double),n,pFile)){
		fputs ("Read error mu\n",stderr); exit (1);
	}
	fclose(pFile);

	pFile = fopen((path+"a.bin").c_str(), "rb");
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	if(n!=fread(a,sizeof(double),n,pFile)){
		fputs ("Read error a\n",stderr); exit (1);
	}
	fclose(pFile);
	
	pFile = fopen((path+"b.bin").c_str(), "rb");
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	if(n!=fread(b,sizeof(double),n,pFile)){
		fputs ("Read error b\n",stderr); exit (1);
	}
	fclose(pFile);
		
	for(i=0;i<n;i++){
		a[i] -= mu[i];
		b[i] -= mu[i];
	}

	if (chol == 0 && ind_p==1){ 
		// only a subregtion is considererd.
		if(verbose==1){
			cout << "Calculating reordering..." << endl;
		}	
		pFile = fopen((path+"ind.bin").c_str(), "rb");
		if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
		if(n!=fread(ind,sizeof(int),n,pFile)){
			fputs ("Read error ind\n",stderr); exit (1);
		}
		fclose(pFile);
		
		// use CAMD to find ordering
		double Control [CAMD_CONTROL], Info [CAMD_INFO];
		camd_l_defaults(Control);
		camd_order(n,(int*)(*Qt).p,(int*)(*Qt).i,reo,Control, Info, ind);
		
		double * a_sort = new double[n];
		double * b_sort = new double[n];
		vector< map<int,double> > R;
				
		vec_permute(n,a,a_sort,reo);
		vec_permute(n,b,b_sort,reo);
				
		(&cm)->nmethods = 1;
		(&cm)->method[0].ordering =1;
		(&cm)->postorder = 0;
					
		cholmod_factor* Rt = cholmod_analyze_p(Qt,reo,NULL,0,&cm);
		if(verbose==1){
			cout << "Calculating Cholesky factor..." << endl;
		}
		cholmod_factorize(Qt,Rt,&cm);	
		cholmod_sparse * Rtr = cholmod_factor_to_sparse(Rt,&cm);
		cholmod_free_factor(&Rt,&cm);
		convert_from_ccs(&R, &Rtr,&cm);
		cholmod_free_sparse(&Rtr,&cm);
				
		if(verbose==1){
			cout << "Calculating integral..." << endl;
		}
		shapeInt(&R, a_sort, b_sort, n_itr, n, 1-alpha, Pv, Ev,max_size,n_threads);
		result[0] = Pv[n-max_size]; //extract the probability  for the 
		result[1] = Ev[n-max_size]; //final node considered.
		delete[] a_sort;
		delete[] b_sort;

	} else if (chol == 0 && ind_p == 0 ) { 
		//entire region is considered, reorder for sparsity
		if(verbose==1){
			cout << "Calculating Cholesky factor..." << endl;
		}
		cholmod_factor * Rt = cholmod_analyze(Qt,&cm);
		cholmod_factorize(Qt,Rt,&cm);
		cholmod_sparse * Rtr = cholmod_factor_to_sparse(Rt,&cm);
			
		for (i=0; i<n; i++) {
			reo[i] = ((int*)(*Rt).Perm)[i];
		}
		for (i = 0; i < n; i++ ){
			ireo[ reo[i] ] = i;
		}
		cholmod_free_factor(&Rt,&cm);
		convert_from_ccs(&R, &Rtr,&cm);
		cholmod_free_sparse(&Rtr,&cm);	
			
		double * a_sort = new double[n];
		double * b_sort = new double[n];
				
		vec_permute(n,a,a_sort,reo);
		vec_permute(n,b,b_sort,reo);
		if(verbose==1){
			cout << "Calculating integral..." << endl;
		}
		shapeInt(&R, a_sort, b_sort, n_itr, n, 1-alpha, Pv, Ev,max_size,n_threads);
		
		result[0] = Pv[0];
		result[1] = Ev[0];
		delete[] a_sort;
		delete[] b_sort;
			
	} else { // cholesky factor is provided for entire region
		if(verbose==1){
			cout << "Calculating integral..." << endl;
		}
		shapeInt(&Q, a, b, n_itr, n, 1-alpha, Pv, Ev,max_size,n_threads);
		result[0] = Pv[0];
		result[1] = Ev[0];
	}

	cholmod_free_sparse(&Qt,&cm);
	cholmod_finish(&cm);
	
	if(verbose==1){
		cout << "Done! Writing data..." << endl;	
	}
	
	pFile = fopen((path+"results.bin").c_str(), "wb");
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	if(2!=fwrite(result, sizeof(double), 2, pFile)){
		fputs ("Write error b\n",stderr); exit (1);
	}
	fclose(pFile);

	delete[] result;
	delete[] mu;
	delete[] a;
	delete[] b;
	delete[] Pv;
	delete[] Ev;
	delete[] reo;
	delete[] ireo;
	delete[] ind;
	
    return 0;
}

