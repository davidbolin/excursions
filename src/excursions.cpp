/* excursions.cpp
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

//#define TIMING
#include "utils.h"
#include "integration.h"

#ifdef TIMING
	#include<ctime>
#endif

using namespace std;

int main (int argc, char * const argv[]) {
	/*
	 alpha: 	probability
	 u: 		level
	 mu: 	posterior mean E(X(s))
	 vars: 	posterior variances V(X(s))
	 rho = 	marginal probabilities P(X(s)>u)
	 n_itr = number of MC iterations
	 type:	0 = positive excursion '>', 1= negative excursion '<', 
	 		2= contour uncertainty, 3= contour map
	 reo:	reordering
	 Q:		precision matrix or its cholesky factor
	 chol:	1=Q is cholesky factor, 0= Q is precision
	 n:		number of nodes in field
	 rho_p	1 if rho provided, 2 if QC method is used.
	 reo_p	1 if reo provided
	 ind_p	1 if a subset of indices are provided
	 */
	
	/*	 
	 TODO: implement latent Gaussian (numerical integration)
	 This can be handeled in R and Matlab to begin with; however, we can take
	 advantage of the sparsity structure of Q to not have to calculate the
	 cholesky factor from scratch for each parameter estimate. Thus, we should 
	 problably implement this in C at some point.
	 
	 TODO: implement two-parameter family (this can probably wait)
	 
	 TODO: implement program for covariance matrix (this can probably wait)
	 */
	
	
#ifdef TIMING
	clock_t start = clock();
	const double ticks_per_ms = static_cast<double>(CLOCKS_PER_SEC)/1000;
#endif
	
	string path;
	
	if (argc > 1) //input path to temp directory
    {
		path = argv[1];
    }
	else //if no input is used, take normal system default.
    {
		path =  "/tmp/";
    }
	
	FILE * pFile;
	int i, n_itr, type, chol, rho_p, reo_p, vars_p;
	int ind_p, nz, max_size, verbose, n_threads;
	int m_size[2];

	double alpha,u, sigma;
	vector< map<int,double> > Q;
	vector< map<int,double> > R;
	cholmod_sparse * Qt;
	
#ifdef TIMING
	double time_read_data = 0;
	double time_comp_rho  = 0;	
	double time_comp_reo  = 0;
	double time_comp_perm = 0;
	double time_comp_chol = 0;
	double time_comp_int  = 0;
	double time_write_data= 0;
	double init_time = static_cast<double>(clock()-start)  / ticks_per_ms;
	start = clock();
#endif
	
	/* Read parameters */
	pFile = fopen((path + "initdata.bin").c_str(), "rb");
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	
	if(1 != fread(&alpha,sizeof(double),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1 != fread(&u,sizeof(double),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1 != fread(&n_itr,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1 !=  fread(&type,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1 !=  fread(&chol,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1 !=  fread(&rho_p,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1 !=  fread(&reo_p,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1 !=  fread(&vars_p,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1 !=  fread(&ind_p,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1 !=  fread(&max_size,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1 !=  fread(&verbose,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}
	if(1 != fread(&n_threads,sizeof(int),1,pFile)){
		fputs ("Read error initdata\n",stderr); exit (1);
	}	
	fclose (pFile);
	
	/* Read precision matrix */
	readSparseMatrix(&Q,m_size,(char*) (path + "precI.bin").c_str(), (char*) (path + "precV.bin").c_str());
	
	int n = m_size[0];
	nz = m_size[1];
	
	double * vars = new double[n];
	double * mu = new double[n]; 
	double * rho = new double[n];
	double * rho_sort = new double[n]; 
	double * rho_u = new double[n];
	double * rho_ngu = new double[n];
	double * rho_l = new double[n];
	double * a = new double[n];
	double * b = new double[n];
	double * Pv = new double[n];
	double * Ev = new double[n];
	double * rho_ng = new double[n];
	double * uv = new double[n];
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
	/* Read mean */
	pFile = fopen((path+"mu.bin").c_str(), "rb");
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	if( n!= fread(mu,sizeof(double),n,pFile)){
		fputs ("Read error mu\n",stderr); exit (1);
	}
	fclose(pFile);

	if(verbose==1){
		if (vars_p==1) {
			cout << "Reading variances..." << endl;	
		} else {
			cout << "Calculating variances..." << endl;
		}
	}

	/* Read or calculate variances */
	if (vars_p==1) { // read variances
		pFile = fopen((path+"vars.bin").c_str(), "rb");
		if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
		if(n!=fread(vars,sizeof(double),n,pFile)){
			fputs ("Read error vars\n",stderr); exit (1);
		}
		fclose(pFile);
	} else { // calculate variances
		if (chol==1) {
			Qinv((int*)(*Qt).i, (int*)(*Qt).p, (double*)(*Qt).x,vars,n);
		} else { //calculate cholesky factor and then variances
			cholmod_factor * Rt = cholmod_analyze(Qt,&cm);
			cholmod_factorize(Qt,Rt,&cm);
			cholmod_sparse * Rtr = cholmod_factor_to_sparse(Rt,&cm);
			double * vars_reo = new double[n];
			Qinv((int*)(*Rtr).i, (int*)(*Rtr).p, (double*)(*Rtr).x,vars_reo,n);
			for (i=0; i<n; i++) {
				vars[((int*)(*Rt).Perm)[i]] = vars_reo[i];
			}
			cholmod_free_factor(&Rt,&cm);
			delete[] vars_reo;
		}
	}
	
	#ifdef TIMING
	time_read_data += static_cast<double>(clock()-start) / ticks_per_ms;
	start = clock();
	#endif
	
	if(verbose==1){
		if (rho_p==1) {
			cout << "Reading marginal excursion probabilities..." << endl;	
		} else if (rho_p==0){
			cout << "Calculating marginal excursion probabilities..." << endl;
		} else {
			cout << "Reading marginal excursion probabilities and calculating Gaussian quantiles..." << endl;
		}
	}

	/* Read or calculate marginal excursion probabilities */
	if (type ==2) { // contour uncertainty
		if (rho_p==1) {
			pFile = fopen((path+"rho.bin").c_str(), "rb");
			if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
			if (n!= fread(rho_u,sizeof(double),n,pFile)){
				fputs ("Read error rho\n",stderr); exit (1);
			}
			fclose(pFile);	
			for (i=0; i<n; i++) {
				rho_l[i] = 1-rho_u[i];
				rho[i] = max(rho_u[i],rho_l[i]);
			}
		} else if (rho_p == 0) {
			for (i=0;i<n;++i) {
				sigma = sqrt(vars[i]);
				rho_l[i] = gsl_cdf_gaussian_Q(mu[i]-u, sigma);
				rho_u[i] = 1-rho_l[i];
				rho[i] = max(rho_u[i],rho_l[i]);
			}
		} else if (rho_p == 2) { // QC method
			//Read non-gaussian probabilities
			pFile = fopen((path+"rho.bin").c_str(), "rb");
			if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
			if(n != fread(rho_ngu,sizeof(double),n, pFile)){
				fputs ("Read error rho_ngu\n",stderr); exit (1);
			}
			fclose(pFile);
				
			//calculate gaussian
			for (i=0;i<n;++i) {
				sigma = sqrt(vars[i]);
				rho_l[i] = gsl_cdf_gaussian_Q(mu[i]-u, sigma);
				rho_u[i] = 1-rho_l[i];
				rho[i] = max(rho_u[i],rho_l[i]);
				rho_ng[i] = max(rho_ngu[i],1-rho_ngu[i]);
				
			}
		}
	} else if (type == 3){ // contour function
		pFile = fopen((path+"rho.bin").c_str(), "rb");
		if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
		if(n != fread(rho,sizeof(double),n,pFile)){
			fputs ("Read error rho\n",stderr); exit (1);
		}
		fclose(pFile);	
	} else {
		if (rho_p==1) {
			pFile = fopen((path+"rho.bin").c_str(), "rb");
			if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
			if (n != fread(rho,sizeof(double),n,pFile)){
				fputs ("Read error rho\n",stderr); exit (1);
			}
			fclose(pFile);		
		} else if (rho_p == 0){
			if (type == 0) { // positive excursion
				for (i=0;i<n;i++) {
					sigma = sqrt(vars[i]);
					rho[i] = gsl_cdf_gaussian_P(mu[i]-u, sigma);
				}
			} else { // negative excursion
				for (i=0;i<n;++i) {
					sigma = sqrt(vars[i]);
					rho[i] = gsl_cdf_gaussian_Q(mu[i]-u, sigma);
				}
			} 
		} else if (rho_p == 2) {
			//read non-Gaussian probabilities
			pFile = fopen((path+"rho.bin").c_str(), "rb");
			if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
			if (n != fread(rho_ng,sizeof(double),n,pFile)){
				fputs ("Read error rho_ng\n",stderr); exit (1);
			}
			fclose(pFile);	
			//calculate Gaussian
			for (i=0;i<n;i++) {
				sigma = sqrt(vars[i]);
				rho[i] = (type == 0) ? gsl_cdf_gaussian_P(mu[i]-u, sigma) : gsl_cdf_gaussian_Q(mu[i]-u, sigma);
			}
		}
	}
	
	if (ind_p==1){
		/*
		 * If a reordering is not provided and only a subset of indices are 
		 * to be considered, set marginal probabilities to zero for all indices
		 * that are not considered.
		*/
		pFile = fopen((path+"ind.bin").c_str(), "rb");
		if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
		if (n != fread(ind,sizeof(int),n,pFile)){
			fputs ("Read error ind\n",stderr); exit (1);
		}
		fclose(pFile);
		for(i=0; i<n; i++){
			if(ind[i]==0){
				rho[i]= -1.0;
				rho_ng[i] = -1.0;
			} 
		}
	} 
	
	
#ifdef TIMING
	time_comp_rho += static_cast<double>(clock()-start) / ticks_per_ms;
	start = clock();
#endif
	
	if(verbose==1){
		if (reo_p==1) {
			cout << "Reading permutation..." << endl;	
		} else {
			cout << "Calculating permutation..." << endl;
		}
	}

	/* Read or calculate permutation */	
	if (reo_p==1) {
		pFile = fopen((path+"reo.bin").c_str(), "rb");
		if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
		if (n!= fread(reo,sizeof(int),n,pFile)){
			fputs ("Read error reo\n",stderr); exit (1);
		}
		fclose(pFile);
	} else {
		/*
		 * Simple sorting without sparsity optimization:
		 */
		if (rho_p ==2){ // if QC is used, sort non-gaussian quantiles
			//sort_vector(n, rho_ng, rho_sort,reo);
			sort_vector(n, rho_ng, rho_sort,reo);
		} else {
			sort_vector(n, rho, rho_sort,reo);
		}
		for (i = 0; i < n; i++ ){
			ireo[ reo[i] ] = i;
		}
				
		/*
		 * Reordering using bounds and CAMD if alpha < 0.5:
		 * Set alpha > 0.5 if interested in F(s).
		 */
		
		//if ( (type==2 && alpha < 0.5) || (type != 2 && alpha < 1) ) {
		if ( ind_p || alpha < 1 ) {
		
			int * cindr = new int[n];
			int * cind = new int[n];
			
			
			int k=0;
			i = n-1;
			/*
			 * If alpha > 0.5, use upper and lower bounds, else only 
			 * upper bound: Using a lower bound does not work for F(s).
			 */
			/*
			if (0){//(alpha < 0.5) {
				while (rho_sort[i] > 1-alpha/(i+1) && i>=0){ 
					cindr[i] = k;
					i--;
				}
				if (i<n-1) {
					k=1;
				}
			} 
			*/
			// Add nodes in order down to the lower bound
			while (rho_sort[i] > 1 - alpha &&  i>=0) {
				cindr[i] = k; 
				i--;
				k++;
			}

			if(i>0){
				// Reorder nodes below the lower bound for sparsity
				while (i>=0) {
					cindr[i] = k;
					i--;
				}
				
				// change back to original ordering before calling camd:			
				for (i=0; i<n; i++) {
					cind[i] = k - cindr[ireo[i]];
				}
				
				// use CAMD to find ordering
				//double Control [CAMD_CONTROL], Info [CAMD_INFO];
				
				//camd_l_defaults(Control);
				//camd_defaults(Control);
				int * reo2 = new int[n];
				if(verbose==1){
					cout << "call CAMD...";
				}
				//int status = camd_order(n,(int*)(*Qt).p,(int*)(*Qt).i,reo2,Control, Info, cind);
				int status = camd_order(n,(int*)(*Qt).p,(int*)(*Qt).i,reo2,(double *) NULL, (double *) NULL, cind);
				if(verbose==1){
					cout << " done." << endl;
				}
				memcpy(reo, reo2, n*sizeof(int));
				
			}
			delete[] cind;
			delete[] cindr;	
		}
	}	
	
	if (ind_p==1){
		for(i=0; i<n; i++){
			if(ind[i]==0){
				rho[i]= 0.0; //change from -1 to 0 to avoid nan
				rho_ng[i] = 0.0;
			} 
		}
	} 
	
#ifdef TIMING
	time_comp_reo += static_cast<double>(clock()-start) / ticks_per_ms;
	start = clock();
#endif
	
	/* 
	 * Set limit vectors for integral 
	 */	
	if(verbose==1){
		cout << "Set limit vectors..." << endl;	
	}
	if(type != 3){	
		for (i=0; i<n; i++) {	// set uv depending on if QC method is used.
			if(rho_p == 2){
				if (type==1) { //negative
					uv[i] = sqrt(vars[i])*gsl_cdf_ugaussian_Pinv(rho_ng[i]);
				} else if (type==0){ //positive or contour
					uv[i] = sqrt(vars[i])*gsl_cdf_ugaussian_Qinv(rho_ng[i]);
				} else if (type==2){ //contour
					uv[i] = sqrt(vars[i])*gsl_cdf_ugaussian_Qinv(rho_ngu[i]);
				}
			} else {
			uv[i] = u-mu[i];
			}
		}
	}
	if (type == 2) { //contour region
		if(rho_p ==2){
			for (i=0;i<n;i++) {
			a[i] =(rho_ngu[i]>0.5)?uv[i]:-numeric_limits<double>::infinity();
			b[i] = (rho_ngu[i]>0.5)?numeric_limits<double>::infinity() : uv[i];
			}
		} else {
			for (i=0;i<n;i++) {
			a[i] = (rho_u[i]>0.5) ? uv[i] : -numeric_limits<double>::infinity();
			b[i] = (rho_u[i]>0.5) ? numeric_limits<double>::infinity() : uv[i];
			}
		}

	} else if (type == 0) { //positive excursions
		for (i=0;i<n;i++) {
			a[i] = uv[i];
			b[i] = numeric_limits<double>::infinity();
		}
	
	} else if (type == 1){ //negative excursions
		for (i=0;i<n;i++) {
			a[i] = -numeric_limits<double>::infinity();
			b[i] = uv[i];
		}
	} else if (type == 3){
		//read limit vectors from file
		pFile = fopen((path+"a.bin").c_str(), "rb");
		if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
		if (n != fread(a,sizeof(double),n,pFile)){
			fputs ("Read error a\n",stderr); exit (1);
		}
		fclose(pFile);
		pFile = fopen((path+"b.bin").c_str(), "rb");
		if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
		if (n != fread(b,sizeof(double),n,pFile)){
			fputs ("Read error b\n",stderr); exit (1);
		}
		fclose(pFile);
	}
	
	/*
	 * Calculate excursion function and first Cholesky factor if needed.
	 */	
	if(chol == 0){//(reo_p==0 && chol ==0){ 
	//reorder, calculate Choleksy factor, and then calculate excursion function
		
		if(verbose==1){
			cout << "Reordering..." << endl;
		}
		double * a_sort = new double[n];
		double * b_sort = new double[n];
		vector< map<int,double> > R;
				
		#ifdef TIMING
		start = clock();
		#endif

		vec_permute(n,a,a_sort,reo);
		vec_permute(n,b,b_sort,reo);
		
		#ifdef TIMING
		time_comp_perm += static_cast<double>(clock()-start) / ticks_per_ms;
		start = clock();
		#endif
		
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
				
		#ifdef TIMING
		time_comp_chol += static_cast<double>(clock()-start) / ticks_per_ms;
		start = clock();
		#endif
		if(verbose==1){
			cout << "Calculating excursion function..." << endl;
		}
		shapeInt(&R, a_sort, b_sort, n_itr, n, 1-alpha, Pv, Ev,max_size,n_threads);
		delete[] a_sort;
		delete[] b_sort;
		#ifdef TIMING
		time_comp_int += static_cast<double>(clock()-start) / ticks_per_ms;
		#endif
		
	} else if (0){//(reo_p == 1 && chol ==0) { 
	// Calculate Cholesky factor and then calculate excursion function
		if(verbose==1){
			cout << "Calculating Cholesky factor..." << endl;
		}
		vector< map<int,double> > R;
		(&cm)->nmethods = 1;
		(&cm)->method[0].ordering =1;
		(&cm)->postorder = 0;
		int * reo2 = new int[n];
		for(i=0;i<n;i++){
			reo2[i] = i;
		}
		cholmod_factor* Rt = cholmod_analyze_p(Qt,reo2,NULL,0,&cm);
		cholmod_factorize(Qt,Rt,&cm);
		cholmod_sparse * Rtr = cholmod_factor_to_sparse(Rt,&cm);
		cholmod_free_factor(&Rt,&cm);
		convert_from_ccs(&R, &Rtr,&cm);
		cholmod_free_sparse(&Rtr,&cm);
		if(verbose==1){
			cout << "Calculating excursion function..." << endl;
		}
		shapeInt(&R, a, b, n_itr, n, 1-alpha, Pv, Ev,max_size,n_threads);
		
	} else { 
	// Calculate excursion function
		if(verbose==1){
			cout << "Calculating excursion function..." << endl;
		}
		shapeInt(&Q, a, b, n_itr, n, 1-alpha, Pv, Ev,max_size,n_threads);
	}
	cholmod_free_sparse(&Qt,&cm);
	cholmod_finish(&cm);
	
if(verbose==1){
	cout << "Done! Writing data..." << endl;	
}
	
#ifdef TIMING
	time_comp_chol += static_cast<double>(clock()-start) / ticks_per_ms;
	start = clock();
#endif
	
	pFile = fopen((path+"results_p.bin").c_str(), "wb");
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	if(n!=fwrite(Pv, sizeof(double), n, pFile)){
	fputs ("Write error b\n",stderr); exit (1);
	}
	fclose(pFile);
	pFile = fopen((path+"results_e.bin").c_str(), "wb");
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	if(n!=fwrite(Ev, sizeof(double), n, pFile)){
	fputs ("Write error b\n",stderr); exit (1);
	}
	fclose(pFile);
	pFile = fopen((path+"results_rho.bin").c_str(), "wb");
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	if(n!=fwrite(rho, sizeof(double), n, pFile)){
	fputs ("Write error b\n",stderr); exit (1);
	}
	fclose(pFile);
	pFile = fopen((path+"results_reo.bin").c_str(), "wb");
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	if(n!=fwrite(reo, sizeof(int), n, pFile)){
	fputs ("Write error b\n",stderr); exit (1);
	}
	fclose(pFile);
	pFile = fopen((path+"vars.bin").c_str(), "wb");
	if (pFile==NULL) {fputs ("File error",stderr); exit (1);}
	if(n!=fwrite(vars, sizeof(double), n, pFile)){
	fputs ("Write error b\n",stderr); exit (1);
	}
	fclose(pFile);
	
#ifdef TIMING
	time_write_data += static_cast<double>(clock()-start) / ticks_per_ms;
	cout << "TIMING:" << endl;
	//cout << "Read data : " << time_read_data << endl;
	//cout << "Calc rho  : " << time_comp_rho << endl;
	cout << "Calc reo  : " << time_comp_reo << endl;
	//cout << "Calc perm : " << time_comp_perm << endl;
	cout << "Calc chol : " << time_comp_chol << endl;
	cout << "Calc int  : " << time_comp_int << endl;
	//cout << "Write data: " << time_write_data << endl;
	cout << "Total time: " << init_time + time_read_data + time_comp_rho + time_comp_reo + time_comp_perm + time_comp_chol + time_comp_int + time_write_data << endl;
#endif
	delete[] vars;
	delete[] mu;
	delete[] rho;
	delete[] rho_sort;
	delete[] rho_u;
	delete[] rho_ngu;
	delete[] rho_l; 
	delete[] a;
	delete[] b;
	delete[] Pv;
	delete[] Ev;
	delete[] rho_ng;
	delete[] uv;
	delete[] reo;
	delete[] ireo;
	delete[] ind;
	
    return 0;
}

