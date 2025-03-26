
/*******************************************************************************************************************
								Linear algebra methods
								(Jacobian Method to compute the eigenvalues 
								and the eignenvectors of a real symmetric matrix)
									Institut Curie
									UMR - 144
									by Victor Racine
									04/15/02 
*******************************************************************************************************************/


#ifndef	__linear_algebra__
		#define __linear_algebra__
		#include <math.h>
		#include <malloc.h>

		// sort eigenvalues
		void eigsrt(double *d, double **v, int n);
		
		double *vector(int nl,int nh);
		double **matrix(int nrl,int nrh,int ncl,int nch);
		void free_vector(double *v,int nl,int nh);
		void free_matrix(double **m,int nrl,int nrh,int ncl,int nch);
		void jacobi(double **a, int n, double *d, double **v, int *nrot);
		void eulerRotation(double *x, double *y, double *z, double phi, double teta, double psi);
		

#endif	//__linear_algebra__
