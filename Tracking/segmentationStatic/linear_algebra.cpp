
/*******************************************************************************************************************
								Linear algebra methods
								(Jacobian Method to compute the eigenvalues 
								and the eignenvectors of a real symmetric matrix)
									Institut Curie
									UMR - 144
									by Victor Racine
									04/15/02 
*******************************************************************************************************************/





#include "linear_algebra.h"

void eigsrt(double *d, double **v, int n)
{
	int k,j,i;
	double p;

	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}

void free_matrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


void free_vector(double *v,int nl,int nh)
{
	free((char*) (v+nl));
}


double **matrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		m[i] -= ncl;
	}
	return m;
}

double *vector(int nl,int nh)
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	return v-nl;
}

void eulerRotation(double *x, double *y, double *z, double phi, double teta, double psi){
	double rotX, rotY, rotZ;
	//ref http://mathworld.wolfram.com/EulerAngles.html

	rotX=(cos(psi)*cos(phi)-cos(teta)*sin(phi)*sin(psi))**x;	//a11
	rotX+=(cos(psi)*sin(phi)+cos(teta)*cos(phi)*sin(psi))**y;	//a12
	rotX+=sin(psi)*sin(teta)**z;								//a13
	
	rotY=-(sin(psi)*cos(phi)-cos(teta)*sin(phi)*cos(psi))**x;	//a21
	rotY+=-(sin(psi)*sin(phi)+cos(teta)*cos(phi)*cos(psi))**y;	//a22
	rotY+=cos(psi)*sin(teta)**z;								//a23
	
	rotZ=sin(teta)*sin(phi)**x;									//a31
	rotZ+=-sin(teta)*cos(phi)**y;								//a32
	rotZ+=cos(teta)**z;											//a33

	*x=rotX;
	*y=rotY;
	*z=rotZ;

}






#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);



void jacobi(double **a, int n, double *d, double **v, int *nrot){

		/* void jacobi(float **a, int n, float *d, float **v, int *nrot)
		** 
		** return the eigenvalues and the eignenvectors of a real symmetric matrix
		** 
		** a: input real symmetric matrix
		** n: input matrix size (a is square)
		** d: output vector of eigenvalues
		** v: output matrix whose columns contains the normalized  eigenvectors
		** nrot: number of Jacobi rotations 
		**
		** From numerical recipes in C p.467. The procedure is modified to keep 
		** the input matrix a unchange. The procedure operate on a copy a1 of matrix 
		*/

	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b=vector(1,n);
	z=vector(1,n);
	for (ip=1;ip<=n;ip++) {
		for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
		v[ip][ip]=1.0;
	}
	for (ip=1;ip<=n;ip++) {
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0.0;
	}
	*nrot=0;
	for (i=1;i<=50;i++) {
		sm=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++)
				sm += (double)fabs(a[ip][iq]);
		}
		if (sm == 0.0) {
			free_vector(z,1,n);
			free_vector(b,1,n);
			return;
		}
		if (i < 4)
			tresh=(double)(0.2*sm/(n*n));
		else
			tresh=0.0;
		for (ip=1;ip<=n-1;ip++) {
			for (iq=ip+1;iq<=n;iq++) {
				g=(double)(100.0*fabs(a[ip][iq]));
				if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
					&& (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
					a[ip][iq]=0.0;
				else if (fabs(a[ip][iq]) > tresh) {
					h=d[iq]-d[ip];
					if ((double)(fabs(h)+g) == (double)fabs(h))
						t=(a[ip][iq])/h;
					else {
						theta=(double)(0.5*h/(a[ip][iq]));
						t=(double)(1.0/(fabs(theta)+sqrt(1.0+theta*theta)));
						if (theta < 0.0) t = -t;
					}
					c=(double)(1.0/sqrt(1+t*t));
					s=t*c;
					tau=(double)(s/(1.0+c));
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0.0;
					for (j=1;j<=ip-1;j++) {
						ROTATE(a,j,ip,j,iq)
					}
					for (j=ip+1;j<=iq-1;j++) {
						ROTATE(a,ip,j,j,iq)
					}
					for (j=iq+1;j<=n;j++) {
						ROTATE(a,ip,j,iq,j)
					}
					for (j=1;j<=n;j++) {
						ROTATE(v,j,ip,j,iq)
					}
					++(*nrot);
				}
			}
		}
		for (ip=1;ip<=n;ip++) {
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0.0;
		}
	}
}
#undef ROTATE
