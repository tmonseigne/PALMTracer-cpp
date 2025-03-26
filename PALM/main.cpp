#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl\gsl_randist.h>
#include <gsl\gsl_rng.h>
#include <gsl\gsl_statistics.h>

using namespace std;
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#define  ARMA_DONT_USE_WRAPPER
#include <armadillo>
using namespace arma;

/* number of data points to fit */
#define N 201

/* number of fit coefficients */
#define NCOEFFS 12

/* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
#define NBREAK (NCOEFFS - 2)

class splitstring : public string
{
	vector<string> flds;

public:
	splitstring(char* s) : string(s) { };
	vector<string>& split(char delim, int rep = 0);
};

// split: receives a char delimiter; returns a vector of strings
// By default ignores repeated delimiters, unless argument rep == 1.
vector<string>& splitstring::split(char delim, int rep)
{
	if (!flds.empty()) flds.clear(); // empty vector if necessary
	string work = data();
	string buf = "";
	int i = 0;
	while (i < work.length()) {
		if (work[i] != delim) buf += work[i];
		else if (rep == 1) {
			flds.push_back(buf);
			buf = "";
		}
		else if (buf.length() > 0) {
			flds.push_back(buf);
			buf = "";
		}
		i++;
	}
	if (!buf.empty()) flds.push_back(buf);
	return flds;
}

double evalFunction(double x, int nbPoints, vec x_sp, double* a, double* b, double* c, double* d)
{
	double y = -1.0;
	for (int i = 0; i < nbPoints - 1; i++) { if ((x < x_sp(i + 1)) && (x >= x_sp(i))) { y = a[i] * pow(x - x_sp(i), 3) + b[i] * pow(x - x_sp(i), 2) + c[i] * (x - x_sp(i)) + d[i]; } }

	return y;
}

void computeBsplines(int n, gsl_vector* x, gsl_vector* y, int* nbPoints, vec* x_sp, vec* y_sp)
{
	const size_t ncoeffs = NCOEFFS;
	const size_t nbreak = NBREAK;
	size_t i, j;
	gsl_bspline_workspace* bw;
	gsl_vector* B;
	double dy;
	gsl_rng* r;
	gsl_vector *c, *w;

	gsl_matrix *X, *cov;
	gsl_multifit_linear_workspace* mw;
	double chisq, Rsq, dof, tss;

	gsl_rng_env_setup();
	r = gsl_rng_alloc(gsl_rng_default);

	/* allocate a cubic bspline workspace (k = 4) */
	bw = gsl_bspline_alloc(4, nbreak);
	B = gsl_vector_alloc(ncoeffs);

	x = gsl_vector_alloc(n);
	y = gsl_vector_alloc(n);
	X = gsl_matrix_alloc(n, ncoeffs);
	c = gsl_vector_alloc(ncoeffs);
	w = gsl_vector_alloc(n);
	cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
	mw = gsl_multifit_linear_alloc(n, ncoeffs);

	ifstream f("Calibration.txt");
	string line;
	getline(f, line); // Header

	/* load data */
	double xMin = +9999;
	double xMax = -9999;
	for (i = 0; i < n; ++i) {
		getline(f, line);
		splitstring s((char*)line.c_str());
		std::vector<std::string> tokens = s.split('\t', 2);
		double xi = atof(tokens[0].c_str());

		if (xi < xMin) xMin = xi;
		if (xi > xMax) xMax = xi;

		double yi = atof(tokens[1].c_str());

		gsl_vector_set(x, i, xi);
		gsl_vector_set(y, i, yi);
		gsl_vector_set(w, i, 1.0);
	}

	/* use uniform breakpoints on [xMin, xMax] */
	gsl_bspline_knots_uniform(xMin, xMax, bw);

	/* construct the fit matrix X */
	for (i = 0; i < n; ++i) {
		double xi = gsl_vector_get(x, i);

		/* compute B_j(xi) for all j */
		gsl_bspline_eval(xi, B, bw);

		/* fill in row i of X */
		for (j = 0; j < ncoeffs; ++j) {
			double Bj = gsl_vector_get(B, j);
			gsl_matrix_set(X, i, j, Bj);
		}
	}

	/* do the fit */
	gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);

	dof = n - ncoeffs;
	tss = gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
	Rsq = 1.0 - chisq / tss;

	*nbPoints = (int)floor((xMax - xMin) / 0.1) + 1;
	*x_sp = vec(*nbPoints);
	*y_sp = vec(*nbPoints);
	/* output the smoothed curve */
	{
		double xi, yi, yerr;

		FILE* f = fopen("Calibration_Smoothed.txt", "w");
		int ind = 0;
		for (xi = xMin; xi < xMax; xi += 0.1) {
			gsl_bspline_eval(xi, B, bw);
			gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
			(*x_sp)(ind) = xi;
			(*y_sp)(ind) = yi;
			fprintf(f, "%f\t%f\n", (*x_sp)(ind), (*y_sp)(ind));
			ind += 1;
		}
		fclose(f);
	}

	gsl_rng_free(r);
	gsl_bspline_free(bw);
	gsl_vector_free(B);
	gsl_matrix_free(X);
	gsl_vector_free(c);
	gsl_vector_free(w);
	gsl_matrix_free(cov);
	gsl_multifit_linear_free(mw);
}

int main(void)
{
	const size_t n = N;
	gsl_vector *x, *y;
	x = gsl_vector_alloc(n);
	y = gsl_vector_alloc(n);
	vec x_sp, y_sp;
	int nbPoints;

	computeBsplines(n, x, y, &nbPoints, &x_sp, &y_sp);
	fprintf(stderr, "nbPoints: %d\n", nbPoints);
	double h = 0.1;

	mat A = mat(nbPoints - 2, nbPoints - 2, fill::zeros);
	vec M = vec(nbPoints - 2, fill::zeros);
	vec B = vec(nbPoints - 2, fill::zeros);

	for (int i = 0; i < nbPoints - 2; i++) {
		B(i) = (y_sp(i) - 2 * y_sp(i + 1) + y_sp(i + 2)) * 6 / (h * h);
		A(i, i) = 4;

		if (i < nbPoints - 3) {
			A(i, i + 1) = 1;
			A(i + 1, i) = 1;
		}
	}
	M = solve(A, B);

	double* a = (double*)malloc((nbPoints - 1) * sizeof(double));
	double* b = (double*)malloc((nbPoints - 1) * sizeof(double));
	double* c = (double*)malloc((nbPoints - 1) * sizeof(double));
	double* d = (double*)malloc((nbPoints - 1) * sizeof(double));
	double* m = (double*)malloc(nbPoints * sizeof(double));

	m[0] = 0.0;
	m[nbPoints - 1] = 0.0;
	for (int i = 1; i < nbPoints - 1; i++) { m[i] = M(i - 1); }

	for (int i = 0; i < nbPoints - 1; i++) {
		a[i] = (m[i + 1] - m[i]) / (6 * h);
		b[i] = m[i] / 2;
		c[i] = (y_sp(i + 1) - y_sp(i)) / h - ((m[i + 1] + 2 * m[i]) / 6) * h;
		d[i] = y_sp(i);
	}

	FILE* f = fopen("Results.txt", "w");
	for (double xi = x_sp(0); xi <= x_sp(nbPoints - 1); xi += 0.01) fprintf(f, "%f\t%f\n", xi, evalFunction(xi, nbPoints, x_sp, a, b, c, d));
	fclose(f);

	return 0;
} /* main() */
