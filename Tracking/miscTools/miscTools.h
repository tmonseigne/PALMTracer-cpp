



/*******************************************************************************************************************
								Miscellaneous tools 
*******************************************************************************************************************/


#ifndef __MISC_TOOLD__
	#define __MISC_TOOLD__

	#include <stdio.h>
	#include <errno.h>
	#include <math.h>
	#include "../array_library_static/allDataLib.h"

	void printErrno(FILE *file, int errno);
	double tenTo(double a);
	unsigned long twoToThe(unsigned long val);
	wchar_t* stringReplace(wchar_t s, wchar_t d, wchar_t* st);
	int computeMSD(DblArray &mSD,DblArray &positions);
	double getAutoThreshold(FltArray &im);
	void fitPow(double *x, double *y, int ndata, double &pow, double &k);
	double distance(double X1, double Y1, double X2, double Y2);
	double distance(double X1, double Y1, double Z1, double X2, double Y2, double Z2);
	double distance2(double X1, double Y1, double X2, double Y2);
	double distance2(double X1, double Y1, double Z1, double X2, double Y2, double Z2);
	int printRegionsW(wchar_t* filename,  USrtArray &mainSegmentation, long nbStructures);
#endif //__MISC_TOOLD__ 