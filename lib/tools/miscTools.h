/*******************************************************************************************************************
								Miscellaneous tools 
*******************************************************************************************************************/


#ifndef __MISC_TOOLD__
#define __MISC_TOOLD__

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include "allDataLib.h"

void printErrno(FILE* file, int errno);
double tenTo(double a);
unsigned long twoToThe(unsigned long val);
wchar_t* stringReplace(wchar_t s, wchar_t d, wchar_t* st);
double getAutoThreshold(FltArray& im);
double distance(double X1, double Y1, double X2, double Y2);
double distance(double X1, double Y1, double Z1, double X2, double Y2, double Z2);
double distance2(double X1, double Y1, double X2, double Y2);
double distance2(double X1, double Y1, double Z1, double X2, double Y2, double Z2);
#endif //__MISC_TOOLD__ 
