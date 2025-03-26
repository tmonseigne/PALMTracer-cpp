#include <STDIO.H>
#include "../print_debug.h"
extern FILE *debugStreamALS;

#include "allComplexNumber.h"

void test(void){
	ComplexChr a,b;

	a.re=1;
	a.im=-8;

	a*=2;
	b=(const char)2;

	openDebugStreamALS();
	fprintf(debugStreamALS, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\t%d %d %d %d \n", a.re, a.im, b.re, b.im); 
	closeDebugStreamALS();


}