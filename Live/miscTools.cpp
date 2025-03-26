


#include "stdafx.h"
#include "miscTools.h"


double tenTo(double a){
	return(exp(a*log(10.)));
}

wchar_t* stringReplace(wchar_t s, wchar_t d, wchar_t* st){
	wchar_t*ptrSt=st;
	while(*ptrSt!='\0'){
		if(*ptrSt==s)
			*ptrSt=d;
		ptrSt++;
	}
	return(st);
}


unsigned long twoToThe(unsigned long val){
	if(val==0) return(1);
	return(((unsigned long)1) << val);
}


double distance(double X1, double Y1, double X2, double Y2){
	return(sqrt(distance2(X1,Y1,X2,Y2)));
}
double distance(double X1, double Y1, double Z1, double X2, double Y2, double Z2){
	return(sqrt(distance2(X1,Y1,Z1,X2,Y2,Z2)));
}

double distance2(double X1, double Y1, double X2, double Y2){
	double res;
	res=(X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2);
	if(res<0)
		return(0.);
	return(res);
}
double distance2(double X1, double Y1, double Z1, double X2, double Y2, double Z2){
	double res;
	res=(X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2);
	if(res<0)
		return(0.);
	return(res);
} 


