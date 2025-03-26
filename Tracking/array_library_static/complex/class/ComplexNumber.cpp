


/*******************************************************************************************************************
								Pactage for compex numbers
									Institut Curie
									UMR - 144
									by Victor Racine
									27/06/2002 
*******************************************************************************************************************/

#ifdef		_COMPLEX_NUMBER_

//#include "ComplexNumber.h"

#ifndef ComplexNumber
#error ComplexNumber
#endif



ComplexNumber::ComplexNumber(){
	;
}

ComplexNumber::ComplexNumber(DataTypeComplex _im, DataTypeComplex _re){
	im=_im;
	re=_re;
}

ComplexNumber::~ComplexNumber(){
}

ComplexNumber& ComplexNumber::operator=(char &dataTypeComplex2){
	im=(DataTypeComplex)0;
	re=(DataTypeComplex)dataTypeComplex2;
	return(*this);
}

ComplexNumber& ComplexNumber::operator=(const char &dataTypeComplex2){
	im=(DataTypeComplex)0;
	re=(DataTypeComplex)dataTypeComplex2;
	return(*this);
}

/*#define DataTypeComplex2  char
#include "complexOperators2.cpp"
#undef DataTypeComplex2*/

#define DataTypeComplex2  unsigned char
#include "complexOperators2.cpp"
#undef DataTypeComplex2

#define DataTypeComplex2  short
#include "complexOperators2.cpp"
#undef DataTypeComplex2

#define DataTypeComplex2  unsigned short
#include "complexOperators2.cpp"
#undef DataTypeComplex2

#define DataTypeComplex2  long
#include "complexOperators2.cpp"
#undef DataTypeComplex2

#define DataTypeComplex2  unsigned long
#include "complexOperators2.cpp"
#undef DataTypeComplex2

#define DataTypeComplex2  float
#include "complexOperators2.cpp"
#undef DataTypeComplex2

#define DataTypeComplex2  double
#include "complexOperators2.cpp"
#undef DataTypeComplex2

#define DataTypeComplex2  bool
#include "complexOperators2.cpp"
#undef DataTypeComplex2


#define ComplexNumber2  ComplexChr
#include "complexOperators.cpp"
#undef ComplexNumber2
#define ComplexNumber2  ComplexSrt
#include "complexOperators.cpp"
#undef ComplexNumber2
#define ComplexNumber2  ComplexLng
#include "complexOperators.cpp"
#undef ComplexNumber2
#define ComplexNumber2  ComplexFlt
#include "complexOperators.cpp"
#undef ComplexNumber2
#define ComplexNumber2  ComplexDbl
#include "complexOperators.cpp"
#undef ComplexNumber2


ComplexNumber& ComplexNumber::operator+=(const ComplexNumber &complexNumber){
	im+=complexNumber.im;
	re+=complexNumber.re;
	return(*this);
}

ComplexNumber& ComplexNumber::operator-=(const ComplexNumber &complexNumber){
	im-=complexNumber.im;
	re-=complexNumber.re;
	return(*this);
}

ComplexNumber& ComplexNumber::operator*=(const ComplexNumber &complexNumber){
	DataTypeComplex lastIm=im;
	im=complexNumber.im*re+complexNumber.re*im;
	re=-complexNumber.im*lastIm+complexNumber.re*re;
	return(*this);
}

ComplexNumber& ComplexNumber::operator/=(const ComplexNumber &complexNumber){
	DataTypeComplex lastIm=im;
	DataTypeComplex denom=complexNumber.im*complexNumber.im+complexNumber.re*complexNumber.re;
	im=(im*complexNumber.re-complexNumber.im*re)/denom;
	re=(lastIm*complexNumber.im+re*complexNumber.re)/denom;
	return(*this);
}


ComplexNumber& ComplexNumber::operator+=(const DataTypeComplex &dataTypeComplex){
	re+=dataTypeComplex;
	return(*this);
}

ComplexNumber& ComplexNumber::operator-=(const DataTypeComplex &dataTypeComplex){
	re-=dataTypeComplex;
	return(*this);
}

ComplexNumber& ComplexNumber::operator*=(const DataTypeComplex &dataTypeComplex){
	re*=dataTypeComplex;
	im*=dataTypeComplex;
	return(*this);
}



ComplexNumber& ComplexNumber::operator/=(const DataTypeComplex &dataTypeComplex){
#ifdef DEBUG
	if(dataTypeComplex==(DataTypeComplex)0.)
		exit(DIV_ZERO);
#endif
	re/=dataTypeComplex;
	im/=dataTypeComplex;
	return(*this);
}

/*ComplexNumber& ComplexNumber::operator=(const DataTypeComplex &dataTypeComplex){
	re=dataTypeComplex;
	im=(DataTypeComplex)0.;
	return(*this);
}*/





DataTypeComplex	ComplexNumber::getAngle(void){
	DataTypeComplex	angle;

	angle=(DataTypeComplex)atan2(im,re);
	if(re>0)
		angle+=(DataTypeComplex)PI;
	return(angle);
}

DataTypeComplex	ComplexNumber::getAbs(void){
	return((DataTypeComplex)sqrt(im*im+re*re));
}


bool operator>(ComplexNumber &complexNumber1, ComplexNumber &complexNumber2){
	return(complexNumber1.getAbs()>complexNumber2.getAbs());
}
bool operator<(ComplexNumber &complexNumber1, ComplexNumber &complexNumber2){
	return(complexNumber1.getAbs()<complexNumber2.getAbs());
}


ostream& operator<< ( ostream& os, ComplexNumber &complexNumber )
{
   os << complexNumber.re; 
   if (complexNumber.im>0) os << "+";
   os << complexNumber.im << "i";
   return os;
}

#endif //		_ARRAY_TYPE_

