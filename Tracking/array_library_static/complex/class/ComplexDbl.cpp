

/*******************************************************************************************************************
								Pactage for general complex real 64bit
									Institut Curie
									UMR - 144
									by Victor Racine
									sept, 23rd 2003
*******************************************************************************************************************/


#include "allComplexNumber.h"



#define	_COMPLEX_NUMBER_
#define ComplexNumber	ComplexDbl
#define DataTypeComplex	double
#include "ComplexNumber.cpp"
#undef _COMPLEX_NUMBER_
#undef ComplexNumber
#undef DataTypeComplex

 