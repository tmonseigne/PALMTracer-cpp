

/*******************************************************************************************************************
								Pactage for general complex real 32bit
									Institut Curie
									UMR - 144
									by Victor Racine
									sept, 23rd 2003
*******************************************************************************************************************/


#include "allComplexNumber.h"



#define	_COMPLEX_NUMBER_
#define ComplexNumber	ComplexFlt
#define DataTypeComplex	float
#include "ComplexNumber.cpp"
#undef _COMPLEX_NUMBER_
#undef ComplexNumber
#undef DataTypeComplex
 