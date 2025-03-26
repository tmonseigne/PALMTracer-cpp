


/*******************************************************************************************************************
								Pactage for general complex floating point 64bit array
									Institut Curie
									UMR - 144
									by Victor Racine
									sept, 23rd 2003
*******************************************************************************************************************/



#include "allComplexDataLib.h"


#define	_ARRAY_TYPE_
#define DataType	ComplexDbl
#define ArrayType	CxDblArray
#define COMPLEX
#define DataTypeSimple	double
#include "../real array/ArrayType.cpp"
#undef DataTypeSimple
#undef COMPLEX
#undef DataType
#undef ArrayType
#undef _ARRAY_TYPE_

