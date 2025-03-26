


/*******************************************************************************************************************
								Pactage for general complex integer 32bit array
									Institut Curie
									UMR - 144
									by Victor Racine
									sept, 23rd 2003
*******************************************************************************************************************/



#include "allComplexDataLib.h"


#define	_ARRAY_TYPE_
#define DataType	ComplexLng
#define ArrayType	CxLngArray
#define COMPLEX
#define DataTypeSimple	long
#include "../real array/ArrayType.cpp"
#undef DataTypeSimple
#undef COMPLEX
#undef DataType
#undef ArrayType
#undef _ARRAY_TYPE_

