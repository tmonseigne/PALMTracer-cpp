#include "stdafx.h"

/*******************************************************************************************************************
								Pactage for general floating point 64bit array
*******************************************************************************************************************/


#include "allRealDataLib.h"


#define	_ARRAY_TYPE_
#define DataType	double
#define ArrayType	DblArray
#include "ArrayType.cpp"
#undef DataType
#undef ArrayType
#undef _ARRAY_TYPE_
