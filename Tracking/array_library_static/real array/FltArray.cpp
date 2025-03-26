#include "stdafx.h"

/*******************************************************************************************************************
								Pactage for general floating point array use
*******************************************************************************************************************/

#include "allRealDataLib.h"


#define	_ARRAY_TYPE_
#define DataType	float
#define ArrayType	FltArray
#include "ArrayType.cpp"
#undef DataType
#undef ArrayType
#undef _ARRAY_TYPE_
