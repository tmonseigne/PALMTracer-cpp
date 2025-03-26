

#include "stdafx.h"

#include "allRealDataLib.h"


#define	_ARRAY_TYPE_
#define DataType	float
#define ArrayType	FltArray
#include "ArrayType.cpp"
#pragma message("after ArrayType")
#undef DataType
#undef ArrayType
#undef _ARRAY_TYPE_
