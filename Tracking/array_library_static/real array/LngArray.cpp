#include "stdafx.h"


/*******************************************************************************************************************
								Pactage for general interger 32bit array
*******************************************************************************************************************/


#include "allRealDataLib.h"


#define	_ARRAY_TYPE_
#define DataType	long
#define ArrayType	LngArray
#define INTEGER
#define MAXVAL 2147483647
#define MINVAL -2147483647
#include "ArrayType.cpp"
#undef DataType
#undef ArrayType
#undef _ARRAY_TYPE_
#undef INTEGER
#undef MAXVAL
#undef MINVAL
