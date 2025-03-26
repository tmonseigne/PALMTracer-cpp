#include "stdafx.h"


/*******************************************************************************************************************
								Pactage for general unsigned interger 32bit array
*******************************************************************************************************************/


#include "allRealDataLib.h"


#define	_ARRAY_TYPE_
#define DataType	unsigned long
#define ArrayType	ULngArray
#define UNSIGNED
#define INTEGER
#define MAXVAL 4294967291
#define MINVAL 0
#include "ArrayType.cpp"
#undef DataType
#undef ArrayType
#undef UNSIGNED
#undef _ARRAY_TYPE_
#undef INTEGER
#undef MAXVAL
#undef MINVAL
