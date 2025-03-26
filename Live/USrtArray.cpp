#include "stdafx.h"
#include "allRealDataLib.h"


#define	_ARRAY_TYPE_
#define DataType	unsigned short
#define ArrayType	USrtArray
#define UNSIGNED
#define INTEGER
#define MAXVAL 65535
#define MINVAL 0
#include "ArrayType.cpp"
#undef UNSIGNED
#undef DataType
#undef ArrayType
#undef _ARRAY_TYPE_
#undef INTEGER
#undef MAXVAL
#undef MINVAL
