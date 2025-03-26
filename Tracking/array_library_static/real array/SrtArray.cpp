/*******************************************************************************************************************
								Pactage for general interger 16bit array
*******************************************************************************************************************/


#include "allRealDataLib.h"

#define	_ARRAY_TYPE_
#define DataType	short
#define ArrayType	SrtArray
#define INTEGER
#define MAXVAL 32767
#define MINVAL -32768
#include "ArrayType.cpp"
#undef DataType
#undef ArrayType
#undef _ARRAY_TYPE_
#undef INTEGER
#undef MAXVAL
#undef MINVAL
