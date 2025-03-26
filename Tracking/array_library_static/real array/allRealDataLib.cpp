#include "stdafx.h"


/*******************************************************************************************************************
								Load all real data types
*******************************************************************************************************************/
/*

class FltArray;

#include "allRealDataLib.h"



#ifdef		_CHR_ARRAY_

	#define	_ARRAY_TYPE_
	#define DataType	char
	#define ArrayType	ChrArray
	#define NB_TYPE	1
	#include "ArrayType.cpp"
	#undef _ARRAY_TYPE_
	#undef DataType
	#undef ArrayType
	#undef NB_TYPE
#endif

#ifdef		_UCHR_ARRAY_
	#define	_ARRAY_TYPE_
	#define UNSIGNED	
	#define NB_TYPE	2
	#define DataType	unsigned char
	#define ArrayType	UChrArray
	#include "ArrayType.cpp"
	#undef UNSIGNED
	#undef _ARRAY_TYPE_
	#undef DataType
	#undef ArrayType
	#undef NB_TYPE

#endif



#ifdef		_SRT_ARRAY_
	#define	_ARRAY_TYPE_
	#define NB_TYPE	3
	#define DataType	short
	#define ArrayType	SrtArray
	#include "ArrayType.cpp"
	#undef DataType
	#undef ArrayType
	#undef _ARRAY_TYPE_
	#undef NB_TYPE

#endif


#ifdef		_USRT_ARRAY_
	#define	_ARRAY_TYPE_
	#define NB_TYPE	4
	#define DataType	unsigned short
	#define ArrayType	USrtArray
	#define UNSIGNED
	#include "ArrayType.cpp"
	#undef UNSIGNED
	#undef DataType
	#undef ArrayType
	#undef _ARRAY_TYPE_
	#undef NB_TYPE

#endif

#ifdef		_LNG_ARRAY_
	#define	_ARRAY_TYPE_
	#define NB_TYPE	5
	#define DataType	long
	#define ArrayType	LngArray
	#include "ArrayType.cpp"
	#undef DataType
	#undef ArrayType
	#undef _ARRAY_TYPE_
	#undef NB_TYPE

#endif

#ifdef		_ULNG_ARRAY_
	#define	_ARRAY_TYPE_
	#define NB_TYPE	6
	#define DataType	unsigned long
	#define ArrayType	ULngArray
	#define UNSIGNED
	#include "ArrayType.cpp"
	#undef DataType
	#undef ArrayType
	#undef UNSIGNED
	#undef _ARRAY_TYPE_
	#undef NB_TYPE

#endif


#ifdef		_FLT_ARRAY_
	#define	_ARRAY_TYPE_
	#define NB_TYPE	7
	#define DataType	float
	#define ArrayType	FltArray
	#include "ArrayType.cpp"
	#undef DataType
	#undef ArrayType
	#undef _ARRAY_TYPE_
	#undef NB_TYPE
#endif


#ifdef		_DBL_ARRAY_

	#define	_ARRAY_TYPE_
	#define NB_TYPE	8
	#define DataType	double
	#define ArrayType	DblArray
	#include "ArrayType.cpp"
	#undef DataType
	#undef ArrayType
	#undef _ARRAY_TYPE_
	#undef NB_TYPE

#endif
*/


/*

#ifdef _COPY_CHR_
	#define _COPY_TYPE_
	#define DestData	char
	#define DestType	ChrArray

	#define SrcType	ChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	UChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	SrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	USrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	LngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	ULngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	FltArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	DblArray
	#include "transTypes.cpp"
	#undef SrcType

	#undef DestType
	#undef DestData
	#undef _COPY_TYPE_
#endif*/
/*

#ifdef _COPY_UCHR_
	#define _COPY_TYPE_
	#define DestData	unsigned char
	#define DestType	UChrArray

	#define SrcType	ChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	UChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	SrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	USrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	LngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	ULngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	FltArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	DblArray
	#include "transTypes.cpp"
	#undef SrcType

	#undef DestType
	#undef DestData
	#undef _COPY_TYPE_
#endif


#ifdef _COPY_SRT_
	#define _COPY_TYPE_
	#define DestData	short
	#define DestType	SrtArray

	#define SrcType	ChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	UChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	SrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	USrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	LngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	ULngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	FltArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	DblArray
	#include "transTypes.cpp"
	#undef SrcType

	#undef DestType
	#undef DestData
	#undef _COPY_TYPE_
#endif

#ifdef _COPY_USRT_
	#define _COPY_TYPE_
	#define DestData	unsigned short
	#define DestType	USrtArray

	#define SrcType	ChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	UChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	SrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	USrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	LngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	ULngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	FltArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	DblArray
	#include "transTypes.cpp"
	#undef SrcType

	#undef DestType
	#undef DestData
	#undef _COPY_TYPE_
#endif

#ifdef _COPY_LNG_
	#define _COPY_TYPE_
	#define DestData	long
	#define DestType	LngArray

	#define SrcType	ChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	UChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	SrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	USrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	LngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	ULngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	FltArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	DblArray
	#include "transTypes.cpp"
	#undef SrcType

	#undef DestType
	#undef DestData
	#undef _COPY_TYPE_
#endif


#ifdef _COPY_ULNG_
	#define _COPY_TYPE_
	#define DestData	unsigned long
	#define DestType	ULngArray

	#define SrcType	ChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	UChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	SrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	USrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	LngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	ULngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	FltArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	DblArray
	#include "transTypes.cpp"
	#undef SrcType

	#undef DestType
	#undef DestData
	#undef _COPY_TYPE_
#endif


#ifdef _COPY_FLT_
	#define _COPY_TYPE_
	#define DestData	float
	#define DestType	FltArray

	#define SrcType	ChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	UChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	SrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	USrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	LngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	ULngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	FltArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	DblArray
	#include "transTypes.cpp"
	#undef SrcType

	#undef DestType
	#undef DestData
	#undef _COPY_TYPE_
#endif


#ifdef _COPY_DBL_
	#define _COPY_TYPE_
	#define DestData	double
	#define DestType	DblArray

	#define SrcType	ChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	UChrArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	SrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	USrtArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	LngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	ULngArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	FltArray
	#include "transTypes.cpp"
	#undef SrcType
	#define SrcType	DblArray
	#include "transTypes.cpp"
	#undef SrcType

	#undef DestType
	#undef DestData
	#undef _COPY_TYPE_
#endif


*/
