


/*******************************************************************************************************************
								Pactage for general interger 8bit array
*******************************************************************************************************************/


#ifndef		_COMPLEX_CHR_
#define		_COMPLEX_CHR_ 

	#define	_ARRAY_TYPE_
		#undef DataType
		#define DataType	char
			#define ArrayType	ChrArray
				#define INTEGER
					#define NB_TYPE	1
						#define TYPE_NAME L"ChrArray"
							#include "ArrayType.h"
						#undef TYPE_NAME
					#undef NB_TYPE
				#undef INTEGER
			#undef ArrayType
		#undef DataType
	#undef _ARRAY_TYPE_

#endif 
