/*******************************************************************************************************************
								Control of allocation of arrayType
*******************************************************************************************************************/

#include "../generalParam/generalParam.h"


#ifdef	DEBUG_AT


#ifndef CONTROL_ALLOC_ARRAY_TYPE
#define CONTROL_ALLOC_ARRAY_TYPE


#include <stdio.h>
#include "print_debug.h"

class PtrListArrayType{
public:
	PtrListArrayType(){};
	PtrListArrayType(void* _ptr);
	void* ptr;
	PtrListArrayType* next;
};


class ControlAllocArrayType{
private:
	long nbArrayAllocated;
	PtrListArrayType* headPtr;
public:
	ControlAllocArrayType();
	~ControlAllocArrayType();
	int allocNewArray(void *_ptr);
	int deallocArray(void *_ptr);
};

#endif //CONTROL_ALLOC_ARRAY_TYPE
#endif //DEBUG_AT
