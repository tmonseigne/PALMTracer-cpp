#include "stdafx.h"

/*******************************************************************************************************************
								Control of allocation of arrayType
*******************************************************************************************************************/
#include "controlAlloc.h"
#ifdef	DEBUG_AT

ControlAllocArrayType* controlAllocArrayTypePublic=NULL;
extern FILE *debugStreamALS;

PtrListArrayType::PtrListArrayType(void* _ptr){
	ptr=_ptr;
	next=NULL;
}

ControlAllocArrayType::ControlAllocArrayType(){
	nbArrayAllocated=0;
	headPtr=NULL;
}

ControlAllocArrayType::~ControlAllocArrayType(){
	openDebugStreamALS();
	fprintf(debugStreamALS, "------Control of allocations------\n"); 
	if(nbArrayAllocated>0){
		fprintf(debugStreamALS, "---%d array are not deallocated\n", nbArrayAllocated);
		PtrListArrayType* ptr;
		for(ptr=headPtr;ptr!=NULL;ptr=ptr->next){
			fprintf(debugStreamALS, "---array not desallocated :0x%p\n", ptr->ptr);
		}
	}else
		fprintf(debugStreamALS, "---all arrays are desallocated----\n");
	fprintf(debugStreamALS, "------End of the control of allocations------\n"); 
	closeDebugStreamALS();
}


int ControlAllocArrayType::allocNewArray(void *_ptr){
	bool errorFinded=false;
	PtrListArrayType* ptr;
	for(ptr=headPtr;ptr!=NULL;ptr=ptr->next)
		if(_ptr==ptr->ptr){
			openDebugStreamALS();
			fprintf(debugStreamALS, "------Control of allocations------\n"); 
			fprintf(debugStreamALS, "----error the array 0x%p has already been allocated---\n");
			fprintf(debugStreamALS, "------End of the control of allocations------\n");
			closeDebugStreamALS();
			errorFinded=true;
		}
	ptr=headPtr;
	headPtr=new PtrListArrayType(_ptr);
	headPtr->next=ptr;
	nbArrayAllocated++;
	if(errorFinded)
		return(ERR_CONTROL_ALLOC_1);
	else
		return(RET_OK);
}

int ControlAllocArrayType::deallocArray(void *_ptr){
	PtrListArrayType* ptr;
	PtrListArrayType* ptrBefore;
	PtrListArrayType* ptrBeforeFinded=NULL;
	for(ptr=headPtr, ptrBefore=NULL;(ptr!=NULL)&&(ptrBeforeFinded==NULL);ptr=ptr->next){
		if(_ptr==ptr->ptr)
			ptrBeforeFinded=ptrBefore;
		ptrBefore=ptr;
	}
	if((ptrBeforeFinded!=NULL)||(_ptr==headPtr->ptr)){
		if(_ptr==headPtr->ptr){
			ptr=headPtr;
			headPtr=headPtr->next;
		}else{
			ptr=ptrBeforeFinded->next;
			ptrBeforeFinded->next=ptr->next;
		}
		delete ptr;
		nbArrayAllocated--;
	}else{
		openDebugStreamALS();
		fprintf(debugStreamALS, "------Control of allocations------\n"); 
		fprintf(debugStreamALS, "----error the array 0x%p is not allocated, deallocation is impossible---\n", _ptr);
		fprintf(debugStreamALS, "------End of the control of allocations------\n");
		closeDebugStreamALS();
		return(ERR_CONTROL_ALLOC_2);
	}
	return(RET_OK);
}
#endif //DEBUG_AT
