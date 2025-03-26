#include "segmentation.h"
int __nbInfoStrucAlloc;
int __nbInfoStrucDealloc;

Segment::Segment()
{
	__nbInfoStrucAlloc = 0;
	__nbInfoStrucDealloc = 0;
}

Segment::~Segment() { }


void Segment::inc(unsigned long& nSeg, unsigned long** ptr, long& sizePtr)
{
	unsigned long* newPtr;
	int i;
	if (nSeg == 65535) return;
	nSeg++;
	if ((long)nSeg < sizePtr) return;


	newPtr = (unsigned long*)malloc(sizeof(long) * sizePtr);
	//backup of data
	for (i = 0; i < (long)sizePtr; i++) newPtr[i] = (*ptr)[i];

	//reallocation of ptr
	sizePtr += MAX_TAB;
	(*ptr) = (unsigned long*)realloc((*ptr), sizeof(long) * sizePtr);

	//copy back the data
	for (i = 0; i < ((long)sizePtr - MAX_TAB); i++) (*ptr)[i] = newPtr[i];
	for (i = sizePtr - MAX_TAB; i < (long)sizePtr; i++) (*ptr)[i] = i;
	free(newPtr);
}

/*************************************************************************************************************************************
locateStringInFile try to find the string st in the stream "file". This file have to be opened.
*************************************************************************************************************************************/
/*int Segment::locateStringInFile(wchar_t* st, FILE* file){
	bool stingFounded=false;
	wchar_t line[1024];
	wchar_t *res;
	long offs;
	if(fseek(file,-30,SEEK_END)!=0)	//change file position
		return(S_READ_ERR);
	while(fgetws(line, 256, file)!=NULL){  //reading a new line
		if((res=wcsstr(line,st))!=NULL){	//if st exists in "line"
			offs=ftell(file)-wcslen(line)+res-line-1;
			return(offs);
		}
	}
	return(S_READ_ERR);
}*/

unsigned long Segment::getNbStructs() { return (nbInfoStruct); }
