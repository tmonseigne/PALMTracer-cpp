#include "segmentation.h"

Segment::Segment() { }

Segment::~Segment() { }


void Segment::inc(unsigned long& nSeg, unsigned long** ptr, long& sizePtr)
{
	int i;
	if (nSeg == 65535) { return; }
	nSeg++;
	if ((long)nSeg < sizePtr) { return; }


	unsigned long* newPtr = (unsigned long*)malloc(sizeof(long) * sizePtr);
	//backup of data
	for (i = 0; i < (long)sizePtr; i++) { newPtr[i] = (*ptr)[i]; }

	//reallocation of ptr
	sizePtr += MAX_TAB;
	*ptr = (unsigned long*)realloc(*ptr, sizeof(long) * sizePtr);

	//copy back the data
	for (i = 0; i < (long)sizePtr - MAX_TAB; i++) { (*ptr)[i] = newPtr[i]; }
	for (i = sizePtr - MAX_TAB; i < (long)sizePtr; i++) { (*ptr)[i] = i; }
	free(newPtr);
}

unsigned long Segment::getNbStructs() { return nbInfoStruct; }
