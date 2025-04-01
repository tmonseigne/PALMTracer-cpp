#pragma once

#include <io.h>
#include <fcntl.h>
#include <string.h>


#include "tools/lists.h"
//#include "../ReadWriteSTKNew/ReadWriteSTK.h"

#define MAX_TAB	10000

#define SIZE_CROSS	4


class Segment
{
protected:
	unsigned long nbInfoStruct;

	static void inc(unsigned long& nSeg, unsigned long** ptr, long& sizePtr);

public:
	Segment();
	virtual ~Segment();
	unsigned long getNbStructs() const;
};
