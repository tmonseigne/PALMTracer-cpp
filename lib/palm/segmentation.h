#pragma once

#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <malloc.h>
#include <string.h>
#include <time.h>


#include "tools/allDataLib.h"
#include "tools/lists.h"
//#include "../ReadWriteSTKNew/ReadWriteSTK.h"
#include "tools/generalParam.h"

#define MAX_TAB	10000

#define SIZE_CROSS	4


class Segment
{
protected:
	unsigned long nbInfoStruct;

	void inc(unsigned long& nSeg, unsigned long** ptr, long& sizePtr);

public:
	Segment();
	virtual ~Segment();
	unsigned long getNbStructs();
};
