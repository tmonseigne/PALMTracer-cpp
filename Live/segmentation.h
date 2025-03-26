#ifndef __segmentation__
#define __segmentation__


#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <malloc.h>
#include <string.h>
#include <time.h>


#include "tools/allDataLib.h"
#include "lists.h"
//#include "../ReadWriteSTKNew/ReadWriteSTK.h"
#include "generalParam.h"

#define MAX_TAB	10000

#define SIZE_CROSS	4


class Segment
{
protected:
	//infoStruct2D *infoStruct;
	unsigned long nbInfoStruct;

	//int alloc(unsigned long nbStruct);
	void inc(unsigned long& nSeg, unsigned long** ptr, long& sizePtr);
	int locateStringInFile(wchar_t* st, FILE* file);


	//virtual int filterStructureWS(FltArray &potential, USrtArray &seg, float thresholdRatio)=0;
	//virtual int getLocalMax(FltArray &image, USrtArray &seg, float threshold)=0;
	//virtual int watershed(FltArray &potential, USrtArray &seg, float threshold)=0;

	//virtual int filterStructureWS(USrtArray &potential, USrtArray &seg, float thresholdRatio)=0;
	//virtual int getLocalMax(USrtArray &image, USrtArray &seg, float threshold)=0;
	//virtual int watershed(USrtArray &potential, USrtArray &seg, float threshold)=0;

public:
	Segment();
	virtual ~Segment();
	//virtual unsigned long  filterData(USrtArray &seg, double* param, unsigned short* useParam)=0;
	//virtual int alloc(unsigned long nbStruct)=0;
	//virtual unsigned short segmentation(FltArray &data, USrtArray &seg, double val, bool diaglonal, bool fillHoles)=0;
	////virtual unsigned short segmentation(USrtArray &data, USrtArray &seg, double val, bool diaglonal, bool fillHoles)=0;
	//virtual int computeStatistics(FltArray &data, USrtArray &seg, bool doRealloc=true)=0;

	//virtual void** getTabInfo()=0;
	unsigned long getNbStructs();

	//virtual int watershedIntelligent(FltArray &potential, USrtArray &seg, float threshold, float thresholdRatio)=0;
};

#endif
