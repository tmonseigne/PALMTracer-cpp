
/*******************************************************************************************************************
								Data segmentation for image (2D & 3D signal) (16bit only)
									Institut Curie
									UMR - 144
									by Victor Racine
									04/11/02 
*******************************************************************************************************************/

#ifndef __segmentation__
#define __segmentation__


#include <io.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <malloc.h>
#include <string.h>
#include <time.h>


#include "../array_library_static/allDataLib.h"
#include "../miscTools/lists.h" 
#include "../ReadWriteSTKNew/ReadWriteSTK.h"
#include "../generalParam\generalParam.h"
#define MAX_TAB	10000

#define SIZE_CROSS	4


class Segment{
protected:
	//infoStruct2D *infoStruct;
	unsigned long nbInfoStruct;

	//int alloc(unsigned long nbStruct);
	void inc( unsigned long &nSeg, unsigned long **ptr, long &sizePtr );
	int locateStringInFile(wchar_t* st, FILE* file);
	void printHTMLCell(wchar_t* st, FILE* file, long countRaw, long countCol);
	void printFloatHTMLCell(double f, FILE* file);
	void printIntHTMLCell(long i, FILE* file, long count=1);

	virtual int filterStructureWS(FltArray &potential, USrtArray &seg, float thresholdRatio)=0;
	virtual int getLocalMax(FltArray &image, USrtArray &seg, float threshold)=0;
	virtual int watershed(FltArray &potential, USrtArray &seg, float threshold)=0;

	virtual int filterStructureWS(USrtArray &potential, USrtArray &seg, float thresholdRatio)=0;
	virtual int getLocalMax(USrtArray &image, USrtArray &seg, float threshold)=0;
	virtual int watershed(USrtArray &potential, USrtArray &seg, float threshold)=0;

public:
	Segment();
	virtual ~Segment();
	virtual unsigned long  filterData(USrtArray &seg, double* param, unsigned short* useParam)=0;
	virtual int alloc(unsigned long nbStruct)=0;
	virtual unsigned short segmentation(FltArray &data, USrtArray &seg, double val, bool diaglonal, bool fillHoles)=0;
	virtual unsigned short segmentation(USrtArray &data, USrtArray &seg, double val, bool diaglonal, bool fillHoles)=0;
	virtual int computeStatistics(FltArray &data, USrtArray &seg)=0;
	virtual int computeStatistics(USrtArray &data, USrtArray &seg)=0;
	virtual void drawSymbolOnStructures(USrtArray &seg, unsigned short typeOfSign)=0;
	virtual int printStatistics(wchar_t* fileName, unsigned short resultType, unsigned long time, unsigned long nbTime,
		unsigned short atEndOfFile,
		unsigned long cptPoints, unsigned short sizeX,
		unsigned short sizeY, unsigned short sizeStk)=0;
	//virtual unsigned short readStatistics(wchar_t* fileName)=0;
	virtual void** getTabInfo()=0;

	virtual int watershedIntelligent(FltArray &potential, USrtArray &seg, float threshold, float thresholdRatio)=0;
	virtual int watershedIntelligent(USrtArray &potential, USrtArray &seg, float threshold, float thresholdRatio)=0;
};

#endif
