
/*******************************************************************************************************************
								Data segmentation for 3D image or image stack  (3D signal) (16bit only)
									Institut Curie
									UMR - 144
									by Victor Racine
									04/15/02 
*******************************************************************************************************************/

#ifndef __segmentation_3D__
#define __segmentation_3D__

#include "segmentation.h"
#include "linear_algebra.h"
#include "info_3D.h"

extern FILE *debugStreamALS;

class Segment3D:public Segment{
private:
	InfoStruct3D* headInfo;
	InfoStruct3D** tabInfo;
	//unsigned long nbInfoStruct;

	int alloc(unsigned long nbStruct);
	void desalloc();

	int makeLinks(unsigned long nbStruct);
	void getSize(long &sizeDataX, long &sizeDataY, long &sizeDataZ, USrtArray &tab);

	int filterStructureWS(FltArray &potential, USrtArray &seg, float thresholdRatio);
	int getLocalMax(FltArray &image, USrtArray &seg, float threshold);
	int watershed(FltArray &potential, USrtArray &seg, float threshold);

	int filterStructureWS(USrtArray &potential, USrtArray &seg, float thresholdRatio);
	int getLocalMax(USrtArray &image, USrtArray &seg, float threshold);
	int watershed(USrtArray &potential, USrtArray &seg, float threshold);

	unsigned short segmentationPositiveValues(USrtArray &data, USrtArray &seg, double val, bool diagonal);
	unsigned short segmentationPositiveValues(USrtArray &data, USrtArray &seg, double val, bool diagonal, USrtArray &segHoles);
	unsigned short segmentationNegativeValues(USrtArray &data, USrtArray &seg, double val, bool diagonal);

	unsigned short segmentationPositiveValues(FltArray &data, USrtArray &seg, double val, bool diagonal);
	unsigned short segmentationPositiveValues(FltArray &data, USrtArray &seg, double val, bool diagonal, USrtArray &segHoles);
	unsigned short segmentationNegativeValues(FltArray &data, USrtArray &seg, double val, bool diagonal);
public:
	Segment3D();
	~Segment3D();
	unsigned long filterData(USrtArray &seg, double* param, unsigned short* useParam);
	unsigned short segmentation(FltArray &data, USrtArray &seg, double val, bool diagonal, bool fillHoles);
	unsigned short segmentation(USrtArray &data, USrtArray &seg, double val, bool diagonal, bool fillHoles);
	int computeStatistics(FltArray &data, USrtArray &seg);
	int computeStatistics(USrtArray &data, USrtArray &seg);
	void drawSymbolOnStructures(USrtArray &seg, unsigned short typeOfSign);
	int printStatistics(wchar_t* fileName, unsigned short resultType, unsigned long time, unsigned long nbTime, unsigned short atEndOfFile,
		unsigned long cptPoints, unsigned short sizeX,
		unsigned short sizeY, unsigned short sizeStk);
	//unsigned short readStatistics(wchar_t* fileName);
	void** getTabInfo();


	int watershedIntelligent(FltArray &potential, USrtArray &seg, float threshold, float thresholdRatio);
	int watershedIntelligent(USrtArray &potential, USrtArray &seg, float threshold, float thresholdRatio);
	void mixIndex(USrtArray &seg);
};




#endif //end __segmentation_3D__
