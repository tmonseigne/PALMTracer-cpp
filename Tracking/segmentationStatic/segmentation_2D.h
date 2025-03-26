
/*******************************************************************************************************************
								Data segmentation for image (2D signal) (16bit only)
									Institut Curie
									UMR - 144
									by Victor Racine
									04/11/02 
*******************************************************************************************************************/

#ifndef __segmentation_2D__
#define __segmentation_2D__

#include "segmentation.h"
#include "info_2D.h"


class Segment2D:public Segment{
private:
	InfoStruct2D* headInfo;
	InfoStruct2D** tabInfo;
	//infoStruct2D *infoStruct;
	//unsigned long nbInfoStruct;

	int alloc(unsigned long nbStruct);
	void desalloc();

	int makeLinks(unsigned long nbStruct);

	int filterStructureWS(FltArray &potential, USrtArray &seg, float thresholdRatio);
	int getLocalMax(FltArray &image, USrtArray &seg, float threshold);
	int watershed(FltArray &potential, USrtArray &seg, float threshold);

	int filterStructureWS(USrtArray &potential, USrtArray &seg, float thresholdRatio);
	int getLocalMax(USrtArray &image, USrtArray &seg, float threshold);
	int watershed(USrtArray &potential, USrtArray &seg, float threshold);
	int computeMaxLength(USrtArray &seg, long maxThreshold);
	unsigned short segmentationPositiveValues(USrtArray &data, USrtArray &seg, double val, bool diagonal);
	unsigned short segmentationPositiveValues(USrtArray &data, USrtArray &seg, double val, bool diagonal, USrtArray &segHoles);
	unsigned short segmentationNegativeValues(USrtArray &data, USrtArray &seg, double val, bool diagonal);

	unsigned short segmentationPositiveValues(FltArray &data, USrtArray &seg, double val, bool diagonal);
	unsigned short segmentationPositiveValues(FltArray &data, USrtArray &seg, double val, bool diagonal, USrtArray &segHoles);
	unsigned short segmentationNegativeValues(FltArray &data, USrtArray &seg, double val, bool diagonal);
	
public:
	Segment2D();
	~Segment2D();
	unsigned long  filterData(USrtArray &seg, double* param, unsigned short* useParam);
	unsigned short segmentation(FltArray &data, USrtArray &seg, double val, bool diagonal, bool fillHoles);
	unsigned short segmentation(USrtArray &data, USrtArray &seg, double val, bool diagonal, bool fillHoles);
	
	
	int computeStatistics(FltArray &data, USrtArray &seg);
	int computeStatistics(USrtArray &data, USrtArray &seg);
	void drawSymbolOnStructures(USrtArray &seg, unsigned short typeOfSign);
	int printStatistics(wchar_t* fileName, unsigned short resultType, unsigned long time, unsigned long nbTime, unsigned short atEndOfFile,
		unsigned long cptPoints, unsigned short sizeX, unsigned short sizeY, unsigned short sizeStk);
	//unsigned short readStatistics(wchar_t* fileName);
	void** getTabInfo();

	
	int watershedIntelligent(FltArray &potential, USrtArray &seg, float threshold, float thresholdRatio);
	int watershedIntelligent(USrtArray &potential, USrtArray &seg, float threshold, float thresholdRatio);
	void mixIndex(USrtArray &seg);
};


#endif //end __segmentation_2D__


