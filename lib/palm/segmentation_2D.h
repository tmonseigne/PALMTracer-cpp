#ifndef __segmentation_2D__
#define __segmentation_2D__

#include "Windows.h"
#include "segmentation.h"
#include "info_2D.h"


class Segment2D : public Segment
{
private:
	InfoStruct2D* headInfo;
	InfoStruct2D** tabInfo;


	int alloc(unsigned long nbStruct);
	void desalloc();

	int filterStructureWS(FltArray& potential, USrtArray& seg, float thresholdRatio);
	int getLocalMax(FltArray& image, USrtArray& seg, float threshold);
	int watershed(FltArray& potential, USrtArray& seg, float threshold);

	int filterStructureWS(USrtArray& potential, USrtArray& seg, float thresholdRatio);
	int getLocalMax(USrtArray& image, USrtArray& seg, float threshold);
	int watershed(USrtArray& potential, USrtArray& seg, float threshold);

	unsigned short segmentationPositiveValues(FltArray& data, USrtArray& seg, double val, bool diagonal);
	unsigned short segmentationPositiveValues(FltArray& data, USrtArray& seg, double val, bool diagonal, USrtArray& segHoles);
	unsigned short segmentationNegativeValues(FltArray& data, USrtArray& seg, double val, bool diagonal);

	float tic();

public:
	Segment2D();
	Segment2D(unsigned long _nbInfoStruct);
	~Segment2D() override;

	unsigned short segmentation(FltArray& data, USrtArray& seg, double val, bool diagonal, bool fillHoles);

	//int computeStatistics(FltArray &data, USrtArray &seg, bool doRealoc=true);
	int computeStatistics(FltArray& data, USrtArray& seg, float* waveData, unsigned short gaussFit, bool doRealoc = true);

	unsigned long filterData(USrtArray& seg, double* param, unsigned short* useParam);
	void** getTabInfo();

	int watershedIntelligent(FltArray& potential, USrtArray& seg, float threshold, float thresholdRatio);
	void mixIndex(USrtArray& seg);
};


#endif //end __segmentation_2D__
