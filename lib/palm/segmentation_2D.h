#pragma once

#include "Windows.h"
#include "segmentation.h"
#include "info_2D.h"
#include "tools/allDataLib.h"


class Segment2D : public Segment
{
private:
	InfoStruct2D* headInfo;
	InfoStruct2D** tabInfo;


	int alloc(unsigned long nbStruct);
	void desalloc();

	int filterStructureWS(FltArray& potential, USrtArray& seg, float thresholdRatio);
	static int getLocalMax(FltArray& image, USrtArray& seg, float threshold);
	int watershed(FltArray& potential, USrtArray& seg, float threshold) const;

	int filterStructureWS(USrtArray& potential, USrtArray& seg, float thresholdRatio);
	static int getLocalMax(USrtArray& image, USrtArray& seg, float threshold);
	int watershed(USrtArray& potential, USrtArray& seg, float threshold) const;

	static float tic();

public:
	Segment2D();
	Segment2D(unsigned long _nbInfoStruct);
	~Segment2D() override;

	int computeStatistics(FltArray& data, USrtArray& seg, float* waveData, unsigned short gaussFit, bool doRealoc = true);

	unsigned long filterData(USrtArray& seg, double* param, unsigned short* useParam);
	void** getTabInfo() const;

	int watershedIntelligent(FltArray& potential, USrtArray& seg, float threshold, float thresholdRatio);
	void mixIndex(USrtArray& seg) const;
};
