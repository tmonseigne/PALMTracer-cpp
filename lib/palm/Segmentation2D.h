#pragma once

#include "Windows.h"
#include "Segmentation.h"
#include "info_2D.h"
#include "tools/allDataLib.h"


class Segmentation2D : public Segmentation
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
	Segmentation2D();
	Segmentation2D(unsigned long _nbInfoStruct);
	~Segmentation2D() override;

	int computeStatistics(FltArray& data, USrtArray& seg, float* waveData, unsigned short gaussFit, bool doRealoc = true);

	unsigned long filterData(USrtArray& seg, double* param, unsigned short* useParam);
	void** getTabInfo() const;

	int watershedIntelligent(FltArray& potential, USrtArray& seg, float threshold, float thresholdRatio);
	void mixIndex(USrtArray& seg) const;
};
