#pragma once

#include <vector>
#include "info_2D.h"
#include "Segmentation.h"


class Segmentation2D final : public Segmentation
{
public:
	Segmentation2D(const size_t height, const size_t width)
		: height_(height), width_(width), size_(height * width), headInfo_(nullptr) { nbInfoStruct_ = 0; }

	~Segmentation2D() override { desalloc(); }

	int ComputeStatistics(const std::vector<double>& data, const std::vector<unsigned short>& seg, const std::vector<double>& waveData, size_t gaussFit, bool doRealoc = true);

	size_t FilterData(std::vector<unsigned short>& seg, const double* param, const unsigned short* useParam);
	std::vector<InfoStruct2D*> GetTabInfo() const { return tabInfo_; }

	int WatershedIntelligent(const std::vector<double>& potential, std::vector<unsigned short>& seg, double threshold, double thresholdRatio);

private:
	size_t height_, width_, size_;
	InfoStruct2D* headInfo_;
	std::vector<InfoStruct2D*> tabInfo_;

	int alloc(size_t nbStruct);
	void desalloc();

	int filterStructureWS(const std::vector<double>& potential, std::vector<unsigned short>& seg, double thresholdRatio);
	unsigned short getLocalMax(const std::vector<double>& potential, std::vector<unsigned short>& seg, double threshold) const;
	size_t watershed(const std::vector<double>& potential, std::vector<unsigned short>& seg, double threshold) const;
};
