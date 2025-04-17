#pragma once

#include <vector>
#include "Info2D.h"


class Segmentation2D
{
public:
	Segmentation2D(const size_t height, const size_t width) : height_(height), width_(width), size_(height * width) { }

	~Segmentation2D() { desalloc(); }

	void ComputeStatistics(const std::vector<double>& data, const std::vector<unsigned short>& seg, const std::vector<double>& waveData, size_t gaussFit);

	void FilterData(std::vector<unsigned short>& seg, const double* param, const unsigned short* useParam);
	std::vector<Info2D> GetInfos() const { return infos_; }

	void WatershedIntelligent(const std::vector<double>& potential, std::vector<unsigned short>& seg, double threshold, double thresholdRatio);

private:
	size_t height_, width_, size_;
	std::vector<Info2D> infos_;

	void alloc(size_t nbStruct);
	void desalloc();

	void filterStructureWS(const std::vector<double>& potential, std::vector<unsigned short>& seg, double thresholdRatio);
	void getLocalMax(const std::vector<double>& potential, std::vector<unsigned short>& seg, double threshold);
	void watershed(const std::vector<double>& potential, std::vector<unsigned short>& seg, double threshold) const;
};
