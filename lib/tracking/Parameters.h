#pragma once

class Parameters
{
public:
	Parameters() = default;
	~Parameters() = default;

	// Static
	static constexpr double DIST_MAX = 20;
	static constexpr size_t N_ITERATIONS = 300;
	static constexpr size_t N_SEG_PARAMS = 13;
	static constexpr size_t N_DIM = 2;
	static constexpr size_t N_TEMPORAL_WINDOWS = 1;

	// Allocation
	int nNodeAllocated_;
	int nNodeDeallocated_;
	int nChainNeighborAllocated_;
	int nChainNeighborDeallocated_;

	// IN/Out
	double* points_;
	double* tracks_;
	size_t nPoints_, nTracks_;
	size_t nTimes_;
	size_t nTraj_;

	// Params
	double maxSpeed_;
	double dZdx_;
	size_t model_; //the number of the used model
	double decrease_;
	size_t minLife_;
	bool allowFusionAndFission_;
	double costBirthDeath_;

	// Computes
	double stdIntensity_, stdMean_, stdDistance_, covMeanIntensity_;
	double maxX_, maxY_, maxZ_;
};
