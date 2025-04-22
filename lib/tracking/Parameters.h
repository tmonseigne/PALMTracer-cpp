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
	static constexpr size_t N_TEMPORAL_WINDOWS = 1;
	static constexpr size_t MODEL = 2; //the number of the used model

	// Params
	double maxSpeed_;
	bool allowFusionAndFission_;
	double costBirthDeath_;

	// Computes
	double stdIntensity_, stdMean_, stdDistance_, covMeanIntensity_;
	double maxX_, maxY_, maxZ_;
};
