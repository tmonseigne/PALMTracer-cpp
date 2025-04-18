/*******************************************************************************************************************
									TimeNodesDistanceIntensity Class
									It is derived from the class TimeNodes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/
#pragma once

#include "tracking/TimeNodes.h"


class TimeNodesDistanceIntensity : public TimeNodes
{
protected:
	Node* allocNode() override;

public:
	TimeNodesDistanceIntensity(int time, int nNodes, Parameters* params);
	void ComputeConstantSigmas(double& stdIntensityComputing, long& countStdInt, double& stdMeanComputing, long& countStdMean, double& stdDistanceComputing,
	                           long& countStdDist, bool sigmaClipping, long* hist = nullptr) const;
	double MakeAllInitializationConnections() const;
	double TryRandomChange(double limit) const;

	double GetLikenessTime(long lastTime) const;
};
