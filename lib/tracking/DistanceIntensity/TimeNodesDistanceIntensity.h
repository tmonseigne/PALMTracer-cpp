/*******************************************************************************************************************
									TimeNodesDistanceIntensity Class
									It is derived from the class TimeNodes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/
#pragma once

#include "tracking/TimeNodes.h"


class TimeNodesDistanceIntensity final : public TimeNodes
{
protected:
	Node* allocNode() override;

public:
	TimeNodesDistanceIntensity(int time, int nNodes, Parameters* params);
	void ComputeConstantSigmas(double& stdIntensity, long& countStdInt, double& stdMean, long& countStdMean, double& stdDistance,
	                           long& countStdDist, bool sigmaClipping, long* hist = nullptr) const;
	double MakeAllInitializationConnections() const;
	double TryRandomChange(double limit) const;

	double GetLikenessTime(long lastTime) const;
};
