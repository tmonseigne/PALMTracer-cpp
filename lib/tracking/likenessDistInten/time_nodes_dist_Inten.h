/*******************************************************************************************************************
									TimeNodesDistanceIntensity Class
									It is derived from the class TimeNodes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/
#pragma once

#include "tracking/time_nodes.h"


class TimeNodesDistanceIntensity : public TimeNodes
{
protected:
	Node* allocNode() override;

public:
	TimeNodesDistanceIntensity(int _time, int _nbNodes, PublicParameters* _publicParameters);
	void computeConstantSigmas(double& stdIntensityComputing, long& countStdInt, double& stdMeanComputing, long& countStdMean, double& stdDistanceComputing, long& countStdDist, bool sigmaClipping, long* hist = NULL);
	double makeAllInitializationConnections();
	double tryRandomChange(double limit);

	double getLikenessTime(long lastTime);
};
