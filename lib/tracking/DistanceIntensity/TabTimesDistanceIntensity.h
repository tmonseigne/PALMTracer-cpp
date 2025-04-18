/*******************************************************************************************************************
									TabTimesDistanceIntensity Class
									It is derived from the class Tracking
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/
#pragma once

#include "tracking/Tracking.h"

class TabTimesDistanceIntensity final : public Tracking
{
protected:
	TimeNodes* allocTimeNodes(long time, long nNodes) override;
	void updateConstantes() const;
	void computeConstantes() const;

public:
	TabTimesDistanceIntensity(Parameters* params) : Tracking(params) { }
	double GetLikenessGlobal() const;
	int RunModel() override;
};
