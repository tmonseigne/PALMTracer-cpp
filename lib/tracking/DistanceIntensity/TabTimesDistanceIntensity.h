/*******************************************************************************************************************
									TabTimesDistanceIntensity Class
									It is derived from the class TabTimes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/
#pragma once

#include "tracking/TabTimes.h"

class TabTimesDistanceIntensity final : public TabTimes
{
protected:
	TimeNodes* allocTimeNodes(long time, long nNodes) override;
	void updateConstantes() const;
	void computeConstantes() const;

public:
	TabTimesDistanceIntensity(Parameters* params) : TabTimes(params) { }
	double GetLikenessGlobal() const;
	int RunModel() override;
};
