/*******************************************************************************************************************
									TabTimesDistanceIntensity Class
									It is derived from the class TabTimes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/
#pragma once

#include "tracking/tab_times.h"

class TabTimesDistanceIntensity : public TabTimes
{
protected:
	TimeNodes* allocTimeNodes(long time, long nbNodes) override;
	void updateConstantes() const;
	void computeConstantes() const;

public:
	TabTimesDistanceIntensity(PublicParameters* publicParameters);
	double getLikenessGlobal() const;
	int runModel() override;
};
