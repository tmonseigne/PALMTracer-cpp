/*******************************************************************************************************************
									TabTimesDistanceIntensity Class
									It is derived from the class TabTimes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/
#ifndef __tab_times_distance_intensity__
#define __tab_times_distance_intensity__

#include "time_nodes_dist_inten.h"
#include "tracking/tab_times.h"

class TabTimesDistanceIntensity : public TabTimes
{
protected:
	TimeNodes* allocTimeNodes(long time, long nbNodes) override;
	void updateConstantes();
	void computeConstantes();

public:
	TabTimesDistanceIntensity(PublicParameters* publicParameters);
	double getLikenessGlobal();
	int runModel() override;
};


#endif	//__tab_times_distance_intensity__
