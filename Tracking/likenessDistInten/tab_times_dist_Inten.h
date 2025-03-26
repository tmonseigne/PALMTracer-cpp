/*******************************************************************************************************************
									TabTimesDistanceIntensity Class
									It is derived from the class TabTimes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/
#ifndef __tab_times_distance_intensity__
#define __tab_times_distance_intensity__

#include "../kernel/tab_times.h"
#include "time_nodes_dist_inten.h"

class TabTimesDistanceIntensity : public TabTimes
{
protected:
	TimeNodes* allocTimeNodes(long time, long nbNodes);
	void allocNode();
	void updateConstantes(void);
	void computeConstantes(void);

public:
	TabTimesDistanceIntensity(PublicParameters* publicParameters);
	~TabTimesDistanceIntensity();
	double getLikenessGlobal();
	int runModel();
};


#endif	//__tab_times_distance_intensity__
