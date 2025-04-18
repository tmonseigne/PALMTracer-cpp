/*******************************************************************************************************************
									TabTimesCloserNeighbors Class
									It is derived from the class TabTimes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/
#pragma once
#include "tracking/tab_times.h"

class TabTimesCloserNeighbors : public TabTimes
{
protected:
	TimeNodes* allocTimeNodes(long time, long nbNodes) override;

public:
	TabTimesCloserNeighbors(Parameters* publicParameters);

	int runModel() override;
};
