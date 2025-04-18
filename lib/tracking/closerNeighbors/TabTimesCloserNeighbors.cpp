/*******************************************************************************************************************
									TabTimesCloserNeighbors Class
									It is derived from the class TabTimes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/


#include "TabTimesCloserNeighbors.h"
#include "TimeNodesCloserNeighbors.h"

/********************************************************************************************************************
Constructuor() of TabTimesCloserNeighbors class.
********************************************************************************************************************/
TabTimesCloserNeighbors::TabTimesCloserNeighbors(Parameters* params): TabTimes(params) { parameters_ = params; }


/********************************************************************************************************************
allocTabTimes(int _nbTimes) alloc an array of _nbTimes long 
********************************************************************************************************************/
TimeNodes* TabTimesCloserNeighbors::allocTimeNodes(const long time, const long nbNodes) { return new TimeNodesCloserNeighbors(time, nbNodes, parameters_); }

/********************************************************************************************************************
runModel() is the core of the application, it links the structures according the selected model
********************************************************************************************************************/
int TabTimesCloserNeighbors::RunModel()
{
	int ret;
	for (int time = 0; time < nbTimes_ - 1; time++) { if ((ret = ((TimeNodesCloserNeighbors*)timeNodes_[time])->MakeAllConnections()) < 0) { return ret; } }
	return RET_OK;
}
