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
TabTimesCloserNeighbors::TabTimesCloserNeighbors(Parameters* _publicParameters): TabTimes(_publicParameters) { parameters_ = _publicParameters; }


/********************************************************************************************************************
allocTabTimes(int _nbTimes) alloc an array of _nbTimes long 
********************************************************************************************************************/
TimeNodes* TabTimesCloserNeighbors::allocTimeNodes(const long _time, const long _nbNodes) { return new TimeNodesCloserNeighbors(_time, _nbNodes, parameters_); }

/********************************************************************************************************************
runModel() is the core of the application, it links the structures according the selected model
********************************************************************************************************************/
int TabTimesCloserNeighbors::RunModel()
{
	int ret;
	for (int _time = 0; _time < nbTimes_ - 1; _time++) { if ((ret = ((TimeNodesCloserNeighbors*)timeNodes_[_time])->makeAllConnections()) < 0) { return ret; } }
	return RET_OK;
}
