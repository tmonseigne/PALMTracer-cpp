/*******************************************************************************************************************
									TabTimesCloserNeighbors Class
									It is derived from the class TabTimes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/


#include "tab_times_closer_neighbors.h"
#include "time_nodes_closer_neighbors.h"

/********************************************************************************************************************
Constructuor() of TabTimesCloserNeighbors class.
********************************************************************************************************************/
TabTimesCloserNeighbors::TabTimesCloserNeighbors(PublicParameters* _publicParameters): TabTimes(_publicParameters) { publicParameters = _publicParameters; }


/********************************************************************************************************************
allocTabTimes(int _nbTimes) alloc an array of _nbTimes long 
********************************************************************************************************************/
TimeNodes* TabTimesCloserNeighbors::allocTimeNodes(const long _time, const long _nbNodes) { return new TimeNodesCloserNeighbors(_time, _nbNodes, publicParameters); }

/********************************************************************************************************************
runModel() is the core of the application, it links the structures according the selected model
********************************************************************************************************************/
int TabTimesCloserNeighbors::runModel()
{
	int ret;
	for (int _time = 0; _time < nbTimes - 1; _time++) { if ((ret = ((TimeNodesCloserNeighbors*)timeNodes[_time])->makeAllConnections()) < 0) { return ret; } }
	return RET_OK;
}
