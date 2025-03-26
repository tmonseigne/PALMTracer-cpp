


/*******************************************************************************************************************
									TabTimesCloserNeighbors Class
									It is derived from the class TabTimes
								This class find the closer neighbors for each structure.
*******************************************************************************************************************/
#ifndef __tab_times_closer_neighbors__
		#define __tab_times_closer_neighbors__
		#
		#include "../kernel/tab_times.h"
		#include "time_nodes_closer_neighbors.h"
	
	class TabTimesCloserNeighbors:public TabTimes{
	protected:
		TimeNodes* allocTimeNodes(long time, long nbNodes);
		void allocNode();
	public:
		TabTimesCloserNeighbors(PublicParameters *publicParameters);
		~TabTimesCloserNeighbors();
		
		int runModel();
		
	};


#endif	//__tab_times_closer_neighbors__