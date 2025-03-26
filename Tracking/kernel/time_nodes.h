/*******************************************************************************************************************
									TimeNodes Class
							Array containing all the Node s of a particular instant
*******************************************************************************************************************/


#ifndef __time_nodes__
#define __time_nodes__


#include "general.h"
#include "nodes.h"


class TimeNodes
{
protected:
	int time;
	virtual Node* allocNode(void) =0;

public:
	Node** tabNodes;
	TimeNodes(int _time, int _nbNodes, PublicParameters* _publicParameters);
	~TimeNodes();

	PublicParameters* publicParameters;
	Node* operator()(int i);
	int nbNodes;
	int linkNeighborsPresent();
	int linkNeighborsBeforeAndAfter(TimeNodes& futureTimeNodes);
	int linkNeighborsBlinkBeforeAndAfter(TimeNodes& futureTimeNodes);
	int sortByIntensityAllNeighbors();

	void checkLinksHomogeneity(wchar_t* st);
};


#endif //__time_nodes__
