/*******************************************************************************************************************
									TabTimes Class
							Array containing all the TimeNodes for all instants
*******************************************************************************************************************/

/********************************************************************************************************************
Class TabTimes contains several TabTimes objets for each time (nbTimes).
********************************************************************************************************************/


#include "tab_times.h"

/********************************************************************************************************************
Constructuor(int _nbTimes) of TabTimes class. 
Create nbTimes objets TimeNodes.
********************************************************************************************************************/
TabTimes::TabTimes(PublicParameters* _publicParameters)
{
	publicParameters = _publicParameters;
	publicParameters->nbNodeAllocated = 0;
	publicParameters->nbNodeDeallocated = 0;
	publicParameters->nbChainNeighborDeallocated = 0;
	publicParameters->nbChainNeighborAllocated = 0;
	timeNodes = NULL;
	nbTimes = 0;
}


/********************************************************************************************************************
Destructuor() of TabTimes class. 
Delete all objets TimeNodes.
********************************************************************************************************************/
TabTimes::~TabTimes()
{
	for (int i = 0; i < nbTimes; i++) { delete timeNodes[i]; }
	delete [] timeNodes;
}


/********************************************************************************************************************
allocTabTimes() alloc the timeNodes array of nbTimes long 
********************************************************************************************************************/
void TabTimes::allocTabTimes() { timeNodes = new TimeNodes*[nbTimes]; }


/********************************************************************************************************************
initialisation() initialisation of the data in memory
********************************************************************************************************************/
int TabTimes::initialisation()
{
	int ret;
	if (publicParameters->statFile[0] != L'\0') { if ((ret = readStatFile()) < 0) { return ret; } }
	else if ((ret = readPointsList()) < 0) { return ret; }

	if ((ret = linkNeighborsPresent()) < 0) { return ret; }
	if ((ret = linkNeighborsBeforeAndAfter()) < 0) { return ret; }
	if ((ret = sortByIntensityAllNeighbors()) < 0) { return ret; }

	if (publicParameters->nbTemporalWindows > 1) { if ((ret = linkNeighborsBlinkBeforeAndAfter()) < 0) { return ret; } }

	checkLinksHomogeneity(L"initalisation");
	srand((unsigned)time(NULL));
	return RET_OK;
}


/********************************************************************************************************************
printConnectionsASCII() print the result of the connections in a ASCII file as Daniel Choquet has defined:
Colonne 1         numero de la trace
Colonne 2         numero de l'image, ou temps (commencant à 1 SVP)
Colonne 3         X
Colonne 4         Y
(Colonne 5        Z)
Colonne 5 ou 6         nombre de pixels de l'objet
Colonne 6 ou 7         intensité intégrée
********************************************************************************************************************/
int TabTimes::printConnectionsASCII()
{
	FILE* stream;
	if ((stream = _wfopen(publicParameters->fileResASCII, L"w")) == NULL) {
		wprintf(L"The file %ls can not be opened\n", publicParameters->fileResASCII);
		return T_ERR_OPENNING_FILE;
	}
	for (int t = 0; t < nbTimes; t++) { const TimeNodes* _localTimeNodes = timeNodes[t];
		for (int n = 0; n < _localTimeNodes->nbNodes; n++) {
			Node* node = _localTimeNodes->tabNodes[n];
			if (node->connectionPast == NULL && node->trackedStructureIndex > 0) {
				for (const Node* nodeConnection = node; nodeConnection != NULL; nodeConnection = nodeConnection->connectionFuture) {
					fprintf(stream, "%d\t", node->trackedStructureIndex);
					fprintf(stream, "%d\t", nodeConnection->time + 1);
					fprintf(stream, "%f\t", nodeConnection->info->getX());
					fprintf(stream, "%f\t", nodeConnection->info->getY());
					if (publicParameters->nbDim == 3) {
						fprintf(stream, "%f\t", nodeConnection->info->getZ());
						fprintf(stream, "%d\t", nodeConnection->info->getVolume());
					}
					else { fprintf(stream, "%d\t", nodeConnection->info->getSurface()); }
					fprintf(stream, "%f\n", nodeConnection->info->getIntensity());
				}
			}
		}
	}
	fclose(stream);
	return RET_OK;
}

/********************************************************************************************************************
printConnectionsASCII_trc() print the result of the connections in a pointer as Daniel Choquet has defined:
Colonne 1         numero de la trace
Colonne 2         numero de l'image, ou temps (commencant à 1 SVP)
Colonne 3         X
Colonne 4         Y
(Colonne 5        Z)
Colonne 5 ou 6         nombre de pixels de l'objet
Colonne 6 ou 7         intensité intégrée
********************************************************************************************************************/
int TabTimes::printConnectionsASCII_trc()
{
	unsigned long index = 0;
	for (int t = 0; t < nbTimes; t++) { const TimeNodes* _localTimeNodes = timeNodes[t];
		for (int n = 0; n < _localTimeNodes->nbNodes; n++) {
			Node* node = _localTimeNodes->tabNodes[n];
			if (node->connectionPast == NULL && node->trackedStructureIndex > 0) {
				for (const Node* nodeConnection = node; nodeConnection != NULL; nodeConnection = nodeConnection->connectionFuture) {
					publicParameters->trc[index++] = node->trackedStructureIndex;
					publicParameters->trc[index++] = nodeConnection->time + 1;
					publicParameters->trc[index++] = nodeConnection->info->getX();
					publicParameters->trc[index++] = nodeConnection->info->getY();

					if (publicParameters->nbDim == 3) {
						publicParameters->trc[index++] = nodeConnection->info->getZ();
						publicParameters->trc[index++] = nodeConnection->info->getVolume();
					}
					else { publicParameters->trc[index++] = nodeConnection->info->getSurface(); }

					publicParameters->trc[index++] = nodeConnection->info->getIntensity();

					publicParameters->trc[index++] = nodeConnection->info->getZ();

					publicParameters->trc[index++] = nodeConnection->info->getColor();
					publicParameters->trc[index++] = nodeConnection->info->getMolId();

					// The number of trajectory
					publicParameters->nbTraj = node->trackedStructureIndex;
				}
			}
		}
	}
	return RET_OK;
}

/********************************************************************************************************************
readStatFile() reads the file specified in publicParameters and allocates the corresponding TimeNodes
********************************************************************************************************************/
int TabTimes::readStatFile()
{
	//unsigned long nbTimes;
	unsigned long nbNodes;
	int fileHandle;


	if ((fileHandle = _wopen(publicParameters->statFile, O_RDONLY | O_BINARY, S_IREAD)) != -1) {
		int count = read(fileHandle, &publicParameters->nbDim, sizeof(long)); //reading the number of dimentions
		if (count != sizeof(long)) { return T_ERR_CORRUPT_FILE; }
		if (publicParameters->nbDim != 2 && publicParameters->nbDim != 3) { return T_ERR_NBDIM; }


		count = read(fileHandle, &nbTimes, sizeof(long)); //reading the number of times
		if (count != sizeof(long)) { return T_ERR_CORRUPT_FILE; }


		allocTabTimes(); //Allocation of tabTimes

		long incre = 0;

		publicParameters->maxX = 0;
		publicParameters->maxY = 0;
		publicParameters->maxZ = 0;

		for (int time_ = 0; time_ < (long)nbTimes; time_++) {
			count = read(fileHandle, &nbNodes, sizeof(long)); //reading the number of nodes
			if (count != sizeof(long)) { return T_ERR_CORRUPT_FILE; }


			timeNodes[time_] = allocTimeNodes(time_, (int)nbNodes); //Allocation of timeNodes
			//new TimeNodes(time_,(int)nbNodes);
			TimeNodes* localTimeNodes = timeNodes[time_];

			if (localTimeNodes == NULL) { return T_ERR_ALLOC; }


			for (long node = 0; node < localTimeNodes->nbNodes; node++) {
				count = read(fileHandle, (*localTimeNodes)(node)->info->getPtrInfoData(), (*localTimeNodes)(node)->info->sizeofInfoData()); //reading one node

				if ((*localTimeNodes)(node)->info->getIntensity() < 0) { wprintf(L"%d : <0\n", ++incre); }

				if (count != (long)(*localTimeNodes)(node)->info->sizeofInfoData()) { return T_ERR_CORRUPT_FILE; }

				if (publicParameters->maxX < (*localTimeNodes)(node)->info->getX()) { publicParameters->maxX = (*localTimeNodes)(node)->info->getX(); }
				if (publicParameters->maxY < (*localTimeNodes)(node)->info->getY()) { publicParameters->maxY = (*localTimeNodes)(node)->info->getY(); }
				if (publicParameters->nbDim == 3) { if (publicParameters->maxZ < (*localTimeNodes)(node)->info->getZ()) { publicParameters->maxZ = (*localTimeNodes)(node)->info->getZ(); } }
			}

			//Sorting the nodes by the intensity for one time
			int indexMax;
			for (int i = 0; i < timeNodes[time_]->nbNodes - 1; i++) {
				double maxIntensity = -1.;
				for (int j = i; j < localTimeNodes->nbNodes; j++) {
					if (maxIntensity < (*localTimeNodes)(j)->info->getIntensity()) {
						maxIntensity = (*localTimeNodes)(j)->info->getIntensity();
						indexMax = j;
					}
				}
				Node* nodeMax = localTimeNodes->tabNodes[indexMax];
				localTimeNodes->tabNodes[indexMax] = localTimeNodes->tabNodes[i];
				localTimeNodes->tabNodes[i] = nodeMax;
			}
		}
		close(fileHandle);
	}
	else { return T_ERR_OPENNING_FILE; }

	return RET_OK;
}

/********************************************************************************************************************
readPointsList() reads the 'pointsList' specified in publicParameters and allocates the corresponding TimeNodes
********************************************************************************************************************/
int TabTimes::readPointsList()
{
	//unsigned long nbTimes;
	unsigned long index = 0;
	int i;

	nbTimes = publicParameters->nbTimes;

	allocTabTimes(); //Allocation of tabTimes

	publicParameters->maxX = 0;
	publicParameters->maxY = 0;
	publicParameters->maxZ = 0;

	for (int time_ = 0; time_ < (long)nbTimes; time_++) {
		// Computing the number of nodes by frame
		unsigned long nbNodes = 0;
		for (i = index; ; /*i+=4*/ i += publicParameters->NbSegParams) {
			if ((int)publicParameters->pointsList[i + 9] != -1 && (int)publicParameters->pointsList[i + 9] != 0) { nbNodes++; } // surface
			else { break; }
		}


		timeNodes[time_] = allocTimeNodes(time_, (int)nbNodes); //Allocation of timeNodes
		TimeNodes* localTimeNodes = timeNodes[time_];

		if (localTimeNodes == NULL) { return T_ERR_ALLOC; }


		for (long node = 0; node < localTimeNodes->nbNodes; node++) {
			// Just for the case of '2D'
			InfoNode* str = (*localTimeNodes)(node)->info;
			// version correcte
			str->setX(publicParameters->pointsList[index + node * publicParameters->NbSegParams + 3]);
			str->setY(publicParameters->pointsList[index + node * publicParameters->NbSegParams + 4]);
			str->setZ(publicParameters->pointsList[index + node * publicParameters->NbSegParams + 10]);

			str->setColor(publicParameters->pointsList[index + node * publicParameters->NbSegParams + 11]);
			str->setMolId(publicParameters->pointsList[index + node * publicParameters->NbSegParams + 12]);

			str->setIndex(node + 1);

			str->setIntensity(publicParameters->pointsList[index + node * publicParameters->NbSegParams + 8]);

			str->setTime(time_ + 1);
			str->setSurface((unsigned long)publicParameters->pointsList[index + node * publicParameters->NbSegParams + 9]);
			str->setMean(publicParameters->pointsList[index + node * publicParameters->NbSegParams + 8] / (unsigned long)publicParameters->pointsList[index + node * publicParameters->NbSegParams + 9]);

			if (publicParameters->maxX < (*localTimeNodes)(node)->info->getX()) { publicParameters->maxX = (*localTimeNodes)(node)->info->getX(); }
			if (publicParameters->maxY < (*localTimeNodes)(node)->info->getY()) { publicParameters->maxY = (*localTimeNodes)(node)->info->getY(); }
			if (publicParameters->nbDim == 3) { if (publicParameters->maxZ < (*localTimeNodes)(node)->info->getZ()) { publicParameters->maxZ = (*localTimeNodes)(node)->info->getZ(); } }
		}

		//Sorting the nodes by the intensity for one time
		int indexMax;
		for (i = 0; i < timeNodes[time_]->nbNodes - 1; i++) {
			double maxIntensity = -1.;
			for (int j = i; j < localTimeNodes->nbNodes; j++) {
				if (maxIntensity < (*localTimeNodes)(j)->info->getIntensity()) {
					maxIntensity = (*localTimeNodes)(j)->info->getIntensity();
					indexMax = j;
				}
			}
			Node* nodeMax = localTimeNodes->tabNodes[indexMax];
			localTimeNodes->tabNodes[indexMax] = localTimeNodes->tabNodes[i];
			localTimeNodes->tabNodes[i] = nodeMax;
		}

		index += (nbNodes + 1) * publicParameters->NbSegParams; //(A vérifier)
	}
	return RET_OK;
}

/********************************************************************************************************************
operator()(int i) to access to all TimeNodes objets. 
********************************************************************************************************************/
TimeNodes* TabTimes::operator()(int i) { return timeNodes[i]; }

/********************************************************************************************************************
operator[](int i) to access to pointors to all TimeNodes objets. 
********************************************************************************************************************/
TimeNodes** TabTimes::operator[](int i) { return timeNodes + i; }

/********************************************************************************************************************
int linkNeighborsPresent() to link of the neighbors in the presents of all nodes of all times
********************************************************************************************************************/
int TabTimes::linkNeighborsPresent()
{
	int ret;
	for (long time_ = 0; time_ < nbTimes; time_++) { if ((ret = timeNodes[time_]->linkNeighborsPresent()) < 0) { return ret; } }
	return RET_OK;
}


/********************************************************************************************************************
int linkNeighborsBeforeAndAfter() to link of the neighbors in the past and in the future of all nodes of all times
********************************************************************************************************************/
int TabTimes::linkNeighborsBeforeAndAfter()
{
	int ret;
	for (long time_ = 0; time_ < nbTimes - 1; time_++) { if ((ret = timeNodes[time_]->linkNeighborsBeforeAndAfter(*timeNodes[time_ + 1])) < 0) { return ret; } }
	return RET_OK;
}

/********************************************************************************************************************
int linkNeighborsBlinkBeforeAndAfter() to link of the neighbors in the past and in the future of all nodes of 
all times for blinking purpose
********************************************************************************************************************/
int TabTimes::linkNeighborsBlinkBeforeAndAfter()
{
	int ret;
	for (long time_ = 0; time_ < nbTimes; time_++) { for (long dt = 2; dt <= publicParameters->nbTemporalWindows; dt++) { if (time_ + dt < nbTimes) { if ((ret = timeNodes[time_]->linkNeighborsBlinkBeforeAndAfter(*timeNodes[time_ + dt])) < 0) { return ret; } } } }
	return RET_OK;
}

/********************************************************************************************************************
int sortByIntensityAllNeighbors() sort by intensity of the node neighbors of each time 
********************************************************************************************************************/
int TabTimes::sortByIntensityAllNeighbors()
{
	int ret;
	for (long time_ = 0; time_ < nbTimes; time_++) { if ((ret = timeNodes[time_]->sortByIntensityAllNeighbors()) < 0) { return ret; } }
	return RET_OK;
}


/********************************************************************************************************************
void setColorToTrackedStructures() set a color index to each tracked structures
********************************************************************************************************************/
void TabTimes::setColorToTrackedStructures()
{
	Node* nodePastFission;
	Node* nodePastFission2;
	int newIndex = 1;
	for (long time_ = 0; time_ < nbTimes; time_++) {
		for (long indexNode = 0; indexNode < timeNodes[time_]->nbNodes; indexNode++) {
			Node* node = timeNodes[time_]->tabNodes[indexNode];
			if (node->trackedStructureIndex != 0) {
				if (node->connectionPast == NULL) {
					if (node->connectionBlinkPast == NULL) {
						if (node->connectionFusion1Past != NULL) {
							if (node->likeness(node->connectionFusion1Past) > node->likeness(node->connectionFusion2Past)) {
								node->colorIndex[0] = node->connectionFusion1Past->colorIndex[0];
								node->colorIndex[1] = node->connectionFusion2Past->colorIndex[0];
							}
							else {
								node->colorIndex[0] = node->connectionFusion2Past->colorIndex[0];
								node->colorIndex[1] = node->connectionFusion1Past->colorIndex[0];
							}
						}
						else { //node->connectionFusion1Past==NULL
							if (node->connectionFissionPast != NULL) {
								Node* node2 = node == node->connectionFissionPast->connectionFission1Future ? node->connectionFissionPast->connectionFission2Future : node->connectionFissionPast->connectionFission1Future;
								for (nodePastFission = node->connectionFissionPast; nodePastFission != NULL; nodePastFission = nodePastFission->connectionPast) { nodePastFission2 = nodePastFission; }
								if (nodePastFission2->connectionFusion1Past != NULL) {
									Node* nodePastFission3 = nodePastFission2->connectionFusion1Past;
									Node* nodePastFission4 = nodePastFission2->connectionFusion2Past;
									if (node->likeness(nodePastFission3) + node2->likeness(nodePastFission4) > node2->likeness(nodePastFission3) + node->likeness(nodePastFission4)) {
										node->colorIndex[0] = nodePastFission3->colorIndex[0];
										node->colorIndex[1] = nodePastFission3->colorIndex[0];
										node2->colorIndex[0] = nodePastFission4->colorIndex[0];
										node2->colorIndex[1] = nodePastFission4->colorIndex[0];
									}
									else {
										node2->colorIndex[0] = nodePastFission3->colorIndex[0];
										node2->colorIndex[1] = nodePastFission3->colorIndex[0];
										node->colorIndex[0] = nodePastFission4->colorIndex[0];
										node->colorIndex[1] = nodePastFission4->colorIndex[0];
									}
								}
								else {
									if (node->connectionFissionPast->likeness(node) > node->connectionFissionPast->likeness(node2)) {
										node->colorIndex[0] = node->connectionFissionPast->colorIndex[0];
										node->colorIndex[1] = node->connectionFissionPast->colorIndex[0];
									}
									else {
										node->colorIndex[0] = node->connectionFissionPast->colorIndex[1];
										node->colorIndex[1] = node->connectionFissionPast->colorIndex[1];
									}
								}
							}
							else { //node->connectionFissionPast==NULL  and node->connectionFusion1Past==NULL
								if (node->connectionFission1Future != NULL) {
									node->colorIndex[0] = newIndex;
									newIndex++;
									node->colorIndex[1] = newIndex;
									newIndex++;
								}
								else {
									node->colorIndex[0] = newIndex;
									node->colorIndex[1] = newIndex;
									newIndex++;
								}
							}
						}
					}
					else {
						node->colorIndex[0] = node->connectionBlinkPast->colorIndex[0];
						node->colorIndex[1] = node->connectionBlinkPast->colorIndex[0];
					}
				}
				else { //node->connectionPast!=NULL
					node->colorIndex[0] = node->connectionPast->colorIndex[0];
					node->colorIndex[1] = node->connectionPast->colorIndex[0];
				}
				if (node->connectionFission1Future != NULL) {
					for (nodePastFission = node; nodePastFission != NULL; nodePastFission = nodePastFission->connectionPast) { nodePastFission2 = nodePastFission; }
					if (nodePastFission2->connectionFusion1Past != NULL) {
						node->colorIndex[0] = nodePastFission2->connectionFusion1Past->colorIndex[0];
						node->colorIndex[1] = nodePastFission2->connectionFusion2Past->colorIndex[0];
					}
					else {
						node->colorIndex[1] = newIndex;
						newIndex++;
					}
				}
			}
			else { //node->trackedStructureIndex==0
				node->colorIndex[0] = 0;
				node->colorIndex[1] = 0;
			}
		}
	}
}

void TabTimes::setIndexToTrackedStructures()
{
	unsigned long newIndex = 1;
	long nbSuccessiveNodes;
	Node* node;
	Node* oldNode;
	for (long time_ = 0; time_ < nbTimes; time_++) {
		for (long indexNode = 0; indexNode < timeNodes[time_]->nbNodes; indexNode++) {
			if (timeNodes[time_]->tabNodes[indexNode]->connectionPast == NULL) {
				for (node = timeNodes[time_]->tabNodes[indexNode], oldNode = NULL, nbSuccessiveNodes = 0; node != NULL && nbSuccessiveNodes < publicParameters->nbMinSuccessiveStructures; node = node->connectionFuture, nbSuccessiveNodes++) {
					oldNode = node;
					node->trackedStructureIndex = 0;
				}
				bool isNewIndex = false;
				if (nbSuccessiveNodes >= publicParameters->nbMinSuccessiveStructures) //it was == 
				{
					isNewIndex = true;
				}
				if (timeNodes[time_]->tabNodes[indexNode]->connectionFissionPast != NULL) { isNewIndex = true; }
				if (timeNodes[time_]->tabNodes[indexNode]->connectionFusion1Past != NULL) { isNewIndex = true; }
				if (timeNodes[time_]->tabNodes[indexNode]->connectionBlinkPast != NULL) { isNewIndex = true; }
				if (oldNode != NULL) {
					if (oldNode->connectionFission1Future != NULL) { isNewIndex = true; }
					if (oldNode->connectionFusionFuture != NULL) { isNewIndex = true; }
					if (oldNode->connectionBlinkFuture != NULL) { isNewIndex = true; }
				}

				if (isNewIndex) {
					timeNodes[time_]->tabNodes[indexNode]->trackedStructureIndex = newIndex;
					newIndex++;
				}
			}
			else { timeNodes[time_]->tabNodes[indexNode]->trackedStructureIndex = timeNodes[time_]->tabNodes[indexNode]->connectionPast->trackedStructureIndex; }
		}
	}
}

/********************************************************************************************************************
void drawSegment(UStrArray &image, long time, Node* node1, Node*node2)  In the 2D image "image" at the time "time", 
draw a segment between the nodes "node1" and "node2"
********************************************************************************************************************/
void TabTimes::drawSegment(USrtArray& image, long time_, Node* node1, Node* node2)
{
	const double dist = node1->distance(node2);
	if (dist > 1e-8) {
		for (double t = 0; t < dist; t++) { const long x = (long)(node2->info->getX() * (dist - t) / dist + node1->info->getX() * t / dist);
			const long y = (long)(node2->info->getY() * (dist - t) / dist + node1->info->getY() * t / dist);
			image.data[time_ * image.size(1) * image.size(2) + image.size(1) * y + x] = (unsigned short)node1->colorIndex[0];
		}
	}
}


/********************************************************************************************************************
int giveResult() just call the function to display the result
********************************************************************************************************************/
int TabTimes::giveResult()
{
	int ret;

	setIndexToTrackedStructures();
	setColorToTrackedStructures();

	if (publicParameters->fileResASCII[0] != L'\0') {
		if ((ret = printConnectionsASCII()) < 0) { return ret; }
		if ((ret = printConnectionsASCII_trc()) < 0) { return ret; } //(A enlever)
	}
	else if ((ret = printConnectionsASCII_trc()) < 0) { return ret; }

	return RET_OK;
}

/********************************************************************************************************************
void checkLinksHomogeneity() just call the function to display the result
********************************************************************************************************************/
void TabTimes::checkLinksHomogeneity(const wchar_t* st) { for (long time_ = 0; time_ < nbTimes; time_++) { timeNodes[time_]->checkLinksHomogeneity(st); } }
