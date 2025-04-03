#include "Segmentation2D.h"
#include "tools/NoError.h"

#define MAX_ELEMENTS 10000000
#define EPSILON 0.0000001

#define SEG(i,j)	seg.data[(i)*sizeDataY+(j)]
#define SEGSTRUCT(i,j)	segStructure.data[(i)*sizeDataY+(j)]
#define SEGHoles(i,j)	segHoles.data[(i)*sizeDataY+(j)]
#define DATA(i,j)	data.data[(i)*sizeDataY+(j)]
#define TEST_ALLOC(ptr) if((ptr)==NULL){printf("Error of allocation at line %d, file %s\n", __LINE__, __FILE__);return(ALLOC_ERR);}

inline long min(const long a, const long b) { return a <= b ? a : b; }
inline long min(const long a, const long b, const long c, const long d) { return min(min(a, b), min(c, d)); }
inline long max(const long a, const long b) { return a >= b ? a : b; }
inline double min(const double a, const double b) { return a <= b ? a : b; }
inline double max(const double a, const double b) { return a >= b ? a : b; }

void** Segmentation2D::getTabInfo() const { return (void**)tabInfo; }

void Segmentation2D::desalloc()
{
	InfoStruct2D* ptr = headInfo;
	while (ptr != NULL) {
		headInfo = ptr;
		ptr = ptr->next_;
		delete headInfo;
	}
	if (tabInfo != NULL) { delete [] tabInfo; }
	tabInfo = NULL;
	nbInfoStruct = 0;
	headInfo = NULL;
}

Segmentation2D::Segmentation2D(const unsigned long _nbInfoStruct)
{
	headInfo = NULL;
	nbInfoStruct = _nbInfoStruct;
	tabInfo = NULL;
}

Segmentation2D::Segmentation2D()
{
	headInfo = NULL;
	nbInfoStruct = 0;
	tabInfo = NULL;
}

Segmentation2D::~Segmentation2D() { desalloc(); }


int Segmentation2D::alloc(const unsigned long nbStruct)
{
	desalloc();

	tabInfo = new InfoStruct2D*[nbStruct + 1];
	if (tabInfo == NULL) { return S_ALLOC_ERR; }


	if (nbStruct > 0) {
		headInfo = new InfoStruct2D(1);

		tabInfo[0] = NULL;
		tabInfo[1] = headInfo;
		InfoStruct2D* ptr = headInfo;
		for (unsigned long ind = 2; ind <= nbStruct; ind++) {
			ptr->next_ = new InfoStruct2D(ind);
			tabInfo[ind] = ptr->next_;
			if (ptr->next_ == NULL) { return S_ALLOC_ERR; }
			ptr = ptr->next_;
		}
	}

	nbInfoStruct = nbStruct;
	return RET_OK;
}


unsigned long Segmentation2D::filterData(USrtArray& seg, double* param, unsigned short* useParam)
{
	const unsigned short sizeDataY = seg.size(1);

	InfoStruct2D* ptr = headInfo;
	InfoStruct2D* ptrMax;
	InfoStruct2D* ptrOld = NULL;
	unsigned long ind = 0;
	long i;
	bool unselected;
	float valMax;
	const long nbInfoStructOld = nbInfoStruct;
	long oldIndex;

	char* tabBound = NULL;
	if (useParam[12] == 1) {
		tabBound = new char[nbInfoStruct + 1];
		for (i = 0; i < (long)(nbInfoStruct + 1); i++) { tabBound[i] = 0; }
		for (i = 0; i < seg.size(0); i++) {
			if (seg.data[i * sizeDataY] != 0) { tabBound[seg.data[i * sizeDataY]] = 1; }
			if (seg.data[i * sizeDataY + seg.size(1) - 1] != 0) { tabBound[seg.data[i * sizeDataY + seg.size(1) - 1]] = 1; }
		}
		for (i = 1; i < seg.size(1) - 1; i++) {
			if (seg.data[i] != 0) { tabBound[seg.data[i]] = 1; }
			if (seg.data[(seg.size(0) - 1) * sizeDataY + i] != 0) { tabBound[seg.data[(seg.size(0) - 1) * sizeDataY + i]] = 1; }
		}
	}

	while (ptr != NULL) {
		ind++;
		unselected = false;
		if (useParam[0] == 1 && ptr->data_.surf <= param[0]) { unselected = true; }                       //minSurf
		else if (useParam[1] == 1 && ptr->data_.surf >= param[1]) { unselected = true; }                  //maxSurf
		else if (useParam[2] == 1 && ptr->data_.mean * ptr->data_.surf <= param[2]) { unselected = true; } //minInten
		else if (useParam[3] == 1 && ptr->data_.mean * ptr->data_.surf >= param[3]) { unselected = true; } //maxInten
		else if (useParam[4] == 1 && ptr->data_.max.val <= param[4]) { unselected = true; }               //minMaxi
		else if (useParam[5] == 1 && ptr->data_.max.val >= param[5]) { unselected = true; }               //maxMaxi
		else if (useParam[7] == 1) {
			if (min(ptr->data_.sigma.x, ptr->data_.sigma.y) / max(ptr->data_.sigma.x, ptr->data_.sigma.y) >= param[7]) { unselected = true; }
		} //minMorphParam
		else if (useParam[8] == 1) {
			if (min(ptr->data_.sigma.x, ptr->data_.sigma.y) / max(ptr->data_.sigma.x, ptr->data_.sigma.y) <= param[8]) { unselected = true; }
		}                                                                              //maxMorphParam
		else if (useParam[12] == 1) { if (tabBound[ind] == 1) { unselected = true; } } //no struct on boundaries
		if (unselected) {
			if (ptrOld == NULL) {
				headInfo = ptr->next_;
				delete ptr;
				ptr = headInfo;
			}
			else {
				ptrOld->next_ = ptr->next_;
				delete ptr;
				ptr = ptrOld->next_;
			}
			tabInfo[ind] = NULL;
			nbInfoStruct--;
		}
		else {
			ptrOld = ptr;
			ptr = ptr->next_;
		}
	}
	if (useParam[12] == 1) { delete [] tabBound; }


	if (useParam[6] == 1) {
		ptr = headInfo;
		ind = 0;
		ptrOld = NULL;
		while (ptr != NULL) {
			ind++;
			unselected = false;
			for (i = 1; i <= nbInfoStructOld && !unselected; i++) {
				if (tabInfo[i] != NULL) {
					const double d = sqrt(
						(double)((ptr->data_.max.x - tabInfo[i]->data_.max.x) * (ptr->data_.max.x - tabInfo[i]->data_.max.x) + (
							         ptr->data_.max.y - tabInfo[i]->data_.max.y) * (ptr->data_.max.y - tabInfo[i]->data_.max.y)));
					if (d > 1 && d < param[6] && tabInfo[i]->data_.max.val > ptr->data_.max.val) { unselected = true; }
				}
			}

			if (unselected) {
				if (ptrOld == NULL) {
					headInfo = ptr->next_;
					oldIndex = ptr->data_.index;
					delete ptr;
					ptr = headInfo;
				}
				else {
					ptrOld->next_ = ptr->next_;
					oldIndex = ptr->data_.index;
					delete ptr;
					ptr = ptrOld->next_;
				}
				tabInfo[oldIndex] = NULL;
				nbInfoStruct--;
			}
			else {
				ptrOld = ptr;
				ptr = ptr->next_;
			}
		}
	}


	for (long indexFilter = 9; indexFilter <= 11; indexFilter++) {
		if (useParam[indexFilter] == 1 && nbInfoStruct > param[indexFilter] && (long)param[indexFilter] > 0) {
			//N first structure with the biggest surface (indexFilter=9), the biggest surface (indexFilter=10) or 
			//the biggest total intensity (indexFilter=11)

			ptrOld = NULL;
			for (i = 0; i < (long)param[indexFilter]; i++) {
				if (ptrOld == NULL) { ptr = headInfo; }
				else { ptr = ptrOld->next_; }

				switch (indexFilter) {
					case 9: valMax = (float)ptr->data_.surf;
						break;
					case 10: valMax = (float)(ptr->data_.surf * ptr->data_.mean);
						break;
					case 11:
					default: valMax = (float)ptr->data_.max.val;
						break;
				}

				ptrMax = NULL;

				while (ptr->next_ != NULL) {
					switch (indexFilter) {
						case 9: if (valMax <= ptr->next_->data_.surf) {
								valMax = (float)ptr->next_->data_.surf;
								ptrMax = ptr;
							}
							break;
						case 10: if (valMax <= (long)ptr->next_->data_.surf * ptr->next_->data_.mean) {
								valMax = (float)((long)ptr->next_->data_.surf * ptr->next_->data_.mean);
								ptrMax = ptr;
							}
							break;
						case 11: if (valMax <= (float)ptr->next_->data_.max.val) {
								valMax = (float)ptr->next_->data_.max.val;
								ptrMax = ptr;
							}
							break;
						default: ;
					}
					ptr = ptr->next_;
				}
				if (ptrMax != NULL) {
					if (ptrOld != NULL) {
						ptr = ptrOld->next_;
						ptrOld->next_ = ptrMax->next_;
						ptrMax->next_ = ptrMax->next_->next_;
						ptrOld->next_->next_ = ptr;
						ptrOld = ptrOld->next_;
					}
					else {
						ptr = headInfo;
						headInfo = ptrMax->next_;
						ptrMax->next_ = headInfo->next_;
						headInfo->next_ = ptr;
						ptrOld = headInfo;
					}
				}
				else if (ptrOld == NULL) { ptrOld = headInfo; }
				else { ptrOld = ptrOld->next_; }
			}
			ptr = ptrOld->next_;
			ptrOld->next_ = NULL;
			while (ptr != NULL) {
				ptrMax = ptr->next_;
				tabInfo[ptr->data_.index] = NULL;
				delete ptr;
				nbInfoStruct--;
				ptr = ptrMax;
			}
		}
	}


	ptr = headInfo;
	for (ind = 1; ind <= nbInfoStruct; ind++) {
		ptr->data_.index = ind;
		ptr = ptr->next_;
	}


	for (i = 0; i < (long)seg.nbPoints; i++) {
		if (seg.data[i] != 0) {
			if (tabInfo[seg.data[i]] == NULL) { seg.data[i] = 0; }
			else { seg.data[i] = (unsigned short)tabInfo[seg.data[i]]->data_.index; }
		}
	}


	delete [] tabInfo;
	tabInfo = new InfoStruct2D*[nbInfoStruct + 1];
	tabInfo[0] = NULL;
	tabInfo[1] = headInfo;
	//tabInfo[1]->data.index=1;
	ptr = headInfo;
	for (ind = 2; ind <= nbInfoStruct; ind++) {
		tabInfo[ind] = ptr->next_;
		//tabInfo[ind]->data.index=ind;
		ptr = ptr->next_;
	}

	return nbInfoStruct;
}


float Segmentation2D::tic()
{
	SYSTEMTIME st_time;
	GetSystemTime(&st_time);
	return (float)(st_time.wMilliseconds + 1000 * (st_time.wSecond + 60 * (st_time.wMinute + 60 * st_time.wHour)));
}

int Segmentation2D::computeStatistics(FltArray& data, USrtArray& seg, float* waveData, const unsigned short gaussFit, const bool doRealoc)
{
	double val;
	unsigned long ind;
	unsigned long cptPoints = 0;
	int ret;

	const unsigned short sizeDataX = data.size(0);
	const unsigned short sizeDataY = data.size(1);


	if (doRealoc) { if ((ret = alloc(nbInfoStruct)) < 0) { return ret; } }


	//initalisation of array

	for (ind = 1; ind <= nbInfoStruct; ind++) {
		tabInfo[ind]->data_.peri = 0;
		tabInfo[ind]->data_.surf = 0;
		tabInfo[ind]->data_.morpho = 0.;
		tabInfo[ind]->data_.max.x = 0;
		tabInfo[ind]->data_.max.y = 0;
		tabInfo[ind]->data_.max.val = 0.;
		tabInfo[ind]->data_.min.x = 0;
		tabInfo[ind]->data_.min.y = 0;
		tabInfo[ind]->data_.min.val = 0;
		tabInfo[ind]->data_.stdev = 0.;
		tabInfo[ind]->data_.moment.x = 0.;
		tabInfo[ind]->data_.moment.y = 0.;
		tabInfo[ind]->data_.momentSurf.x = 0.;
		tabInfo[ind]->data_.momentSurf.y = 0.;
		tabInfo[ind]->data_.mean = 0.;
		tabInfo[ind]->data_.intensity = 0.;
		tabInfo[ind]->data_.var.xx = 0.;
		tabInfo[ind]->data_.var.yy = 0.;
		tabInfo[ind]->data_.var.xy = 0.;
		tabInfo[ind]->data_.sigma.x = 0.;
		tabInfo[ind]->data_.sigma.y = 0.;
		tabInfo[ind]->data_.angle = 0.;
		tabInfo[ind]->data_.boundingBox.xMax = 0;
		tabInfo[ind]->data_.boundingBox.yMax = 0;
		tabInfo[ind]->data_.boundingBox.xMin = sizeDataX;
		tabInfo[ind]->data_.boundingBox.yMin = sizeDataY;
		tabInfo[ind]->data_.nbSubstructures = 0;
	}

	for (unsigned short i = 0; i < sizeDataX; i++) {
		for (unsigned short j = 0; j < sizeDataY; j++) {
			if (seg.data[i * sizeDataY + j] > 0) {
				cptPoints++;
				val = data.data[i * sizeDataY + j];
				ind = seg.data[i * sizeDataY + j];


				tabInfo[ind]->data_.surf++;
				tabInfo[ind]->data_.moment.x += val * (double)i;
				tabInfo[ind]->data_.moment.y += val * (double)j;
				tabInfo[ind]->data_.momentSurf.x += i;
				tabInfo[ind]->data_.momentSurf.y += j;
				tabInfo[ind]->data_.mean += val;
				tabInfo[ind]->data_.var.xx += val * i * i;
				tabInfo[ind]->data_.var.yy += val * j * j;
				tabInfo[ind]->data_.var.xy += val * i * j;
				tabInfo[ind]->data_.stdev += val * val;

				if (tabInfo[ind]->data_.boundingBox.xMax < i) { tabInfo[ind]->data_.boundingBox.xMax = i; }
				if (tabInfo[ind]->data_.boundingBox.yMax < j) { tabInfo[ind]->data_.boundingBox.yMax = j; }
				if (tabInfo[ind]->data_.boundingBox.xMin > i) { tabInfo[ind]->data_.boundingBox.xMin = i; }
				if (tabInfo[ind]->data_.boundingBox.yMin > j) { tabInfo[ind]->data_.boundingBox.yMin = j; }


				if (val > tabInfo[ind]->data_.max.val || tabInfo[ind]->data_.surf == 1) {
					tabInfo[ind]->data_.max.val = val;
					tabInfo[ind]->data_.max.x = i;
					tabInfo[ind]->data_.max.y = j;
				}

				if (val < tabInfo[ind]->data_.min.val || tabInfo[ind]->data_.surf == 1) {
					tabInfo[ind]->data_.min.val = val;
					tabInfo[ind]->data_.min.x = i;
					tabInfo[ind]->data_.min.y = j;
				}

				if (i == 0) { tabInfo[ind]->data_.peri++; }
				if (i == sizeDataX - 1) { tabInfo[ind]->data_.peri++; }
				if (j == 0) { tabInfo[ind]->data_.peri++; }
				if (j == sizeDataY - 1) { tabInfo[ind]->data_.peri++; }
				if (i != 0 && i != sizeDataX - 1 && j != 0 && j != sizeDataY - 1) {
					if (seg.data[(i - 1) * sizeDataY + j] != seg.data[i * sizeDataY + j]) { tabInfo[ind]->data_.peri++; }
					if (seg.data[(i + 1) * sizeDataY + j] != seg.data[i * sizeDataY + j]) { tabInfo[ind]->data_.peri++; }
					if (seg.data[i * sizeDataY + j - 1] != seg.data[i * sizeDataY + j]) { tabInfo[ind]->data_.peri++; }
					if (seg.data[i * sizeDataY + j + 1] != seg.data[i * sizeDataY + j]) { tabInfo[ind]->data_.peri++; }
				}
			}
		}
	}


	for (ind = 1; ind <= nbInfoStruct; ind++) {
		if (tabInfo[ind]->data_.surf > 0) {
			if (gaussFit == 5) {
				/*************************************/
				tabInfo[ind]->data_.moment.x = 0.;
				tabInfo[ind]->data_.moment.y = 0.;
				double mean = 0.;

				for (int ii = -2; ii <= 2; ii++) // sizeDataX
				{
					for (int jj = -2; jj <= 2; jj++) // sizeDataY
					{
						const int refx = (int)tabInfo[ind]->data_.max.x;
						const int refy = (int)tabInfo[ind]->data_.max.y;

						if (ii + refx >= 0 && ii + refx < sizeDataX && jj + refy < sizeDataY && jj + refy >= 0) {
							val = (double)waveData[(ii + refx) * sizeDataY + (jj + refy)];
							tabInfo[ind]->data_.moment.x += val * (double)(ii + refx);
							tabInfo[ind]->data_.moment.y += val * (double)(jj + refy);
							mean += val;
						}
					}
				}

				tabInfo[ind]->data_.moment.x /= mean;
				tabInfo[ind]->data_.moment.y /= mean;
				tabInfo[ind]->data_.mean /= tabInfo[ind]->data_.surf;
				tabInfo[ind]->data_.intensity = tabInfo[ind]->data_.mean;

				/***********************************/
			}
			else {
				tabInfo[ind]->data_.moment.x /= tabInfo[ind]->data_.mean;
				tabInfo[ind]->data_.moment.y /= tabInfo[ind]->data_.mean;

				tabInfo[ind]->data_.momentSurf.x /= tabInfo[ind]->data_.surf;
				tabInfo[ind]->data_.momentSurf.y /= tabInfo[ind]->data_.surf;
				tabInfo[ind]->data_.var.xx /= tabInfo[ind]->data_.mean;
				tabInfo[ind]->data_.var.yy /= tabInfo[ind]->data_.mean;
				tabInfo[ind]->data_.var.xy /= tabInfo[ind]->data_.mean;
				tabInfo[ind]->data_.var.xx -= tabInfo[ind]->data_.moment.x * tabInfo[ind]->data_.moment.x;
				tabInfo[ind]->data_.var.yy -= tabInfo[ind]->data_.moment.y * tabInfo[ind]->data_.moment.y;
				tabInfo[ind]->data_.var.xy -= tabInfo[ind]->data_.moment.x * tabInfo[ind]->data_.moment.y;
				tabInfo[ind]->data_.intensity = tabInfo[ind]->data_.mean;
				tabInfo[ind]->data_.stdev -= tabInfo[ind]->data_.mean;
				tabInfo[ind]->data_.mean /= tabInfo[ind]->data_.surf;
				tabInfo[ind]->data_.stdev /= tabInfo[ind]->data_.surf;
				tabInfo[ind]->data_.stdev = sqrt(tabInfo[ind]->data_.stdev);
			}
		}
	}


	return cptPoints;
}


void Segmentation2D::mixIndex(USrtArray& seg) const
{
	USrtArray* mixer = new USrtArray(nbInfoStruct + 1);
	(*mixer)(0) = 0;
	srand(111);

	long i;
	for (i = 1; i <= (long)nbInfoStruct; i++) { (*mixer)(i) = (unsigned short)i; }
	for (i = 1; i <= (long)nbInfoStruct; i++) {
		const long ind = rand() % nbInfoStruct + 1;
		const long val = (*mixer)(ind);
		(*mixer)(ind) = (*mixer)(i);
		(*mixer)(i) = (unsigned short)val;
	}


	for (i = 0; i < (long)seg.nbPoints; i++) { seg[i] = (*mixer)(seg.data[i]); }
	delete mixer;
}

/********************************************************************************************************************
int Segmentation2D::filterStructureWS(...) filter the structure of the watershed.
FltArray &potential: input image
USrtArray &localMax: local maximum
float thresholdRatio: value to distinguish if 2 strucutres in contact have to be considered as a single structure (
							oversegmentation)
********************************************************************************************************************/
int Segmentation2D::filterStructureWS(FltArray& potential, USrtArray& seg, const float thresholdRatio)
{
	wsFusion2D* ptrStr = (wsFusion2D*)malloc((nbInfoStruct + 1) * sizeof(wsFusion2D));
	long i, j;

	const long sizeDataY = potential.size(1);
	const long sizeDataX = potential.size(0);
	long label;
	bool change = true;
	while (change) {
		const long oldNbInfoStruct = nbInfoStruct;

		for (i = 1; i <= (long)nbInfoStruct; i++) {
			ptrStr[i].maxPotentialWithNeighbor = -1e10;
			ptrStr[i].minPotentialOnTheStructure = +1e10;
			ptrStr[i].maxPotentialOnTheStructure = -1e10;
			ptrStr[i].structureInContact = 0;
			ptrStr[i].kept = true;
			ptrStr[i].locked = false;
			ptrStr[i].structureToBecaume = i;
			ptrStr[i].newIndex = 0;
		}


		for (i = 0; i < sizeDataX; i++) {
			for (j = 0; j < sizeDataY; j++) {
				if (seg.data[i * sizeDataY + j] > 0) {
					if (ptrStr[seg.data[i * sizeDataY + j]].minPotentialOnTheStructure > potential.data[i * sizeDataY + j]) {
						ptrStr[seg.data[i * sizeDataY + j]].minPotentialOnTheStructure = potential.data[i * sizeDataY + j];
					}
					if (ptrStr[seg.data[i * sizeDataY + j]].maxPotentialOnTheStructure < potential.data[i * sizeDataY + j]) {
						ptrStr[seg.data[i * sizeDataY + j]].x = i;
						ptrStr[seg.data[i * sizeDataY + j]].y = j;
						ptrStr[seg.data[i * sizeDataY + j]].maxPotentialOnTheStructure = potential.data[i * sizeDataY + j];
					}
					for (long ii = __max(i - 1, 0); ii <= __min(i + 1, sizeDataX - 1); ii++) {
						for (long jj = __max(j - 1, 0); jj <= __min(j + 1, sizeDataY - 1); jj++) {
							if (ii != i || jj != j) {
								if (seg.data[ii * sizeDataY + jj] != 0) {
									if (seg.data[i * sizeDataY + j] != seg.data[ii * sizeDataY + jj]) {
										if (potential.data[i * sizeDataY + j] > ptrStr[seg.data[i * sizeDataY + j]].maxPotentialWithNeighbor) {
											ptrStr[seg.data[i * sizeDataY + j]].maxPotentialWithNeighbor = potential.data[i * sizeDataY + j];
											ptrStr[seg.data[i * sizeDataY + j]].structureInContact = seg.data[ii * sizeDataY + jj];
										}
									}
								}
							}
						}
					}
				}
			}
		}


		change = false;
		for (i = 1; i <= oldNbInfoStruct; i++) {
			if (ptrStr[i].structureInContact > 0) {
				if (ptrStr[i].locked == false && ptrStr[ptrStr[i].structureInContact].locked == false) {
					const double ratio = (ptrStr[i].maxPotentialOnTheStructure - ptrStr[i].maxPotentialWithNeighbor) / (
						                     ptrStr[i].maxPotentialOnTheStructure - ptrStr[i].minPotentialOnTheStructure);

					if (ratio < thresholdRatio) {
						if (ptrStr[i].maxPotentialOnTheStructure < ptrStr[ptrStr[i].structureInContact].maxPotentialOnTheStructure) {
							ptrStr[i].structureToBecaume = ptrStr[i].structureInContact;
							nbInfoStruct--;
							ptrStr[ptrStr[i].structureInContact].locked = true;
							ptrStr[i].locked = true;
							change = true;
						}
					}
				}
			}
		}

		label = 0;

		for (i = 1; i <= oldNbInfoStruct; i++) {
			if (ptrStr[i].structureToBecaume == i) {
				label++;
				ptrStr[i].newIndex = label;
			}
		}


		for (i = 0; i < sizeDataX; i++) {
			for (j = 0; j < sizeDataY; j++) {
				if (seg.data[i * sizeDataY + j] != 0) {
					seg.data[i * sizeDataY + j] = (unsigned short)ptrStr[ptrStr[seg.data[i * sizeDataY + j]].structureToBecaume].newIndex;
				}
			}
		}
	}

	seg = (unsigned short)0;

	label = 0;
	if (nbInfoStruct == 0) { return 0; }
	for (i = 1; i <= (long)nbInfoStruct; i++) {
		if (ptrStr[i].structureToBecaume == i) {
			label++;
			seg.data[ptrStr[i].x * sizeDataY + ptrStr[i].y] = (unsigned short)label;
		}
	}


	free(ptrStr);

	return label - 1;
}

/********************************************************************************************************************
int Segmentation2D::getLocalMax(...) compute local maximum of the image an mark them in the image seg.
FltArray &image: input image
USrtArray &localMax: local maximum
float threshold: a point to be considered as a local maximum have to be > to threshold.
********************************************************************************************************************/
int Segmentation2D::getLocalMax(FltArray& potential, USrtArray& seg, const float threshold)
{
	seg = (unsigned short)0;
	long ii, jj;
	long label = 1;
	const long sizeDataY = potential.size(1);
	const long sizeDataX = potential.size(0);


	for (long i = 0; i < sizeDataX; i++) {
		for (long j = 0; j < sizeDataY; j++) {
			if (potential.data[i * sizeDataY + j] > threshold) {
				bool err = false;
				//test if it is a new maximum from the vicinity
				for (ii = __max(i - 1, 0); ii <= __min(i + 1, sizeDataX - 1) && !err; ii++) {
					for (jj = __max(j - 1, 0); jj <= __min(j + 1, sizeDataY - 1) && !err; jj++) {
						if (ii != i || jj != j) { if (potential.data[i * sizeDataY + j] < potential.data[ii * sizeDataY + jj]) { err = true; } }
					}
				}
				if (!err) { //new maximum found
					unsigned short maxi = 0;
					err = false;
					for (ii = __max(i - 1, 0); ii <= __min(i + 1, sizeDataX - 1) && !err; ii++) {
						for (jj = __max(j - 1, 0); jj <= __min(j + 1, sizeDataY - 1) && !err; jj++) {
							if (ii != i || jj != j) {
								if (maxi < seg.data[ii * sizeDataY + jj]) {
									maxi = seg.data[ii * sizeDataY + jj];
									err = true;
								}
							}
						}
					}
					if (maxi > 0) { seg.data[i * sizeDataY + j] = maxi; }
					else {
						seg.data[i * sizeDataY + j] = (unsigned short)label;
						label++;
					}
				}
			}
		}
	}

	return label - 1;
}


/********************************************************************************************************************
int Segmentation2D::watershed(...) compute watershed segmentation insied the reagion upper the threshold value
FltArray &potential: defines the potential function the makes the watershed segmented area growing
USrtArray &seg: it contains the result of the watershed segmentation.
********************************************************************************************************************/
int Segmentation2D::watershed(FltArray& potential, USrtArray& seg, const float threshold) const
{
	float cle;

	coor2D p;
	coor2D p1;


	const long sizeDataY = potential.size(1);
	const long sizeDataX = potential.size(0);

	Heap<coor2D> tas(sizeDataY * sizeDataX);

	// Initialisation des marqueurs
	for (p.y = 0; p.y < sizeDataY; p.y++) {
		for (p.x = 0; p.x < sizeDataX; p.x++) {
			bool vois = false;

			if (seg.data[p.x * sizeDataY + p.y] != 0) {
				for (p1.x = p.x - 1; p1.x <= p.x + 1 && !vois; p1.x++) {
					for (p1.y = p.y - 1; p1.y <= p.y + 1 && !vois; p1.y++) {
						if (p1.x != p.x || p1.y != p.y) {
							if (p1.x >= 0) {
								if (p1.x < potential.size(0)) {
									if (p1.y >= 0) { if (p1.y < potential.size(1)) { if (seg.data[p1.x * sizeDataY + p1.y] == 0) { vois = true; } } }
								}
							}
						}
					}
				}

				// On empile les points frontieres.
				if (vois) {
					cle = -potential.data[p1.x * sizeDataY + p1.y];
					tas.Push(p, cle);
				}
			}
		}
	}


	// Croissance prioritaire.
	int it = 0;
	while (!tas.Empty()) {
		it = it + 1;
		p = tas.Pop();
		for (p1.x = p.x - 1; p1.x <= p.x + 1; p1.x++) {
			for (p1.y = p.y - 1; p1.y <= p.y + 1; p1.y++) {
				if (p1.x != p.x || p1.y != p.y) {
					if (p1.x >= 0) {
						if (p1.x < potential.size(0)) {
							if (p1.y >= 0) {
								if (p1.y < potential.size(1)) {
									if (potential.data[p1.x * sizeDataY + p1.y] > threshold) {
										if (seg.data[p1.x * sizeDataY + p1.y] == 0) {
											// On empile ce nouveau avec une priorité legèrement supérieure
											// a sa vrai valeur de potentiel.
											cle = (float)(-potential.data[p1.x * sizeDataY + p1.y] + (float)it * EPSILON);
											tas.Push(p1, cle);
											seg.data[p1.x * sizeDataY + p1.y] = seg.data[p.x * sizeDataY + p.y];
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	return nbInfoStruct;
}


/********************************************************************************************************************
int Segmentation2D::watershedIntelligent(...) compute watershed segmentation insied the reagion upper the threshold value.
It tries to minimize the oversegmentation.
FltArray &potential: defines the potential function the makes the watershed segmented area growing
USrtArray &seg: it contains the result of the watershed segmentation.
float threshold: define the threshold value to distinguish the segmented structures
float thresholdRatio: value 0<= X <=1. Makes a fusion to minimize the oversegmentation
						if X=0 simple watershed
						if X=1 no watershed
return value: number of detected structures
********************************************************************************************************************/
int Segmentation2D::watershedIntelligent(FltArray& potential, USrtArray& seg, const float threshold, const float thresholdRatio)
{
	int ret;

	//first, find the local maximum to set the different germinates
	if ((nbInfoStruct = getLocalMax(potential, seg, threshold)) < 0) { return nbInfoStruct; }

	//first watershed
	if ((ret = watershed(potential, seg, threshold)) < 0) { return ret; }

	if (thresholdRatio == 0) {
		mixIndex(seg);
		return nbInfoStruct;
	}

	//filter the germinates to minimize the oversegmentation
	if ((ret = filterStructureWS(potential, seg, thresholdRatio)) < 0) { return ret; }

	//first watershed
	if ((ret = watershed(potential, seg, threshold)) < 0) { return ret; }
	mixIndex(seg);
	return nbInfoStruct;
}

#undef MAX_ELEMENTS
#undef EPSILON
#undef DATA
#undef SEGHoles
#undef SEG
