#include "segmentation_2D.h"

#define SEG(i,j)	seg.data[(i)*sizeDataY+(j)]
#define SEGSTRUCT(i,j)	segStructure.data[(i)*sizeDataY+(j)]
#define SEGHoles(i,j)	segHoles.data[(i)*sizeDataY+(j)]
#define DATA(i,j)	data.data[(i)*sizeDataY+(j)]
#define TEST_ALLOC(ptr) if((ptr)==NULL){printf("Error of allocation at line %d, file %s\n", __LINE__, __FILE__);return(ALLOC_ERR);}

inline long min(long a, long b) { return (a <= b ? a : b); }
inline long min(long a, long b, long c, long d) { return (min(min(a, b), min(c, d))); }


inline long max(long a, long b) { return (a >= b ? a : b); }

inline double min(double a, double b) { return (a <= b ? a : b); }

inline double max(double a, double b) { return (a >= b ? a : b); }

void** Segment2D::getTabInfo() { return ((void**)tabInfo); }

void Segment2D::desalloc()
{
	InfoStruct2D* ptr;
	ptr = headInfo;
	while (ptr != NULL) {
		headInfo = ptr;
		ptr = ptr->next;
		delete headInfo;
	}
	if (tabInfo != NULL) delete [] tabInfo;
	tabInfo = NULL;
	nbInfoStruct = 0;
	headInfo = NULL;
}

Segment2D::Segment2D(unsigned long _nbInfoStruct)
{
	headInfo = NULL;
	nbInfoStruct = _nbInfoStruct;
	tabInfo = NULL;
}

Segment2D::Segment2D()
{
	headInfo = NULL;
	nbInfoStruct = 0;
	tabInfo = NULL;
}

Segment2D::~Segment2D() { desalloc(); }


int Segment2D::alloc(unsigned long nbStruct)
{
	InfoStruct2D* ptr;
	unsigned long ind;

	desalloc();

	tabInfo = new InfoStruct2D*[nbStruct + 1];
	if (tabInfo == NULL) return (S_ALLOC_ERR);


	if (nbStruct > 0) {
		headInfo = new InfoStruct2D(1);

		tabInfo[0] = NULL;
		tabInfo[1] = headInfo;
		ptr = headInfo;
		for (ind = 2; ind <= nbStruct; ind++) {
			ptr->next = new InfoStruct2D(ind);
			tabInfo[ind] = ptr->next;
			if (ptr->next == NULL) return (S_ALLOC_ERR);
			ptr = ptr->next;
		}
	}

	nbInfoStruct = nbStruct;
	return (RET_OK);
}


unsigned long Segment2D::filterData(USrtArray& seg, double* param, unsigned short* useParam)
{
	unsigned short sizeDataY = seg.size(1);

	InfoStruct2D* ptr = headInfo;
	InfoStruct2D* ptrMax;
	InfoStruct2D* ptrOld = NULL;
	unsigned long ind = 0;
	long i, indexFilter;
	bool unselected;
	float valMax;
	long nbInfoStructOld = nbInfoStruct;
	long oldIndex;

	char* tabBound = NULL;
	if (useParam[12] == 1) {
		tabBound = new char[nbInfoStruct + 1];
		for (i = 0; i < (long)(nbInfoStruct + 1); i++) tabBound[i] = 0;
		for (i = 0; i < seg.size(0); i++) {
			if (seg.data[i * sizeDataY] != 0) tabBound[seg.data[i * sizeDataY]] = 1;
			if (seg.data[i * sizeDataY + seg.size(1) - 1] != 0) tabBound[seg.data[i * sizeDataY + seg.size(1) - 1]] = 1;
		}
		for (i = 1; i < seg.size(1) - 1; i++) {
			if (seg.data[i] != 0) tabBound[seg.data[i]] = 1;
			if (seg.data[(seg.size(0) - 1) * sizeDataY + i] != 0) tabBound[seg.data[(seg.size(0) - 1) * sizeDataY + i]] = 1;
		}
	}

	while (ptr != NULL) {
		ind++;
		unselected = false;
		if ((useParam[0] == 1) && (ptr->data.surf <= param[0])) unselected = true;                       //minSurf
		else if ((useParam[1] == 1) && (ptr->data.surf >= param[1])) unselected = true;                  //maxSurf
		else if ((useParam[2] == 1) && (ptr->data.mean * ptr->data.surf <= param[2])) unselected = true; //minInten
		else if ((useParam[3] == 1) && (ptr->data.mean * ptr->data.surf >= param[3])) unselected = true; //maxInten
		else if ((useParam[4] == 1) && (ptr->data.max.val <= param[4])) unselected = true;               //minMaxi
		else if ((useParam[5] == 1) && (ptr->data.max.val >= param[5])) unselected = true;               //maxMaxi
		else if (useParam[7] == 1) {
			if (min(ptr->data.sigma.x, ptr->data.sigma.y) / max(ptr->data.sigma.x, ptr->data.sigma.y) >= param[7]) unselected = true; //minMorphParam
		}
		else if (useParam[8] == 1) {
			if (min(ptr->data.sigma.x, ptr->data.sigma.y) / max(ptr->data.sigma.x, ptr->data.sigma.y) <= param[8]) unselected = true; //maxMorphParam
		}
		else if (useParam[12] == 1) if (tabBound[ind] == 1) unselected = true; //no struct on boundaries
		if (unselected) {
			if (ptrOld == NULL) {
				headInfo = ptr->next;
				delete ptr;
				ptr = headInfo;
			}
			else {
				ptrOld->next = ptr->next;
				delete ptr;
				ptr = ptrOld->next;
			}
			tabInfo[ind] = NULL;
			nbInfoStruct--;
		}
		else {
			ptrOld = ptr;
			ptr = ptr->next;
		}
	}
	if (useParam[12] == 1) delete [] tabBound;


	if (useParam[6] == 1) {
		ptr = headInfo;
		ind = 0;
		ptrOld = NULL;
		while (ptr != NULL) {
			ind++;
			unselected = false;
			for (i = 1; (i <= nbInfoStructOld) && !unselected; i++)
				if (tabInfo[i] != NULL) {
					double d = sqrt((double)((ptr->data.max.x - tabInfo[i]->data.max.x) * (ptr->data.max.x - tabInfo[i]->data.max.x) + (ptr->data.max.y - tabInfo[i]->data.max.y) * (ptr->data.max.y - tabInfo[i]->data.max.y)));
					if ((d > 1) && (d < param[6]) && (tabInfo[i]->data.max.val > ptr->data.max.val)) unselected = true;
				}

			if (unselected) {
				if (ptrOld == NULL) {
					headInfo = ptr->next;
					oldIndex = ptr->data.index;
					delete ptr;
					ptr = headInfo;
				}
				else {
					ptrOld->next = ptr->next;
					oldIndex = ptr->data.index;
					delete ptr;
					ptr = ptrOld->next;
				}
				tabInfo[oldIndex] = NULL;
				nbInfoStruct--;
			}
			else {
				ptrOld = ptr;
				ptr = ptr->next;
			}
		}
	}


	for (indexFilter = 9; indexFilter <= 11; indexFilter++)
		if ((useParam[indexFilter] == 1) && (nbInfoStruct > param[indexFilter]) && ((long)param[indexFilter] > 0)) {
			//N first structure with the biggest surface (indexFilter=9), the biggest surface (indexFilter=10) or 
			//the biggest total intensity (indexFilter=11)

			ptrOld = NULL;
			for (i = 0; i < (long)param[indexFilter]; i++) {
				if (ptrOld == NULL) ptr = headInfo;
				else ptr = ptrOld->next;

				switch (indexFilter) {
					case 9: valMax = (float)ptr->data.surf;
						break;
					case 10: valMax = (float)(ptr->data.surf * ptr->data.mean);
						break;
					case 11: valMax = (float)(ptr->data.max.val);
						break;
				}

				ptrMax = NULL;

				while (ptr->next != NULL) {
					switch (indexFilter) {
						case 9: if (valMax <= ptr->next->data.surf) {
								valMax = (float)ptr->next->data.surf;
								ptrMax = ptr;
							}
							break;
						case 10: if (valMax <= (long)ptr->next->data.surf * ptr->next->data.mean) {
								valMax = (float)((long)ptr->next->data.surf * ptr->next->data.mean);
								ptrMax = ptr;
							}
							break;
						case 11: if (valMax <= (float)ptr->next->data.max.val) {
								valMax = (float)(ptr->next->data.max.val);
								ptrMax = ptr;
							}
							break;
					}
					ptr = ptr->next;
				}
				if (ptrMax != NULL) {
					if (ptrOld != NULL) {
						ptr = ptrOld->next;
						ptrOld->next = ptrMax->next;
						ptrMax->next = ptrMax->next->next;
						ptrOld->next->next = ptr;
						ptrOld = ptrOld->next;
					}
					else {
						ptr = headInfo;
						headInfo = ptrMax->next;
						ptrMax->next = headInfo->next;
						headInfo->next = ptr;
						ptrOld = headInfo;
					}
				}
				else if (ptrOld == NULL) ptrOld = headInfo;
				else ptrOld = ptrOld->next;
			}
			ptr = ptrOld->next;
			ptrOld->next = NULL;
			while (ptr != NULL) {
				ptrMax = ptr->next;
				tabInfo[ptr->data.index] = NULL;
				delete ptr;
				nbInfoStruct--;
				ptr = ptrMax;
			}
		}


	ptr = headInfo;
	for (ind = 1; ind <= nbInfoStruct; ind++) {
		ptr->data.index = ind;
		ptr = ptr->next;
	}


	for (i = 0; i < (long)seg.nbPoints; i++)
		if (seg.data[i] != 0) {
			if (tabInfo[seg.data[i]] == NULL) seg.data[i] = 0;
			else seg.data[i] = (unsigned short)tabInfo[seg.data[i]]->data.index;
		}


	delete [] tabInfo;
	tabInfo = new InfoStruct2D*[nbInfoStruct + 1];
	tabInfo[0] = NULL;
	tabInfo[1] = headInfo;
	//tabInfo[1]->data.index=1;
	ptr = headInfo;
	for (ind = 2; ind <= nbInfoStruct; ind++) {
		tabInfo[ind] = ptr->next;
		//tabInfo[ind]->data.index=ind;
		ptr = ptr->next;
	}

	return (nbInfoStruct);
}


float Segment2D::tic(void)
{
	SYSTEMTIME st_time;
	GetSystemTime(&st_time);
	return ((float)(st_time.wMilliseconds + 1000 * (st_time.wSecond + 60 * (st_time.wMinute + 60 * st_time.wHour))));
}

int Segment2D::computeStatistics(FltArray& data, USrtArray& seg, float* waveData, unsigned short gaussFit, bool doRealoc)
{
	unsigned short i, j;
	double val;
	unsigned long ind;
	unsigned long cptPoints = 0;
	int ret;

	unsigned short sizeDataX = data.size(0);
	unsigned short sizeDataY = data.size(1);


	if (doRealoc) if ((ret = alloc(nbInfoStruct)) < 0) return (ret);


	//initalisation of array

	for (ind = 1; ind <= nbInfoStruct; ind++) {
		tabInfo[ind]->data.peri = 0;
		tabInfo[ind]->data.surf = 0;
		tabInfo[ind]->data.morpho = 0.;
		tabInfo[ind]->data.max.x = 0;
		tabInfo[ind]->data.max.y = 0;
		tabInfo[ind]->data.max.val = 0.;
		tabInfo[ind]->data.min.x = 0;
		tabInfo[ind]->data.min.y = 0;
		tabInfo[ind]->data.min.val = 0;
		tabInfo[ind]->data.stdev = 0.;
		tabInfo[ind]->data.moment.x = 0.;
		tabInfo[ind]->data.moment.y = 0.;
		tabInfo[ind]->data.momentSurf.x = 0.;
		tabInfo[ind]->data.momentSurf.y = 0.;
		tabInfo[ind]->data.mean = 0.;
		tabInfo[ind]->data.intensity = 0.;
		tabInfo[ind]->data.var.xx = 0.;
		tabInfo[ind]->data.var.yy = 0.;
		tabInfo[ind]->data.var.xy = 0.;
		tabInfo[ind]->data.sigma.x = 0.;
		tabInfo[ind]->data.sigma.y = 0.;
		tabInfo[ind]->data.angle = 0.;
		tabInfo[ind]->data.boundingBox.xMax = 0;
		tabInfo[ind]->data.boundingBox.yMax = 0;
		tabInfo[ind]->data.boundingBox.xMin = sizeDataX;
		tabInfo[ind]->data.boundingBox.yMin = sizeDataY;
		tabInfo[ind]->data.nbSubstructures = 0;
	}

	for (i = 0; i < sizeDataX; i++)
		for (j = 0; j < sizeDataY; j++)
			if (seg.data[i * sizeDataY + j] > 0) {
				cptPoints++;
				val = data.data[i * sizeDataY + j];
				ind = seg.data[i * sizeDataY + j];


				tabInfo[ind]->data.surf++;
				tabInfo[ind]->data.moment.x += val * (double)i;
				tabInfo[ind]->data.moment.y += val * (double)j;
				tabInfo[ind]->data.momentSurf.x += i;
				tabInfo[ind]->data.momentSurf.y += j;
				tabInfo[ind]->data.mean += val;
				tabInfo[ind]->data.var.xx += val * i * i;
				tabInfo[ind]->data.var.yy += val * j * j;
				tabInfo[ind]->data.var.xy += val * i * j;
				tabInfo[ind]->data.stdev += val * val;

				if (tabInfo[ind]->data.boundingBox.xMax < i) tabInfo[ind]->data.boundingBox.xMax = i;
				if (tabInfo[ind]->data.boundingBox.yMax < j) tabInfo[ind]->data.boundingBox.yMax = j;
				if (tabInfo[ind]->data.boundingBox.xMin > i) tabInfo[ind]->data.boundingBox.xMin = i;
				if (tabInfo[ind]->data.boundingBox.yMin > j) tabInfo[ind]->data.boundingBox.yMin = j;


				if ((val > tabInfo[ind]->data.max.val) || (tabInfo[ind]->data.surf == 1)) {
					tabInfo[ind]->data.max.val = val;
					tabInfo[ind]->data.max.x = i;
					tabInfo[ind]->data.max.y = j;
				}

				if ((val < tabInfo[ind]->data.min.val) || (tabInfo[ind]->data.surf == 1)) {
					tabInfo[ind]->data.min.val = val;
					tabInfo[ind]->data.min.x = i;
					tabInfo[ind]->data.min.y = j;
				}

				if (i == 0) tabInfo[ind]->data.peri++;
				if (i == sizeDataX - 1) tabInfo[ind]->data.peri++;
				if (j == 0) tabInfo[ind]->data.peri++;
				if (j == sizeDataY - 1) tabInfo[ind]->data.peri++;
				if ((i != 0) && (i != sizeDataX - 1) && (j != 0) && (j != sizeDataY - 1)) {
					if (seg.data[(i - 1) * sizeDataY + j] != seg.data[i * sizeDataY + j]) tabInfo[ind]->data.peri++;
					if (seg.data[(i + 1) * sizeDataY + j] != seg.data[i * sizeDataY + j]) tabInfo[ind]->data.peri++;
					if (seg.data[i * sizeDataY + j - 1] != seg.data[i * sizeDataY + j]) tabInfo[ind]->data.peri++;
					if (seg.data[i * sizeDataY + j + 1] != seg.data[i * sizeDataY + j]) tabInfo[ind]->data.peri++;
				}
			}


	for (ind = 1; ind <= nbInfoStruct; ind++) {
		if (tabInfo[ind]->data.surf > 0) {
			if (gaussFit == 5) {
				/*************************************/
				tabInfo[ind]->data.moment.x = 0.;
				tabInfo[ind]->data.moment.y = 0.;
				double mean = 0.;

				for (int ii = -2; ii <= 2; ii++) // sizeDataX
				{
					for (int jj = -2; jj <= 2; jj++) // sizeDataY
					{
						int refx = (int)tabInfo[ind]->data.max.x;
						int refy = (int)tabInfo[ind]->data.max.y;

						if (((ii + refx) >= 0) && ((ii + refx) < sizeDataX) && ((jj + refy) < sizeDataY) && ((jj + refy) >= 0)) {
							val = (double)waveData[(ii + refx) * sizeDataY + (jj + refy)];
							tabInfo[ind]->data.moment.x += val * (double)(ii + refx);
							tabInfo[ind]->data.moment.y += val * (double)(jj + refy);
							mean += val;
						}
					}
				}

				tabInfo[ind]->data.moment.x /= mean;
				tabInfo[ind]->data.moment.y /= mean;
				tabInfo[ind]->data.mean /= tabInfo[ind]->data.surf;
				tabInfo[ind]->data.intensity = tabInfo[ind]->data.mean;

				/***********************************/
			}
			else {
				tabInfo[ind]->data.moment.x /= tabInfo[ind]->data.mean;
				tabInfo[ind]->data.moment.y /= tabInfo[ind]->data.mean;

				tabInfo[ind]->data.momentSurf.x /= tabInfo[ind]->data.surf;
				tabInfo[ind]->data.momentSurf.y /= tabInfo[ind]->data.surf;
				tabInfo[ind]->data.var.xx /= tabInfo[ind]->data.mean;
				tabInfo[ind]->data.var.yy /= tabInfo[ind]->data.mean;
				tabInfo[ind]->data.var.xy /= tabInfo[ind]->data.mean;
				tabInfo[ind]->data.var.xx -= tabInfo[ind]->data.moment.x * tabInfo[ind]->data.moment.x;
				tabInfo[ind]->data.var.yy -= tabInfo[ind]->data.moment.y * tabInfo[ind]->data.moment.y;
				tabInfo[ind]->data.var.xy -= tabInfo[ind]->data.moment.x * tabInfo[ind]->data.moment.y;
				tabInfo[ind]->data.intensity = tabInfo[ind]->data.mean;
				tabInfo[ind]->data.stdev -= tabInfo[ind]->data.mean;
				tabInfo[ind]->data.mean /= tabInfo[ind]->data.surf;
				tabInfo[ind]->data.stdev /= tabInfo[ind]->data.surf;
				tabInfo[ind]->data.stdev = sqrt(tabInfo[ind]->data.stdev);
			}
		}
	}


	return (cptPoints);
}


unsigned short Segment2D::segmentation(FltArray& data, USrtArray& seg, double val, bool diagonal, bool fillHoles)
{
	if (!fillHoles) return (segmentationPositiveValues(data, seg, val, diagonal));
	int ret;
	long i;
	long sizeDataX = data.size(0);
	long sizeDataY = data.size(1);

	USrtArray* segHoles = new USrtArray;
	if ((ret = segHoles->alloc(sizeDataX, sizeDataY)) < 0) return (ret);

	ret = segmentationNegativeValues(data, *segHoles, val, !diagonal);

	unsigned long* ptr;
	ptr = (unsigned long*)malloc(sizeof(long) * (ret + 1));
	for (i = 0; i <= ret; i++) ptr[i] = i;
	for (i = 0; i < sizeDataY; i++) {
		ptr[(*segHoles)(0, i)] = 0;
		ptr[(*segHoles)(sizeDataX - 1, i)] = 0;
	}
	for (i = 0; i < sizeDataX; i++) {
		ptr[(*segHoles)(i, 0)] = 0;
		ptr[(*segHoles)(i, sizeDataY - 1)] = 0;
	}
	for (i = 0; i < sizeDataX * sizeDataY; i++) segHoles->data[i] = (unsigned short)ptr[segHoles->data[i]];
	ret = segmentationPositiveValues(data, seg, val, diagonal, *segHoles);
	delete segHoles;

	// test to remove the memory leak detected by ccpcheck
	//20090206
	free(ptr);
	return (ret);
}


unsigned short Segment2D::segmentationNegativeValues(FltArray& data, USrtArray& seg, double val, bool diagonal)
{
	/**
		segmentation2D segment the 2D signal "data" which is lower than "val". It gives a map of segmentation "seg". 
		A label is set to each structure in the map seg. 
	**/
	unsigned long* ptr;

	unsigned long* ptr2;
	long i, j;
	unsigned long nSeg;
	int structFinded;
	int mini;
	//int max;
	long sizePtr = MAX_TAB;
	long sizeDataX = data.size(0);
	long sizeDataY = data.size(1);


	//allocation of ptr
	//ptr is an array of index to enumerate different structure
	//if ptr[i]=j then structure i and structure j are same.
	ptr = (unsigned long*)malloc(sizeof(long) * MAX_TAB);
	//initialisation of ptr
	for (i = 0; i < MAX_TAB; i++) ptr[i] = (unsigned short)i;
	//nSeg indicate the number of the next unused structure
	nSeg = 1;


	//initialisation of seg
	seg = (unsigned short)0;

	//first pixel
	if (fabs(data.data[0]) <= val) {
		seg.data[0] = (unsigned short)nSeg;
		inc(nSeg, &ptr, sizePtr);
	}

	//first line
	for (j = 1; j < sizeDataY; j++)
		if (fabs(data.data[j]) <= val) {
			if (seg.data[j - 1] == 0) {
				seg.data[j] = (unsigned short)nSeg;
				inc(nSeg, &ptr, sizePtr);
			}
			else seg.data[j] = seg.data[j - 1];
		}

	//the rest
	for (i = 1; i < sizeDataX; i++) {
		if (DATA(i, 0) <= val) {
			if (SEG(i-1, 0) != 0)
				SEG(i, 0) = SEG(i-1, 0);
			else if (SEG(i-1, 1) != 0)
				SEG(i, 0) = SEG(i-1, 1);
			else {
				SEG(i, 0) = (unsigned short)nSeg;
				inc(nSeg, &ptr, sizePtr);
			}
		}
		for (j = 1; j < sizeDataY; j++) {
			structFinded = 0;
			if (DATA(i, j) <= val) {
				if (SEG(i, j-1) != 0) {
					structFinded = SEG(i, j-1);
					SEG(i, j) = SEG(i, j-1);
				}

				if (diagonal)
					if (SEG(i-1, j-1) != 0) {
						if (structFinded == 0) {
							SEG(i, j) = SEG(i-1, j-1);
							structFinded = SEG(i-1, j-1);
						}
						else if (structFinded != SEG(i-1, j-1)) {
							mini = min(structFinded,SEG(i-1, j-1), ptr[structFinded], ptr[SEG(i-1, j-1)]);
							if ((long)ptr[SEG(i-1, j-1)] <= mini)
								for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i-1, j-1)) if (ptr[SEG(i-1, j-1)] == ptr[k]) ptr[k] = mini;
							if ((long)ptr[SEG(i, j)] <= mini)
								for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i, j)) if (ptr[SEG(i, j)] == ptr[k]) ptr[k] = mini;
							ptr[SEG(i-1, j-1)] = mini;
							ptr[structFinded] = mini;
							structFinded = mini;
						}
					}

				if (SEG(i-1, j) != 0) {
					if (structFinded == 0) {
						SEG(i, j) = SEG(i-1, j);
						structFinded = SEG(i-1, j);
					}
					else if (structFinded != SEG(i-1, j)) {
						mini = min(structFinded,SEG(i-1, j), ptr[structFinded], ptr[SEG(i-1, j)]);
						if ((long)ptr[SEG(i-1, j)] > mini)
							for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i-1, j)) if (ptr[SEG(i-1, j)] == ptr[k]) ptr[k] = mini;
						if ((long)ptr[SEG(i, j)] > mini)
							for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i, j)) if (ptr[SEG(i, j)] == ptr[k]) ptr[k] = mini;
						ptr[SEG(i-1, j)] = mini;
						ptr[structFinded] = mini;
						structFinded = mini;
					}
				}

				if (diagonal)
					if (j + 1 < sizeDataY)
						if (SEG(i-1, j+1) != 0) {
							if (structFinded == 0) {
								SEG(i, j) = SEG(i-1, j+1);
								structFinded = SEG(i-1, j+1);
							}
							else if (structFinded != SEG(i-1, j+1)) {
								mini = min(structFinded,SEG(i-1, j+1), ptr[structFinded], ptr[SEG(i-1, j+1)]);
								if ((long)ptr[SEG(i-1, j+1)] > mini)
									for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i-1, j+1)) if (ptr[SEG(i-1, j+1)] == ptr[k]) ptr[k] = mini;
								if ((long)ptr[SEG(i, j)] > mini)
									for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i, j)) if (ptr[SEG(i, j)] == ptr[k]) ptr[k] = mini;
								ptr[SEG(i-1, j+1)] = mini;
								ptr[structFinded] = mini;
								structFinded = mini;
							}
						}

				if (structFinded == 0) {
					SEG(i, j) = (unsigned short)nSeg;
					inc(nSeg, &ptr, sizePtr);
				}
			} // end DATA(i,j)<=val
		}     // end for j
	}         //end for i


	//assembling structures


	ptr2 = (unsigned long*)malloc(sizeof(long) * nSeg);
	for (i = 0; i < (long)nSeg; i++) ptr2[i] = 0;
	structFinded = 1;
	for (i = 1; i < (long)nSeg; i++)
		if (ptr2[ptr[i]] == 0) {
			ptr2[ptr[i]] = structFinded;
			structFinded++;
		}


	for (i = 0; i < sizeDataX * sizeDataY; i++) seg.data[i] = (unsigned short)ptr2[ptr[seg.data[i]]];

	free(ptr);
	free(ptr2);
	nbInfoStruct = structFinded - 1;
	mixIndex(seg);

	return (structFinded - 1);
}


unsigned short Segment2D::segmentationPositiveValues(FltArray& data, USrtArray& seg, double val, bool diagonal)
{
	/**
		segmentation2D segment the 2D signal "data" which is upper than "val". It gives a map of segmentation "seg". 
		A label is set to each structure in the map seg. 
	**/
	unsigned long* ptr;

	unsigned long* ptr2;
	long i, j;
	unsigned long nSeg;
	int structFinded;
	int mini;
	//int max;
	long sizePtr = MAX_TAB;
	long sizeDataX = data.size(0);
	long sizeDataY = data.size(1);


	//allocation of ptr
	//ptr is an array of index to enumerate different structure
	//if ptr[i]=j then structure i and structure j are same.
	ptr = (unsigned long*)malloc(sizeof(long) * MAX_TAB);
	//initialisation of ptr
	for (i = 0; i < MAX_TAB; i++) ptr[i] = (unsigned short)i;
	//nSeg indicate the number of the next unused structure
	nSeg = 1;


	//initialisation of seg
	seg = (unsigned short)0;

	//first pixel
	if (fabs(data.data[0]) > val) {
		seg.data[0] = (unsigned short)nSeg;
		inc(nSeg, &ptr, sizePtr);
	}

	//first line
	for (j = 1; j < sizeDataY; j++)
		if (fabs(data.data[j]) > val) {
			if (seg.data[j - 1] == 0) {
				seg.data[j] = (unsigned short)nSeg;
				inc(nSeg, &ptr, sizePtr);
			}
			else seg.data[j] = seg.data[j - 1];
		}

	//the rest
	for (i = 1; i < sizeDataX; i++) {
		if (DATA(i, 0) > val) {
			if (SEG(i-1, 0) != 0)
				SEG(i, 0) = SEG(i-1, 0);
			else if (SEG(i-1, 1) != 0)
				SEG(i, 0) = SEG(i-1, 1);
			else {
				SEG(i, 0) = (unsigned short)nSeg;
				inc(nSeg, &ptr, sizePtr);
			}
		}
		for (j = 1; j < sizeDataY; j++) {
			structFinded = 0;
			if (DATA(i, j) > val) {
				if (SEG(i, j-1) != 0) {
					structFinded = SEG(i, j-1);
					SEG(i, j) = SEG(i, j-1);
				}

				if (diagonal)
					if (SEG(i-1, j-1) != 0) {
						if (structFinded == 0) {
							SEG(i, j) = SEG(i-1, j-1);
							structFinded = SEG(i-1, j-1);
						}
						else if (structFinded != SEG(i-1, j-1)) {
							mini = min(structFinded,SEG(i-1, j-1), ptr[structFinded], ptr[SEG(i-1, j-1)]);
							if ((long)ptr[SEG(i-1, j-1)] > mini)
								for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i-1, j-1)) if (ptr[SEG(i-1, j-1)] == ptr[k]) ptr[k] = mini;
							if ((long)ptr[SEG(i, j)] > mini)
								for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i, j)) if (ptr[SEG(i, j)] == ptr[k]) ptr[k] = mini;
							ptr[SEG(i-1, j-1)] = mini;
							ptr[structFinded] = mini;
							structFinded = mini;
						}
					}

				if (SEG(i-1, j) != 0) {
					if (structFinded == 0) {
						SEG(i, j) = SEG(i-1, j);
						structFinded = SEG(i-1, j);
					}
					else if (structFinded != SEG(i-1, j)) {
						mini = min(structFinded,SEG(i-1, j), ptr[structFinded], ptr[SEG(i-1, j)]);
						if ((long)ptr[SEG(i-1, j)] > mini)
							for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i-1, j)) if (ptr[SEG(i-1, j)] == ptr[k]) ptr[k] = mini;
						if ((long)ptr[SEG(i, j)] > mini)
							for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i, j)) if (ptr[SEG(i, j)] == ptr[k]) ptr[k] = mini;
						ptr[SEG(i-1, j)] = mini;
						ptr[structFinded] = mini;
						structFinded = mini;
					}
				}

				if (diagonal)
					if (j + 1 < sizeDataY)
						if (SEG(i-1, j+1) != 0) {
							if (structFinded == 0) {
								SEG(i, j) = SEG(i-1, j+1);
								structFinded = SEG(i-1, j+1);
							}
							else if (structFinded != SEG(i-1, j+1)) {
								mini = min(structFinded,SEG(i-1, j+1), ptr[structFinded], ptr[SEG(i-1, j+1)]);
								if ((long)ptr[SEG(i-1, j+1)] > mini)
									for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i-1, j+1)) if (ptr[SEG(i-1, j+1)] == ptr[k]) ptr[k] = mini;
								if ((long)ptr[SEG(i, j)] > mini)
									for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i, j)) if (ptr[SEG(i, j)] == ptr[k]) ptr[k] = mini;
								ptr[SEG(i-1, j+1)] = mini;
								ptr[structFinded] = mini;
								structFinded = mini;
							}
						}

				if (structFinded == 0) {
					SEG(i, j) = (unsigned short)nSeg;
					inc(nSeg, &ptr, sizePtr);
				}
			} // end DATA(i,j)>val
		}     // end for j
	}         //end for i


	//assembling structures


	ptr2 = (unsigned long*)malloc(sizeof(long) * nSeg);
	for (i = 0; i < (long)nSeg; i++) ptr2[i] = 0;
	structFinded = 1;
	for (i = 1; i < (long)nSeg; i++)
		if (ptr2[ptr[i]] == 0) {
			ptr2[ptr[i]] = structFinded;
			structFinded++;
		}


	for (i = 0; i < sizeDataX * sizeDataY; i++) seg.data[i] = (unsigned short)ptr2[ptr[seg.data[i]]];

	free(ptr);
	free(ptr2);
	nbInfoStruct = structFinded - 1;
	mixIndex(seg);

	return (structFinded - 1);
}

unsigned short Segment2D::segmentationPositiveValues(FltArray& data, USrtArray& seg, double val, bool diagonal, USrtArray& segHoles)
{
	/**
		segmentation2D segment the 2D signal "data" which is upper than "val". It gives a map of segmentation "seg". 
		A label is set to each structure in the map seg. 
	**/
	unsigned long* ptr;

	unsigned long* ptr2;
	long i, j;
	unsigned long nSeg;
	int structFinded;
	int mini;
	//int max;
	long sizePtr = MAX_TAB;
	long sizeDataX = data.size(0);
	long sizeDataY = data.size(1);


	//allocation of ptr
	//ptr is an array of index to enumerate different structure
	//if ptr[i]=j then structure i and structure j are same.
	ptr = (unsigned long*)malloc(sizeof(long) * MAX_TAB);
	//initialisation of ptr
	for (i = 0; i < MAX_TAB; i++) ptr[i] = (unsigned short)i;
	//nSeg indicate the number of the next unused structure
	nSeg = 1;


	//initialisation of seg
	seg = (unsigned short)0;

	//first pixel
	if ((fabs(data.data[0]) > val) || (segHoles.data[0] != 0)) {
		seg.data[0] = (unsigned short)nSeg;
		inc(nSeg, &ptr, sizePtr);
	}

	//first line
	for (j = 1; j < sizeDataY; j++)
		if ((fabs(data.data[j]) > val) || (segHoles.data[j] != 0)) {
			if (seg.data[j - 1] == 0) {
				seg.data[j] = (unsigned short)nSeg;
				inc(nSeg, &ptr, sizePtr);
			}
			else seg.data[j] = seg.data[j - 1];
		}

	//the rest
	for (i = 1; i < sizeDataX; i++) {
		if ((DATA(i, 0) > val) || (SEGHoles(i, 0) != 0)) {
			if (SEG(i-1, 0) != 0)
				SEG(i, 0) = SEG(i-1, 0);
			else if (SEG(i-1, 1) != 0)
				SEG(i, 0) = SEG(i-1, 1);
			else {
				SEG(i, 0) = (unsigned short)nSeg;
				inc(nSeg, &ptr, sizePtr);
			}
		}
		for (j = 1; j < sizeDataY; j++) {
			structFinded = 0;
			if ((DATA(i, j) > val) || (SEGHoles(i, j) != 0)) {
				if (SEG(i, j-1) != 0) {
					structFinded = SEG(i, j-1);
					SEG(i, j) = SEG(i, j-1);
				}

				if (diagonal)
					if (SEG(i-1, j-1) != 0) {
						if (structFinded == 0) {
							SEG(i, j) = SEG(i-1, j-1);
							structFinded = SEG(i-1, j-1);
						}
						else if (structFinded != SEG(i-1, j-1)) {
							mini = min(structFinded,SEG(i-1, j-1), ptr[structFinded], ptr[SEG(i-1, j-1)]);
							if ((long)ptr[SEG(i-1, j-1)] > mini)
								for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i-1, j-1)) if (ptr[SEG(i-1, j-1)] == ptr[k]) ptr[k] = mini;
							if ((long)ptr[SEG(i, j)] > mini)
								for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i, j)) if (ptr[SEG(i, j)] == ptr[k]) ptr[k] = mini;
							ptr[SEG(i-1, j-1)] = mini;
							ptr[structFinded] = mini;
							structFinded = mini;
						}
					}

				if (SEG(i-1, j) != 0) {
					if (structFinded == 0) {
						SEG(i, j) = SEG(i-1, j);
						structFinded = SEG(i-1, j);
					}
					else if (structFinded != SEG(i-1, j)) {
						mini = min(structFinded,SEG(i-1, j), ptr[structFinded], ptr[SEG(i-1, j)]);
						if ((long)ptr[SEG(i-1, j)] > mini)
							for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i-1, j)) if (ptr[SEG(i-1, j)] == ptr[k]) ptr[k] = mini;
						if ((long)ptr[SEG(i, j)] > mini)
							for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i, j)) if (ptr[SEG(i, j)] == ptr[k]) ptr[k] = mini;
						ptr[SEG(i-1, j)] = mini;
						ptr[structFinded] = mini;
						structFinded = mini;
					}
				}

				if (diagonal)
					if (j + 1 < sizeDataY)
						if (SEG(i-1, j+1) != 0) {
							if (structFinded == 0) {
								SEG(i, j) = SEG(i-1, j+1);
								structFinded = SEG(i-1, j+1);
							}
							else if (structFinded != SEG(i-1, j+1)) {
								mini = min(structFinded,SEG(i-1, j+1), ptr[structFinded], ptr[SEG(i-1, j+1)]);
								if ((long)ptr[SEG(i-1, j+1)] > mini)
									for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i-1, j+1)) if (ptr[SEG(i-1, j+1)] == ptr[k]) ptr[k] = mini;
								if ((long)ptr[SEG(i, j)] > mini)
									for (long k = mini; k < (long)nSeg; k++) if (k != SEG(i, j)) if (ptr[SEG(i, j)] == ptr[k]) ptr[k] = mini;
								ptr[SEG(i-1, j+1)] = mini;
								ptr[structFinded] = mini;
								structFinded = mini;
							}
						}

				if (structFinded == 0) {
					SEG(i, j) = (unsigned short)nSeg;
					inc(nSeg, &ptr, sizePtr);
				}
			} // end DATA(i,j)>val
		}     // end for j
	}         //end for i


	//assembling structures


	ptr2 = (unsigned long*)malloc(sizeof(long) * nSeg);
	for (i = 0; i < (long)nSeg; i++) ptr2[i] = 0;
	structFinded = 1;
	for (i = 1; i < (long)nSeg; i++)
		if (ptr2[ptr[i]] == 0) {
			ptr2[ptr[i]] = structFinded;
			structFinded++;
		}


	for (i = 0; i < sizeDataX * sizeDataY; i++) seg.data[i] = (unsigned short)ptr2[ptr[seg.data[i]]];

	free(ptr);
	free(ptr2);
	nbInfoStruct = structFinded - 1;
	mixIndex(seg);

	return (structFinded - 1);
}

void Segment2D::mixIndex(USrtArray& seg)
{
	USrtArray* mixer = new USrtArray(nbInfoStruct + 1);
	(*mixer)(0) = 0;
	srand(111);

	long ind, val, i;
	for (i = 1; i <= (long)nbInfoStruct; i++) { (*mixer)(i) = (unsigned short)i; }
	for (i = 1; i <= (long)nbInfoStruct; i++) {
		ind = rand() % nbInfoStruct + 1;
		val = (*mixer)(ind);
		(*mixer)(ind) = (*mixer)(i);
		(*mixer)(i) = (unsigned short)val;
	}


	for (i = 0; i < (long)seg.nbPoints; i++) { seg[i] = (*mixer)(seg.data[i]); }
	delete mixer;
}

#undef DATA
#undef SEGHoles
#undef SEG
