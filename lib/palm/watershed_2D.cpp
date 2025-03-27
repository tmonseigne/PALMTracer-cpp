#include "segmentation_2D.h"

#define MAX_ELEMENTS 10000000
#define EPSILON 0.0000001


/********************************************************************************************************************
int Segment2D::filterStructureWS(...) filter the structure of the watershed. 
FltArray &potential: input image
USrtArray &localMax: local maximum
float thresholdRatio: value to distinguish if 2 strucutres in contact have to be considered as a single structure (
							oversegmentation)
********************************************************************************************************************/
int Segment2D::filterStructureWS(FltArray& potential, USrtArray& seg, float thresholdRatio)
{
	wsFusion2D* ptrStr;
	ptrStr = (wsFusion2D*)malloc((nbInfoStruct + 1) * sizeof(wsFusion2D));
	long i, j, ii, jj;

	long sizeDataY = potential.size(1);
	long sizeDataX = potential.size(0);
	long oldNbInfoStruct;
	double ratio;
	long label;
	bool change = true;
	while (change) {
		oldNbInfoStruct = nbInfoStruct;

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


		for (i = 0; i < sizeDataX; i++)
			for (j = 0; j < sizeDataY; j++)
				if (seg.data[i * sizeDataY + j] > 0) {
					if (ptrStr[seg.data[i * sizeDataY + j]].minPotentialOnTheStructure > potential.data[i * sizeDataY + j]) ptrStr[seg.data[i * sizeDataY + j]].minPotentialOnTheStructure = potential.data[i * sizeDataY + j];
					if (ptrStr[seg.data[i * sizeDataY + j]].maxPotentialOnTheStructure < potential.data[i * sizeDataY + j]) {
						ptrStr[seg.data[i * sizeDataY + j]].x = i;
						ptrStr[seg.data[i * sizeDataY + j]].y = j;
						ptrStr[seg.data[i * sizeDataY + j]].maxPotentialOnTheStructure = potential.data[i * sizeDataY + j];
					}
					for (ii = __max(i-1, 0); ii <= __min(i+1, sizeDataX-1); ii++)
						for (jj = __max(j-1, 0); jj <= __min(j+1, sizeDataY-1); jj++)
							if ((ii != i) || (jj != j))
								if (seg.data[ii * sizeDataY + jj] != 0)
									if (seg.data[i * sizeDataY + j] != seg.data[ii * sizeDataY + jj])
										if (potential.data[i * sizeDataY + j] > ptrStr[seg.data[i * sizeDataY + j]].maxPotentialWithNeighbor) {
											ptrStr[seg.data[i * sizeDataY + j]].maxPotentialWithNeighbor = potential.data[i * sizeDataY + j];
											ptrStr[seg.data[i * sizeDataY + j]].structureInContact = seg.data[ii * sizeDataY + jj];
										}
				}


		change = false;
		for (i = 1; i <= oldNbInfoStruct; i++) {
			if (ptrStr[i].structureInContact > 0)
				if ((ptrStr[i].locked == false) && (ptrStr[ptrStr[i].structureInContact].locked == false)) {
					ratio = (ptrStr[i].maxPotentialOnTheStructure - ptrStr[i].maxPotentialWithNeighbor) / (ptrStr[i].maxPotentialOnTheStructure - ptrStr[i].minPotentialOnTheStructure);

					if (ratio < thresholdRatio)
						if (ptrStr[i].maxPotentialOnTheStructure < ptrStr[ptrStr[i].structureInContact].maxPotentialOnTheStructure) {
							ptrStr[i].structureToBecaume = ptrStr[i].structureInContact;
							nbInfoStruct--;
							ptrStr[ptrStr[i].structureInContact].locked = true;
							ptrStr[i].locked = true;
							change = true;
						}
				}
		}

		label = 0;

		for (i = 1; i <= oldNbInfoStruct; i++)
			if (ptrStr[i].structureToBecaume == i) {
				label++;
				ptrStr[i].newIndex = label;
			}


		for (i = 0; i < sizeDataX; i++) for (j = 0; j < sizeDataY; j++) { if (seg.data[i * sizeDataY + j] != 0) { seg.data[i * sizeDataY + j] = (unsigned short)ptrStr[ptrStr[seg.data[i * sizeDataY + j]].structureToBecaume].newIndex; } }
	}

	seg = (unsigned short)0;

	label = 0;
	if (nbInfoStruct == 0) return (0);
	for (i = 1; i <= (long)nbInfoStruct; i++) {
		if (ptrStr[i].structureToBecaume == i) {
			label++;
			seg.data[ptrStr[i].x * sizeDataY + ptrStr[i].y] = (unsigned short)label;
		}
	}


	free(ptrStr);

	return (label - 1);
}

/********************************************************************************************************************
int Segment2D::getLocalMax(...) compute local maximum of the image an mark them in the image seg. 
FltArray &image: input image
USrtArray &localMax: local maximum
float threshold: a point to be considered as a local maximum have to be > to threshold.
********************************************************************************************************************/
int Segment2D::getLocalMax(FltArray& image, USrtArray& seg, float threshold)
{
	seg = (unsigned short)0;
	long i, j, ii, jj;
	long label = 1;
	bool err;
	long sizeDataY = image.size(1);
	long sizeDataX = image.size(0);


	for (i = 0; i < sizeDataX; i++)
		for (j = 0; j < sizeDataY; j++)
			if (image.data[i * sizeDataY + j] > threshold) {
				err = false;
				//test if it is a new maximum from the vicinity
				for (ii = __max(i-1, 0); (ii <= __min(i+1, sizeDataX-1)) && !err; ii++) for (jj = __max(j-1, 0); (jj <= __min(j+1, sizeDataY-1)) && !err; jj++) { if ((ii != i) || (jj != j)) if (image.data[i * sizeDataY + j] < image.data[ii * sizeDataY + jj]) err = true; }
				if (!err) { //new maximum found
					unsigned short maxi = 0;
					err = false;
					for (ii = __max(i-1, 0); (ii <= __min(i+1, sizeDataX-1)) && !err; ii++)
						for (jj = __max(j-1, 0); (jj <= __min(j+1, sizeDataY-1)) && !err; jj++)
							if ((ii != i) || (jj != j))
								if (maxi < seg.data[ii * sizeDataY + jj]) {
									maxi = seg.data[ii * sizeDataY + jj];
									err = true;
								}
					if (maxi > 0) seg.data[i * sizeDataY + j] = maxi;
					else {
						seg.data[i * sizeDataY + j] = (unsigned short)label;
						label++;
					}
				}
			}

	return (label - 1);
}


/********************************************************************************************************************
int Segment2D::watershed(...) compute watershed segmentation insied the reagion upper the threshold value
FltArray &potential: defines the potential function the makes the watershed segmented area growing
USrtArray &seg: it contains the result of the watershed segmentation. 
********************************************************************************************************************/
int Segment2D::watershed(FltArray& potential, USrtArray& seg, float threshold)
{
	float cle;

	coor2D p;
	coor2D p1;


	long sizeDataY = potential.size(1);
	long sizeDataX = potential.size(0);

	Heap<coor2D> tas(sizeDataY * sizeDataX);

	// Initialisation des marqueurs
	for (p.y = 0; p.y < sizeDataY; p.y++)
		for (p.x = 0; p.x < sizeDataX; p.x++) {
			bool vois = false;

			if (seg.data[p.x * sizeDataY + p.y] != 0) {
				for (p1.x = p.x - 1; (p1.x <= p.x + 1) && !vois; p1.x++)
					for (p1.y = p.y - 1; (p1.y <= p.y + 1) && !vois; p1.y++)
						if ((p1.x != p.x) || (p1.y != p.y))
							if (p1.x >= 0)
								if (p1.x < potential.size(0)) if (p1.y >= 0) if (p1.y < potential.size(1)) if (seg.data[p1.x * sizeDataY + p1.y] == 0) vois = true;

				// On empile les points frontieres.
				if (vois) {
					cle = -potential.data[p1.x * sizeDataY + p1.y];
					tas.Push(p, cle);
				}
			}
		}


	// Croissance prioritaire.
	int it = 0;
	while (!tas.Empty()) {
		it = it + 1;
		p = tas.Pop();
		for (p1.x = p.x - 1; p1.x <= p.x + 1; p1.x++)
			for (p1.y = p.y - 1; p1.y <= p.y + 1; p1.y++) {
				if ((p1.x != p.x) || (p1.y != p.y))
					if (p1.x >= 0)
						if (p1.x < potential.size(0))
							if (p1.y >= 0)
								if (p1.y < potential.size(1))
									if (potential.data[p1.x * sizeDataY + p1.y] > threshold)
										if (seg.data[p1.x * sizeDataY + p1.y] == 0) {
											// On empile ce nouveau avec une priorité legèrement supérieure
											// a sa vrai valeur de potentiel.
											cle = (float)(-potential.data[p1.x * sizeDataY + p1.y] + (float)it * EPSILON);
											tas.Push(p1, cle);
											seg.data[p1.x * sizeDataY + p1.y] = seg.data[p.x * sizeDataY + p.y];
										}
			}
	}

	return (nbInfoStruct);
}


/********************************************************************************************************************
int Segment2D::watershedIntelligent(...) compute watershed segmentation insied the reagion upper the threshold value. 
It tries to minimize the oversegmentation.
FltArray &potential: defines the potential function the makes the watershed segmented area growing
USrtArray &seg: it contains the result of the watershed segmentation. 
float threshold: define the threshold value to distinguish the segmented structures
float thresholdRatio: value 0<= X <=1. Makes a fusion to minimize the oversegmentation
						if X=0 simple watershed
						if X=1 no watershed
return value: number of detected structures
********************************************************************************************************************/
int Segment2D::watershedIntelligent(FltArray& potential, USrtArray& seg, float threshold, float thresholdRatio)
{
	int ret;

	//first, find the local maximum to set the different germinates
	if ((nbInfoStruct = getLocalMax(potential, seg, threshold)) < 0) return (nbInfoStruct);

	//first watershed
	if ((ret = watershed(potential, seg, threshold)) < 0) return (ret);

	if (thresholdRatio == 0) {
		mixIndex(seg);
		return (nbInfoStruct);
	}

	//filter the germinates to minimize the oversegmentation
	if ((ret = filterStructureWS(potential, seg, thresholdRatio)) < 0) return (ret);

	//first watershed
	if ((ret = watershed(potential, seg, threshold)) < 0) return (ret);
	mixIndex(seg);
	return (nbInfoStruct);
}

/********************************************************************************************************************
int Segment2D::filterStructureWS(...) filter the structure of the watershed. 
USrtArray &potential: input image
USrtArray &localMax: local maximum
float thresholdRatio: value to distinguish if 2 strucutres in contact have to be considered as a single structure (
							oversegmentation)
********************************************************************************************************************/
int Segment2D::filterStructureWS(USrtArray& potential, USrtArray& seg, float thresholdRatio)
{
	wsFusion2D* ptrStr;
	ptrStr = (wsFusion2D*)malloc((nbInfoStruct + 1) * sizeof(wsFusion2D));
	long i, j, ii, jj;

	long sizeDataY = potential.size(1);
	long sizeDataX = potential.size(0);
	long oldNbInfoStruct;
	double ratio;
	long label;
	bool change = true;
	while (change) {
		oldNbInfoStruct = nbInfoStruct;


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


		for (i = 0; i < sizeDataX; i++)
			for (j = 0; j < sizeDataY; j++)
				if (seg.data[i * sizeDataY + j] > 0) {
					if (ptrStr[seg.data[i * sizeDataY + j]].minPotentialOnTheStructure > potential.data[i * sizeDataY + j]) ptrStr[seg.data[i * sizeDataY + j]].minPotentialOnTheStructure = potential.data[i * sizeDataY + j];
					if (ptrStr[seg.data[i * sizeDataY + j]].maxPotentialOnTheStructure < potential.data[i * sizeDataY + j]) {
						ptrStr[seg.data[i * sizeDataY + j]].x = i;
						ptrStr[seg.data[i * sizeDataY + j]].y = j;
						ptrStr[seg.data[i * sizeDataY + j]].maxPotentialOnTheStructure = potential.data[i * sizeDataY + j];
					}
					for (ii = __max(i-1, 0); ii <= __min(i+1, sizeDataX-1); ii++)
						for (jj = __max(j-1, 0); jj <= __min(j+1, sizeDataY-1); jj++)
							if ((ii != i) || (jj != j))
								if (seg.data[ii * sizeDataY + jj] != 0)
									if (seg.data[i * sizeDataY + j] != seg.data[ii * sizeDataY + jj])
										if (potential.data[i * sizeDataY + j] > ptrStr[seg.data[i * sizeDataY + j]].maxPotentialWithNeighbor) {
											ptrStr[seg.data[i * sizeDataY + j]].maxPotentialWithNeighbor = potential.data[i * sizeDataY + j];
											ptrStr[seg.data[i * sizeDataY + j]].structureInContact = seg.data[ii * sizeDataY + jj];
										}
				}


		change = false;
		for (i = 1; i <= oldNbInfoStruct; i++) {
			if (ptrStr[i].structureInContact > 0)
				if ((ptrStr[i].locked == false) && (ptrStr[ptrStr[i].structureInContact].locked == false)) {
					ratio = (ptrStr[i].maxPotentialOnTheStructure - ptrStr[i].maxPotentialWithNeighbor) / (ptrStr[i].maxPotentialOnTheStructure - ptrStr[i].minPotentialOnTheStructure);

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

		label = 0;

		for (i = 1; i <= oldNbInfoStruct; i++)
			if (ptrStr[i].structureToBecaume == i) {
				label++;
				ptrStr[i].newIndex = label;
			}


		for (i = 0; i < sizeDataX; i++) for (j = 0; j < sizeDataY; j++) { if (seg.data[i * sizeDataY + j] != 0) { seg.data[i * sizeDataY + j] = (unsigned short)ptrStr[ptrStr[seg.data[i * sizeDataY + j]].structureToBecaume].newIndex; } }
	}

	seg = (unsigned short)0;

	label = 0;
	for (i = 1; i <= (long)nbInfoStruct; i++) {
		if (ptrStr[i].structureToBecaume == i) {
			label++;
			seg.data[ptrStr[i].x * sizeDataY + ptrStr[i].y] = (unsigned short)label;
		}
	}


	free(ptrStr);

	return (label - 1);
}

/********************************************************************************************************************
int Segment2D::getLocalMax(...) compute local maximum of the image an mark them in the image seg. 
USrtArray &image: input image
USrtArray &localMax: local maximum
float threshold: a point to be considered as a local maximum have to be > to threshold.
********************************************************************************************************************/
int Segment2D::getLocalMax(USrtArray& image, USrtArray& seg, float threshold)
{
	seg = (unsigned short)0;
	long i, j, ii, jj;
	long label = 1;
	bool err;
	long sizeDataY = image.size(1);
	long sizeDataX = image.size(0);


	for (i = 0; i < sizeDataX; i++)
		for (j = 0; j < sizeDataY; j++)
			if (image.data[i * sizeDataY + j] > threshold) {
				err = false;
				//test if it is a new maximum from the vicinity
				for (ii = __max(i-1, 0); (ii <= __min(i+1, sizeDataX-1)) && !err; ii++) for (jj = __max(j-1, 0); (jj <= __min(j+1, sizeDataY-1)) && !err; jj++) { if ((ii != i) || (jj != j)) if ((image.data[i * sizeDataY + j] < image.data[ii * sizeDataY + jj]) || (image.data[ii * sizeDataY + jj] <= threshold)) err = true; }
				if (!err) { //new maximum found. 
					unsigned short maxi = 0;
					err = false;
					for (ii = __max(i-1, 0); (ii <= __min(i+1, sizeDataX-1)) && !err; ii++)
						for (jj = __max(j-1, 0); (jj <= __min(j+1, sizeDataY-1)) && !err; jj++)
							if ((ii != i) || (jj != j))
								if (maxi < seg.data[ii * sizeDataY + jj]) {
									maxi = seg.data[ii * sizeDataY + jj];
									err = true;
								}
					if (maxi > 0) seg.data[i * sizeDataY + j] = maxi;
					else {
						seg.data[i * sizeDataY + j] = (unsigned short)label;
						label++;
					}
				}
			}

	return (label - 1);
}


/********************************************************************************************************************
int Segment2D::watershed(...) compute watershed segmentation insied the reagion upper the threshold value
USrtArray &potential: defines the potential function the makes the watershed segmented area growing
USrtArray &seg: it contains the result of the watershed segmentation. 
********************************************************************************************************************/
int Segment2D::watershed(USrtArray& potential, USrtArray& seg, float threshold)
{
	float cle;

	coor2D p;
	coor2D p1;

	long sizeDataY = potential.size(1);
	long sizeDataX = potential.size(0);

	Heap<coor2D> tas(sizeDataY * sizeDataX);


	// Initialisation des marqueurs
	for (p.y = 0; p.y < sizeDataY; p.y++)
		for (p.x = 0; p.x < sizeDataX; p.x++) {
			bool vois = false;

			if (seg.data[p.x * sizeDataY + p.y] != 0) {
				for (p1.x = p.x - 1; (p1.x <= p.x + 1) && !vois; p1.x++)
					for (p1.y = p.y - 1; (p1.y <= p.y + 1) && !vois; p1.y++)
						if ((p1.x != p.x) || (p1.y != p.y))
							if (p1.x >= 0)
								if (p1.x < potential.size(0)) if (p1.y >= 0) if (p1.y < potential.size(1)) if (seg.data[p1.x * sizeDataY + p1.y] == 0) vois = true;

				// On empile les points frontieres.
				if (vois) {
					cle = (float)(-potential.data[p1.x * sizeDataY + p1.y]);
					tas.Push(p, cle);
				}
			}
		}


	// Croissance prioritaire.
	int it = 0;
	while (!tas.Empty()) {
		it = it + 1;
		p = tas.Pop();
		for (p1.x = p.x - 1; p1.x <= p.x + 1; p1.x++)
			for (p1.y = p.y - 1; p1.y <= p.y + 1; p1.y++) {
				if ((p1.x != p.x) || (p1.y != p.y))
					if (p1.x >= 0)
						if (p1.x < potential.size(0))
							if (p1.y >= 0)
								if (p1.y < potential.size(1))
									if (potential.data[p1.x * sizeDataY + p1.y] > threshold)
										if (seg.data[p1.x * sizeDataY + p1.y] == 0) {
											// On empile ce nouveau avec une priorité legèrement supérieure
											// a sa vrai valeur de potentiel.
											cle = (float)(-potential.data[p1.x * sizeDataY + p1.y] + (float)it * EPSILON);
											tas.Push(p1, cle);
											seg.data[p1.x * sizeDataY + p1.y] = seg.data[p.x * sizeDataY + p.y];
										}
			}
	}

	return (nbInfoStruct);
}


#undef MAX_ELEMENTS
#undef EPSILON
