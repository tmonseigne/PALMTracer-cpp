


/*******************************************************************************************************************
								WaterShed segmentation (3D signal) 
									Institut Curie
									UMR - 144
									by Victor Racine
									2005 01 05
*******************************************************************************************************************/

#include "segmentation_3D.h"

extern FILE *debugStreamALS;
#define MAX_ELEMENTS 10000000
#define EPSILON 0.0000001

/********************************************************************************************************************
int Segment3D::filterStructureWS(...) filter the structure of the watershed. 
FltArray &potential: input image
USrtArray &localMax: local maximum
float thresholdRatio: value to distinguish if 2 strucutres in contact have to be considered as a single structure (
							oversegmentation)
********************************************************************************************************************/
int Segment3D::filterStructureWS(FltArray &potential, USrtArray &seg, float thresholdRatio){
	
	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "->-> Segment3D::filterStructureWS(FltArray 0x%p &potential(%d, %d, %d), USrtArray 0x%p &seg(%d, %d, %d), float %f threshold)  <-<-\n",		
		&potential, potential.size(0), potential.size(1), potential.size(2), &seg, seg.size(0), seg.size(1), seg.size(2), thresholdRatio);
	closeDebugStreamALS();
	#endif

	wsFusion3D* ptrStr;
	ptrStr=(wsFusion3D*)malloc((nbInfoStruct+1)*sizeof(wsFusion3D));
	long i,j,k,ii, jj,kk;

	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);

	long oldNbInfoStruct;
	double ratio;
	long label;
	bool change=true;
	while(change){

		oldNbInfoStruct=nbInfoStruct;

		for(i=1;i<=nbInfoStruct;i++){
			ptrStr[i].maxPotentialWithNeighbor=-1e10;
			ptrStr[i].minPotentialOnTheStructure=+1e10;
			ptrStr[i].maxPotentialOnTheStructure=-1e10;
			ptrStr[i].structureInContact=0;
			ptrStr[i].kept=true;
			ptrStr[i].locked=false;
			ptrStr[i].structureToBecaume=i;
			ptrStr[i].newIndex=0;
		}

		
		for(i=0;i<sizeDataX;i++)
		for(j=0;j<sizeDataY;j++)
		for(k=0;k<sizeDataZ;k++)
			if(seg.data[(i*sizeDataY+j)*sizeDataZ+k]>0){
				if(ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].minPotentialOnTheStructure>potential.data[(i*sizeDataY+j)*sizeDataZ+k])
					ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].minPotentialOnTheStructure=potential.data[(i*sizeDataY+j)*sizeDataZ+k];
				if(ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].maxPotentialOnTheStructure<potential.data[(i*sizeDataY+j)*sizeDataZ+k]){
					ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].x=i;
					ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].y=j;
					ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].z=k;
					ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].maxPotentialOnTheStructure=potential.data[(i*sizeDataY+j)*sizeDataZ+k];
				}
				for(ii=__max(i-1,0);ii<=__min(i+1,sizeDataX-1);ii++)
				for(jj=__max(j-1,0);jj<=__min(j+1,sizeDataY-1);jj++)
				for(kk=__max(k-1,0);kk<=__min(k+1,sizeDataZ-1);kk++)
					if((ii!=i)||(jj!=j)||(kk!=k))
						if(seg.data[(ii*sizeDataY+jj)*sizeDataZ+kk]!=0)
							if(seg.data[(i*sizeDataY+j)*sizeDataZ+k]!=seg.data[(ii*sizeDataY+jj)*sizeDataZ+kk])
								if(potential.data[(i*sizeDataY+j)*sizeDataZ+k]>ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].maxPotentialWithNeighbor){
									ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].maxPotentialWithNeighbor=potential.data[(i*sizeDataY+j)*sizeDataZ+k];
									ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].structureInContact=seg.data[(ii*sizeDataY+jj)*sizeDataZ+kk];
								}
			} 

	
		change=false;
		for(i=1;i<=oldNbInfoStruct;i++){
			if(ptrStr[i].structureInContact>0)
				if((ptrStr[i].locked==false)&&(ptrStr[ptrStr[i].structureInContact].locked==false)){
					ratio=(ptrStr[i].maxPotentialOnTheStructure-ptrStr[i].maxPotentialWithNeighbor)/
						(ptrStr[i].maxPotentialOnTheStructure-ptrStr[i].minPotentialOnTheStructure);

					if(ratio<thresholdRatio)
						if(ptrStr[i].maxPotentialOnTheStructure<ptrStr[ptrStr[i].structureInContact].maxPotentialOnTheStructure){
							ptrStr[i].structureToBecaume=ptrStr[i].structureInContact;
							nbInfoStruct--;
							ptrStr[ptrStr[i].structureInContact].locked=true;
							ptrStr[i].locked=true;
							change=true;

						}
				}
		}

		label=0;

		for(i=1;i<=oldNbInfoStruct;i++)
			if(ptrStr[i].structureToBecaume==i){
				label++;
				ptrStr[i].newIndex=label;
			}


		for(i=0;i<sizeDataX*sizeDataY*sizeDataZ;i++){
			if(seg.data[i]!=0){
				seg.data[i]=ptrStr[ptrStr[seg.data[i]].structureToBecaume].newIndex;
			}
		}

	}

	seg=(unsigned short)0;

	label=0;
	for(i=1;i<=nbInfoStruct;i++){
		if(ptrStr[i].structureToBecaume==i){
			label++;
			seg.data[(ptrStr[i].x*sizeDataY+ptrStr[i].y)*sizeDataZ+ptrStr[i].z]=label;
		}
	}



	free(ptrStr);

	return(label-1);
}

/********************************************************************************************************************
int Segment3D::getLocalMax(...) compute local maximum of the image an mark them in the image seg. 
FltArray &image: input image
USrtArray &localMax: local maximum
float threshold: a point to be considered as a local maximum have to be > to threshold.
********************************************************************************************************************/
int Segment3D::getLocalMax(FltArray &image, USrtArray &seg, float threshold){
	

	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "->-> Segment3D::getLocalMax(FltArray 0x%p &potential(%d, %d, %d), USrtArray 0x%p &seg(%d, %d, %d), float %f threshold)  <-<-\n",
		&image, image.size(0), image.size(1), image.size(2), &seg, seg.size(0), seg.size(1), seg.size(2),threshold);
	closeDebugStreamALS();
	#endif

	seg=(unsigned short)0;
	long i,j,k, ii,jj,kk;
	long label=1;
	bool err;
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);

	for(i=0;i<sizeDataX;i++)
	for(j=0;j<sizeDataY;j++)
	for(k=0;k<sizeDataZ;k++)
		if(image.data[(i*sizeDataY+j)*sizeDataZ+k]>threshold){
			err=false;
			//test if it is a new maximum from the vicinity
			for(ii=__max(i-1,0);(ii<=__min(i+1,sizeDataX-1))&&!err;ii++)
			for(jj=__max(j-1,0);(jj<=__min(j+1,sizeDataY-1))&&!err;jj++)
			for(kk=__max(k-1,0);(kk<=__min(k+1,sizeDataZ-1))&&!err;kk++){
				if((ii!=i)||(jj!=j)||(kk!=k))
					if(image.data[(i*sizeDataY+j)*sizeDataZ+k]<image.data[(ii*sizeDataY+jj)*sizeDataZ+kk])
						err=true;
			}
			if(!err){ //new maximum found
				unsigned short maxi=0;
				err=false;
				for(ii=__max(i-1,0);(ii<=__min(i+1,sizeDataX-1))&&!err;ii++)
				for(jj=__max(j-1,0);(jj<=__min(j+1,sizeDataY-1))&&!err;jj++)
				for(kk=__max(k-1,0);(kk<=__min(k+1,sizeDataZ-1))&&!err;kk++)
					if((ii!=i)||(jj!=j)||(kk!=k))
						if(maxi<seg.data[(ii*sizeDataY+jj)*sizeDataZ+kk]){
							maxi=seg.data[(ii*sizeDataY+jj)*sizeDataZ+kk];
							err=true;
						}
				if(maxi>0)
					seg.data[(i*sizeDataY+j)*sizeDataZ+k]=maxi;
				else{
					seg.data[(i*sizeDataY+j)*sizeDataZ+k]=(unsigned short)label;
					label++;
				}
			}

		}
	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "->-> Segment3D::getLocalMax End  <-<-\n");
	closeDebugStreamALS();
	#endif

	return(label-1);

}






/********************************************************************************************************************
int Segment3D::watershed(...) compute watershed segmentation insied the reagion upper the threshold value
FltArray &potential: defines the potential function the makes the watershed segmented area growing
USrtArray &seg: it contains the result of the watershed segmentation. 
********************************************************************************************************************/
int Segment3D::watershed(FltArray &potential, USrtArray &seg, float threshold){

	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "->-> Segment3D::watershed(FltArray 0x%p &potential(%d, %d, %d), USrtArray 0x%p &seg(%d, %d, %d), float %f threshold)  <-<-\n",
		&potential, potential.size(0), potential.size(1), potential.size(2), &seg, seg.size(0), seg.size(1), seg.size(2),threshold);
	closeDebugStreamALS();
	#endif

	float cle;

	coor3D p;
	coor3D p1;

	
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);

	Heap <coor3D> tas(sizeDataX*sizeDataY*sizeDataZ);

	// Initialisation des marqueurs
	for (p.z=0;p.z<sizeDataZ;p.z++)
	for (p.y=0;p.y<sizeDataY;p.y++)
	for (p.x=0;p.x<sizeDataX;p.x++){

		bool vois=false;

		if (seg.data[(p.x*sizeDataY+p.y)*sizeDataZ+p.z]!=0){

			for (p1.x=p.x-1;(p1.x<=p.x+1)&&!vois;p1.x++)
			for (p1.y=p.y-1;(p1.y<=p.y+1)&&!vois;p1.y++)
			for (p1.z=p.z-1;(p1.z<=p.z+1)&&!vois;p1.z++)
				if((p1.x!=p.x)||(p1.y!=p.y)||(p1.z!=p.z))
				if(p1.x>=0)	if(p1.x<sizeDataX)
				if(p1.y>=0)	if(p1.y<sizeDataY)
				if(p1.z>=0)	if(p1.z<sizeDataZ)
				if(seg.data[p1.x*sizeDataY+p1.y]==0)
					vois=true;

			// On empile les points frontieres.
			if (vois) {
				cle=-potential.data[(p1.x*sizeDataY+p1.y)*sizeDataZ+p1.z];
				tas.Push(p,cle);
			}
		}
	}


	// Croissance prioritaire.
	int it=0;
	while (!tas.Empty()) {
		it=it+1;
		p=tas.Pop();
		for (p1.x=p.x-1;p1.x<=p.x+1;p1.x++)
		for (p1.y=p.y-1;p1.y<=p.y+1;p1.y++)
		for (p1.z=p.z-1;p1.z<=p.z+1;p1.z++){
			if((p1.x!=p.x)||(p1.y!=p.y)||(p1.z!=p.z))
			if(p1.x>=0)	if(p1.x<sizeDataX)
			if(p1.y>=0)	if(p1.y<sizeDataY)
			if(p1.z>=0)	if(p1.z<sizeDataZ)
			if(potential.data[(p1.x*sizeDataY+p1.y)*sizeDataZ+p1.z]>threshold)
			if(seg.data[(p1.x*sizeDataY+p1.y)*sizeDataZ+p1.z]==0){
			// On empile ce nouveau avec une priorité legèrement supérieure
			// a sa vrai valeur de potentiel.
				cle=-potential.data[(p1.x*sizeDataY+p1.y)*sizeDataZ+p1.z]+(float)it*EPSILON;
				tas.Push(p1,cle);
				seg.data[(p1.x*sizeDataY+p1.y)*sizeDataZ+p1.z]=seg.data[(p.x*sizeDataY+p.y)*sizeDataZ+p.z];
			}
		}
	}


	return(nbInfoStruct);

}



/********************************************************************************************************************
int Segment3D::watershedIntelligent(...) compute watershed segmentation insied the reagion upper the threshold value. 
It tries to minimize the oversegmentation.
FltArray &potential: defines the potential function the makes the watershed segmented area growing
USrtArray &seg: it contains the result of the watershed segmentation. 
float threshold: define the threshold value to distinguish the segmented structures
float thresholdRatio: value 0<= X <=1. Makes a fusion to minimize the oversegmentation
						if X=0 simple watershed
						if X=1 no watershed
return value: number of detected structures
********************************************************************************************************************/
int Segment3D::watershedIntelligent(FltArray &potential, USrtArray &seg, float threshold, float thresholdRatio){

	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "->-> Segment3D::watershedIntelligent(FltArray 0x%p &potential(%d, %d, %d), USrtArray 0x%p &seg(%d, %d, %d), float %f threshold, float %f thresholdRatio)  <-<-\n",
		&potential, potential.size(0), potential.size(1), potential.size(2), &seg, seg.size(0), seg.size(1), seg.size(2),threshold ,thresholdRatio);
	closeDebugStreamALS();
	#endif
	int ret;

	//first, find the local maximum to set the different germinates
	if((nbInfoStruct=getLocalMax(potential, seg, threshold))<0)
		return(nbInfoStruct);

	
	if((ret=watershed(potential, seg, threshold))<0) //first watershed
		return(ret);

	if(thresholdRatio==0){
		mixIndex(seg);
		return(nbInfoStruct);
	}



	if((ret=filterStructureWS(potential, seg, thresholdRatio))<0) //filter the germinates to minimize the oversegmentation
		return(ret);

	
	if((ret=watershed(potential, seg, threshold))<0) //first watershed
		return(ret);

	mixIndex(seg);

	return(nbInfoStruct);

}


/********************************************************************************************************************
int Segment3D::filterStructureWS(...) filter the structure of the watershed. 
USrtArray &potential: input image
USrtArray &localMax: local maximum
float thresholdRatio: value to distinguish if 2 strucutres in contact have to be considered as a single structure (
							oversegmentation)
********************************************************************************************************************/
int Segment3D::filterStructureWS(USrtArray &potential, USrtArray &seg, float thresholdRatio){
	
	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "->-> Segment3D::filterStructureWS(USrtArray 0x%p &potential(%d, %d, %d), USrtArray 0x%p &seg(%d, %d, %d), float %f threshold)  <-<-\n",		
		&potential, potential.size(0), potential.size(1), potential.size(2), &seg, seg.size(0), seg.size(1), seg.size(2), thresholdRatio);
	closeDebugStreamALS();
	#endif

	wsFusion3D* ptrStr;
	ptrStr=(wsFusion3D*)malloc((nbInfoStruct+1)*sizeof(wsFusion3D));
	long i,j,k,ii, jj,kk;

	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);

	long oldNbInfoStruct;
	double ratio;
	long label;
	bool change=true;
	while(change){

		oldNbInfoStruct=nbInfoStruct;

		for(i=1;i<=nbInfoStruct;i++){
			ptrStr[i].maxPotentialWithNeighbor=-1e10;
			ptrStr[i].minPotentialOnTheStructure=+1e10;
			ptrStr[i].maxPotentialOnTheStructure=-1e10;
			ptrStr[i].structureInContact=0;
			ptrStr[i].kept=true;
			ptrStr[i].locked=false;
			ptrStr[i].structureToBecaume=i;
			ptrStr[i].newIndex=0;
		}

		
		for(i=0;i<sizeDataX;i++)
		for(j=0;j<sizeDataY;j++)
		for(k=0;k<sizeDataZ;k++)
			if(seg.data[(i*sizeDataY+j)*sizeDataZ+k]>0){
				if(ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].minPotentialOnTheStructure>potential.data[(i*sizeDataY+j)*sizeDataZ+k])
					ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].minPotentialOnTheStructure=potential.data[(i*sizeDataY+j)*sizeDataZ+k];
				if(ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].maxPotentialOnTheStructure<potential.data[(i*sizeDataY+j)*sizeDataZ+k]){
					ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].x=i;
					ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].y=j;
					ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].z=k;
					ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].maxPotentialOnTheStructure=potential.data[(i*sizeDataY+j)*sizeDataZ+k];
				}
				for(ii=__max(i-1,0);ii<=__min(i+1,sizeDataX-1);ii++)
				for(jj=__max(j-1,0);jj<=__min(j+1,sizeDataY-1);jj++)
				for(kk=__max(k-1,0);kk<=__min(k+1,sizeDataZ-1);kk++)
					if((ii!=i)||(jj!=j)||(kk!=k))
						if(seg.data[(ii*sizeDataY+jj)*sizeDataZ+kk]!=0)
							if(seg.data[(i*sizeDataY+j)*sizeDataZ+k]!=seg.data[(ii*sizeDataY+jj)*sizeDataZ+kk])
								if(potential.data[(i*sizeDataY+j)*sizeDataZ+k]>ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].maxPotentialWithNeighbor){
									ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].maxPotentialWithNeighbor=potential.data[(i*sizeDataY+j)*sizeDataZ+k];
									ptrStr[seg.data[(i*sizeDataY+j)*sizeDataZ+k]].structureInContact=seg.data[(ii*sizeDataY+jj)*sizeDataZ+kk];
								}
			} 

	
		change=false;
		for(i=1;i<=oldNbInfoStruct;i++){
			if(ptrStr[i].structureInContact>0)
				if((ptrStr[i].locked==false)&&(ptrStr[ptrStr[i].structureInContact].locked==false)){
					ratio=(ptrStr[i].maxPotentialOnTheStructure-ptrStr[i].maxPotentialWithNeighbor)/
						(ptrStr[i].maxPotentialOnTheStructure-ptrStr[i].minPotentialOnTheStructure);

					if(ratio<thresholdRatio)
						if(ptrStr[i].maxPotentialOnTheStructure<ptrStr[ptrStr[i].structureInContact].maxPotentialOnTheStructure){
							ptrStr[i].structureToBecaume=ptrStr[i].structureInContact;
							nbInfoStruct--;
							ptrStr[ptrStr[i].structureInContact].locked=true;
							ptrStr[i].locked=true;
							change=true;

						}
				}
		}

		label=0;

		for(i=1;i<=oldNbInfoStruct;i++)
			if(ptrStr[i].structureToBecaume==i){
				label++;
				ptrStr[i].newIndex=label;
			}


		for(i=0;i<sizeDataX*sizeDataY*sizeDataZ;i++){
			if(seg.data[i]!=0){
				seg.data[i]=ptrStr[ptrStr[seg.data[i]].structureToBecaume].newIndex;
			}
		}

	}

	seg=(unsigned short)0;

	label=0;
	for(i=1;i<=nbInfoStruct;i++){
		if(ptrStr[i].structureToBecaume==i){
			label++;
			seg.data[(ptrStr[i].x*sizeDataY+ptrStr[i].y)*sizeDataZ+ptrStr[i].z]=label;
		}
	}



	free(ptrStr);

	return(label-1);
}

/********************************************************************************************************************
int Segment3D::getLocalMax(...) compute local maximum of the image an mark them in the image seg. 
USrtArray &image: input image
USrtArray &localMax: local maximum
float threshold: a point to be considered as a local maximum have to be > to threshold.
********************************************************************************************************************/
int Segment3D::getLocalMax(USrtArray &image, USrtArray &seg, float threshold){
	

	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "->-> Segment3D::getLocalMax(USrtArray 0x%p &potential(%d, %d, %d), USrtArray 0x%p &seg(%d, %d, %d), float %f threshold)  <-<-\n",
		&image, image.size(0), image.size(1), image.size(2), &seg, seg.size(0), seg.size(1), seg.size(2),threshold);
	closeDebugStreamALS();
	#endif

	seg=(unsigned short)0;
	long i,j,k, ii,jj,kk;
	long label=1;
	bool err;
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);

	for(i=0;i<sizeDataX;i++)
	for(j=0;j<sizeDataY;j++)
	for(k=0;k<sizeDataZ;k++)
		if(image.data[(i*sizeDataY+j)*sizeDataZ+k]>threshold){
			err=false;
			//test if it is a new maximum from the vicinity
			for(ii=__max(i-1,0);(ii<=__min(i+1,sizeDataX-1))&&!err;ii++)
			for(jj=__max(j-1,0);(jj<=__min(j+1,sizeDataY-1))&&!err;jj++)
			for(kk=__max(k-1,0);(kk<=__min(k+1,sizeDataZ-1))&&!err;kk++){
				if((ii!=i)||(jj!=j)||(kk!=k))
					if(image.data[(i*sizeDataY+j)*sizeDataZ+k]<image.data[(ii*sizeDataY+jj)*sizeDataZ+kk])
						err=true;
			}
			if(!err){ //new maximum found
				unsigned short maxi=0;
				err=false;
				for(ii=__max(i-1,0);(ii<=__min(i+1,sizeDataX-1))&&!err;ii++)
				for(jj=__max(j-1,0);(jj<=__min(j+1,sizeDataY-1))&&!err;jj++)
				for(kk=__max(k-1,0);(kk<=__min(k+1,sizeDataZ-1))&&!err;kk++)
					if((ii!=i)||(jj!=j)||(kk!=k))
						if(maxi<seg.data[(ii*sizeDataY+jj)*sizeDataZ+kk]){
							maxi=seg.data[(ii*sizeDataY+jj)*sizeDataZ+kk];
							err=true;
						}
				if(maxi>0)
					seg.data[(i*sizeDataY+j)*sizeDataZ+k]=maxi;
				else{
					seg.data[(i*sizeDataY+j)*sizeDataZ+k]=(unsigned short)label;
					label++;
				}
			}

		}
	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "->-> Segment3D::getLocalMax End  <-<-\n");
	closeDebugStreamALS();
	#endif

	return(label-1);

}






/********************************************************************************************************************
int Segment3D::watershed(...) compute watershed segmentation insied the reagion upper the threshold value
USrtArray &potential: defines the potential function the makes the watershed segmented area growing
USrtArray &seg: it contains the result of the watershed segmentation. 
********************************************************************************************************************/
int Segment3D::watershed(USrtArray &potential, USrtArray &seg, float threshold){

	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "->-> Segment3D::watershed(USrtArray 0x%p &potential(%d, %d, %d), USrtArray 0x%p &seg(%d, %d, %d), float %f threshold)  <-<-\n",
		&potential, potential.size(0), potential.size(1), potential.size(2), &seg, seg.size(0), seg.size(1), seg.size(2),threshold);
	closeDebugStreamALS();
	#endif

	float cle;

	coor3D p;
	coor3D p1;

	

	
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);
	
	Heap <coor3D> tas(sizeDataX*sizeDataY*sizeDataZ);

	// Initialisation des marqueurs
	for (p.z=0;p.z<sizeDataZ;p.z++)
	for (p.y=0;p.y<sizeDataY;p.y++)
	for (p.x=0;p.x<sizeDataX;p.x++){

		bool vois=false;

		if (seg.data[(p.x*sizeDataY+p.y)*sizeDataZ+p.z]!=0){

			for (p1.x=p.x-1;(p1.x<=p.x+1)&&!vois;p1.x++)
			for (p1.y=p.y-1;(p1.y<=p.y+1)&&!vois;p1.y++)
			for (p1.z=p.z-1;(p1.z<=p.z+1)&&!vois;p1.z++)
				if((p1.x!=p.x)||(p1.y!=p.y)||(p1.z!=p.z))
				if(p1.x>=0)	if(p1.x<sizeDataX)
				if(p1.y>=0)	if(p1.y<sizeDataY)
				if(p1.z>=0)	if(p1.z<sizeDataZ)
				if(seg.data[p1.x*sizeDataY+p1.y]==0)
					vois=true;

			// On empile les points frontieres.
			if (vois) {
				cle=-potential.data[(p1.x*sizeDataY+p1.y)*sizeDataZ+p1.z];
				tas.Push(p,cle);
			}
		}
	}


	// Croissance prioritaire.
	int it=0;
	while (!tas.Empty()) {
		it=it+1;
		p=tas.Pop();
		for (p1.x=p.x-1;p1.x<=p.x+1;p1.x++)
		for (p1.y=p.y-1;p1.y<=p.y+1;p1.y++)
		for (p1.z=p.z-1;p1.z<=p.z+1;p1.z++){
			if((p1.x!=p.x)||(p1.y!=p.y)||(p1.z!=p.z))
			if(p1.x>=0)	if(p1.x<sizeDataX)
			if(p1.y>=0)	if(p1.y<sizeDataY)
			if(p1.z>=0)	if(p1.z<sizeDataZ)
			if(potential.data[(p1.x*sizeDataY+p1.y)*sizeDataZ+p1.z]>threshold)
			if(seg.data[(p1.x*sizeDataY+p1.y)*sizeDataZ+p1.z]==0){
			// On empile ce nouveau avec une priorité legèrement supérieure
			// a sa vrai valeur de potentiel.
				cle=-potential.data[(p1.x*sizeDataY+p1.y)*sizeDataZ+p1.z]+(float)it*EPSILON;
				tas.Push(p1,cle);
				seg.data[(p1.x*sizeDataY+p1.y)*sizeDataZ+p1.z]=seg.data[(p.x*sizeDataY+p.y)*sizeDataZ+p.z];
			}
		}
	}


	return(nbInfoStruct);

}



/********************************************************************************************************************
int Segment3D::watershedIntelligent(...) compute watershed segmentation insied the reagion upper the threshold value. 
It tries to minimize the oversegmentation.
USrtArray &potential: defines the potential function the makes the watershed segmented area growing
USrtArray &seg: it contains the result of the watershed segmentation. 
float threshold: define the threshold value to distinguish the segmented structures
float thresholdRatio: value 0<= X <=1. Makes a fusion to minimize the oversegmentation
						if X=0 simple watershed
						if X=1 no watershed
return value: number of detected structures
********************************************************************************************************************/
int Segment3D::watershedIntelligent(USrtArray &potential, USrtArray &seg, float threshold, float thresholdRatio){

	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "->-> Segment3D::watershedIntelligent(USrtArray 0x%p &potential(%d, %d, %d), USrtArray 0x%p &seg(%d, %d, %d), float %f threshold, float %f thresholdRatio)  <-<-\n",
		&potential, potential.size(0), potential.size(1), potential.size(2), &seg, seg.size(0), seg.size(1), seg.size(2),threshold ,thresholdRatio);
	closeDebugStreamALS();
	#endif
	int ret;

	//first, find the local maximum to set the different germinates
	if((nbInfoStruct=getLocalMax(potential, seg, threshold))<0)
		return(nbInfoStruct);

	
	if((ret=watershed(potential, seg, threshold))<0) //first watershed
		return(ret);

	if(thresholdRatio==0){
		mixIndex(seg);
		return(nbInfoStruct);
	}



	if((ret=filterStructureWS(potential, seg, thresholdRatio))<0) //filter the germinates to minimize the oversegmentation
		return(ret);

	
	if((ret=watershed(potential, seg, threshold))<0) //first watershed
		return(ret);
	mixIndex(seg);


	return(nbInfoStruct);

}


#undef MAX_ELEMENTS
#undef EPSILON