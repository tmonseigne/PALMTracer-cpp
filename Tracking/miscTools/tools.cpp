
#include "stdafx.h"
#include "tools.h"

#define min(a,b) ((((a))<((b)))?((a)):((b)))
#define max(a,b) ((((a))>((b)))?((a)):((b)))

inline long turnClockwiseX(const long x, const long y){
	if(x!=0)
		return(0);
	else
		if(y==1)
			return(-1);
		else
			return(1);
}

inline long turnClockwiseY(const long x, const long y){
	if(y!=0)
		return(0);
	else
		if(x==1)
			return(1);
		else
			return(-1);
}

inline void turnClockwise(long & x, long & y){
	if(y==0){
		if(x==1){
			x=0;
			y=1;
		}else{
			x=0;
			y=-1;
		}
	}else if(y==1){
		x=-1;
		y=0;
	}else{
		x=1;
		y=0;
	}
}
inline void turnAntiClockwise(long & x, long & y){
	if(y==0){
		if(x==1){
			x=0;
			y=-1;
		}else{
			x=0;
			y=1;
		}
	}else if(y==1){
		x=1;
		y=0;
	}else{
		x=-1;
		y=0;
	}
}

#define ISOUT(x,y) (isOutside((mainSegmentation), (val), (x), (y), (xSize), (ySize)))

bool isOutside(USrtArray &mainSegmentation, const long structure, const long x, const long y, const long xSize, const long ySize){
	if(x<0)
		return(true);
	if(y<0)
		return(true);
	if(x>=xSize)
		return(true);
	if(y>=ySize)
		return(true);
	return(mainSegmentation.data[x*ySize+y]!=structure);
}

int printRegions(char* filename,  USrtArray &mainSegmentation, long nbStructures){
	FILE* file;
	int ret;
	unsigned short val;
	long posX, posY;
	long dirX, dirY;
	long dirX0, dirY0;
	long nbPoints;
	long xSize=mainSegmentation.size(0);
	long ySize=mainSegmentation.size(1);

	USrtArray usedBorder;
	ret=usedBorder.alloc(xSize, ySize); if(ret<0) return(ret);
	usedBorder=(unsigned short)0;


	SrtArray isAlreadyWritten(nbStructures+1,3);
	isAlreadyWritten=(short)-1;
	ret=fopen_s(&file, filename, "w");
	if(ret<0) return(ret);
	for(long x=0;x<xSize;x++)
	for(long y=0;y<ySize;y++){
		val=mainSegmentation.data[x*ySize+y];
		if(val>0)
			if(usedBorder.data[x*ySize+y]==0){ // new region to write
				
				if(x==0){
					dirX=-1;
					dirY=0;
				}else if(x==xSize){
					dirX=1;
					dirY=0;
				}else if(y==0){
					dirX=0;
					dirY=-1;
				}else if(y==ySize){
					dirX=0;
					dirY=1;
				}else if(val!=mainSegmentation.data[(x-1)*ySize+y]){
					dirX=-1;
					dirY=0;
				}else if(val!=mainSegmentation.data[(x)*ySize+y-1]){
					dirX=0;
					dirY=-1;
				
				}else{
					continue;
					//return(ERR_ALGORITHMIC); // algo error!!!
				}
				isAlreadyWritten(val,0)=1;
				isAlreadyWritten(val,1)=(unsigned short)x;
				isAlreadyWritten(val,2)=(unsigned short)y;
				posX=x;
				posY=y;
				dirX0=dirX;
				dirY0=dirY;
				usedBorder.data[posX*ySize+posY]=1;
				nbPoints=0;
				do{ // first round of the structure to count the points and update isAlreadyWritten(val,[1,2])
					if(ISOUT(posX+turnClockwiseX(dirX,dirY), posY+turnClockwiseY(dirX,dirY)))
						turnClockwise(dirX,dirY);
					else if(ISOUT(posX+turnClockwiseX(dirX,dirY)+dirX,posY+turnClockwiseY(dirX,dirY)+dirY)){
						posX+=turnClockwiseX(dirX,dirY);
						posY+=turnClockwiseY(dirX,dirY);
						nbPoints++;
						if(posX<isAlreadyWritten(val,1))
							isAlreadyWritten(val,1)=(unsigned short)posX;
						if(posY<isAlreadyWritten(val,2))
							isAlreadyWritten(val,2)=(unsigned short)posY;
					}else{
						posX+=turnClockwiseX(dirX,dirY)+dirX;
						posY+=turnClockwiseY(dirX,dirY)+dirY;
						nbPoints++;
						turnAntiClockwise(dirX,dirY);
						if(posX<isAlreadyWritten(val,1))
							isAlreadyWritten(val,1)=(unsigned short)posX;
						if(posY<isAlreadyWritten(val,2))
							isAlreadyWritten(val,2)=(unsigned short)posY;
					}
				}while((posX!=x)||(posY!=y)||(dirX!=dirX0)||(dirY!=dirY0));
				if(nbPoints>1){
					ret=fwprintf_s(file, L"0 3, 1 65535, 2 %d %d, 3 0 0, 4 0, 5 1, 6 %d", isAlreadyWritten(val,2), isAlreadyWritten(val,1), nbPoints);
					do{ // second round of the structure to write the points in the textfile
						if(ISOUT(posX+turnClockwiseX(dirX,dirY), posY+turnClockwiseY(dirX,dirY)))
							turnClockwise(dirX,dirY);
						else if(ISOUT(posX+turnClockwiseX(dirX,dirY)+dirX,posY+turnClockwiseY(dirX,dirY)+dirY)){
							posX+=turnClockwiseX(dirX,dirY);
							posY+=turnClockwiseY(dirX,dirY);
							usedBorder.data[posX*ySize+posY]=1;
							ret=fwprintf_s(file, L" %d %d", posY, posX);
						}else{
							posX+=turnClockwiseX(dirX,dirY)+dirX;
							posY+=turnClockwiseY(dirX,dirY)+dirY;
							turnAntiClockwise(dirX,dirY);
							usedBorder.data[posX*ySize+posY]=1;
							ret=fwprintf_s(file, L" %d %d", posY, posX);
						}
					}while((posX!=x)||(posY!=y)||(dirX!=dirX0)||(dirY!=dirY0));
					ret=fwprintf_s(file, L", 7 %d\n", val);
				}
			}
	}
	fclose(file);
	return(RET_OK);
}

int computeColocalization(USrtArray **allObjects, USrtArray &mainSeg, HCSData	*results, long nbStructures, long nbChannels){
	for(int ind=1;ind<nbStructures;ind++)
		for(int i=0;i<nbChannels;i++)
		for(int j=0;j<nbChannels;j++)
			results[ind].chColocalization[i][j]=0;

	for(int i=0;i<nbChannels;i++){
		for(int j=i+1;j<nbChannels;j++){
			for(int x=0;x<(long)allObjects[0]->nbPoints;x++){
				if(allObjects[i]->data[x]!=0)
					if(allObjects[j]->data[x]!=0)
						if(mainSeg.data[x]!=0)
							results[allObjects[i]->data[x]].chColocalization[i][j]++;
			}
		}
	}
	for(int ind=1;ind<nbStructures;ind++){
		for(int i=1;i<nbChannels;i++){
			for(int j=0;j<i;j++){
				results[ind].chColocalization[i][j]=results[ind].chColocalization[j][i];
			}
		}
	}
	for(int ind=1;ind<nbStructures;ind++){
		for(int i=0;i<nbChannels;i++)
		for(int j=0;j<nbChannels;j++)
			if(i!=j){
				if(results[ind].chObjects[i].surf>0)
					results[ind].chColocalization[i][j]/=results[ind].chObjects[i].surf;
			}
	}
	return(RET_OK);
}
int printListFeatures(wchar_t* generalOutputFilename){
	FILE* file;
	int ret;
	int cpt=0;
	ret=_wfopen_s(&file, generalOutputFilename, L"w");
	if(ret<0) return(ret);
	ret=fwprintf_s(file, L"%d\n", sizeof(HCSData)/sizeof(double));
	if(ret<0) return(ret);
	
	if((ret=fwprintf_s(file, L"index\n")) <0) return(ret); cpt++;
	if((ret=fwprintf_s(file, L"selected\n")) <0) return(ret); cpt++;
	if((ret=fwprintf_s(file, L"seg.nbPixels\n")) <0) return(ret); cpt++;
	if((ret=fwprintf_s(file, L"seg.center.x\n")) <0) return(ret); cpt++;
	if((ret=fwprintf_s(file, L"seg.center.y\n")) <0) return(ret); cpt++;
	if((ret=fwprintf_s(file, L"seg.perimeter\n")) <0) return(ret); cpt++;
	if((ret=fwprintf_s(file, L"seg.perimeterWithoutNeighbor\n")) <0) return(ret); cpt++;
	if((ret=fwprintf_s(file, L"seg.boundingBox.xMin\n")) <0) return(ret); cpt++;
	if((ret=fwprintf_s(file, L"seg.boundingBox.xMax\n")) <0) return(ret); cpt++;
	if((ret=fwprintf_s(file, L"seg.boundingBox.yMin\n")) <0) return(ret); cpt++;
	if((ret=fwprintf_s(file, L"seg.boundingBox.yMax\n")) <0) return(ret); cpt++;
	for(int i=1;i<=4;i++){
		if((ret=fwprintf_s(file, L"ch[%d].Int.valBackground\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Int.valMax\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Int.valMin\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Int.mean\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Int.totalIntensity\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Int.stdev\n", i)) <0) return(ret); cpt++;
		for(int z=0;z<36;z++){
			if((ret=fwprintf_s(file, L"ch[%d].Int.zernikeMoments[%d]\n", i, z)) <0) return(ret); cpt++;
		}
	}
	for(int i=1;i<=4;i++){
		if((ret=fwprintf_s(file, L"ch[%d].Obj.mean\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Obj.valMorpho\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Obj.peri\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Obj.surf\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Obj.distFromObjectCenterToSegCenter\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Obj.distFromObjectCenterToSegCenterSurf\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Obj.sigma.x\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Obj.sigma.y\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Obj.angle\n", i)) <0) return(ret); cpt++;
		if((ret=fwprintf_s(file, L"ch[%d].Obj.nbSubstructures\n", i)) <0) return(ret); cpt++;
	}
	for(int i=1;i<=4;i++){
		for(int j=1;j<=4;j++){
			if((ret=fwprintf_s(file, L"coloc[%d][%d]\n", i,j)) <0) return(ret); cpt++;
		}
	}
	
	fclose(file);
	if(cpt!=sizeof(HCSData)/sizeof(double)){
		printf("Error in the file of features (%S). sizeof(HCSData)/sizeof(double)=%d, cpt=%d\n", generalOutputFilename, sizeof(HCSData)/sizeof(double), cpt);
		return(-2);
	}

	

	return(RET_OK);
}
int printResults(char* filename, HCSData	*results, long nbStructures, long newNbStructs){
	int ret;
	int fh;
	long val;
	ret=_sopen_s(&fh, filename, _O_CREAT | _O_BINARY | _O_WRONLY, _SH_DENYRW, _S_IREAD | _S_IWRITE);
	if(ret!=0) return(ERR_OPEN_STAT);
	val=0;
	ret=_write(fh, &val, sizeof(long));
	if(ret!=sizeof(long)){
		printf("Error writing the first number in the file %s\n", filename);
		return(ERR_WRITE_STAT);
	}
	ret=_write(fh, &newNbStructs, sizeof(long));
	if(ret!=sizeof(long)){
		printf("Error writing the second number in the file %s\n", filename);
		return(ERR_WRITE_STAT);
	}
	for(long ind=1; ind<=nbStructures;ind++){
		if(results[ind].selected==1){
			ret=_write(fh, results+ind, sizeof(HCSData));
			if(ret!=sizeof(HCSData)){
				printf("Error writing the data field (%d/%d) in the file %s\n", ind , newNbStructs, filename);
				return(ERR_WRITE_STAT);
			}
		}
	}

	ret=_lseek(fh, 0, SEEK_SET);
	if(ret!=0){
		printf("Error writing the validation of the data file (%s)\n", filename);
		return(ERR_WRITE_STAT);
	}
 
	val=1;
	ret=_write(fh, &val, sizeof(long));
	if(ret!=sizeof(long)){
		printf("Error writing the validation number in the file %s\n", filename);
		return(ERR_WRITE_STAT);
	}

	_close(fh);

	return(RET_OK);
}


int printResultsASCII(char* filename, HCSData	*results, long nbStructures, long newNbStructs){
	int ret;

	//long val;
	FILE* file;
	ret=fopen_s(&file, filename, "wt");
	if(ret!=0) return(ERR_OPEN_STAT);
	
	ret=fprintf_s(file, "%d\n", newNbStructs);
	//ret=_fwrite(file, &newNbStructs, sizeof(long));
	if(ret<=0){
		printf("Error writing the second number in the file %s\n", filename);
		return(ERR_WRITE_STAT);
	}
	for(long ind=1; ind<=nbStructures;ind++){
		if(results[ind].selected==1){
			for(int w=0; w<sizeof(HCSData)/sizeof(double);w++){
				if(w!=sizeof(HCSData)/sizeof(double)-1)
					ret=fprintf_s(file, "%f\t", ((double*)(results+ind))[w]);
				else
					ret=fprintf_s(file, "%f\n", ((double*)(results+ind))[w]);
				if(ret<=0){
					printf("Error writing the second number in the file %s\n", filename);
					return(ERR_WRITE_STAT);
				}
			}
		}
	}

	fclose(file);

	return(RET_OK);
}


int renumberResults(HCSData	*results, USrtArray& mainSeg, long nbStructures){

	int newIndex;
	newIndex=1;
	for(long ind=1;ind<=nbStructures;ind++){
		if(results[ind].selected==1.){
			results[ind].index=newIndex;
			newIndex++;
		}else
			results[ind].index=0;
	}
	for(unsigned long x=0;x<mainSeg.nbPoints;x++){
		if(mainSeg.data[x]>0)
			mainSeg.data[x]=(unsigned short)results[mainSeg.data[x]].index;
	}
	return(newIndex-1);
}

int downSample(USrtArray& im, USrtArray& downSampledIm,long downSampleRange){
	int ret;
	long downSampleX=(long)ceil((double)im.size(0)/(double)downSampleRange);
	long downSampleY=(long)ceil((double)im.size(1)/(double)downSampleRange);
	ret=downSampledIm.alloc(downSampleX, downSampleY);
	TEST_RET;
	long sum;
	long saturated;
	saturated=0;
	for(long dsI=0; dsI<downSampleX; dsI++)
	for(long dsJ=0; dsJ<downSampleY; dsJ++){
		sum=0;
		for(long i=0; i<downSampleRange; i++)
		for(long j=0; j<downSampleRange; j++)
			if((dsI*downSampleRange+i<im.size(0))&&(dsJ*downSampleRange+j<im.size(1)))
				sum+=im(dsI*downSampleRange+i,dsJ*downSampleRange+j);
		if(sum>downSampledIm.maxVal){
			downSampledIm(dsI, dsJ)=(unsigned short)downSampledIm.maxVal;
			saturated++;
		}else{
			downSampledIm(dsI, dsJ)=(unsigned short)sum;
		}
	}
	if(saturated>0)
		printf("Saturation detected during downsample: %d\n", saturated);

	return(RET_OK);
}



int computeBackground(USrtArray &mainSeg, USrtArray &image, FltArray &imageOutF, HCSData* results, long nbStruct, long halfWindowsSize, long minNbPoints, long channel){

	int ret;
	for(long i=1;i<=nbStruct;i++){
		if(results[i].selected==1.){
			long nbPoints=0;
			double totalIntensity=0.;

			for(long x=max(0, (long)results[i].seg.center.x-halfWindowsSize);x<min(mainSeg.size(0), (long)results[i].seg.center.x+halfWindowsSize+1);x++)
			for(long y=max(0, (long)results[i].seg.center.y-halfWindowsSize);y<min(mainSeg.size(1), (long)results[i].seg.center.y+halfWindowsSize+1);y++){
				if(mainSeg.data[x*mainSeg.size(1)+y]==0){
					nbPoints++;
					totalIntensity+=image.data[x*mainSeg.size(1)+y];
				}
			}
			if(nbPoints>=minNbPoints){
				results[i].chIntensities[channel].valBackground=totalIntensity/(double)nbPoints;
			}else{
				results[i].selected=0;
				results[i].chIntensities[channel].valBackground=-1000.;
			}
		}
	}

	ret=imageOutF.alloc(image.size(0), image.size(1));
	TEST_RET;
	double val;
	for(long x=0;x<mainSeg.size(0);x++)
	for(long y=0;y<mainSeg.size(1);y++){
		val=0.;
		if(mainSeg.data[x*mainSeg.size(1)+y]!=0){
			if(results[mainSeg.data[x*mainSeg.size(1)+y]].selected!=0){
				val=(double)image.data[x*mainSeg.size(1)+y]-results[mainSeg.data[x*mainSeg.size(1)+y]].chIntensities[channel].valBackground;
			}
		}
		imageOutF.data[x*mainSeg.size(1)+y]=(float)val;
	}


	return(RET_OK);
}


// gives the distance between the points 1 and 2
double d(double x1, double y1, double x2, double y2){
	return(sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)));
}

// gives the distance square between the points 1 and 2
double d2(double x1, double y1, double x2, double y2){
	return((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}
int computeMainSegmentation(USrtArray &mainSegmentation, HCSData* results, const long nbStruct, const double distMin){
	double distMin2=distMin*distMin;
	long distMinL=(long)ceil(distMin);
	for(long x=0;x<mainSegmentation.size(0);x++)
	for(long y=0;y<mainSegmentation.size(1);y++){
		unsigned long indexCloser;
		double		  distCloser2;
		double		  distLocal2;
		indexCloser=0;
		for(long i=1;i<=nbStruct;i++){
			if((abs(x-results[i].seg.center.x)<=distMinL)&&(abs(y-results[i].seg.center.y)<=distMinL)){
				distLocal2=d2(x,y,results[i].seg.center.x, results[i].seg.center.y);
				if(distLocal2<=distMin2){
					if(indexCloser==0){
						distCloser2=distLocal2;
						indexCloser=i;
					}else 
						if(distLocal2<=distCloser2){
							distCloser2=distLocal2;
							indexCloser=i;
						}
				}
			}
		}
		mainSegmentation.data[x*mainSegmentation.size(1)+y]=(unsigned short)indexCloser;
	}
	for(long i=1;i<=nbStruct;i++){
		results[i].seg.nbPixels=0;
		results[i].seg.perimeter=0;
		results[i].seg.perimeterWithoutNeighbor=0;
		results[i].seg.boundingBox.xMin=mainSegmentation.size(0);
		results[i].seg.boundingBox.xMax=0;
		results[i].seg.boundingBox.yMin=mainSegmentation.size(1);
		results[i].seg.boundingBox.yMax=0;
		
	}
	for(long x=0;x<mainSegmentation.size(0);x++)
	for(long y=0;y<mainSegmentation.size(1);y++){
		long objectIndex=mainSegmentation.data[x*mainSegmentation.size(1)+y];
		bool isPerimeter=false;
		bool isPerimeterWithoutNeighbor=false;
		if(objectIndex>0){
			results[objectIndex].seg.nbPixels++;
			if(x<results[objectIndex].seg.boundingBox.xMin)
				results[objectIndex].seg.boundingBox.xMin=x;
			if(x>results[objectIndex].seg.boundingBox.xMax)
				results[objectIndex].seg.boundingBox.xMax=x;
			if(y<results[objectIndex].seg.boundingBox.yMin)
				results[objectIndex].seg.boundingBox.yMin=y;
			if(y>results[objectIndex].seg.boundingBox.yMax)
				results[objectIndex].seg.boundingBox.yMax=y;
			if(x-1>=0)
				if(mainSegmentation.data[(x-1)*mainSegmentation.size(1)+y]!=objectIndex){
					isPerimeter=true;
					if(mainSegmentation.data[(x-1)*mainSegmentation.size(1)+y]==0)
						isPerimeterWithoutNeighbor=true;
				}
			if(x+1<mainSegmentation.size(0))
				if(mainSegmentation.data[(x+1)*mainSegmentation.size(1)+y]!=objectIndex){
					isPerimeter=true;
					if(mainSegmentation.data[(x+1)*mainSegmentation.size(1)+y]==0)
						isPerimeterWithoutNeighbor=true;
				}
			if(y-1>=0)
				if(mainSegmentation.data[x*mainSegmentation.size(1)+y-1]!=objectIndex){
					isPerimeter=true;
					if(mainSegmentation.data[x*mainSegmentation.size(1)+y-1]==0)
						isPerimeterWithoutNeighbor=true;
				}
			if(y+1<mainSegmentation.size(1))
				if(mainSegmentation.data[x*mainSegmentation.size(1)+y+1]!=objectIndex){
					isPerimeter=true;
					if(mainSegmentation.data[x*mainSegmentation.size(1)+y+1]==0)
						isPerimeterWithoutNeighbor=true;
				}
			if(isPerimeter){
					results[objectIndex].seg.perimeter++;
				if(isPerimeterWithoutNeighbor)
					results[objectIndex].seg.perimeterWithoutNeighbor++;
			}
		}
		
	}

	return(RET_OK);
}


int filterMainSegmentation(USrtArray &mainSegmentation, HCSData* results, const long nbStruct, const double pctMinBorderNoCells){
	
	for(long i=1;i<=nbStruct;i++){
		if(results[i].seg.perimeterWithoutNeighbor/results[i].seg.perimeter<pctMinBorderNoCells)
			results[i].selected=0.;
	}
	for(long x=0;x<mainSegmentation.size(0);x++)
	for(long y=0;y<mainSegmentation.size(1);y++){
		long objectIndex=mainSegmentation.data[x*mainSegmentation.size(1)+y];
		if(objectIndex>0)
			if(results[objectIndex].selected==0)
				mainSegmentation.data[x*mainSegmentation.size(1)+y]=0;
	}
	return(RET_OK);
}


int waveletSeg(USrtArray &channel, USrtArray &channelSeg, const unsigned long noWavelet, const float thresholdVal, const float minVol, const float minInt, Segment2D **segmentHCSPub){
	int ret;
	double structFiltParam[13]=           {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
	unsigned short structFiltParamUse[13]={0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 };
	unsigned short planesToUse[5]={0,0,0,0,0};
	unsigned short operators[5]=  {0,0,0,0,0};
	double threshold[5]=          {0,0,0,0,0};
	
	planesToUse[noWavelet]=2;
	threshold[noWavelet]=thresholdVal;

	//structFiltParamUse[12]=1; // no border
	if(minVol>0.){
		structFiltParamUse[0]=1; // minimum volume
		structFiltParam[0]=250;
	}
	if(minInt>0.){
		structFiltParamUse[2]=1; // minimum intensity
		structFiltParam[2]=minInt; 
	}
	if(segmentHCSPub!=NULL){
		ret=waveletSegmentation2DHCS( 
							//nucleusDownSample->data, output->data, 
							channel.data, channelSeg.data, 
							//nucleusDownSample->size(0),nucleusDownSample->size(1), 
							channel.size(0),channel.size(1),
							0/*unsigned short nbAngle*/, 
							0/*unsigned short anisoTransforType*/, 
							/*unsigned short**/ planesToUse, 
							/*unsigned short**/ operators, 
							/*double**/ threshold, 
							2/*unsigned short renormalize*/, 
							0/*double renorm1*/, 
							0/*double renorm2*/, 
							3/*unsigned short visualization*/, 
							1/*unsigned short segmentResult*/, 
							1/*unsigned short useDiagonalForSeg*/, 
							2/*unsigned short computeStatistics*/, 
							0/*unsigned short realProj*/, 
							1/*unsigned short structFilt*/, 
							/*double**/ structFiltParam, 
							/*unsigned short**/ structFiltParamUse, 
							1/*unsigned short watershed*/, 
							0.4/*double watershedThresholdRatio*/, 
							0/*short autoThresholdModel*/, segmentHCSPub);
	}else{
		ret=waveletSegmentation2DHCS( 
							//nucleusDownSample->data, output->data, 
							channel.data, channelSeg.data, 
							//nucleusDownSample->size(0),nucleusDownSample->size(1), 
							channel.size(0),channel.size(1),
							0/*unsigned short nbAngle*/, 
							0/*unsigned short anisoTransforType*/, 
							/*unsigned short**/ planesToUse, 
							/*unsigned short**/ operators, 
							/*double**/ threshold, 
							2/*unsigned short renormalize*/, 
							0/*double renorm1*/, 
							0/*double renorm2*/, 
							2/*unsigned short visualization*/, 
							0/*unsigned short segmentResult*/, 
							1/*unsigned short useDiagonalForSeg*/, 
							0/*unsigned short computeStatistics*/, 
							0/*unsigned short realProj*/, 
							0/*unsigned short structFilt*/, 
							/*double**/ structFiltParam, 
							/*unsigned short**/ structFiltParamUse, 
							0/*unsigned short watershed*/, 
							0.4/*double watershedThresholdRatio*/, 
							0/*short autoThresholdModel*/, segmentHCSPub);
	}
	TEST_RET;
	return(RET_OK);
}
