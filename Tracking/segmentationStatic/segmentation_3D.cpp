
/*******************************************************************************************************************
								Data segmentation for 3D image or image stack  (3D signal) (16bit only)
									Institut Curie
									UMR - 144
									by Victor Racine
									04/15/02 
*******************************************************************************************************************/

#include "segmentation_3D.h"

#define SEG(i,j,k)	seg.data[( (i)*sizeDataY+(j) )*sizeDataZ+(k)]
#define SEGHoles(i,j,k)	seg.data[( (i)*sizeDataY+(j) )*sizeDataZ+(k)]

#define DATA(i,j,k)	data.data[( (i)*sizeDataY+(j) )*sizeDataZ+(k)]

inline long min(long a, long b){
	return(a<=b ? a:b);
}
inline long min(long a, long b, long c, long d){
	return(min(min(a,b),min(c,d)));
}

inline double min(double a, double b){
	return(a<=b ? a:b);
}

inline double min(double a, double b, double c){
	return(min(min(a,b),c));
}

inline long max(long a, long b){
	return(a>=b ? a:b);
}



inline double max(double a, double b){
	return(a>=b ? a:b);
}

inline double max(double a, double b, double c){
	return(max(max(a,b),c));
}


void Segment3D::getSize(long &sizeDataX, long &sizeDataY, long &sizeDataZ, USrtArray &tab){
	sizeDataX=tab.size(2);
	sizeDataY=tab.size(0);
	sizeDataZ=tab.size(1);
}

void** Segment3D::getTabInfo(){
	return((void**)tabInfo);
}

void Segment3D::desalloc(){
	InfoStruct3D* ptr;
	ptr=headInfo;
	while(ptr!=NULL){
		headInfo=ptr;
		ptr=ptr->next;
		delete headInfo;
	}
	if(tabInfo!=NULL)
		delete [] tabInfo;
	tabInfo=NULL;
	nbInfoStruct=0;
	headInfo=NULL;
}

Segment3D::Segment3D(){

	headInfo=NULL;
	nbInfoStruct=0;
	tabInfo=NULL;
}

Segment3D::~Segment3D(){

	desalloc();
}


int Segment3D::alloc(unsigned long nbStruct){
	InfoStruct3D* ptr;
	unsigned long ind;
	desalloc();
	tabInfo=new InfoStruct3D*[nbStruct+1];
	if(tabInfo==NULL) return(S_ALLOC_ERR);

	if(nbStruct>0){
		headInfo=new InfoStruct3D(1);
		tabInfo[0]=NULL;
		tabInfo[1]=headInfo;
		ptr=headInfo;
		for(ind=2;ind<=nbStruct;ind++){
			ptr->next=new InfoStruct3D(ind);
			tabInfo[ind]=ptr->next;
			if(ptr->next==NULL) return(S_ALLOC_ERR);
			ptr=ptr->next;
		}
	}
	nbInfoStruct=nbStruct;
	return(RET_OK);
}

unsigned long Segment3D::filterData(USrtArray &seg, double* param, unsigned short* useParam){

	InfoStruct3D* ptr=headInfo;
	InfoStruct3D* ptrMax;
	InfoStruct3D* ptrOld=NULL;
	unsigned long ind=0;
	long i,j,indexFilter;
	bool unselected;
	float valMax;
	long nbInfoStructOld=nbInfoStruct;
	long oldIndex;

	unsigned short sizeDataX=seg.size(0);
	unsigned short sizeDataY=seg.size(1);
	unsigned short sizeDataZ=seg.size(2);

	char *tabBound=NULL;
	if(useParam[12]==1){
		tabBound=new char[nbInfoStruct+1];
		for(i=0;i<(long)(nbInfoStruct+1);i++)
			tabBound[i]=0;
		for(i=0;i<seg.size(0);i++)
		for(j=0;j<seg.size(1);j++){
			if(SEG(i,j,0)!=0)
				tabBound[SEG(i,j,0)]=1;
			if(SEG(i,j,seg.size(2)-1)!=0)
				tabBound[SEG(i,j,seg.size(2)-1)]=1;
		}
		for(i=0;i<seg.size(1);i++)
		for(j=0;j<seg.size(2);j++){
			if(SEG(0,i,j)!=0)
				tabBound[SEG(0,i,j)]=1;
			if(SEG(seg.size(0)-1,i,j)!=0)
				tabBound[SEG(seg.size(0)-1,i,j)]=1;
		}
		for(i=1;i<seg.size(2);i++)
		for(j=0;j<seg.size(0);j++){
			if(SEG(j,0,i)!=0)
				tabBound[SEG(j,0,i)]=1;
			if(SEG(j,seg.size(1)-1,i)!=0)
				tabBound[SEG(j,seg.size(1)-1,i)]=1;
		}
	}

	while(ptr!=NULL){
		ind++;
		unselected=false;
		if((useParam[0]==1)&&(ptr->data.volume<=param[0])) unselected=true; //minVol
		else if((useParam[1]==1)&&(ptr->data.volume>=param[1])) unselected=true; //maxSurf
		else if((useParam[2]==1)&&(ptr->data.mean*ptr->data.volume<=param[2])) unselected=true; //minInten
		else if((useParam[3]==1)&&(ptr->data.mean*ptr->data.volume>=param[3])) unselected=true; //maxInten
		else if((useParam[4]==1)&&(ptr->data.max.val<=param[4])) unselected=true; //minMaxi
		else if((useParam[5]==1)&&(ptr->data.max.val>=param[5])) unselected=true; //maxMaxi
		else if((useParam[7]==1)&&(fabs(ptr->data.det)>1e-10)){
			if(min(ptr->data.sigma.x,ptr->data.sigma.y,ptr->data.sigma.z)
				/max(ptr->data.sigma.x,ptr->data.sigma.y,ptr->data.sigma.z)>=param[7]) unselected=true; //minMorphParam
		}else if((useParam[8]==1)&&((fabs(ptr->data.det)>1e-10))){
			if(min(ptr->data.sigma.x,ptr->data.sigma.y,ptr->data.sigma.z)
				/max(ptr->data.sigma.x,ptr->data.sigma.y,ptr->data.sigma.z)<=param[8]) unselected=true; //maxMorphParam
		}else if(useParam[12]==1)if(tabBound[ind]==1) unselected=true; //no struct on boundaries
		if(unselected){
			if(ptrOld==NULL){
				headInfo=ptr->next;
				delete ptr;
				ptr=headInfo;
			}else{
				ptrOld->next=ptr->next;
				delete ptr;
				ptr=ptrOld->next;
			}
			tabInfo[ind]=NULL;
			nbInfoStruct--;
		}else{
			ptrOld=ptr;
			ptr=ptr->next;
		}
	}
	if(useParam[12]==1)
		delete [] tabBound;

	if(useParam[6]==1){
		ptr=headInfo;
		ind=0;
		ptrOld=NULL;
		while(ptr!=NULL){
			ind++;
			unselected=false;
			for(i=1;(i<=nbInfoStructOld)&&!unselected;i++)
				if(tabInfo[i]!=NULL){
					double d=sqrt( (ptr->data.max.x-tabInfo[i]->data.max.x)*(ptr->data.max.x-tabInfo[i]->data.max.x) 
						+ (ptr->data.max.y-tabInfo[i]->data.max.y)*(ptr->data.max.y-tabInfo[i]->data.max.y)
						+ (ptr->data.max.z-tabInfo[i]->data.max.z)*(ptr->data.max.z-tabInfo[i]->data.max.z));
					if((d>1)&&(d<param[6])&&(tabInfo[i]->data.max.val>ptr->data.max.val))
						unselected=true;
				}
			if(unselected){
				if(ptrOld==NULL){
					headInfo=ptr->next;
					oldIndex=ptr->data.index;
					delete ptr;
					ptr=headInfo;
				}else{
					ptrOld->next=ptr->next;
					oldIndex=ptr->data.index;
					delete ptr;
					ptr=ptrOld->next;
				}
				tabInfo[oldIndex]=NULL;
				nbInfoStruct--;
			}else{
				ptrOld=ptr;
				ptr=ptr->next;
			}
		}
	}

	for(indexFilter=9;indexFilter<=11;indexFilter++)
		if( (useParam[indexFilter]==1) && (nbInfoStruct>param[indexFilter]) && ((long)param[indexFilter]>0) ){	
					//N first structure with the biggest surface (indexFilter=7), the biggest surface (indexFilter=8) or 
						//the biggest total intensity (indexFilter=9)
			ptrOld=NULL;
			for(i=0;i<(long)param[indexFilter];i++){
				if(ptrOld==NULL)
					ptr=headInfo;
				else
					ptr=ptrOld->next;

				switch(indexFilter){
				case 9: valMax=(float)ptr->data.volume; break;
				case 10: valMax=(float)(ptr->data.volume*ptr->data.mean); break;
				case 11: valMax=(float)(ptr->data.max.val); break;
				}

				ptrMax=NULL;

				while(ptr->next!=NULL){
					switch(indexFilter){
					case 9: 
						if(valMax<=ptr->next->data.volume){
							valMax=(float)ptr->next->data.volume;
							ptrMax=ptr;
						}
						break;
					case 10:
						if(valMax<=(long)ptr->next->data.volume*ptr->next->data.mean){
							valMax=(float)((long)ptr->next->data.volume*ptr->next->data.mean);
							ptrMax=ptr;
						}
						break;
					case 11:
						if(valMax<=(float)ptr->next->data.max.val){
							valMax=(float)(ptr->next->data.max.val);
							ptrMax=ptr;
						}
						break;
					}
					ptr=ptr->next;
				}
				if(ptrMax!=NULL){
					if(ptrOld!=NULL){
						ptr=ptrOld->next;
						ptrOld->next=ptrMax->next;
						ptrMax->next=ptrMax->next->next;
						ptrOld->next->next=ptr;
						ptrOld=ptrOld->next;
					}else{
						ptr=headInfo;
						headInfo=ptrMax->next;
						ptrMax->next=headInfo->next;
						headInfo->next=ptr;
						ptrOld=headInfo;
					}
				}else
					if(ptrOld==NULL)
						ptrOld=headInfo;
					else
						ptrOld=ptrOld->next;
			}
			ptr=ptrOld->next;
			ptrOld->next=NULL;
			while(ptr!=NULL){
				ptrMax=ptr->next;
				tabInfo[ptr->data.index]=NULL;
				delete ptr;
				nbInfoStruct--;
				ptr=ptrMax;
			}
		}
	
	ptr=headInfo;
	for(ind=1;ind<=nbInfoStruct;ind++){
		ptr->data.index=ind;
		ptr=ptr->next;
	}

	for(i=0;i<(long)seg.nbPoints;i++)
		if(seg.data[i]!=0){
			if(tabInfo[seg.data[i]]==NULL)
				seg.data[i]=0;
			else
				seg.data[i]=(unsigned short)tabInfo[seg.data[i]]->data.index;
		}
	delete [] tabInfo;
	tabInfo=new InfoStruct3D*[nbInfoStruct+1];
	tabInfo[0]=NULL;
	tabInfo[1]=headInfo;
	//tabInfo[1]->data.index=1;
	ptr=headInfo;
	for(ind=2;ind<=nbInfoStruct;ind++){
		tabInfo[ind]=ptr->next;
		//tabInfo[ind]->data.index=ind;
		ptr=ptr->next;
	}

	return(nbInfoStruct);
}

void Segment3D::drawSymbolOnStructures(USrtArray &seg, unsigned short typeOfSign){
	unsigned long ind;
	long x,y,z;
	double a,b,c;
	double rotX, rotY, rotZ;
	
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);

	seg=(unsigned short)0;

	if(typeOfSign==4){	//point
		for(ind=1;ind<=nbInfoStruct;ind++){
			seg.data[((unsigned long)tabInfo[ind]->data.moment.x*sizeDataY+(unsigned long)tabInfo[ind]->data.moment.y )*sizeDataZ+(unsigned long)tabInfo[ind]->data.moment.z ]=(unsigned short)ind;		
		}
	}
	if(typeOfSign==5)	//cross
		for(ind=1;ind<=nbInfoStruct;ind++){
			for(x=-SIZE_CROSS;x<SIZE_CROSS+1;x++)
				if((tabInfo[ind]->data.moment.x+x>=0) && (tabInfo[ind]->data.moment.x+x<sizeDataX))
					seg.data[(((unsigned long)tabInfo[ind]->data.moment.x+x)*sizeDataY+(unsigned long)tabInfo[ind]->data.moment.y )*sizeDataZ+
						(unsigned long)tabInfo[ind]->data.moment.z ]=(unsigned short)ind;
			for(y=-SIZE_CROSS;y<SIZE_CROSS+1;y++)
				if((tabInfo[ind]->data.moment.y+y>=0) && (tabInfo[ind]->data.moment.y+y<sizeDataY))
					seg.data[((unsigned long)tabInfo[ind]->data.moment.x*sizeDataY+((unsigned long)tabInfo[ind]->data.moment.y+y) )*sizeDataZ+
						(unsigned long)tabInfo[ind]->data.moment.z ]=(unsigned short)ind;
			for(z=-SIZE_CROSS;z<SIZE_CROSS+1;z++)
				if((tabInfo[ind]->data.moment.z+z>=0) && (tabInfo[ind]->data.moment.z+z<sizeDataZ))
					seg.data[((unsigned long)tabInfo[ind]->data.moment.x*sizeDataY+(unsigned long)tabInfo[ind]->data.moment.y )*sizeDataZ+
						((unsigned long)tabInfo[ind]->data.moment.z+z) ]=(unsigned short)ind;
		}

#define DELTA	0.2
	if(typeOfSign==6)	//ellipse
		for(ind=1;ind<=nbInfoStruct;ind++)
			if(fabs(tabInfo[ind]->data.det>1e-10))		//det is significant
				if((tabInfo[ind]->data.sigma.x>1e-10) && (tabInfo[ind]->data.sigma.y>1e-10) && (tabInfo[ind]->data.sigma.z>1e-10))
					for(a=-tabInfo[ind]->data.sigma.x-DELTA; a<=tabInfo[ind]->data.sigma.x+DELTA; a+=DELTA)
					for(b=-tabInfo[ind]->data.sigma.y-DELTA; b<=tabInfo[ind]->data.sigma.y+DELTA; b+=DELTA)
					for(c=-tabInfo[ind]->data.sigma.z-DELTA; c<=tabInfo[ind]->data.sigma.z+DELTA; c+=DELTA)
						if( a*a/(tabInfo[ind]->data.sigma.x*tabInfo[ind]->data.sigma.x) + b*b/(tabInfo[ind]->data.sigma.y*tabInfo[ind]->data.sigma.y) + 
								c*c/(tabInfo[ind]->data.sigma.z*tabInfo[ind]->data.sigma.z) <= 1. ){
							rotX=a;
							rotY=b;
							rotZ=c;
							eulerRotation(&rotX, &rotY, &rotZ, tabInfo[ind]->data.angle.phi, tabInfo[ind]->data.angle.teta, tabInfo[ind]->data.angle.psi);
							if( ((long)(tabInfo[ind]->data.moment.x+rotX)>=0) && ((long)(tabInfo[ind]->data.moment.x+rotX)<sizeDataX) )
							if( ((long)(tabInfo[ind]->data.moment.y+rotY)>=0) && ((long)(tabInfo[ind]->data.moment.y+rotY)<sizeDataY) )
							if( ((long)(tabInfo[ind]->data.moment.z+rotZ)>=0) && ((long)(tabInfo[ind]->data.moment.z+rotZ)<sizeDataZ) )
								SEG( (long)(tabInfo[ind]->data.moment.x+rotX), (long)(tabInfo[ind]->data.moment.y+rotY), 
									(long)(tabInfo[ind]->data.moment.z+rotZ) )=(unsigned short)ind;
						}
#undef DELTA
}



int Segment3D::printStatistics(wchar_t *fileName, unsigned short resultType, unsigned long time, unsigned long nbTime, unsigned short atEndOfFile,
							   unsigned long cptPoints, unsigned short sizeX, unsigned short sizeY, unsigned short sizeStk){
	unsigned long ind;
	FILE *stream;
	int fileHandle;
	long nbDim=3;
	long ret;

	switch(resultType){
	case 1:	//Binary file
		if((atEndOfFile==0) || (time==1)){
			if( (fileHandle = _wopen( fileName, O_RDWR | O_CREAT | O_BINARY | _O_TRUNC, S_IWRITE )) == -1 ) return(S_OPEN_FILE_ERR);
			if(write(fileHandle, &nbDim, sizeof(long))!=sizeof(long)) return(S_WRTIE_ERR);
			if(write(fileHandle, &nbTime, sizeof(long))!=sizeof(long)) return(S_WRTIE_ERR);
		}else
			if( (fileHandle = _wopen( fileName, O_RDWR | O_CREAT | O_BINARY | _O_APPEND, S_IWRITE )) == -1 ) return(S_OPEN_FILE_ERR);

		write(fileHandle, &nbInfoStruct, sizeof(long));
		for(ind=1;ind<=nbInfoStruct;ind++){
			tabInfo[ind]->data.time=time;
			if(write( fileHandle, &(tabInfo[ind]->data), sizeof(InfoData3D))!=sizeof(InfoData3D)) return(S_WRTIE_ERR);
		}
		close( fileHandle );
		break;
	case 2:	//ASCII file
		if((atEndOfFile==0) || (time==1)){
			if((stream=_wfopen( fileName, L"w" ))==NULL){ 
				return(S_OPEN_FILE_ERR);
			}
			fprintf(stream,"3D stat, total of %d times\n",nbTime);
		}else
			if((stream=_wfopen( fileName, L"a+" ))==NULL){ 
				return(S_OPEN_FILE_ERR);
			}

		fprintf(stream,"3D stat at time %d, pourcentage of significant value:\t%f\t",time, (double)cptPoints/((double)(sizeX*sizeY*sizeStk)));
		fprintf(stream,"NbStructuresDetected: %d\n",nbInfoStruct);

		fprintf(stream,"Index\tMean\tVolume\tIntensity\tVolumicSurface\tMorpho\tmaxX\tmaxY\tmaxZ\tvalMax\tCentroidX\tCentroidY\tCentroidZ\tBinaryCentroidX\tBinaryCentroidY\tBinaryCentroidZ\tSigmaX\tSigmaY\tSigmaZ\tAnglePhiRad\tAngleThetaRad\tAnglePsiRad\n");


		for(ind=1;ind<=nbInfoStruct;ind++){
			fprintf(stream,"%d\t%f\t%d\t%f\t%d\t%f\t",ind,	tabInfo[ind]->data.mean, 
				tabInfo[ind]->data.volume, tabInfo[ind]->data.intensity ,tabInfo[ind]->data.surfVol, tabInfo[ind]->data.morpho);
			fprintf(stream,"%d\t%d\t%d\t%f\t", tabInfo[ind]->data.max.x,tabInfo[ind]->data.max.y,tabInfo[ind]->data.max.z,tabInfo[ind]->data.max.val);
			fprintf(stream,"%f\t%f\t%f\t",tabInfo[ind]->data.moment.x,tabInfo[ind]->data.moment.y,tabInfo[ind]->data.moment.z);
			fprintf(stream,"%f\t%f\t%f\t",tabInfo[ind]->data.momentSurf.x,tabInfo[ind]->data.momentSurf.y,
				tabInfo[ind]->data.momentSurf.z);
			if(fabs(tabInfo[ind]->data.det)>1e-10){
				fprintf(stream,"%f\t%f\t%f\t", tabInfo[ind]->data.sigma.x,tabInfo[ind]->data.sigma.y,tabInfo[ind]->data.sigma.z);
				fprintf(stream,"%f\t",tabInfo[ind]->data.angle.phi);
				fprintf(stream,"%f\t",tabInfo[ind]->data.angle.teta);
				fprintf(stream,"%f\t",tabInfo[ind]->data.angle.psi);
			}
			fprintf(stream,"\n");

		}
		fclose( stream );
	break;
	case 3: //HTML file
		if((atEndOfFile==0) || (time==1)){
			if((stream=_wfopen( fileName, L"w" ))==NULL){ 
				return(S_OPEN_FILE_ERR);
			}
			//HTML 4.01
			fprintf(stream,"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n");

			//header
			fprintf(stream,"<html>\n<head>\n");
			fprintf(stream,"\t<META NAME=\"Author\" LANG=\"en\" CONTENT=\"Victor RACINE\">\n");
			fprintf(stream,"\t<META NAME=\"Publisher\" CONTENT=\"MIA_C++.dll\">\n");
			fprintf(stream,"\t<META NAME=\"Reply-to\" CONTENT=\"victor.racine@curie.fr (Victor RACINE)\">\n");
			fprintf(stream,"\t<META NAME=\"Indentifier-URL\" CONTENT=\"http://materiels.curie.fr/imagerie/\">\n");
			fprintf(stream,"\t<title>\n\t\tMIA statistics\n\t</title>\n\n");
			fprintf(stream,"</head>\n<body>\n\n");


			//infos
			fprintf(stream,"<p>This file has been generated by <i>MIA_C++.dll</i>.<br>\n");
			//swapping of axis to be compatible with MM spec
			fprintf(stream,"3D statistics of an image %dx%dx%d, total of %d times.</p>\n",sizeY, sizeX,sizeStk, nbTime);

			//the begining of the table
			fprintf(stream,"<table>\n");
			
			//first line
			fprintf(stream,"\t<tr align=\"center\" bgcolor=\"#222288\">\n");
			printHTMLCell(L"Time", stream,2,1);
			printHTMLCell(L"Index", stream,2,1);
			printHTMLCell(L"Center of mass", stream,1,3);
			printHTMLCell(L"Volume center of mass", stream,1,3);
			printHTMLCell(L"Mean", stream,2,1);
			printHTMLCell(L"Volume", stream,2,1);
			printHTMLCell(L"Intensity", stream,2,1);
			printHTMLCell(L"Surface volume", stream,2,1);
			printHTMLCell(L"Shape factor", stream,2,1);
			printHTMLCell(L"Max", stream,1,4);
			printHTMLCell(L"Ellipse axis", stream,1,3);
			printHTMLCell(L"Angle (rad)", stream,1,3);
			fprintf(stream,"\t</tr>\n");

			//second line
			fprintf(stream,"\t<tr  align=\"center\" bgcolor=\"#222288\">\n");
			printHTMLCell(L"X", stream,1,1);
			printHTMLCell(L"Y", stream,1,1);
			printHTMLCell(L"Z", stream,1,1);
			printHTMLCell(L"X", stream,1,1);
			printHTMLCell(L"Y", stream,1,1);
			printHTMLCell(L"Z", stream,1,1);
			printHTMLCell(L"X", stream,1,1);
			printHTMLCell(L"Y", stream,1,1);
			printHTMLCell(L"Z", stream,1,1);
			printHTMLCell(L"Val", stream,1,1);
			printHTMLCell(L"Major", stream,1,1);
			printHTMLCell(L"Middle", stream,1,1);
			printHTMLCell(L"Minor", stream,1,1);
			printHTMLCell(L"Phi", stream,1,1);
			printHTMLCell(L"Teta", stream,1,1);
			printHTMLCell(L"Psi", stream,1,1);
			fprintf(stream,"\t</tr>\n");

		}else{
			if((stream=_wfopen( fileName, L"r+" ))==NULL){ 
				return(S_OPEN_FILE_ERR);
			}
			if((ret=locateStringInFile(L"</table>", stream))<0){ 
				openDebugStreamALS();
				fprintf(debugStreamALS, "Error to find the string in the file in locateStringInFile\n");
				closeDebugStreamALS();
				return(ret);

			}else{
				if(fseek(stream,ret,SEEK_SET)!=0){
					openDebugStreamALS();
					fprintf(debugStreamALS, "Error to move in the file: offset=%d\n", -ret);
					closeDebugStreamALS();
					return(S_READ_ERR);
				}
				
			}
		}

		for(ind=1;ind<=nbInfoStruct;ind++){
			//new line
			if(ind%2)
				fprintf(stream,"\t<tr bgcolor=\"#EEFFDD\">\n");  //yellow
			else
				fprintf(stream,"\t<tr>\n");
			if(ind==1)
				printIntHTMLCell(time, stream, nbInfoStruct);  //Time
			printIntHTMLCell(ind, stream);  //Index
			printFloatHTMLCell(tabInfo[ind]->data.moment.z, stream);  //Center of mass X //swapping of axis to be compatible with MM spec
			printFloatHTMLCell(tabInfo[ind]->data.moment.y, stream);  //Center of mass Y //swapping of axis to be compatible with MM spec
			printFloatHTMLCell(tabInfo[ind]->data.moment.x, stream);  //Center of mass Z //swapping of axis to be compatible with MM spec
			printFloatHTMLCell(tabInfo[ind]->data.momentSurf.z, stream);  //Surface center of mass X //swapping of axis to be compatible with MM spec
			printFloatHTMLCell(tabInfo[ind]->data.momentSurf.y, stream);  //Surface center of mass Y //swapping of axis to be compatible with MM spec
			printFloatHTMLCell(tabInfo[ind]->data.momentSurf.x, stream);  //Surface center of mass Z //swapping of axis to be compatible with MM spec
			printFloatHTMLCell(tabInfo[ind]->data.mean, stream);  //Mean
			printIntHTMLCell(tabInfo[ind]->data.volume, stream);  //Volume
			printFloatHTMLCell(tabInfo[ind]->data.intensity, stream);  //Intensity
			printIntHTMLCell(tabInfo[ind]->data.surfVol, stream);  //Surface volume
			printFloatHTMLCell(tabInfo[ind]->data.morpho, stream);  //Shape factor
			printIntHTMLCell(tabInfo[ind]->data.max.z, stream);  //Max X //swapping of axis to be compatible with MM spec
			printIntHTMLCell(tabInfo[ind]->data.max.y, stream);  //Max Y //swapping of axis to be compatible with MM spec
			printIntHTMLCell(tabInfo[ind]->data.max.x, stream);  //Max Z //swapping of axis to be compatible with MM spec
			printFloatHTMLCell(tabInfo[ind]->data.max.val, stream);  //Max Value
			if(fabs(tabInfo[ind]->data.det)>1e-10){
				printFloatHTMLCell(tabInfo[ind]->data.sigma.z, stream);  //Ellipse major axis
				printFloatHTMLCell(tabInfo[ind]->data.sigma.y, stream);  //Ellipse middle axis
				printFloatHTMLCell(tabInfo[ind]->data.sigma.x, stream);  //Ellipse minor axis
				printFloatHTMLCell(tabInfo[ind]->data.angle.phi, stream);  //Angle Phi (rad)
				printFloatHTMLCell(tabInfo[ind]->data.angle.teta, stream);  //Angle Teta (rad)
				printFloatHTMLCell(tabInfo[ind]->data.angle.psi, stream);  //Angle Psi (rad)
			}
			fprintf(stream,"\t</tr>\n");
			
		}
		//foot of file
		fprintf(stream,"</table>\n</body>\n</html>\n");
		fclose( stream );
		
		break;
	}
	return(RET_OK);
}

int Segment3D::computeStatistics(USrtArray &data, USrtArray &seg){

	unsigned short	i,j,k;
	double			val;
	unsigned long	ind;
	unsigned long	cptPoints=0;
	double			powDet;
	int nrot;
	int ret; 
	double	deter;

	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);

	
	//allocation of array
	if((ret=alloc(nbInfoStruct))<0) return(ret);
	//initalisation of array
	for(ind=1;ind<=nbInfoStruct;ind++){
		tabInfo[ind]->data.surfVol = 0;
		tabInfo[ind]->data.volume = 0;
        tabInfo[ind]->data.morpho = 0.;
		tabInfo[ind]->data.max.x = -1;
		tabInfo[ind]->data.max.y = -1;
		tabInfo[ind]->data.max.z = -1;
		tabInfo[ind]->data.max.val = 0.;
		tabInfo[ind]->data.moment.x = 0.;
		tabInfo[ind]->data.moment.y = 0.;
		tabInfo[ind]->data.moment.z = 0.;
		tabInfo[ind]->data.momentSurf.x = 0.;
		tabInfo[ind]->data.momentSurf.y = 0.;
		tabInfo[ind]->data.momentSurf.z = 0.;
		tabInfo[ind]->data.mean = 0.;
		tabInfo[ind]->data.var.xx = 0.;
		tabInfo[ind]->data.var.yy = 0.;
		tabInfo[ind]->data.var.zz = 0.;
		tabInfo[ind]->data.var.xy = 0.;
		tabInfo[ind]->data.var.xz = 0.;
		tabInfo[ind]->data.var.yz = 0.;
		tabInfo[ind]->data.sigma.x = 0.;
		tabInfo[ind]->data.sigma.y = 0.;
		tabInfo[ind]->data.sigma.z = 0.;
	}


	for(i=0;i<sizeDataX;i++)
	for(j=0;j<sizeDataY;j++)
	for(k=0;k<sizeDataZ;k++)
		if(SEG(i,j,k)>0){
			cptPoints++;
			val=DATA(i,j,k);
			ind=SEG(i,j,k);

			tabInfo[ind]->data.volume++;
			tabInfo[ind]->data.moment.x += val*i;
			tabInfo[ind]->data.moment.y += val*j;
			tabInfo[ind]->data.moment.z += val*k;
			tabInfo[ind]->data.momentSurf.x += i;
			tabInfo[ind]->data.momentSurf.y += j;
			tabInfo[ind]->data.momentSurf.z += k;
			tabInfo[ind]->data.mean += val;
			tabInfo[ind]->data.var.xx += val*i*i;
			tabInfo[ind]->data.var.yy += val*j*j;
			tabInfo[ind]->data.var.zz += val*k*k;
			tabInfo[ind]->data.var.xy +=  val*i*j;
			tabInfo[ind]->data.var.xz +=  val*i*k;
			tabInfo[ind]->data.var.yz +=  val*k*j;

			if (fabs(val) > tabInfo[ind]->data.max.val){
				tabInfo[ind]->data.max.val = fabs(val);
				tabInfo[ind]->data.max.x = i;
				tabInfo[ind]->data.max.y = j;
				tabInfo[ind]->data.max.z = k;
			}


			if(i==0) 
				tabInfo[ind]->data.surfVol++;
			else
				if (SEG(i-1,j,k)!=SEG(i,j,k)) tabInfo[ind]->data.surfVol++;

			if(i==sizeDataX-1) 
				tabInfo[ind]->data.surfVol++;
			else
				if (SEG(i+1,j,k)!=SEG(i,j,k)) tabInfo[ind]->data.surfVol++;

			if(j==0) 
				tabInfo[ind]->data.surfVol++;
			else
				if (SEG(i,j-1,k)!=SEG(i,j,k)) tabInfo[ind]->data.surfVol++;

			if(j==sizeDataY-1) 
				tabInfo[ind]->data.surfVol++;   
			else
				if (SEG(i,j+1,k)!=SEG(i,j,k)) tabInfo[ind]->data.surfVol++;

			if(k==0) 
				tabInfo[ind]->data.surfVol++;
			else
				if (SEG(i,j,k-1)!=SEG(i,j,k)) tabInfo[ind]->data.surfVol++;

			if(k==sizeDataZ-1) 
				tabInfo[ind]->data.surfVol++;
			else
				if (SEG(i,j,k+1)!=SEG(i,j,k)) tabInfo[ind]->data.surfVol++;
		}
	

	double **c; //intercorrelation matrix
	double **v; //intercorrelation eigenvectors matrix 
	double *d;  //intercorrelation eigenvalues vector 
	// allocate intercorrelation matrix 
	c=matrix(1,3,1,3);
	// allocate output eigenvectors 
	v=matrix(1, 3, 1, 3);
	// allocate output eigenvalues 
	d=vector(1, 3);


	for(ind=1;ind<=nbInfoStruct;ind++){
		

		tabInfo[ind]->data.moment.x		/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.moment.y		/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.moment.z		/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.momentSurf.x	/= tabInfo[ind]->data.volume;
		tabInfo[ind]->data.momentSurf.y	/= tabInfo[ind]->data.volume;
		tabInfo[ind]->data.momentSurf.z	/= tabInfo[ind]->data.volume;
		tabInfo[ind]->data.var.xx			/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.var.yy			/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.var.zz			/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.var.xy			/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.var.xz			/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.var.yz			/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.var.xx			-= tabInfo[ind]->data.moment.x*tabInfo[ind]->data.moment.x;
		tabInfo[ind]->data.var.yy			-= tabInfo[ind]->data.moment.y*tabInfo[ind]->data.moment.y;
		tabInfo[ind]->data.var.zz			-= tabInfo[ind]->data.moment.z*tabInfo[ind]->data.moment.z;
		tabInfo[ind]->data.var.xy			-= tabInfo[ind]->data.moment.x*tabInfo[ind]->data.moment.y;
		tabInfo[ind]->data.var.xz			-= tabInfo[ind]->data.moment.x*tabInfo[ind]->data.moment.z;
		tabInfo[ind]->data.var.yz			-= tabInfo[ind]->data.moment.y*tabInfo[ind]->data.moment.z;
		tabInfo[ind]->data.mean			/= tabInfo[ind]->data.volume;
		tabInfo[ind]->data.intensity=tabInfo[ind]->data.mean*tabInfo[ind]->data.volume;

		//test if the matrix of variance is a rotation matrix
		if( fabs(tabInfo[ind]->data.det=fabs(tabInfo[ind]->data.var.xx*tabInfo[ind]->data.var.yy*tabInfo[ind]->data.var.zz+2.*tabInfo[ind]->data.var.xy*tabInfo[ind]->data.var.yz*tabInfo[ind]->data.var.xz-tabInfo[ind]->data.var.xz*tabInfo[ind]->data.var.xz*tabInfo[ind]->data.var.yy-tabInfo[ind]->data.var.yz*tabInfo[ind]->data.var.yz*tabInfo[ind]->data.var.xx-tabInfo[ind]->data.var.xy*tabInfo[ind]->data.var.xy*tabInfo[ind]->data.var.zz))>1e-10 ){

			//powDet=1./pow(fabs(tabInfo[ind]->data.det),0.3333333);
			c[1][1]=tabInfo[ind]->data.var.xx;
			c[2][2]=tabInfo[ind]->data.var.yy;
			c[3][3]=tabInfo[ind]->data.var.zz;
			c[1][2]=tabInfo[ind]->data.var.xy;
			c[1][3]=tabInfo[ind]->data.var.xz;
			c[2][1]=tabInfo[ind]->data.var.xy;
			c[3][1]=tabInfo[ind]->data.var.xz;
			c[2][3]=tabInfo[ind]->data.var.yz;
			c[3][2]=tabInfo[ind]->data.var.yz;


			//compute the the eigenvalues and the eignenvectors
			jacobi(c, 3, d, v, &nrot);

			deter=v[1][1]*v[2][2]*v[3][3] + v[1][2]*v[2][3]*v[3][1] + v[1][3]*v[2][1]*v[3][2] 
				- v[3][1]*v[2][2]*v[1][3] - v[3][2]*v[2][3]*v[1][1] - v[3][3]*v[2][1]*v[1][2];
			powDet=1./pow(fabs(deter),0.3333333);


			//acording to the euler angle:
			tabInfo[ind]->data.angle.psi=atan(v[1][3]/v[2][3]); //!!!
			tabInfo[ind]->data.angle.teta=asin(v[1][3]*powDet/sin(tabInfo[ind]->data.angle.psi));//*powDet
			tabInfo[ind]->data.angle.phi=asin(v[3][1]*powDet/sin(tabInfo[ind]->data.angle.teta));//*powDet

			if(sqrt(fabs(d[1]))>sqrt(fabs(d[3]))){ //1>3
				if(sqrt(fabs(d[1]))>sqrt(fabs(d[2]))){ //1>2 & 1>3
					if(sqrt(fabs(d[2]))>sqrt(fabs(d[3]))){ //2>3 => 1>2>3
						tabInfo[ind]->data.sigma.z=sqrt(fabs(d[1]));
						tabInfo[ind]->data.sigma.y=sqrt(fabs(d[2]));
						tabInfo[ind]->data.sigma.x=sqrt(fabs(d[3]));
					}else{//3>2 & 1>3 =>1>3>2
						tabInfo[ind]->data.sigma.z=sqrt(fabs(d[1]));
						tabInfo[ind]->data.sigma.y=sqrt(fabs(d[3]));
						tabInfo[ind]->data.sigma.x=sqrt(fabs(d[2]));
					}
				}else{ //2>1 & 1>3 => 2>1>3
					tabInfo[ind]->data.sigma.z=sqrt(fabs(d[2]));			
					tabInfo[ind]->data.sigma.y=sqrt(fabs(d[1]));
					tabInfo[ind]->data.sigma.x=sqrt(fabs(d[3]));
				}
			}else{ //3>1
				if(sqrt(fabs(d[1]))>sqrt(fabs(d[2]))){ //1>2 => 3>1>2
					tabInfo[ind]->data.sigma.z=sqrt(fabs(d[3]));
					tabInfo[ind]->data.sigma.y=sqrt(fabs(d[1]));
					tabInfo[ind]->data.sigma.x=sqrt(fabs(d[2]));
				}else{ //2>1 & 3>1
					if(sqrt(fabs(d[2]))>sqrt(fabs(d[3]))){ //2>3 => 2>3>1
						tabInfo[ind]->data.sigma.z=sqrt(fabs(d[2]));			
						tabInfo[ind]->data.sigma.y=sqrt(fabs(d[3]));
						tabInfo[ind]->data.sigma.x=sqrt(fabs(d[1]));
					}else{ //3>2 => 3>2>1
						tabInfo[ind]->data.sigma.z=sqrt(fabs(d[3]));			
						tabInfo[ind]->data.sigma.y=sqrt(fabs(d[2]));
						tabInfo[ind]->data.sigma.x=sqrt(fabs(d[1]));
					}
				}
			}
			//tabInfo[ind]->data.sigma.x=sqrt(fabs(d[1]));
			//tabInfo[ind]->data.sigma.y=sqrt(fabs(d[2]));
			//tabInfo[ind]->data.sigma.z=sqrt(fabs(d[3]));
		}

		tabInfo[ind]->data.morpho = tabInfo[ind]->data.sigma.x*tabInfo[ind]->data.sigma.y/tabInfo[ind]->data.sigma.z/tabInfo[ind]->data.sigma.z;
	}
	free_vector(d,1,3);
	free_matrix(c,1,3,1,3);
	free_matrix(v,1,3,1,3);
	
	return(cptPoints);
}


int Segment3D::computeStatistics(FltArray &data, USrtArray &seg){

	unsigned short	i,j,k;
	double			val;
	unsigned long		ind;
	unsigned long	cptPoints=0;
	double			powDet;
	int nrot;
	int ret; 
	double	deter;

	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);

	
	//allocation of array
	if((ret=alloc(nbInfoStruct))<0) return(ret);

	//initalisation of array
	for(ind=1;ind<=nbInfoStruct;ind++){
		tabInfo[ind]->data.surfVol = 0;
		tabInfo[ind]->data.volume = 0;
        tabInfo[ind]->data.morpho = 0.;
		tabInfo[ind]->data.max.x = -1;
		tabInfo[ind]->data.max.y = -1;
		tabInfo[ind]->data.max.z = -1;
		tabInfo[ind]->data.max.val = 0.;
		tabInfo[ind]->data.moment.x = 0.;
		tabInfo[ind]->data.moment.y = 0.;
		tabInfo[ind]->data.moment.z = 0.;
		tabInfo[ind]->data.momentSurf.x = 0.;
		tabInfo[ind]->data.momentSurf.y = 0.;
		tabInfo[ind]->data.momentSurf.z = 0.;
		tabInfo[ind]->data.mean = 0.;
		tabInfo[ind]->data.var.xx = 0.;
		tabInfo[ind]->data.var.yy = 0.;
		tabInfo[ind]->data.var.zz = 0.;
		tabInfo[ind]->data.var.xy = 0.;
		tabInfo[ind]->data.var.xz = 0.;
		tabInfo[ind]->data.var.yz = 0.;
		tabInfo[ind]->data.sigma.x = 0.;
		tabInfo[ind]->data.sigma.y = 0.;
		tabInfo[ind]->data.sigma.z = 0.;
	}

	for(i=0;i<sizeDataX;i++)
	for(j=0;j<sizeDataY;j++)
	for(k=0;k<sizeDataZ;k++)
		if(SEG(i,j,k)>0){
			cptPoints++;
			val=DATA(i,j,k);
			ind=SEG(i,j,k);

			tabInfo[ind]->data.volume++;
			tabInfo[ind]->data.moment.x += val*i;
			tabInfo[ind]->data.moment.y += val*j;
			tabInfo[ind]->data.moment.z += val*k;
			tabInfo[ind]->data.momentSurf.x += i;
			tabInfo[ind]->data.momentSurf.y += j;
			tabInfo[ind]->data.momentSurf.z += k;
			tabInfo[ind]->data.mean += val;
			tabInfo[ind]->data.var.xx += val*i*i;
			tabInfo[ind]->data.var.yy += val*j*j;
			tabInfo[ind]->data.var.zz += val*k*k;
			tabInfo[ind]->data.var.xy +=  val*i*j;
			tabInfo[ind]->data.var.xz +=  val*i*k;
			tabInfo[ind]->data.var.yz +=  val*k*j;

			if (fabs(val) > tabInfo[ind]->data.max.val){
				tabInfo[ind]->data.max.val = fabs(val);
				tabInfo[ind]->data.max.x = i;
				tabInfo[ind]->data.max.y = j;
				tabInfo[ind]->data.max.z = k;
			}

			if(i==0) 
				tabInfo[ind]->data.surfVol++;
			else
				if (SEG(i-1,j,k)!=SEG(i,j,k)) tabInfo[ind]->data.surfVol++;

			if(i==sizeDataX-1) 
				tabInfo[ind]->data.surfVol++;
			else
				if (SEG(i+1,j,k)!=SEG(i,j,k)) tabInfo[ind]->data.surfVol++;

			if(j==0) 
				tabInfo[ind]->data.surfVol++;
			else
				if (SEG(i,j-1,k)!=SEG(i,j,k)) tabInfo[ind]->data.surfVol++;

			if(j==sizeDataY-1) 
				tabInfo[ind]->data.surfVol++;   
			else
				if (SEG(i,j+1,k)!=SEG(i,j,k)) tabInfo[ind]->data.surfVol++;

			if(k==0) 
				tabInfo[ind]->data.surfVol++;
			else
				if (SEG(i,j,k-1)!=SEG(i,j,k)) tabInfo[ind]->data.surfVol++;

			if(k==sizeDataZ-1) 
				tabInfo[ind]->data.surfVol++;
			else
				if (SEG(i,j,k+1)!=SEG(i,j,k)) tabInfo[ind]->data.surfVol++;
		}

	double **c; //intercorrelation matrix
	double **v; //intercorrelation eigenvectors matrix 
	double *d;  //intercorrelation eigenvalues vector 
	// allocate intercorrelation matrix 
	c=matrix(1,3,1,3);
	// allocate output eigenvectors 
	v=matrix(1, 3, 1, 3);
	// allocate output eigenvalues 
	d=vector(1, 3);



	for(ind=1;ind<=nbInfoStruct;ind++){
		if (tabInfo[ind]->data.surfVol>1) 
			tabInfo[ind]->data.morpho = 6. * PI * (double)(tabInfo[ind]->data.volume) / (double)(tabInfo[ind]->data.surfVol) * 
				sqrt(PI/(double)(tabInfo[ind]->data.surfVol));
		else tabInfo[ind]->data.morpho=0.;

		tabInfo[ind]->data.moment.x		/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.moment.y		/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.moment.z		/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.momentSurf.x	/= tabInfo[ind]->data.volume;
		tabInfo[ind]->data.momentSurf.y	/= tabInfo[ind]->data.volume;
		tabInfo[ind]->data.momentSurf.z	/= tabInfo[ind]->data.volume;
		tabInfo[ind]->data.var.xx			/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.var.yy			/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.var.zz			/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.var.xy			/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.var.xz			/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.var.yz			/= tabInfo[ind]->data.mean;
		tabInfo[ind]->data.var.xx			-= tabInfo[ind]->data.moment.x*tabInfo[ind]->data.moment.x;
		tabInfo[ind]->data.var.yy			-= tabInfo[ind]->data.moment.y*tabInfo[ind]->data.moment.y;
		tabInfo[ind]->data.var.zz			-= tabInfo[ind]->data.moment.z*tabInfo[ind]->data.moment.z;
		tabInfo[ind]->data.var.xy			-= tabInfo[ind]->data.moment.x*tabInfo[ind]->data.moment.y;
		tabInfo[ind]->data.var.xz			-= tabInfo[ind]->data.moment.x*tabInfo[ind]->data.moment.z;
		tabInfo[ind]->data.var.yz			-= tabInfo[ind]->data.moment.y*tabInfo[ind]->data.moment.z;
		tabInfo[ind]->data.mean			/= tabInfo[ind]->data.volume;
		tabInfo[ind]->data.intensity=tabInfo[ind]->data.mean*tabInfo[ind]->data.volume;

		//test if the matrix of variance is a rotation matrix
		if( fabs(tabInfo[ind]->data.det=fabs(tabInfo[ind]->data.var.xx*tabInfo[ind]->data.var.yy*tabInfo[ind]->data.var.zz+2.*tabInfo[ind]->data.var.xy*tabInfo[ind]->data.var.yz*tabInfo[ind]->data.var.xz-tabInfo[ind]->data.var.xz*tabInfo[ind]->data.var.xz*tabInfo[ind]->data.var.yy-tabInfo[ind]->data.var.yz*tabInfo[ind]->data.var.yz*tabInfo[ind]->data.var.xx-tabInfo[ind]->data.var.xy*tabInfo[ind]->data.var.xy*tabInfo[ind]->data.var.zz))>1e-10 ){

			//powDet=1./pow(fabs(tabInfo[ind]->data.det),0.3333333);
			c[1][1]=tabInfo[ind]->data.var.xx;
			c[2][2]=tabInfo[ind]->data.var.yy;
			c[3][3]=tabInfo[ind]->data.var.zz;
			c[1][2]=tabInfo[ind]->data.var.xy;
			c[1][3]=tabInfo[ind]->data.var.xz;
			c[2][1]=tabInfo[ind]->data.var.xy;
			c[3][1]=tabInfo[ind]->data.var.xz;
			c[2][3]=tabInfo[ind]->data.var.yz;
			c[3][2]=tabInfo[ind]->data.var.yz;


			//compute the the eigenvalues and the eignenvectors
			jacobi(c, 3, d, v, &nrot);

			deter=v[1][1]*v[2][2]*v[3][3] + v[1][2]*v[2][3]*v[3][1] + v[1][3]*v[2][1]*v[3][2] 
				- v[3][1]*v[2][2]*v[1][3] - v[3][2]*v[2][3]*v[1][1] - v[3][3]*v[2][1]*v[1][2];
			powDet=1./pow(fabs(deter),0.3333333);


			//acording to the euler angle:
			tabInfo[ind]->data.angle.psi=atan(v[1][3]/v[2][3]); //!!!
			tabInfo[ind]->data.angle.teta=asin(v[1][3]*powDet/sin(tabInfo[ind]->data.angle.psi));//*powDet
			tabInfo[ind]->data.angle.phi=asin(v[3][1]*powDet/sin(tabInfo[ind]->data.angle.teta));//*powDet

			if(sqrt(fabs(d[1]))>sqrt(fabs(d[3]))){ //1>3
				if(sqrt(fabs(d[1]))>sqrt(fabs(d[2]))){ //1>2 & 1>3
					if(sqrt(fabs(d[2]))>sqrt(fabs(d[3]))){ //2>3 => 1>2>3
						tabInfo[ind]->data.sigma.z=sqrt(fabs(d[1]));
						tabInfo[ind]->data.sigma.y=sqrt(fabs(d[2]));
						tabInfo[ind]->data.sigma.x=sqrt(fabs(d[3]));
					}else{//3>2 & 1>3 =>1>3>2
						tabInfo[ind]->data.sigma.z=sqrt(fabs(d[1]));
						tabInfo[ind]->data.sigma.y=sqrt(fabs(d[3]));
						tabInfo[ind]->data.sigma.x=sqrt(fabs(d[2]));
					}
				}else{ //2>1 & 1>3 => 2>1>3
					tabInfo[ind]->data.sigma.z=sqrt(fabs(d[2]));			
					tabInfo[ind]->data.sigma.y=sqrt(fabs(d[1]));
					tabInfo[ind]->data.sigma.x=sqrt(fabs(d[3]));
				}
			}else{ //3>1
				if(sqrt(fabs(d[1]))>sqrt(fabs(d[2]))){ //1>2 => 3>1>2
					tabInfo[ind]->data.sigma.z=sqrt(fabs(d[3]));
					tabInfo[ind]->data.sigma.y=sqrt(fabs(d[1]));
					tabInfo[ind]->data.sigma.x=sqrt(fabs(d[2]));
				}else{ //2>1 & 3>1
					if(sqrt(fabs(d[2]))>sqrt(fabs(d[3]))){ //2>3 => 2>3>1
						tabInfo[ind]->data.sigma.z=sqrt(fabs(d[2]));			
						tabInfo[ind]->data.sigma.y=sqrt(fabs(d[3]));
						tabInfo[ind]->data.sigma.x=sqrt(fabs(d[1]));
					}else{ //3>2 => 3>2>1
						tabInfo[ind]->data.sigma.z=sqrt(fabs(d[3]));			
						tabInfo[ind]->data.sigma.y=sqrt(fabs(d[2]));
						tabInfo[ind]->data.sigma.x=sqrt(fabs(d[1]));
					}
				}
			} 
		}
		tabInfo[ind]->data.morpho = tabInfo[ind]->data.sigma.x*tabInfo[ind]->data.sigma.y/tabInfo[ind]->data.sigma.z/tabInfo[ind]->data.sigma.z;


	}
	free_vector(d,1,3);
	free_matrix(c,1,3,1,3);
	free_matrix(v,1,3,1,3);
	return(cptPoints);
}


unsigned short Segment3D::segmentation(USrtArray &data, USrtArray &seg, double val, bool diagonal, bool fillHoles){

	
	if(!fillHoles)
		return(segmentationPositiveValues(data, seg, val, diagonal));
	int ret;
	long i;
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);

	USrtArray *segHoles=new USrtArray;
	if((ret=segHoles->alloc(sizeDataX, sizeDataY, sizeDataZ))<0) return(ret);

	ret=segmentationNegativeValues(data, *segHoles, val, !diagonal);

	unsigned long *ptr;
	ptr=(unsigned long *)malloc(sizeof(long)*(ret+1));
	for(i=0;i<=ret;i++)
		ptr[i]=i;
	for(i=0;i<sizeDataY;i++){
		ptr[(*segHoles)(0,i)]=0;
		ptr[(*segHoles)(sizeDataX-1,i)]=0;
	}
	for(i=0;i<sizeDataX;i++){
		ptr[(*segHoles)(i,0)]=0;
		ptr[(*segHoles)(i,sizeDataY-1)]=0;
	}
	for(i=0;i<sizeDataX*sizeDataY;i++)
		segHoles->data[i]=ptr[segHoles->data[i]];
	ret=segmentationPositiveValues(data, seg, val, diagonal, *segHoles);
	delete segHoles;
	return(ret);
}

unsigned short Segment3D::segmentationNegativeValues(USrtArray &data, USrtArray &seg, double val, bool diagonal){
/**
	segmentation3D segment the 3D signal "data" which is lower than "val". It gives a map of segmentation "seg". 
	A label is set to each structure in the map seg. 
**/

	unsigned long *ptr;
	unsigned long *ptr2;
	long i,j,k;
	unsigned long nSeg;
	int structFinded;
	int mini;
	//int max;
	long sizePtr=MAX_TAB;
	
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);



	//allocation of ptr
	//ptr is an array of index to enumerate different structure
	//if ptr[i]=j then structure i and structure j are same.
	ptr=(unsigned long *)malloc(sizeof(long)*MAX_TAB);
	//initialisation of ptr
	for(i=0;i<MAX_TAB;i++)
		ptr[i]=(unsigned long)i;
	//nSeg indicate the number of the next unused structure
	nSeg=1;


	//initialisation of seg
	seg=(unsigned short)0;


	for(i=0;i<sizeDataX;i++)
	for(j=0;j<sizeDataY;j++)
	for(k=0;k<sizeDataZ;k++)
		if(DATA(i,j,k)<=val){
			structFinded=0;
			//test (i-1,j,k)
			if((i-1>=0))
				if( SEG(i-1,j,k)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i-1,j,k);
						structFinded=SEG(i-1,j,k);
					}else if(structFinded!=SEG(i-1,j,k)){
						mini=min(structFinded,SEG(i-1,j,k),ptr[structFinded], ptr[SEG(i-1,j,k)]);
						ptr[structFinded]=mini;
						ptr[SEG(i-1,j,k)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}
			//test (i,j,k-1)
			if(k-1>=0)
				if( SEG(i,j,k-1)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i,j,k-1);
						structFinded=SEG(i,j,k-1);
					}else if(structFinded!=SEG(i,j,k-1)){
						mini=min(structFinded,SEG(i,j,k-1),ptr[structFinded], ptr[SEG(i,j,k-1)]);
						ptr[structFinded]=mini;
						ptr[SEG(i,j,k-1)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}
			//test (i,j-1,k)
			if(j-1>=0)
				if( SEG(i,j-1,k)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i,j-1,k);
						structFinded=SEG(i,j-1,k);
					}else if(structFinded!=SEG(i,j-1,k)){
						mini=min(structFinded,SEG(i,j-1,k),ptr[structFinded], ptr[SEG(i,j-1,k)]);
						ptr[structFinded]=mini;
						ptr[SEG(i,j-1,k)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}

			if(diagonal){
				//test (i-1,j-1,k-1)
				if((i-1>=0)  && (j-1>=0)  &&  (k-1>=0)){
					
					if( SEG(i-1,j-1,k-1)!=0 ){
						SEG(i,j,k)=SEG(i-1,j-1,k-1);
						structFinded=SEG(i-1,j-1,k-1);
					}
				}
				//test (i-1,j,k-1)
				if((i-1>=0)  &&  (k-1>=0))
					if( SEG(i-1,j,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j,k-1);
							structFinded=SEG(i-1,j,k-1);
						}else if(structFinded!=SEG(i-1,j,k-1)){
							mini=min(structFinded,SEG(i-1,j,k-1),ptr[structFinded], ptr[SEG(i-1,j,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}

				//test (i-1,j+1,k-1)
				if((i-1>=0)  && (j+1<sizeDataY)  &&  (k-1>=0))
					if( SEG(i-1,j+1,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j+1,k-1);
							structFinded=SEG(i-1,j+1,k-1);
						}else if(structFinded!=SEG(i-1,j+1,k-1)){
							mini=min(structFinded,SEG(i-1,j+1,k-1),ptr[structFinded], ptr[SEG(i-1,j+1,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j+1,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}
					
				//test (i-1,j-1,k)
				if((i-1>=0)  && (j-1>=0))
					if( SEG(i-1,j-1,k)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j-1,k);
							structFinded=SEG(i-1,j-1,k);
						}else if(structFinded!=SEG(i-1,j-1,k)){
							mini=min(structFinded,SEG(i-1,j-1,k),ptr[structFinded], ptr[SEG(i-1,j-1,k)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j-1,k)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	

				//test (i-1,j,k+1)
				if((i-1>=0)  &&  (k+1<sizeDataZ))
					if( SEG(i-1,j,k+1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j,k+1);
							structFinded=SEG(i-1,j,k+1);
						}else if(structFinded!=SEG(i-1,j,k+1)){
							mini=min(structFinded,SEG(i-1,j,k+1),ptr[structFinded], ptr[SEG(i-1,j,k+1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j,k+1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	

				//test (i,j-1,k-1)
				if((j-1>=0)  &&  (k-1>=0))
					if( SEG(i,j-1,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i,j-1,k-1);
							structFinded=SEG(i,j-1,k-1);
						}else if(structFinded!=SEG(i,j-1,k-1)){
							mini=min(structFinded,SEG(i,j-1,k-1),ptr[structFinded], ptr[SEG(i,j-1,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i,j-1,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}

				//test (i,j-1,k+1)
				if((j-1>=0)  &&  (k+1<sizeDataZ))
					if( SEG(i,j-1,k+1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i,j-1,k+1);
							structFinded=SEG(i,j-1,k+1);
						}else if(structFinded!=SEG(i,j-1,k+1)){
							mini=min(structFinded,SEG(i,j-1,k+1),ptr[structFinded], ptr[SEG(i,j-1,k+1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i,j-1,k+1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	
			}

			if(structFinded==0){
				SEG(i,j,k)=(unsigned short)nSeg;
				inc(nSeg, &ptr, sizePtr);
			}
		}


	//assembling structures

	for(i=1;i<(long)nSeg;i++)
		if(ptr[i]!=(unsigned long)i)
			ptr[i]=ptr[ptr[i]];

	ptr2=(unsigned long *)malloc(sizeof(long)*nSeg);
	for(i=0;i<(long)nSeg;i++)
		ptr2[i]=0;
	structFinded=1;
	for(i=1;i<(long)nSeg;i++)
		if(ptr2[ptr[i]]==0){
			ptr2[ptr[i]]=structFinded;
			structFinded++;
		}


	for(i=0;i<sizeDataX*sizeDataY*sizeDataZ;i++)
		seg.data[i]=(unsigned short)ptr2[ptr[seg.data[i]]];

	free(ptr);
	free(ptr2);
	nbInfoStruct=structFinded-1;
	mixIndex(seg);
	return(structFinded-1);
}

unsigned short Segment3D::segmentationPositiveValues(USrtArray &data, USrtArray &seg, double val, bool diagonal){
/**
	segmentation3D segment the 3D signal "data" which is upper than "val". It gives a map of segmentation "seg". 
	A label is set to each structure in the map seg. 
**/

	unsigned long *ptr;
	unsigned long *ptr2;
	long i,j,k;
	unsigned long nSeg;
	int structFinded;
	int mini;
	//int max;
	long sizePtr=MAX_TAB;
	
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);

	//allocation of ptr
	//ptr is an array of index to enumerate different structure
	//if ptr[i]=j then structure i and structure j are same.
	ptr=(unsigned long *)malloc(sizeof(long)*MAX_TAB);
	//initialisation of ptr
	for(i=0;i<MAX_TAB;i++)
		ptr[i]=(unsigned long)i;
	//nSeg indicate the number of the next unused structure
	nSeg=1;


	//initialisation of seg
	seg=(unsigned short)0;


	for(i=0;i<sizeDataX;i++)
	for(j=0;j<sizeDataY;j++)
	for(k=0;k<sizeDataZ;k++)
		if(DATA(i,j,k)>val){
			structFinded=0;
			//test (i-1,j,k)
			if((i-1>=0))
				if( SEG(i-1,j,k)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i-1,j,k);
						structFinded=SEG(i-1,j,k);
					}else if(structFinded!=SEG(i-1,j,k)){
						mini=min(structFinded,SEG(i-1,j,k),ptr[structFinded], ptr[SEG(i-1,j,k)]);
						ptr[structFinded]=mini;
						ptr[SEG(i-1,j,k)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}
			//test (i,j,k-1)
			if(k-1>=0)
				if( SEG(i,j,k-1)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i,j,k-1);
						structFinded=SEG(i,j,k-1);
					}else if(structFinded!=SEG(i,j,k-1)){
						mini=min(structFinded,SEG(i,j,k-1),ptr[structFinded], ptr[SEG(i,j,k-1)]);
						ptr[structFinded]=mini;
						ptr[SEG(i,j,k-1)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}
			//test (i,j-1,k)
			if(j-1>=0)
				if( SEG(i,j-1,k)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i,j-1,k);
						structFinded=SEG(i,j-1,k);
					}else if(structFinded!=SEG(i,j-1,k)){
						mini=min(structFinded,SEG(i,j-1,k),ptr[structFinded], ptr[SEG(i,j-1,k)]);
						ptr[structFinded]=mini;
						ptr[SEG(i,j-1,k)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}

			if(diagonal){
				//test (i-1,j-1,k-1)
				if((i-1>=0)  && (j-1>=0)  &&  (k-1>=0)){
					
					if( SEG(i-1,j-1,k-1)!=0 ){
						SEG(i,j,k)=SEG(i-1,j-1,k-1);
						structFinded=SEG(i-1,j-1,k-1);
					}
				}
				//test (i-1,j,k-1)
				if((i-1>=0)  &&  (k-1>=0))
					if( SEG(i-1,j,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j,k-1);
							structFinded=SEG(i-1,j,k-1);
						}else if(structFinded!=SEG(i-1,j,k-1)){
							mini=min(structFinded,SEG(i-1,j,k-1),ptr[structFinded], ptr[SEG(i-1,j,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}

				//test (i-1,j+1,k-1)
				if((i-1>=0)  && (j+1<sizeDataY)  &&  (k-1>=0))
					if( SEG(i-1,j+1,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j+1,k-1);
							structFinded=SEG(i-1,j+1,k-1);
						}else if(structFinded!=SEG(i-1,j+1,k-1)){
							mini=min(structFinded,SEG(i-1,j+1,k-1),ptr[structFinded], ptr[SEG(i-1,j+1,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j+1,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}
					
				//test (i-1,j-1,k)
				if((i-1>=0)  && (j-1>=0))
					if( SEG(i-1,j-1,k)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j-1,k);
							structFinded=SEG(i-1,j-1,k);
						}else if(structFinded!=SEG(i-1,j-1,k)){
							mini=min(structFinded,SEG(i-1,j-1,k),ptr[structFinded], ptr[SEG(i-1,j-1,k)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j-1,k)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	

				//test (i-1,j,k+1)
				if((i-1>=0)  &&  (k+1<sizeDataZ))
					if( SEG(i-1,j,k+1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j,k+1);
							structFinded=SEG(i-1,j,k+1);
						}else if(structFinded!=SEG(i-1,j,k+1)){
							mini=min(structFinded,SEG(i-1,j,k+1),ptr[structFinded], ptr[SEG(i-1,j,k+1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j,k+1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	

				//test (i,j-1,k-1)
				if((j-1>=0)  &&  (k-1>=0))
					if( SEG(i,j-1,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i,j-1,k-1);
							structFinded=SEG(i,j-1,k-1);
						}else if(structFinded!=SEG(i,j-1,k-1)){
							mini=min(structFinded,SEG(i,j-1,k-1),ptr[structFinded], ptr[SEG(i,j-1,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i,j-1,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}

				//test (i,j-1,k+1)
				if((j-1>=0)  &&  (k+1<sizeDataZ))
					if( SEG(i,j-1,k+1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i,j-1,k+1);
							structFinded=SEG(i,j-1,k+1);
						}else if(structFinded!=SEG(i,j-1,k+1)){
							mini=min(structFinded,SEG(i,j-1,k+1),ptr[structFinded], ptr[SEG(i,j-1,k+1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i,j-1,k+1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	
			}

			if(structFinded==0){
				SEG(i,j,k)=(unsigned short)nSeg;
				inc(nSeg, &ptr, sizePtr);
			}
		}


	//assembling structures

	for(i=1;i<(long)nSeg;i++)
		if(ptr[i]!=(unsigned long)i)
			ptr[i]=ptr[ptr[i]];

	ptr2=(unsigned long *)malloc(sizeof(long)*nSeg);
	for(i=0;i<(long)nSeg;i++)
		ptr2[i]=0;
	structFinded=1;
	for(i=1;i<(long)nSeg;i++)
		if(ptr2[ptr[i]]==0){
			ptr2[ptr[i]]=structFinded;
			structFinded++;
		}


	for(i=0;i<sizeDataX*sizeDataY*sizeDataZ;i++)
		seg.data[i]=(unsigned short)ptr2[ptr[seg.data[i]]];

	free(ptr);
	free(ptr2);
	nbInfoStruct=structFinded-1;
	mixIndex(seg);
	return(structFinded-1);
}

unsigned short Segment3D::segmentationPositiveValues(USrtArray &data, USrtArray &seg, double val, bool diagonal, USrtArray &segHoles){
/**
	segmentation3D segment the 3D signal "data" which is upper than "val". It gives a map of segmentation "seg". 
	A label is set to each structure in the map seg. 
**/

		unsigned long *ptr;
	unsigned long *ptr2;
	long i,j,k;
	unsigned long nSeg;
	int structFinded;
	int mini;
	//int max;
	long sizePtr=MAX_TAB;
	
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);



	//allocation of ptr
	//ptr is an array of index to enumerate different structure
	//if ptr[i]=j then structure i and structure j are same.
	ptr=(unsigned long *)malloc(sizeof(long)*MAX_TAB);
	//initialisation of ptr
	for(i=0;i<MAX_TAB;i++)
		ptr[i]=(unsigned long)i;
	//nSeg indicate the number of the next unused structure
	nSeg=1;


	//initialisation of seg
	seg=(unsigned short)0;


	for(i=0;i<sizeDataX;i++)
	for(j=0;j<sizeDataY;j++)
	for(k=0;k<sizeDataZ;k++)
		if((DATA(i,j,k)>val)||(SEGHoles(i,j,k)!=0)){
			structFinded=0;
			//test (i-1,j,k)
			if((i-1>=0))
				if( SEG(i-1,j,k)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i-1,j,k);
						structFinded=SEG(i-1,j,k);
					}else if(structFinded!=SEG(i-1,j,k)){
						mini=min(structFinded,SEG(i-1,j,k),ptr[structFinded], ptr[SEG(i-1,j,k)]);
						ptr[structFinded]=mini;
						ptr[SEG(i-1,j,k)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}
			//test (i,j,k-1)
			if(k-1>=0)
				if( SEG(i,j,k-1)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i,j,k-1);
						structFinded=SEG(i,j,k-1);
					}else if(structFinded!=SEG(i,j,k-1)){
						mini=min(structFinded,SEG(i,j,k-1),ptr[structFinded], ptr[SEG(i,j,k-1)]);
						ptr[structFinded]=mini;
						ptr[SEG(i,j,k-1)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}
			//test (i,j-1,k)
			if(j-1>=0)
				if( SEG(i,j-1,k)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i,j-1,k);
						structFinded=SEG(i,j-1,k);
					}else if(structFinded!=SEG(i,j-1,k)){
						mini=min(structFinded,SEG(i,j-1,k),ptr[structFinded], ptr[SEG(i,j-1,k)]);
						ptr[structFinded]=mini;
						ptr[SEG(i,j-1,k)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}

			if(diagonal){
				//test (i-1,j-1,k-1)
				if((i-1>=0)  && (j-1>=0)  &&  (k-1>=0)){
					
					if( SEG(i-1,j-1,k-1)!=0 ){
						SEG(i,j,k)=SEG(i-1,j-1,k-1);
						structFinded=SEG(i-1,j-1,k-1);
					}
				}
				//test (i-1,j,k-1)
				if((i-1>=0)  &&  (k-1>=0))
					if( SEG(i-1,j,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j,k-1);
							structFinded=SEG(i-1,j,k-1);
						}else if(structFinded!=SEG(i-1,j,k-1)){
							mini=min(structFinded,SEG(i-1,j,k-1),ptr[structFinded], ptr[SEG(i-1,j,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}

				//test (i-1,j+1,k-1)
				if((i-1>=0)  && (j+1<sizeDataY)  &&  (k-1>=0))
					if( SEG(i-1,j+1,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j+1,k-1);
							structFinded=SEG(i-1,j+1,k-1);
						}else if(structFinded!=SEG(i-1,j+1,k-1)){
							mini=min(structFinded,SEG(i-1,j+1,k-1),ptr[structFinded], ptr[SEG(i-1,j+1,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j+1,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}
					
				//test (i-1,j-1,k)
				if((i-1>=0)  && (j-1>=0))
					if( SEG(i-1,j-1,k)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j-1,k);
							structFinded=SEG(i-1,j-1,k);
						}else if(structFinded!=SEG(i-1,j-1,k)){
							mini=min(structFinded,SEG(i-1,j-1,k),ptr[structFinded], ptr[SEG(i-1,j-1,k)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j-1,k)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	

				//test (i-1,j,k+1)
				if((i-1>=0)  &&  (k+1<sizeDataZ))
					if( SEG(i-1,j,k+1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j,k+1);
							structFinded=SEG(i-1,j,k+1);
						}else if(structFinded!=SEG(i-1,j,k+1)){
							mini=min(structFinded,SEG(i-1,j,k+1),ptr[structFinded], ptr[SEG(i-1,j,k+1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j,k+1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	

				//test (i,j-1,k-1)
				if((j-1>=0)  &&  (k-1>=0))
					if( SEG(i,j-1,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i,j-1,k-1);
							structFinded=SEG(i,j-1,k-1);
						}else if(structFinded!=SEG(i,j-1,k-1)){
							mini=min(structFinded,SEG(i,j-1,k-1),ptr[structFinded], ptr[SEG(i,j-1,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i,j-1,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}

				//test (i,j-1,k+1)
				if((j-1>=0)  &&  (k+1<sizeDataZ))
					if( SEG(i,j-1,k+1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i,j-1,k+1);
							structFinded=SEG(i,j-1,k+1);
						}else if(structFinded!=SEG(i,j-1,k+1)){
							mini=min(structFinded,SEG(i,j-1,k+1),ptr[structFinded], ptr[SEG(i,j-1,k+1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i,j-1,k+1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	
			}

			if(structFinded==0){
				SEG(i,j,k)=(unsigned short)nSeg;
				inc(nSeg, &ptr, sizePtr);
			}
		}


	//assembling structures

	for(i=1;i<(long)nSeg;i++)
		if(ptr[i]!=(unsigned long)i)
			ptr[i]=ptr[ptr[i]];

	ptr2=(unsigned long *)malloc(sizeof(long)*nSeg);
	for(i=0;i<(long)nSeg;i++)
		ptr2[i]=0;
	structFinded=1;
	for(i=1;i<(long)nSeg;i++)
		if(ptr2[ptr[i]]==0){
			ptr2[ptr[i]]=structFinded;
			structFinded++;
		}


	for(i=0;i<sizeDataX*sizeDataY*sizeDataZ;i++)
		seg.data[i]=(unsigned short)ptr2[ptr[seg.data[i]]];

	free(ptr);
	free(ptr2);
	nbInfoStruct=structFinded-1;
	mixIndex(seg);
	return(structFinded-1);
}

unsigned short Segment3D::segmentation(FltArray &data, USrtArray &seg, double val, bool diagonal, bool fillHoles){

	
	if(!fillHoles)
		return(segmentationPositiveValues(data, seg, val, diagonal));
	int ret;
	long i;
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);

	USrtArray *segHoles=new USrtArray;
	if((ret=segHoles->alloc(sizeDataX, sizeDataY, sizeDataZ))<0) return(ret);

	ret=segmentationNegativeValues(data, *segHoles, val, !diagonal);

	unsigned long *ptr;
	ptr=(unsigned long *)malloc(sizeof(long)*(ret+1));
	for(i=0;i<=ret;i++)
		ptr[i]=i;
	for(i=0;i<sizeDataY;i++){
		ptr[(*segHoles)(0,i)]=0;
		ptr[(*segHoles)(sizeDataX-1,i)]=0;
	}
	for(i=0;i<sizeDataX;i++){
		ptr[(*segHoles)(i,0)]=0;
		ptr[(*segHoles)(i,sizeDataY-1)]=0;
	}
	for(i=0;i<sizeDataX*sizeDataY;i++)
		segHoles->data[i]=ptr[segHoles->data[i]];
	ret=segmentationPositiveValues(data, seg, val, diagonal, *segHoles);
	delete segHoles;
	return(ret);
}

unsigned short Segment3D::segmentationNegativeValues(FltArray &data, USrtArray &seg, double val, bool diagonal){
/**
	segmentation3D segment the 3D signal "data" which is lower than "val". It gives a map of segmentation "seg". 
	A label is set to each structure in the map seg. 
**/

	unsigned long *ptr;
	unsigned long *ptr2;
	long i,j,k;
	unsigned long nSeg;
	int structFinded;
	int mini;
	//int max;
	long sizePtr=MAX_TAB;
	
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);



	//allocation of ptr
	//ptr is an array of index to enumerate different structure
	//if ptr[i]=j then structure i and structure j are same.
	ptr=(unsigned long *)malloc(sizeof(long)*MAX_TAB);
	//initialisation of ptr
	for(i=0;i<MAX_TAB;i++)
		ptr[i]=(unsigned long)i;
	//nSeg indicate the number of the next unused structure
	nSeg=1;


	//initialisation of seg
	seg=(unsigned short)0;


	for(i=0;i<sizeDataX;i++)
	for(j=0;j<sizeDataY;j++)
	for(k=0;k<sizeDataZ;k++)
		if(DATA(i,j,k)<=val){
			structFinded=0;
			//test (i-1,j,k)
			if((i-1>=0))
				if( SEG(i-1,j,k)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i-1,j,k);
						structFinded=SEG(i-1,j,k);
					}else if(structFinded!=SEG(i-1,j,k)){
						mini=min(structFinded,SEG(i-1,j,k),ptr[structFinded], ptr[SEG(i-1,j,k)]);
						ptr[structFinded]=mini;
						ptr[SEG(i-1,j,k)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}
			//test (i,j,k-1)
			if(k-1>=0)
				if( SEG(i,j,k-1)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i,j,k-1);
						structFinded=SEG(i,j,k-1);
					}else if(structFinded!=SEG(i,j,k-1)){
						mini=min(structFinded,SEG(i,j,k-1),ptr[structFinded], ptr[SEG(i,j,k-1)]);
						ptr[structFinded]=mini;
						ptr[SEG(i,j,k-1)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}
			//test (i,j-1,k)
			if(j-1>=0)
				if( SEG(i,j-1,k)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i,j-1,k);
						structFinded=SEG(i,j-1,k);
					}else if(structFinded!=SEG(i,j-1,k)){
						mini=min(structFinded,SEG(i,j-1,k),ptr[structFinded], ptr[SEG(i,j-1,k)]);
						ptr[structFinded]=mini;
						ptr[SEG(i,j-1,k)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}

			if(diagonal){
				//test (i-1,j-1,k-1)
				if((i-1>=0)  && (j-1>=0)  &&  (k-1>=0)){
					
					if( SEG(i-1,j-1,k-1)!=0 ){
						SEG(i,j,k)=SEG(i-1,j-1,k-1);
						structFinded=SEG(i-1,j-1,k-1);
					}
				}
				//test (i-1,j,k-1)
				if((i-1>=0)  &&  (k-1>=0))
					if( SEG(i-1,j,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j,k-1);
							structFinded=SEG(i-1,j,k-1);
						}else if(structFinded!=SEG(i-1,j,k-1)){
							mini=min(structFinded,SEG(i-1,j,k-1),ptr[structFinded], ptr[SEG(i-1,j,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}

				//test (i-1,j+1,k-1)
				if((i-1>=0)  && (j+1<sizeDataY)  &&  (k-1>=0))
					if( SEG(i-1,j+1,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j+1,k-1);
							structFinded=SEG(i-1,j+1,k-1);
						}else if(structFinded!=SEG(i-1,j+1,k-1)){
							mini=min(structFinded,SEG(i-1,j+1,k-1),ptr[structFinded], ptr[SEG(i-1,j+1,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j+1,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}
					
				//test (i-1,j-1,k)
				if((i-1>=0)  && (j-1>=0))
					if( SEG(i-1,j-1,k)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j-1,k);
							structFinded=SEG(i-1,j-1,k);
						}else if(structFinded!=SEG(i-1,j-1,k)){
							mini=min(structFinded,SEG(i-1,j-1,k),ptr[structFinded], ptr[SEG(i-1,j-1,k)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j-1,k)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	

				//test (i-1,j,k+1)
				if((i-1>=0)  &&  (k+1<sizeDataZ))
					if( SEG(i-1,j,k+1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j,k+1);
							structFinded=SEG(i-1,j,k+1);
						}else if(structFinded!=SEG(i-1,j,k+1)){
							mini=min(structFinded,SEG(i-1,j,k+1),ptr[structFinded], ptr[SEG(i-1,j,k+1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j,k+1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	

				//test (i,j-1,k-1)
				if((j-1>=0)  &&  (k-1>=0))
					if( SEG(i,j-1,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i,j-1,k-1);
							structFinded=SEG(i,j-1,k-1);
						}else if(structFinded!=SEG(i,j-1,k-1)){
							mini=min(structFinded,SEG(i,j-1,k-1),ptr[structFinded], ptr[SEG(i,j-1,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i,j-1,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}

				//test (i,j-1,k+1)
				if((j-1>=0)  &&  (k+1<sizeDataZ))
					if( SEG(i,j-1,k+1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i,j-1,k+1);
							structFinded=SEG(i,j-1,k+1);
						}else if(structFinded!=SEG(i,j-1,k+1)){
							mini=min(structFinded,SEG(i,j-1,k+1),ptr[structFinded], ptr[SEG(i,j-1,k+1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i,j-1,k+1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	
			}

			if(structFinded==0){
				SEG(i,j,k)=(unsigned short)nSeg;
				inc(nSeg, &ptr, sizePtr);
			}
		}


	//assembling structures

	for(i=1;i<(long)nSeg;i++)
		if(ptr[i]!=(unsigned long)i)
			ptr[i]=ptr[ptr[i]];

	ptr2=(unsigned long *)malloc(sizeof(long)*nSeg);
	for(i=0;i<(long)nSeg;i++)
		ptr2[i]=0;
	structFinded=1;
	for(i=1;i<(long)nSeg;i++)
		if(ptr2[ptr[i]]==0){
			ptr2[ptr[i]]=structFinded;
			structFinded++;
		}


	for(i=0;i<sizeDataX*sizeDataY*sizeDataZ;i++)
		seg.data[i]=(unsigned short)ptr2[ptr[seg.data[i]]];

	free(ptr);
	free(ptr2);
	nbInfoStruct=structFinded-1;
	mixIndex(seg);
	return(structFinded-1);
}

unsigned short Segment3D::segmentationPositiveValues(FltArray &data, USrtArray &seg, double val, bool diagonal){
/**
	segmentation3D segment the 3D signal "data" which is upper than "val". It gives a map of segmentation "seg". 
	A label is set to each structure in the map seg. 
**/

	unsigned long *ptr;
	unsigned long *ptr2;
	long i,j,k;
	unsigned long nSeg;
	int structFinded;
	int mini;
	//int max;
	long sizePtr=MAX_TAB;
	
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);

	//allocation of ptr
	//ptr is an array of index to enumerate different structure
	//if ptr[i]=j then structure i and structure j are same.
	ptr=(unsigned long *)malloc(sizeof(long)*MAX_TAB);
	//initialisation of ptr
	for(i=0;i<MAX_TAB;i++)
		ptr[i]=(unsigned long)i;
	//nSeg indicate the number of the next unused structure
	nSeg=1;


	//initialisation of seg
	seg=(unsigned short)0;


	for(i=0;i<sizeDataX;i++)
	for(j=0;j<sizeDataY;j++)
	for(k=0;k<sizeDataZ;k++)
		if(DATA(i,j,k)>val){
			structFinded=0;
			//test (i-1,j,k)
			if((i-1>=0))
				if( SEG(i-1,j,k)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i-1,j,k);
						structFinded=SEG(i-1,j,k);
					}else if(structFinded!=SEG(i-1,j,k)){
						mini=min(structFinded,SEG(i-1,j,k),ptr[structFinded], ptr[SEG(i-1,j,k)]);
						ptr[structFinded]=mini;
						ptr[SEG(i-1,j,k)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}
			//test (i,j,k-1)
			if(k-1>=0)
				if( SEG(i,j,k-1)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i,j,k-1);
						structFinded=SEG(i,j,k-1);
					}else if(structFinded!=SEG(i,j,k-1)){
						mini=min(structFinded,SEG(i,j,k-1),ptr[structFinded], ptr[SEG(i,j,k-1)]);
						ptr[structFinded]=mini;
						ptr[SEG(i,j,k-1)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}
			//test (i,j-1,k)
			if(j-1>=0)
				if( SEG(i,j-1,k)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i,j-1,k);
						structFinded=SEG(i,j-1,k);
					}else if(structFinded!=SEG(i,j-1,k)){
						mini=min(structFinded,SEG(i,j-1,k),ptr[structFinded], ptr[SEG(i,j-1,k)]);
						ptr[structFinded]=mini;
						ptr[SEG(i,j-1,k)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}

			if(diagonal){
				//test (i-1,j-1,k-1)
				if((i-1>=0)  && (j-1>=0)  &&  (k-1>=0)){
					
					if( SEG(i-1,j-1,k-1)!=0 ){
						SEG(i,j,k)=SEG(i-1,j-1,k-1);
						structFinded=SEG(i-1,j-1,k-1);
					}
				}
				//test (i-1,j,k-1)
				if((i-1>=0)  &&  (k-1>=0))
					if( SEG(i-1,j,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j,k-1);
							structFinded=SEG(i-1,j,k-1);
						}else if(structFinded!=SEG(i-1,j,k-1)){
							mini=min(structFinded,SEG(i-1,j,k-1),ptr[structFinded], ptr[SEG(i-1,j,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}

				//test (i-1,j+1,k-1)
				if((i-1>=0)  && (j+1<sizeDataY)  &&  (k-1>=0))
					if( SEG(i-1,j+1,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j+1,k-1);
							structFinded=SEG(i-1,j+1,k-1);
						}else if(structFinded!=SEG(i-1,j+1,k-1)){
							mini=min(structFinded,SEG(i-1,j+1,k-1),ptr[structFinded], ptr[SEG(i-1,j+1,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j+1,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}
					
				//test (i-1,j-1,k)
				if((i-1>=0)  && (j-1>=0))
					if( SEG(i-1,j-1,k)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j-1,k);
							structFinded=SEG(i-1,j-1,k);
						}else if(structFinded!=SEG(i-1,j-1,k)){
							mini=min(structFinded,SEG(i-1,j-1,k),ptr[structFinded], ptr[SEG(i-1,j-1,k)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j-1,k)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	

				//test (i-1,j,k+1)
				if((i-1>=0)  &&  (k+1<sizeDataZ))
					if( SEG(i-1,j,k+1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j,k+1);
							structFinded=SEG(i-1,j,k+1);
						}else if(structFinded!=SEG(i-1,j,k+1)){
							mini=min(structFinded,SEG(i-1,j,k+1),ptr[structFinded], ptr[SEG(i-1,j,k+1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j,k+1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	

				//test (i,j-1,k-1)
				if((j-1>=0)  &&  (k-1>=0))
					if( SEG(i,j-1,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i,j-1,k-1);
							structFinded=SEG(i,j-1,k-1);
						}else if(structFinded!=SEG(i,j-1,k-1)){
							mini=min(structFinded,SEG(i,j-1,k-1),ptr[structFinded], ptr[SEG(i,j-1,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i,j-1,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}

				//test (i,j-1,k+1)
				if((j-1>=0)  &&  (k+1<sizeDataZ))
					if( SEG(i,j-1,k+1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i,j-1,k+1);
							structFinded=SEG(i,j-1,k+1);
						}else if(structFinded!=SEG(i,j-1,k+1)){
							mini=min(structFinded,SEG(i,j-1,k+1),ptr[structFinded], ptr[SEG(i,j-1,k+1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i,j-1,k+1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	
			}

			if(structFinded==0){
				SEG(i,j,k)=(unsigned short)nSeg;
				inc(nSeg, &ptr, sizePtr);
			}
		}


	//assembling structures

	for(i=1;i<(long)nSeg;i++)
		if(ptr[i]!=(unsigned long)i)
			ptr[i]=ptr[ptr[i]];

	ptr2=(unsigned long *)malloc(sizeof(long)*nSeg);
	for(i=0;i<(long)nSeg;i++)
		ptr2[i]=0;
	structFinded=1;
	for(i=1;i<(long)nSeg;i++)
		if(ptr2[ptr[i]]==0){
			ptr2[ptr[i]]=structFinded;
			structFinded++;
		}


	for(i=0;i<sizeDataX*sizeDataY*sizeDataZ;i++)
		seg.data[i]=(unsigned short)ptr2[ptr[seg.data[i]]];

	free(ptr);
	free(ptr2);
	nbInfoStruct=structFinded-1;
	mixIndex(seg);
	return(structFinded-1);
}

unsigned short Segment3D::segmentationPositiveValues(FltArray &data, USrtArray &seg, double val, bool diagonal, USrtArray &segHoles){
/**
	segmentation3D segment the 3D signal "data" which is upper than "val". It gives a map of segmentation "seg". 
	A label is set to each structure in the map seg. 
**/

		unsigned long *ptr;
	unsigned long *ptr2;
	long i,j,k;
	unsigned long nSeg;
	int structFinded;
	int mini;
	//int max;
	long sizePtr=MAX_TAB;
	
	long sizeDataX;
	long sizeDataY;
	long sizeDataZ;
	getSize(sizeDataX, sizeDataY, sizeDataZ, seg);



	//allocation of ptr
	//ptr is an array of index to enumerate different structure
	//if ptr[i]=j then structure i and structure j are same.
	ptr=(unsigned long *)malloc(sizeof(long)*MAX_TAB);
	//initialisation of ptr
	for(i=0;i<MAX_TAB;i++)
		ptr[i]=(unsigned long)i;
	//nSeg indicate the number of the next unused structure
	nSeg=1;


	//initialisation of seg
	seg=(unsigned short)0;


	for(i=0;i<sizeDataX;i++)
	for(j=0;j<sizeDataY;j++)
	for(k=0;k<sizeDataZ;k++)
		if((DATA(i,j,k)>val)||(SEGHoles(i,j,k)!=0)){
			structFinded=0;
			//test (i-1,j,k)
			if((i-1>=0))
				if( SEG(i-1,j,k)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i-1,j,k);
						structFinded=SEG(i-1,j,k);
					}else if(structFinded!=SEG(i-1,j,k)){
						mini=min(structFinded,SEG(i-1,j,k),ptr[structFinded], ptr[SEG(i-1,j,k)]);
						ptr[structFinded]=mini;
						ptr[SEG(i-1,j,k)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}
			//test (i,j,k-1)
			if(k-1>=0)
				if( SEG(i,j,k-1)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i,j,k-1);
						structFinded=SEG(i,j,k-1);
					}else if(structFinded!=SEG(i,j,k-1)){
						mini=min(structFinded,SEG(i,j,k-1),ptr[structFinded], ptr[SEG(i,j,k-1)]);
						ptr[structFinded]=mini;
						ptr[SEG(i,j,k-1)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}
			//test (i,j-1,k)
			if(j-1>=0)
				if( SEG(i,j-1,k)!=0 ){
					if(structFinded==0){
						SEG(i,j,k)=SEG(i,j-1,k);
						structFinded=SEG(i,j-1,k);
					}else if(structFinded!=SEG(i,j-1,k)){
						mini=min(structFinded,SEG(i,j-1,k),ptr[structFinded], ptr[SEG(i,j-1,k)]);
						ptr[structFinded]=mini;
						ptr[SEG(i,j-1,k)]=mini;	
						SEG(i,j,k)=mini;
						structFinded=mini;
					}
				}

			if(diagonal){
				//test (i-1,j-1,k-1)
				if((i-1>=0)  && (j-1>=0)  &&  (k-1>=0)){
					
					if( SEG(i-1,j-1,k-1)!=0 ){
						SEG(i,j,k)=SEG(i-1,j-1,k-1);
						structFinded=SEG(i-1,j-1,k-1);
					}
				}
				//test (i-1,j,k-1)
				if((i-1>=0)  &&  (k-1>=0))
					if( SEG(i-1,j,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j,k-1);
							structFinded=SEG(i-1,j,k-1);
						}else if(structFinded!=SEG(i-1,j,k-1)){
							mini=min(structFinded,SEG(i-1,j,k-1),ptr[structFinded], ptr[SEG(i-1,j,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}

				//test (i-1,j+1,k-1)
				if((i-1>=0)  && (j+1<sizeDataY)  &&  (k-1>=0))
					if( SEG(i-1,j+1,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j+1,k-1);
							structFinded=SEG(i-1,j+1,k-1);
						}else if(structFinded!=SEG(i-1,j+1,k-1)){
							mini=min(structFinded,SEG(i-1,j+1,k-1),ptr[structFinded], ptr[SEG(i-1,j+1,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j+1,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}
					
				//test (i-1,j-1,k)
				if((i-1>=0)  && (j-1>=0))
					if( SEG(i-1,j-1,k)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j-1,k);
							structFinded=SEG(i-1,j-1,k);
						}else if(structFinded!=SEG(i-1,j-1,k)){
							mini=min(structFinded,SEG(i-1,j-1,k),ptr[structFinded], ptr[SEG(i-1,j-1,k)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j-1,k)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	

				//test (i-1,j,k+1)
				if((i-1>=0)  &&  (k+1<sizeDataZ))
					if( SEG(i-1,j,k+1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i-1,j,k+1);
							structFinded=SEG(i-1,j,k+1);
						}else if(structFinded!=SEG(i-1,j,k+1)){
							mini=min(structFinded,SEG(i-1,j,k+1),ptr[structFinded], ptr[SEG(i-1,j,k+1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i-1,j,k+1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	

				//test (i,j-1,k-1)
				if((j-1>=0)  &&  (k-1>=0))
					if( SEG(i,j-1,k-1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i,j-1,k-1);
							structFinded=SEG(i,j-1,k-1);
						}else if(structFinded!=SEG(i,j-1,k-1)){
							mini=min(structFinded,SEG(i,j-1,k-1),ptr[structFinded], ptr[SEG(i,j-1,k-1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i,j-1,k-1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}

				//test (i,j-1,k+1)
				if((j-1>=0)  &&  (k+1<sizeDataZ))
					if( SEG(i,j-1,k+1)!=0 ){
						if(structFinded==0){
							SEG(i,j,k)=SEG(i,j-1,k+1);
							structFinded=SEG(i,j-1,k+1);
						}else if(structFinded!=SEG(i,j-1,k+1)){
							mini=min(structFinded,SEG(i,j-1,k+1),ptr[structFinded], ptr[SEG(i,j-1,k+1)]);
							ptr[structFinded]=mini;
							ptr[SEG(i,j-1,k+1)]=mini;	
							SEG(i,j,k)=mini;
							structFinded=mini;
						}
					}	
			}

			if(structFinded==0){
				SEG(i,j,k)=(unsigned short)nSeg;
				inc(nSeg, &ptr, sizePtr);
			}
		}


	//assembling structures

	for(i=1;i<(long)nSeg;i++)
		if(ptr[i]!=(unsigned long)i)
			ptr[i]=ptr[ptr[i]];

	ptr2=(unsigned long *)malloc(sizeof(long)*nSeg);
	for(i=0;i<(long)nSeg;i++)
		ptr2[i]=0;
	structFinded=1;
	for(i=1;i<(long)nSeg;i++)
		if(ptr2[ptr[i]]==0){
			ptr2[ptr[i]]=structFinded;
			structFinded++;
		}


	for(i=0;i<sizeDataX*sizeDataY*sizeDataZ;i++)
		seg.data[i]=(unsigned short)ptr2[ptr[seg.data[i]]];

	free(ptr);
	free(ptr2);
	nbInfoStruct=structFinded-1;
	mixIndex(seg);
	return(structFinded-1);
}

void  Segment3D::mixIndex(USrtArray &seg){
	
	USrtArray *mixer=new USrtArray(nbInfoStruct+1);
	(*mixer)(0)=0;
	srand( 111 );

	long ind, val,i;
	for(i=1; i<=nbInfoStruct;i++){
		(*mixer)(i)=i;
	}
	for(i=1; i<=nbInfoStruct;i++){
		ind=rand()%nbInfoStruct+1;
		val=(*mixer)(ind);
		(*mixer)(ind)=(*mixer)(i);
		(*mixer)(i)=val;
	}

	
	for(i=0; i<seg.nbPoints;i++){
		seg[i]=(*mixer)(seg.data[i]);	
	}
	delete mixer;
}


#undef DATA
#undef SEG

