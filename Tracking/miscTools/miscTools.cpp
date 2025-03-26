
#include "stdafx.h"


/*******************************************************************************************************************
								Miscellaneous tools
*******************************************************************************************************************/


#include "miscTools.h"
extern FILE *debugStreamALS;

void fitPow(double *x, double *y, int ndata, double &pow, double &k){
//Given a set of data points x[1..ndata], y[1..ndata] fit them to a function y = k.x^pow. Returned are pow and k
int i;
double t,sxoss,sx=0.0,sy=0.0,st2=0.0,ss;
if(x==NULL){
	for(i=0;i<ndata;i++){
		y[i]=log(y[i]);
	}
	pow=0.0;
	for (i=0;i<ndata;i++) { 
		sx += log((double)(i+1));
		sy += y[i];
	}
	ss=ndata;
	sxoss=sx/ss;
	for (i=0;i<ndata;i++) {
		t=log((double)(i+1))-sxoss;
		st2 += t*t;
		pow += t*y[i];
	}
	pow /= st2; 
	k=(sy-sx*(pow))/ss;
}
return;
for(i=0;i<ndata;i++){
	x[i]=log(x[i]);
	y[i]=log(y[i]);
}
pow=0.0;
for (i=0;i<ndata;i++) { 
	sx += x[i];
	sy += y[i];
}
ss=ndata;
sxoss=sx/ss;
for (i=0;i<ndata;i++) {
	t=x[i]-sxoss;
	st2 += t*t;
	pow += t*y[i];
}
pow /= st2; 
k=(sy-sx*(pow))/ss;

}

double getAutoThreshold(FltArray &im){

	double threshold;
	double oldThreshold;
	double m0, m1; //the average of the 2 parts of the histogram
	double count0, count1, count;
	double s0, s1;
	double sigma;


	
	
	long cpt=30;
	threshold=0.;
	oldThreshold=1e10;
	count=0;
	for(long i=0;i<(long)im.nbPoints;i++){
		if(im.data[i]>0.)		{
			count++;
			threshold+=im.data[i];
		}
	}
	threshold/=count;

	while((fabs(oldThreshold-threshold)>1e-12)&&(cpt>0)){
		count0=count1=0.;
		m0=m1=0.;
		for(long i=0;i<(long)im.nbPoints;i++){
			if(im.data[i]>0.)			{
				if(im.data[i]<threshold){
					m0+=im.data[i];
					count0++;
				}else{
					m1+=im.data[i];
					count1++;			
				}
			}
		}
		
		m0/=count0;
		m1/=count1;
		s0=s1=0.;

		long i;
		for(i=0;i<(long)im.nbPoints;i++){
			if(im.data[i]>0.)			{
				if(im.data[i]<threshold){
					s0+=(im.data[i]-m0)*(im.data[i]-m0);
				}else{
					s1+=(im.data[i]-m1)*(im.data[i]-m1);			
				}
			}
		}

		s0=sqrt(s0);
		s1=sqrt(s1);

		sigma=(s0+s1)/count;
		
		oldThreshold=threshold;
		threshold=(m0+m1)/2.+sigma*sigma*log(count0/count1)/(m1-m0);

		cpt--;
	}


	return(threshold);
}

void printErrno(FILE *file, int errno){
	switch(errno){
	case ECHILD:
		fprintf(file,"errno: No spawned processes.\n");
		break;
	case EAGAIN:
		fprintf(file,"errno: No more processes. An attempt to create a new process failed because there are no more process slots, or there is not enough memory, or the maximum nesting level has been reached.\n");
		break;
	case E2BIG:
		fprintf(file,"errno: Argument list too long.\n");
		break;
	case EACCES:
		fprintf(file,"errno: Permission denied. The file's permission setting does not allow the specified access. This error signifies that an attempt was made to access a file (or, in some cases, a directory) in a way that is incompatible with the file's attributes.\n");
		fprintf(file,"For example, the error can occur when an attempt is made to read from a file that is not open, to open an existing read-only file for writing, or to open a directory instead of a file. Under MS-DOS operating system versions 3.0 and later, EACCES may also indicate a locking or sharing violation.\n");
		fprintf(file,"The error can also occur in an attempt to rename a file or directory or to remove an existing directory.\n");
		break;
	case EBADF:
		fprintf(file,"errno: Bad file number. There are two possible causes: 1) The specified file handle is not a valid file-handle value or does not refer to an open file. 2) An attempt was made to write to a file or device opened for read-only access.\n");
		break;
	case EDEADLOCK:
		fprintf(file,"errno: Resource deadlock would occur. The argument to a math function is not in the domain of the function.\n");
		break;
	case EDOM:
		fprintf(file,"errno: Math argument.\n");
		break;
	case EEXIST:
		fprintf(file,"errno: Files exist. An attempt has been made to create a file that already exists. For example, the _O_CREAT and _O_EXCL flags are specified in an _open call, but the named file already exists.\n");
		break;
	case EINVAL:
		fprintf(file,"errno: Invalid argument. An invalid value was given for one of the arguments to a function. For example, the value given for the origin when positioning a file pointer (by means of a call to fseek) is before the beginning of the file.\n");
		break;
	case EMFILE:
		fprintf(file,"errno: Too many open files. No more file handles are available, so no more files can be opened.\n");
		break;
	case ENOENT:
		fprintf(file,"errno: No such file or directory. The specified file or directory does not exist or cannot be found. This message can occur whenever a specified file does not exist or a component of a path does not specify an existing directory.\n");
		break;
	case ENOEXEC:
		fprintf(file,"errno: Exec format error. An attempt was made to execute a file that is not executable or that has an invalid executable-file format.\n");
		break;
	case ENOMEM:
		fprintf(file,"errno: Not enough core. Not enough memory is available for the attempted operator. For example, this message can occur when insufficient memory is available to execute a child process, or when the allocation request in a _getcwd call cannot be satisfied.\n");
		break;
	case ENOSPC:
		fprintf(file,"errno: No space left on device. No more space for writing is available on the device (for example, when the disk is full).\n");
		break;
	case ERANGE:
		fprintf(file,"errno: Result too large. An argument to a math function is too large, resulting in partial or total loss of significance in the result. This error can also occur in other functions when an argument is larger than expected (for example, when the buffer argument to _getcwd is longer than expected).\n");
		break;
	case EXDEV:
		fprintf(file,"errno: Cross-device link. An attempt was made to move a file to a different device (using the rename function).\n");
		break;
	default:
		fprintf(file,"errno: Unknown error.\n");

	}
}

double tenTo(double a){
	return(exp(a*log(10.)));
}

wchar_t* stringReplace(wchar_t s, wchar_t d, wchar_t* st){
	wchar_t*ptrSt=st;
	while(*ptrSt!='\0'){
		if(*ptrSt==s)
			*ptrSt=d;
		ptrSt++;
	}
	return(st);
}


unsigned long twoToThe(unsigned long val){
	if(val==0) return(1);
	return(((unsigned long)1) << val);
}


double distance(double X1, double Y1, double X2, double Y2){
	return(sqrt(distance2(X1,Y1,X2,Y2)));
}
double distance(double X1, double Y1, double Z1, double X2, double Y2, double Z2){
	return(sqrt(distance2(X1,Y1,Z1,X2,Y2,Z2)));
}

double distance2(double X1, double Y1, double X2, double Y2){
	double res;
	res=(X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2);
	if(res<0)
		return(0.);
	return(res);
}
double distance2(double X1, double Y1, double Z1, double X2, double Y2, double Z2){
	double res;
	res=(X1-X2)*(X1-X2)+(Y1-Y2)*(Y1-Y2)+(Z1-Z2)*(Z1-Z2);
	if(res<0)
		return(0.);
	return(res);
} 

/********************************************************************************************************************
computeMSD(DblArray &mSD,DblArray &positions) print the mean square displacement: 
MSD(n)=<r²n>=1/(j-n) * sum(d(i,i+n)², i=0:i-n-1)
It is the definition of Anderson et al. Journal of Cell Science 101, 415-425 (1992).
mSD is the result and positions are the point positions over time.
********************************************************************************************************************/
int computeMSD(DblArray &mSD,DblArray &positions){
	double val;
	if(mSD.nbDim!=1)
		return(ERR_DIM);
	if(mSD.size(0)!=positions.size(0)-2)
		return(ERR_DIM);
	if((positions.size(1)!=2)&&(positions.size(1)!=3))
		return(ERR_DIM);
	long j=positions.size(0);
	if(positions.size(1)==2)
		for(long n=1;n<=mSD.size(0);n++){
			val=0.;
			for(long i=0;i<=j-n-1;i++)
				val+=distance2(positions(i,0),positions(i,1),positions(i+n,0),positions(i+n,1));
			val/=(j-n);
			mSD(n-1)=val;
		}
	else
		for(long n=1;n<=mSD.size(0);n++){
			val=0.;
			for(long i=0;i<=j-n-1;i++)
				val+=distance2(positions(i,0),positions(i,1),positions(i,2),positions(i+n,0),positions(i+n,1),positions(i+n,2));
			val/=(j-n);
			mSD(n-1)=val;
		} 

	return(RET_OK);
}


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


int printRegionsW(wchar_t* filename,  USrtArray &mainSegmentation, long nbStructures){
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
	file=_wfopen( filename, L"w");
	if(file==NULL) return(ret);
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
					ret=fwprintf(file, L"0 3, 1 65535, 2 %d %d, 3 0 0, 4 0, 5 1, 6 %d", isAlreadyWritten(val,2), isAlreadyWritten(val,1), nbPoints);
					do{ // second round of the structure to write the points in the textfile
						if(ISOUT(posX+turnClockwiseX(dirX,dirY), posY+turnClockwiseY(dirX,dirY)))
							turnClockwise(dirX,dirY);
						else if(ISOUT(posX+turnClockwiseX(dirX,dirY)+dirX,posY+turnClockwiseY(dirX,dirY)+dirY)){
							posX+=turnClockwiseX(dirX,dirY);
							posY+=turnClockwiseY(dirX,dirY);
							usedBorder.data[posX*ySize+posY]=1;
							ret=fwprintf(file, L" %d %d", posY, posX);
						}else{
							posX+=turnClockwiseX(dirX,dirY)+dirX;
							posY+=turnClockwiseY(dirX,dirY)+dirY;
							turnAntiClockwise(dirX,dirY);
							usedBorder.data[posX*ySize+posY]=1;
							ret=fwprintf(file, L" %d %d", posY, posX);
						}
					}while((posX!=x)||(posY!=y)||(dirX!=dirX0)||(dirY!=dirY0));
					ret=fwprintf(file, L", 7 %d\n", val);
				}
			}
	}
	fclose(file);
	return(RET_OK);
}

#undef ISOUT