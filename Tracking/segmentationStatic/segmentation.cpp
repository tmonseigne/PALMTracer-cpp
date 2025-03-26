
/*******************************************************************************************************************
								Data segmentation for image (2D & 3D signal) (16bit only)
									Institut Curie
									UMR - 144
									by Victor Racine
									04/11/02 
*******************************************************************************************************************/

#include "segmentation.h"
extern FILE *debugStreamALS;
int __nbInfoStrucAlloc;
int __nbInfoStrucDealloc;

Segment::Segment(){
	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "## Segment() ##\n");
	closeDebugStreamALS();
	#endif
	__nbInfoStrucAlloc=0;
	__nbInfoStrucDealloc=0;
}

Segment::~Segment(){
	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "## ~Segment() ##\n");
	closeDebugStreamALS();
	#endif

	if(__nbInfoStrucAlloc!=__nbInfoStrucDealloc){
		openDebugStreamALS();
		fprintf(debugStreamALS, "!!! Error in allocation of InfoStruct2D/3D: nbInfoStrucAlloc=%d, nbInfoStrucDealloc=%d !!!\n");
		closeDebugStreamALS();
	}
}


void Segment::inc( unsigned long &nSeg, unsigned long **ptr, long &sizePtr ){
	unsigned long *newPtr;
	int i;
	if(nSeg==65535) return;
	nSeg++;
	if((long)nSeg<sizePtr) return;

	#ifdef	DEBUG_AT
	openDebugStreamALS();
	fprintf(debugStreamALS, "-- Segment::inc(nSeg: %d, sizePtr: %d) --\n",sizePtr);
	closeDebugStreamALS();
	#endif
	
	newPtr=(unsigned long *)malloc(sizeof(long)*sizePtr);
	//backup of data
	for(i=0;i<(long)sizePtr;i++)
		newPtr[i]=(*ptr)[i];

	//reallocation of ptr
	sizePtr+=MAX_TAB;
	(*ptr)=(unsigned long *)realloc((*ptr),sizeof(long)*sizePtr);

	//copy back the data
	for(i=0;i<((long)sizePtr-MAX_TAB);i++)
		(*ptr)[i]=newPtr[i];
	for(i=sizePtr-MAX_TAB;i<(long)sizePtr;i++)
		(*ptr)[i]=i;
	free(newPtr);
}

/*************************************************************************************************************************************
locateStringInFile try to find the string st in the stream "file". This file have to be opened.
*************************************************************************************************************************************/
int Segment::locateStringInFile(wchar_t* st, FILE* file){
	bool stingFounded=false;
	wchar_t line[1024];
	wchar_t *res;
	long offs;
	if(fseek(file,-30,SEEK_END)!=0)	//change file position
		return(S_READ_ERR);
	while(fgetws(line, 256, file)!=NULL){  //reading a new line
		if((res=wcsstr(line,st))!=NULL){	//if st exists in "line"
			offs=ftell(file)-wcslen(line)+res-line-1;
			return(offs);
		}
	}
	return(S_READ_ERR);
}




/*************************************************************************************************************************************
printdTMLCell print an array cell in HTML for tde statistic file.
*************************************************************************************************************************************/
void Segment::printHTMLCell(wchar_t* st, FILE* file, long countRaw, long countCol){
	fprintf(file,"\t\t<td");

	if(countRaw!=1)
		fwprintf(file,L" rowspan=\"%d\"", countRaw);
	if(countCol!=1)
		fwprintf(file,L" colspan=\"%d\"", countCol);
	fwprintf(file,L"><font color=\"#ffffff\"><b>%ls</b></font></td>\n", st);
}


/*************************************************************************************************************************************
printFloatdTMLCell print a float array cell in HTML for tde statistic file.
*************************************************************************************************************************************/
void Segment::printFloatHTMLCell(double f, FILE* file){
	fwprintf(file,L"\t\t<td>%f</td>\n", f);
}


/*************************************************************************************************************************************
printIntdTMLCell print a interger array cell in HTML for tde statistic file.
*************************************************************************************************************************************/
void Segment::printIntHTMLCell(long i, FILE* file, long count){
	if(count==1)
		fwprintf(file,L"\t\t<td>%d</td>\n", i);
	else
		fwprintf(file,L"\t\t<td   bgcolor=\"#FFEEEE\" rowspan=\"%d\">%d</td>\n", count, i);
}
