

/*******************************************************************************************************************
								Link class to treat data from segmentation of images 3D (16bit only)
									Institut Curie
									UMR - 144
									by Victor Racine
									april 25th 2003
*******************************************************************************************************************/
#include "info_3D.h"

InfoStruct3D::InfoStruct3D(int index){
	data.index=index;
	next=NULL;
	__nbInfoStrucAlloc++;
}

InfoStruct3D::~InfoStruct3D(){
	__nbInfoStrucDealloc++;
}