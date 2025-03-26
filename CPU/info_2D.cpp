

#include "stdafx.h"
#include "info_2D.h"


InfoStruct2D::InfoStruct2D(int index){
	data.index=index;
	data.zOrderComputed=0;
	next=NULL;
	__nbInfoStrucAlloc++;
}

InfoStruct2D::~InfoStruct2D(){
	__nbInfoStrucDealloc++;
}

