/*******************************************************************************************************************
								Link class to treat data from segmentation of images 2D (16bit only)
*******************************************************************************************************************/
#include "info_2D.h"


InfoStruct2D::InfoStruct2D(int index)
{
	data.index = index;
	data.zOrderComputed = 0;
	next = nullptr;
}

InfoStruct2D::~InfoStruct2D() { }
