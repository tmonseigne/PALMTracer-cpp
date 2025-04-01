/*******************************************************************************************************************
								Link class to treat data from segmentation of images 3D (16bit only)
*******************************************************************************************************************/
#include "info_3D.h"

InfoStruct3D::InfoStruct3D(int index)
{
	data.index = index;
	next = nullptr;
}

InfoStruct3D::~InfoStruct3D() { }
