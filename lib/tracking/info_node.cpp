/*******************************************************************************************************************
									InfoNode Class
								  InfoNode contains the information for a given Node in 2D or in 3D
*******************************************************************************************************************/


#include <cmath>

#include "Info_Node.h"

//extern PublicParameters *publicParameters;

InfoNode::InfoNode(PublicParameters* _publicParameters) { publicParameters = _publicParameters; }

float InfoNode2D::distance(InfoNode* infoNode)
{
	const double dX = getX() - infoNode->getX();
	const double dY = getY() - infoNode->getY();
	const float res = (float)sqrt(dX * dX + dY * dY);
	if (res < 0.) { return 0.; }
	return res;
}

float InfoNode2D::distance(InfoNode* i1, InfoNode* i2)
{
	const double n1n2X = (i1->getX() * i1->getIntensity() + i2->getX() * i2->getIntensity()) / (i1->getIntensity() + i2->getIntensity());
	const double n1n2Y = (i1->getY() * i1->getIntensity() + i2->getY() * i2->getIntensity()) / (i1->getIntensity() + i2->getIntensity());
	const double dX = getX() - n1n2X;
	const double dY = getY() - n1n2Y;
	float res = (float)sqrt(dX * dX + dY * dY);
	if (res < 0.) { return 0.; }
	if (res < 1e-9) { res += 0; }
	return res;
}


double InfoNode2D::getMean() { return info.mean; }
double InfoNode2D::getIntensity() { return info.intensity; }
unsigned long InfoNode2D::getIndex() { return info.index; }
unsigned long InfoNode2D::getTime() { return info.time; }

long InfoNode2D::getSurface() { return info.surf; }

double InfoNode2D::getX() { return info.moment.y; }

double InfoNode2D::getY() { return info.moment.x; }

double InfoNode2D::getZ()
{
	// modified by Adel
	//return(T_ERR_OOP);
	return Z;
}

// Added by Adel
void InfoNode2D::setColor(const double color) { Color = color; }

void InfoNode2D::setMolId(const double molId) { MolId = molId; }

double InfoNode2D::getColor() { return Color; }
double InfoNode2D::getMolId() { return MolId; }
//

unsigned long InfoNode2D::sizeofInfoData() { return sizeof(InfoData2D); }

void* InfoNode2D::getPtrInfoData() { return &info; }

void InfoNode2D::setSurface(const long s) { info.surf = s; }

void InfoNode2D::setX(const double x) { info.moment.y = x; }

void InfoNode2D::setY(const double y) { info.moment.x = y; }

void InfoNode2D::setZ(const double z) { Z = z; }

void InfoNode2D::setIntensity(const double I) { info.intensity = I; }

void InfoNode2D::setIndex(const unsigned long i) { info.index = i; }

void InfoNode2D::setTime(const unsigned long t) { info.time = t; }

void InfoNode2D::setMean(const double m) { info.mean = m; }

//!!!!! invertion of coordinates !!! to restect Amira and UIC convenitons
float InfoNode3D::distance(InfoNode* infoNode)
{
	const double dX = getX() - infoNode->getX();
	const double dY = getY() - infoNode->getY();
	const double dZ = getZ() - infoNode->getZ();
	const float res = (float)sqrt(dX * dX + dY * dY + dZ * dZ);
	if (res < 0.) { return 0.; }
	return res;
}

float InfoNode3D::distance(InfoNode* i1, InfoNode* i2)
{
	const double n1n2X = (i1->getX() * i1->getIntensity() + i2->getX() * i2->getIntensity()) / (i1->getIntensity() + i2->getIntensity());
	const double n1n2Y = (i1->getY() * i1->getIntensity() + i2->getY() * i2->getIntensity()) / (i1->getIntensity() + i2->getIntensity());
	const double n1n2Z = (i1->getZ() * i1->getIntensity() + i2->getZ() * i2->getIntensity()) / (i1->getIntensity() + i2->getIntensity());
	const double dX = getX() - n1n2X;
	const double dY = getY() - n1n2Y;
	const double dZ = getZ() - n1n2Z;
	const float res = (float)sqrt(dX * dX + dY * dY + dZ * dZ);
	if (res < 0.) { return 0.; }
	return res;
}


double InfoNode3D::getMean() { return info.mean; }
double InfoNode3D::getIntensity() { return info.intensity; }
long InfoNode3D::getVolume() { return info.volume; }
unsigned long InfoNode3D::getIndex() { return info.index; }
unsigned long InfoNode3D::getTime() { return info.time; }

unsigned long InfoNode3D::sizeofInfoData() { return sizeof(InfoData3D); }

void* InfoNode3D::getPtrInfoData() { return &info; }


//!!!!! invertion of coordinates !!! to respect Amira and UIC convenitons
double InfoNode3D::getX() { return info.moment.z; }
double InfoNode3D::getY() { return info.moment.y; }
double InfoNode3D::getZ() { return info.moment.x * publicParameters->dZDX; }

void InfoNode3D::setSurface(long s) { }

void InfoNode3D::setX(double x) { }

void InfoNode3D::setY(double y) { }

void InfoNode3D::setZ(double z) { }

void InfoNode3D::setIntensity(double I) { }

void InfoNode3D::setIndex(unsigned long i) { }

void InfoNode3D::setTime(unsigned long t) { }

void InfoNode3D::setMean(double m) { }
