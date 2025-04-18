/*******************************************************************************************************************
									InfoNode Class
								  InfoNode contains the information for a given Node in 2D or in 3D
*******************************************************************************************************************/


#include <cmath>

#include "InfoNode.h"

//extern Parameters *parameters_;

InfoNode::InfoNode(Parameters* params) { parameters_ = params; }

double InfoNode2D::Distance(InfoNode* node)
{
	const double Dx = GetX() - node->GetX();
	const double Dy = GetY() - node->GetY();
	const double Res = sqrt(Dx * Dx + Dy * Dy);
	if (Res < 0.) { return 0.; }
	return Res;
}

double InfoNode2D::Distance(InfoNode* i1, InfoNode* i2)
{
	const double N1N2X = (i1->GetX() * i1->GetIntensity() + i2->GetX() * i2->GetIntensity()) / (i1->GetIntensity() + i2->GetIntensity());
	const double N1N2Y = (i1->GetY() * i1->GetIntensity() + i2->GetY() * i2->GetIntensity()) / (i1->GetIntensity() + i2->GetIntensity());
	const double Dx = GetX() - N1N2X;
	const double Dy = GetY() - N1N2Y;
	double res = sqrt(Dx * Dx + Dy * Dy);
	if (res < 0.) { return 0.; }
	if (res < 1e-9) { res += 0; }
	return res;
}


double InfoNode2D::GetMean() { return info_.mean; }
double InfoNode2D::GetIntensity() { return info_.intensity; }
unsigned long InfoNode2D::GetIndex() { return info_.index; }
long InfoNode2D::GetSurface() { return info_.surf; }
double InfoNode2D::GetX() { return info_.moment.y; }
double InfoNode2D::GetY() { return info_.moment.x; }
double InfoNode2D::GetZ() { return z_; }
void InfoNode2D::SetColor(const double color) { color_ = color; }
void InfoNode2D::SetMolId(const double molId) { molId_ = molId; }
double InfoNode2D::GetColor() { return color_; }
double InfoNode2D::GetMolId() { return molId_; }
unsigned long InfoNode2D::SizeofInfoData() { return sizeof(info_2d); }
void* InfoNode2D::GetPtrInfoData() { return &info_; }
void InfoNode2D::SetSurface(const long s) { info_.surf = s; }
void InfoNode2D::SetX(const double x) { info_.moment.y = x; }
void InfoNode2D::SetY(const double y) { info_.moment.x = y; }
void InfoNode2D::SetZ(const double z) { z_ = z; }
void InfoNode2D::SetIntensity(const double i) { info_.intensity = i; }
void InfoNode2D::SetIndex(const unsigned long i) { info_.index = i; }
void InfoNode2D::SetMean(const double m) { info_.mean = m; }

//!!!!! invertion of coordinates !!! to restect Amira and UIC convenitons
double InfoNode3D::Distance(InfoNode* node)
{
	const double Dx = GetX() - node->GetX();
	const double Dy = GetY() - node->GetY();
	const double Dz = GetZ() - node->GetZ();
	const double Res = sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
	if (Res < 0.) { return 0.; }
	return Res;
}

double InfoNode3D::Distance(InfoNode* i1, InfoNode* i2)
{
	const double N1N2X = (i1->GetX() * i1->GetIntensity() + i2->GetX() * i2->GetIntensity()) / (i1->GetIntensity() + i2->GetIntensity());
	const double N1N2Y = (i1->GetY() * i1->GetIntensity() + i2->GetY() * i2->GetIntensity()) / (i1->GetIntensity() + i2->GetIntensity());
	const double N1N2Z = (i1->GetZ() * i1->GetIntensity() + i2->GetZ() * i2->GetIntensity()) / (i1->GetIntensity() + i2->GetIntensity());
	const double Dx = GetX() - N1N2X;
	const double Dy = GetY() - N1N2Y;
	const double Dz = GetZ() - N1N2Z;
	const double Res = sqrt(Dx * Dx + Dy * Dy + Dz * Dz);
	if (Res < 0.) { return 0.; }
	return Res;
}


double InfoNode3D::GetMean() { return info_.mean; }
double InfoNode3D::GetIntensity() { return info_.intensity; }
long InfoNode3D::GetVolume() { return info_.volume; }
unsigned long InfoNode3D::GetIndex() { return info_.index; }
unsigned long InfoNode3D::SizeofInfoData() { return sizeof(info_3d); }
void* InfoNode3D::GetPtrInfoData() { return &info_; }
double InfoNode3D::GetX() { return info_.moment.z; }
double InfoNode3D::GetY() { return info_.moment.y; }
double InfoNode3D::GetZ() { return info_.moment.x * parameters_->dZdx_; }
void InfoNode3D::SetSurface(long s) { }
void InfoNode3D::SetX(double x) { }
void InfoNode3D::SetY(double y) { }
void InfoNode3D::SetZ(double z) { }
void InfoNode3D::SetIntensity(double i) { }
void InfoNode3D::SetIndex(unsigned long i) { }
void InfoNode3D::SetMean(double m) { }
