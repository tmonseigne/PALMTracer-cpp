/*******************************************************************************************************************
									InfoNode Class
								  InfoNode contains the information for a given Node in 2D or in 3D
*******************************************************************************************************************/


#include <cmath>

#include "InfoNode.h"

//extern Parameters *parameters_;

double InfoNode2D::Distance(InfoNode* node)
{
	const double Dx = GetX() - node->GetX(), Dy = GetY() - node->GetY();
	return sqrt(Dx * Dx + Dy * Dy);
}

double InfoNode2D::Distance(InfoNode* n1, InfoNode* n2)
{
	const double N1N2X = (n1->GetX() * n1->GetIntensity() + n2->GetX() * n2->GetIntensity()) / (n1->GetIntensity() + n2->GetIntensity());
	const double N1N2Y = (n1->GetY() * n1->GetIntensity() + n2->GetY() * n2->GetIntensity()) / (n1->GetIntensity() + n2->GetIntensity());
	const double Dx = GetX() - N1N2X, Dy = GetY() - N1N2Y;
	return sqrt(Dx * Dx + Dy * Dy);
}


double InfoNode2D::GetMean() { return info_.mean; }
double InfoNode2D::GetIntensity() { return info_.intensity; }
long InfoNode2D::GetSurface() { return info_.surf; }
double InfoNode2D::GetX() { return info_.moment.y; }
double InfoNode2D::GetY() { return info_.moment.x; }
double InfoNode2D::GetZ() { return z_; }
void InfoNode2D::SetColor(const double color) { color_ = color; }
void InfoNode2D::SetMolId(const double molId) { molId_ = molId; }
double InfoNode2D::GetColor() { return color_; }
double InfoNode2D::GetMolId() { return molId_; }
void InfoNode2D::SetSurface(const long s) { info_.surf = s; }
void InfoNode2D::SetX(const double x) { info_.moment.y = x; }
void InfoNode2D::SetY(const double y) { info_.moment.x = y; }
void InfoNode2D::SetZ(const double z) { z_ = z; }
void InfoNode2D::SetIntensity(const double i) { info_.intensity = i; }
void InfoNode2D::SetMean(const double m) { info_.mean = m; }
