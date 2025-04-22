/*******************************************************************************************************************
									InfoNode Class
								  InfoNode contains the information for a given Node in 2D or in 3D
*******************************************************************************************************************/
#include "InfoNode.h"
#include <cmath>

double InfoNode::Distance(const InfoNode& node) const
{
	const double Dx = x_ - node.x_, Dy = y_ - node.y_;
	return sqrt(Dx * Dx + Dy * Dy);
}

double InfoNode::Distance(const InfoNode& n1, const InfoNode& n2) const
{
	const double N1N2X = (n1.x_ * n1.intensity_ + n2.x_ * n2.intensity_) / (n1.intensity_ + n2.intensity_);
	const double N1N2Y = (n1.y_ * n1.intensity_ + n2.y_ * n2.intensity_) / (n1.intensity_ + n2.intensity_);
	const double Dx = x_ - N1N2X, Dy = y_ - N1N2Y;
	return sqrt(Dx * Dx + Dy * Dy);
}
