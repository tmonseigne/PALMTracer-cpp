/*******************************************************************************************************************
									InfoNode Class
*******************************************************************************************************************/
#pragma once

class InfoNode
{
public:
	double x_, y_, z_, intensity_, mean_, surf_, color_, molId_;

	InfoNode(): x_(0), y_(0), z_(0), intensity_(0), mean_(0), surf_(0), color_(0), molId_(0) { }
	~InfoNode() = default;

	double Distance(const InfoNode& node) const;
	double Distance(const InfoNode& n1, const InfoNode& n2) const;
};
