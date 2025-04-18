/*******************************************************************************************************************
									InfoNode InfoNode2D InfoNode3D Classes
								  InfoDataXX contains the information for a given Node in 2D or in 3D
*******************************************************************************************************************/
#pragma once

#include "Constants.hpp"
#include "Parameters.h"
#include "palm/Info2D.h"
#include "palm/Info3D.h"

class InfoNode
{
public:
	InfoNode(Parameters* params);
	~InfoNode() = default;

	Parameters* parameters_;
	virtual double Distance(InfoNode* node) { return T_ERR_OOP; }
	virtual double Distance(InfoNode* i1, InfoNode* i2) { return T_ERR_OOP; }
	virtual double GetMean() { return T_ERR_OOP; }
	virtual long GetVolume() { return T_ERR_OOP; }
	virtual long GetSurface() { return T_ERR_OOP; }
	virtual double GetX() { return T_ERR_OOP; }
	virtual double GetY() { return T_ERR_OOP; }
	virtual double GetZ() { return T_ERR_OOP; }
	virtual double GetIntensity() { return T_ERR_OOP; }
	virtual unsigned long GetIndex() { return T_ERR_OOP; }
	virtual unsigned long SizeofInfoData() { return T_ERR_OOP; }
	virtual void* GetPtrInfoData() { return nullptr; }

	virtual void SetSurface(long s) { }
	virtual void SetX(double x) { }
	virtual void SetY(double y) { }
	virtual void SetZ(double z) { }
	virtual void SetIntensity(double i) { }
	virtual void SetIndex(unsigned long i) { }
	virtual void SetMean(double m) { }

	// Added by Adel
	virtual void SetColor(double z) { }
	virtual void SetMolId(double molId) { }
	virtual double GetColor() { return 0; }
	virtual double GetMolId() { return 0; }
};

class InfoNode2D : public InfoNode
{
public:
	InfoNode2D(Parameters* params): InfoNode(params) { }
	double Distance(InfoNode* node) override;
	double Distance(InfoNode* i1, InfoNode* i2) override;
	info_2d info_;
	double GetMean() override;
	double GetX() override;
	double GetY() override;
	double GetZ() override;
	double GetIntensity() override;
	long GetSurface() override;
	unsigned long GetIndex() override;
	unsigned long SizeofInfoData() override;
	void* GetPtrInfoData() override;

	void SetSurface(long s) override;
	void SetX(double x) override;
	void SetY(double y) override;
	void SetZ(double z) override;
	void SetIntensity(double i) override;
	void SetIndex(unsigned long i) override;
	void SetMean(double m) override;

	// added by Adel
	double z_;
	double color_;
	double molId_;
	void SetColor(double color) override;
	void SetMolId(double molId) override;
	double GetColor() override;
	double GetMolId() override;
};


class InfoNode3D : public InfoNode
{
public:
	InfoNode3D(Parameters* params): InfoNode(params) { }
	double Distance(InfoNode* node) override;
	double Distance(InfoNode* i1, InfoNode* i2) override;
	info_3d info_;
	double GetMean() override;
	double GetX() override;
	double GetY() override;
	double GetZ() override;
	double GetIntensity() override;
	long GetVolume() override;
	unsigned long GetIndex() override;
	unsigned long SizeofInfoData() override;
	void* GetPtrInfoData() override;

	void SetSurface(long s) override;
	void SetX(double x) override;
	void SetY(double y) override;
	void SetZ(double z) override;
	void SetIntensity(double i) override;
	void SetIndex(unsigned long i) override;
	void SetMean(double m) override;
};
