/*******************************************************************************************************************
									InfoNode InfoNode2D InfoNode3D Classes
								  InfoDataXX contains the information for a given Node in 2D or in 3D
*******************************************************************************************************************/
#pragma once

#include "general.h"
#include "palm/info_2D.h"
#include "palm/info_3D.h"

class InfoNode
{
public:
	InfoNode(PublicParameters* _publicParameters);
	~InfoNode() { };

	PublicParameters* publicParameters;
	virtual float distance(InfoNode* infoData) { return T_ERR_OOP; };
	virtual float distance(InfoNode* i1, InfoNode* i2) { return T_ERR_OOP; };
	virtual double getMean() { return T_ERR_OOP; };
	virtual long getVolume() { return T_ERR_OOP; };
	virtual long getSurface() { return T_ERR_OOP; };
	virtual double getX() { return T_ERR_OOP; };
	virtual double getY() { return T_ERR_OOP; };
	virtual double getZ() { return T_ERR_OOP; }
	virtual double getIntensity() { return T_ERR_OOP; };
	virtual unsigned long getIndex() { return T_ERR_OOP; };
	virtual unsigned long getTime() { return T_ERR_OOP; };
	virtual unsigned long sizeofInfoData() { return T_ERR_OOP; };
	virtual void* getPtrInfoData() { return NULL; };

	virtual void setSurface(long s) { };
	virtual void setX(double x) { };
	virtual void setY(double y) { };
	virtual void setZ(double z) { };
	virtual void setIntensity(double I) { };
	virtual void setIndex(unsigned long i) { };
	virtual void setTime(unsigned long t) { };
	virtual void setMean(double m) { };

	// Added by Adel
	virtual void setColor(double z) { };
	virtual void setMolId(double molId) { };
	virtual double getColor() { return 0; };
	virtual double getMolId() { return 0; };
};

class InfoNode2D : public InfoNode
{
public:
	InfoNode2D(PublicParameters* _publicParameters): InfoNode(_publicParameters) { };
	float distance(InfoNode* infoData) override;
	float distance(InfoNode* i1, InfoNode* i2) override;
	InfoData2D info;
	double getMean() override;
	double getX() override;
	double getY() override;
	double getZ() override;
	double getIntensity() override;
	long getSurface() override;
	unsigned long getIndex() override;
	unsigned long getTime() override;
	unsigned long sizeofInfoData() override;
	void* getPtrInfoData() override;

	void setSurface(long s) override;
	void setX(double x) override;
	void setY(double y) override;
	void setZ(double z) override;
	void setIntensity(double I) override;
	void setIndex(unsigned long i) override;
	void setTime(unsigned long t) override;
	void setMean(double m) override;

	// added by Adel
	double Z;
	double Color;
	double MolId;
	void setColor(double z) override;
	void setMolId(double molId) override;
	double getColor() override;
	double getMolId() override;
};


class InfoNode3D : public InfoNode
{
public:
	InfoNode3D(PublicParameters* _publicParameters): InfoNode(_publicParameters) { };
	float distance(InfoNode* infoData) override;
	float distance(InfoNode* i1, InfoNode* i2) override;
	InfoData3D info;
	double getMean() override;
	double getX() override;
	double getY() override;
	double getZ() override;
	double getIntensity() override;
	long getVolume() override;
	unsigned long getIndex() override;
	unsigned long getTime() override;
	unsigned long sizeofInfoData() override;
	void* getPtrInfoData() override;

	void setSurface(long s) override;
	void setX(double x) override;
	void setY(double y) override;
	void setZ(double z) override;
	void setIntensity(double I) override;
	void setIndex(unsigned long i) override;
	void setTime(unsigned long t) override;
	void setMean(double m) override;
};
