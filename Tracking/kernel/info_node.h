/*******************************************************************************************************************
									InfoNode InfoNode2D InfoNode3D Classes
								  InfoDataXX contains the information for a given Node in 2D or in 3D
*******************************************************************************************************************/


#ifndef __info_node__
#define __info_node__

#include "general.h"
#include "../segmentationStatic/info_2D.h"
#include "../segmentationStatic/info_3D.h"

class InfoNode
{
public:
	InfoNode(PublicParameters* _publicParameters);
	~InfoNode() { };

	PublicParameters* publicParameters;
	virtual float distance(InfoNode* infoData) { return (T_ERR_OOP); };
	virtual float distance(InfoNode* i1, InfoNode* i2) { return (T_ERR_OOP); };
	virtual double getMean() { return (T_ERR_OOP); };
	virtual long getVolume() { return (T_ERR_OOP); };
	virtual long getSurface() { return (T_ERR_OOP); };
	virtual double getX() { return (T_ERR_OOP); };
	virtual double getY() { return (T_ERR_OOP); };
	virtual double getZ() { return (T_ERR_OOP); }
	virtual double getIntensity() { return (T_ERR_OOP); };
	virtual unsigned long getIndex() { return (T_ERR_OOP); };
	virtual unsigned long getTime() { return (T_ERR_OOP); };
	virtual unsigned long sizeofInfoData() { return (T_ERR_OOP); };
	virtual void* getPtrInfoData() { return (NULL); };

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
	float distance(InfoNode* infoData);
	float distance(InfoNode* i1, InfoNode* i2);
	InfoData2D info;
	double getMean();
	double getX();
	double getY();
	double getZ();
	double getIntensity();
	long getSurface();
	unsigned long getIndex();
	unsigned long getTime();
	unsigned long sizeofInfoData();
	void* getPtrInfoData();

	void setSurface(long s);
	void setX(double x);
	void setY(double y);
	void setZ(double z);
	void setIntensity(double I);
	void setIndex(unsigned long i);
	void setTime(unsigned long t);
	void setMean(double m);

	// added by Adel
	double Z;
	double Color;
	double MolId;
	void setColor(double z);
	void setMolId(double molId);
	double getColor();
	double getMolId();
};


class InfoNode3D : public InfoNode
{
public:
	InfoNode3D(PublicParameters* _publicParameters): InfoNode(_publicParameters) { };
	float distance(InfoNode* infoData);
	float distance(InfoNode* i1, InfoNode* i2);
	InfoData3D info;
	double getMean();
	double getX();
	double getY();
	double getZ();
	double getIntensity();
	long getVolume();
	unsigned long getIndex();
	unsigned long getTime();
	unsigned long sizeofInfoData();
	void* getPtrInfoData();

	void setSurface(long s);
	void setX(double x);
	void setY(double y);
	void setZ(double z);
	void setIntensity(double I);
	void setIndex(unsigned long i);
	void setTime(unsigned long t);
	void setMean(double m);
};

#endif //__info_node__
