#include "aTrous.h"
#include "segmentation_2D.h"
#include "tools/allDataLib.h"

// Cette classe est exportée de PALM.dll
//class PALM_API CPALM {

class CPALM
{
public:
	// TODO : ajoutez ici vos méthodes.

	long CPALMProcessing();
	CPALM(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size);
	~CPALM();
	unsigned long _sizeX;
	unsigned long _sizeY;
	unsigned short* _image;
	double* _pointsList;

private:
	unsigned long nbCalls;

	long _pointsListSize;
	int _nbPointsMaxPerPlane;
	int _nbSegParams;
	unsigned long _waveletNo;
	double _thresholdVal;
	double _volMin;
	double _intMin;
	double _watershedRatio;
	unsigned short _gaussFit;
	float _sigmaX;
	float _sigmaY;
	float _angle;
	unsigned short _size;
	long currentPoint;

	double _structFiltParam[13];
	unsigned short _structFiltParamUse[13];
	unsigned short _planesToUse[10];
	unsigned short _operators[9];
	double _threshold[9];

	ATrous* _aTrousPub;
	USrtArray* _segPub;
	FltArray* _outputPub;
	Segment2D* _segmentPub;
	USrtArray* _input;
	BlArray* _planesToUseB;
	BlArray* _planesToUseB2;
	int _nbPlane;
	int waveletSegmentation2D(unsigned short* outputData);
	int waveletSegmentation2D2(unsigned short* outputData);
};
