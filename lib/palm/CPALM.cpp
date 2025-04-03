// PALM.cpp : définit le point d'entrée pour l'application DLL.
//

#include "CPALM.h"
#include "fit/GaussianFit2D.h"
#include "fit/GaussianFit2D_Sigma.h"
#include "fit/GaussianFit3D.h"
#include "fit/GaussianFit3D_Theta.h"

CPALM::CPALM(unsigned short* image, double* pointsList, const unsigned long pointsListSize, const unsigned long sizeX, const unsigned long sizeY,
             const unsigned long waveletNo, const double thresholdVal, const double watershedRatio, const double volMin, const double intMin,
             const unsigned short gaussFit, const double sigmaX, const double sigmaY, const double angle, const unsigned short roi)
{
	currentPoint = 0;
	_image = image;
	_pointsList = pointsList;
	_pointsListSize = pointsListSize;
	_nbSegParams = 13; /*_nbSegParams = 10;*/
	_nbPointsMaxPerPlane = _pointsListSize / _nbSegParams;
	_sizeX = sizeX;
	_sizeY = sizeY;
	_waveletNo = waveletNo;
	_thresholdVal = thresholdVal;
	_watershedRatio = watershedRatio;
	_aTrousPub = new ATrous();
	_segPub = new USrtArray();
	_outputPub = new FltArray();
	_segmentPub = new Segment2D();
	_volMin = volMin;
	_intMin = intMin;
	_gaussFit = gaussFit;
	_sigmaX = (float)sigmaX;
	_sigmaY = (float)sigmaY;
	_angle = (float)angle;
	_size = roi;

	for (double& i : _structFiltParam) { i = 0; }
	for (unsigned short& i : _structFiltParamUse) { i = 0; }
	for (unsigned short& i : _planesToUse) { i = 0; }
	for (unsigned short& _operator : _operators) { _operator = 0; }
	for (double& i : _threshold) { i = 0; }


	_planesToUse[_waveletNo] = 2;
	_threshold[_waveletNo] = _thresholdVal;

	if (_volMin > 0.) {
		_structFiltParamUse[0] = 1; // minimum volume
		_structFiltParam[0] = _volMin;
	}
	if (_intMin > 0.) {
		_structFiltParamUse[2] = 1; // minimum intensity
		_structFiltParam[2] = _intMin;
	}

	_input = new USrtArray();
	if (_input->alloc(_image, _sizeX, _sizeY) < 0) { return; }
	if (_aTrousPub->wDecomp->alloc(_sizeX, _sizeY, waveletNo + 2) < 0) { return; }
	if (_segPub->realloc(_sizeX, _sizeY) < 0) { return; }
	if (_outputPub->alloc(_sizeX, _sizeY) < 0) { return; }


	_planesToUseB = new BlArray(9);

	_nbPlane = 0;
	for (int i = 0; i < 9; i++) {
		(*_planesToUseB)(i) = _planesToUse[i] != 0;
		if ((*_planesToUseB)(i)) { _nbPlane = i + 2; }
	}
	(*_planesToUseB)(8) = false;
}

long CPALM::CPALMProcessing()
{
	int ret;
	long nbStruct;
	int nbPointsSign;

	if ((ret = _aTrousPub->aTrousDecomposition2D(*_input, *_aTrousPub->wDecomp, *_planesToUseB, _nbPlane, false)) < 0) { return ret; }

	for (long j = 0; j < long(_sizeX * _sizeY); j++) {
		const long offs = _waveletNo * _sizeX * _sizeY;
		if (_aTrousPub->wDecomp->data[offs + j] > _thresholdVal) {
			_outputPub->data[j] = _aTrousPub->wDecomp->data[offs + j];
			_segPub->data[j] = 1;
		}
		else {
			_segPub->data[j] = 0;
			_outputPub->data[j] = 0;
		}
	}

	if (_watershedRatio == 1.) { if ((nbStruct = _segmentPub->segmentation(*_outputPub, *_segPub, 1e-8f, true, false)) < 0) { return nbStruct; } }
	else { if ((nbStruct = _segmentPub->watershedIntelligent(*_outputPub, *_segPub, 1e-8f, (float)_watershedRatio)) < 0) { return nbStruct; } }


	if ((nbPointsSign = _segmentPub->computeStatistics(*_outputPub, *_segPub, _aTrousPub->wDecomp->data + _waveletNo * _sizeX * _sizeY, _gaussFit)) < 0) {
		return nbPointsSign;
	}


	if (_volMin > 0. || _intMin > 0.) { if ((nbStruct = _segmentPub->filterData(*_segPub, _structFiltParam, _structFiltParamUse)) < 0) { return ret; } }

	InfoStruct2D** p = (InfoStruct2D**)_segmentPub->getTabInfo();

	for (int i = 1; i <= min(nbStruct, _nbPointsMaxPerPlane - 2); ++i) {
		switch (_gaussFit) {
			case 0:
			case 5: {
				initGaussianFitParam(p[i]->data.moment.y, p[i]->data.moment.x, 0, 0, 0);
				rss_ = 0;
				fillPointList(p[i]->data);
				break;
			}

			case 1: // Sigma is constant
			{
				// Gaussian fit of 2D
				initGaussianFitParam(p[i]->data.moment.y, p[i]->data.moment.x, _sigmaX, 0, 0);
				GaussianFit2D* fit = new GaussianFit2D();
				fit->Run(_image, _sizeY, _sizeX, x_, y_, intensity_, sigmaX_, offset_, rss_, 10, _size);
				sigmaY_ = sigmaX_;
				fillPointList(p[i]->data);
				break;
			}

			case 2: {
				// Gaussian fit of 2D
				initGaussianFitParam(p[i]->data.moment.y, p[i]->data.moment.x, _sigmaX, 0, 0);
				GaussianFit2D_Sigma* fit = new GaussianFit2D_Sigma();
				fit->Run(_image, _sizeY, _sizeX, x_, y_, intensity_, sigmaX_, offset_, rss_, 10, _size);
				sigmaY_ = sigmaX_;
				fillPointList(p[i]->data);
				break;
			}

			case 3: // Theta is constant
			{
				// Gaussian fit of 3D
				initGaussianFitParam(p[i]->data.moment.y, p[i]->data.moment.x, _sigmaX, _sigmaY, _angle);
				GaussianFit3D* fit = new GaussianFit3D();
				fit->Run(_image, _sizeY, _sizeX, x_, y_, intensity_, sigmaX_, sigmaY_, theta_, offset_, rss_, 10, _size);
				fillPointList(p[i]->data);
				break;
			}

			case 4: {
				// Gaussian fit of 3D
				initGaussianFitParam(p[i]->data.moment.y, p[i]->data.moment.x, _sigmaX, _sigmaY, _angle);
				GaussianFit3D_Theta* fit = new GaussianFit3D_Theta();
				fit->Run(_image, _sizeY, _sizeX, x_, y_, intensity_, sigmaX_, sigmaY_, theta_, offset_, rss_, 20 * 7, _size);
				fillPointList(p[i]->data);
				break;
			}

			default: // Error on Gauss Fit Param
				std::fill_n(_pointsList + currentPoint, 13, -1);
				currentPoint += 13;
				break;
		}
	}

	// Blank Line
	std::fill_n(_pointsList + currentPoint, 13, -1);
	currentPoint += 13;

	return currentPoint;
}


CPALM::~CPALM()
{
	delete _aTrousPub;
	delete _segPub;
	delete _outputPub;
	delete _segmentPub;
	delete _input;
	delete _planesToUseB;
}

void CPALM::initGaussianFitParam(const double x, const double y, const double sigmaX, const double sigmaY, const double theta)
{
	x_ = x;
	y_ = y;
	intensity_ = 0;
	sigmaX_ = sigmaX;
	sigmaY_ = sigmaY;
	theta_ = theta;
	offset_ = 0;
	rss_ = -1;
}

void CPALM::fillPointList(const InfoDataRich2D& data)
{
	_pointsList[currentPoint++] = sigmaX_;
	_pointsList[currentPoint++] = sigmaY_;
	_pointsList[currentPoint++] = theta_;
	_pointsList[currentPoint++] = y_;
	_pointsList[currentPoint++] = x_;
	_pointsList[currentPoint++] = intensity_;
	_pointsList[currentPoint++] = offset_;
	_pointsList[currentPoint++] = rss_;
	_pointsList[currentPoint++] = data.intensity;
	_pointsList[currentPoint++] = data.surf;

	_pointsList[currentPoint++] = 0.; // Z coordinate
	_pointsList[currentPoint++] = 0.; // Spectral or other
	_pointsList[currentPoint++] = 0.; // Other
}
