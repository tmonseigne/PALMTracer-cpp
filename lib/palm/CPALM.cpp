// PALM.cpp : définit le point d'entrée pour l'application DLL.
//

#include "CPALM.h"

#include "fit/GaussianFit2D.h"
#include "fit/GaussianFit2D_Sigma.h"
#include "fit/GaussianFit3D.h"
#include "fit/GaussianFit3D_Theta.h"

CPALM::CPALM(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size)
{
	long ret;
	nbCalls = 0;
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
	_size = size;

	for (int i = 0; i < 13; i++) _structFiltParam[i] = 0;
	for (int i = 0; i < 13; i++) _structFiltParamUse[i] = 0;
	for (int i = 0; i < 10; i++) _planesToUse[i] = 0;
	for (int i = 0; i < 9; i++) _operators[i] = 0;
	for (int i = 0; i < 9; i++) _threshold[i] = 0;


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
	if ((ret = _input->alloc(_image, _sizeX, _sizeY)) < 0) return;
	if ((ret = _aTrousPub->wDecomp->alloc(_sizeX, _sizeY, waveletNo + 2)) < 0) return;
	if ((ret = _segPub->realloc(_sizeX, _sizeY)) < 0) return;
	if ((ret = _outputPub->alloc(_sizeX, _sizeY)) < 0) return;


	_planesToUseB = new BlArray(9);

	_nbPlane = 0;
	for (int i = 0; i < 9; i++) {
		(*_planesToUseB)(i) = (_planesToUse[i] != 0);
		if ((*_planesToUseB)(i)) _nbPlane = i + 2;
	}
	(*_planesToUseB)(8) = false;
}

long CPALM::CPALMProcessing()
{
	//
	//FILE* f = fopen("C:\\debug.txt", "w");
	//fprintf(f, "After CPALMProcessing\n");
	//fclose(f);

	//FILE* f = fopen("C:\\image.txt", "a+");
	//for(int i=0; i<_sizeX*_sizeY; i++)
	//	fprintf(f, "%d\n", _image[i]);
	//fclose(f);

	//FILE* f = fopen("debug.txt", "a+");
	//fprintf(f, "Begin\n");

	int ret;
	nbCalls++;
	long nbStruct;
	int nbPointsSign;

	float begin = clock();

	//FILE* f = fopen("C:\\nbPlane.txt", "a+");
	//fprintf(f, "%d\n", _nbPlane);
	//fclose(f);

	if ((ret = _aTrousPub->aTrousDecomposition2D(*_input, *(_aTrousPub->wDecomp), *_planesToUseB, _nbPlane, false)) < 0) return (ret);

	//float _end = clock();
	//FILE* _fp = fopen("C:\\WaveletTimeCPU.txt", "a+");
	//fprintf(_fp, "Time: %f\n", (_end - begin) / (float)CLOCKS_PER_SEC );
	//fclose(_fp);

	//FILE *f = fopen("C:\\image_CPU.txt", "a+");
	for (long j = 0; j < long(_sizeX * _sizeY); j++) {
		long offs = _waveletNo * _sizeX * _sizeY;
		//fprintf(f,"%f\n", _aTrousPub->wDecomp->data[j]);
		//fprintf(f,"%d\n", _image[j]);
		if (_aTrousPub->wDecomp->data[offs + j] > _thresholdVal) {
			_outputPub->data[j] = _aTrousPub->wDecomp->data[offs + j];
			_segPub->data[j] = 1;
		}
		else {
			_segPub->data[j] = 0;
			_outputPub->data[j] = 0;
		}
	}
	//fclose(f);
	begin = clock();

	if (_watershedRatio == 1.) {
		if ((nbStruct = _segmentPub->segmentation(*_outputPub, *_segPub, 1e-8f, true, false)) < 0) return (nbStruct);
		//fprintf(f, "segmentation\n");
	}
	else {
		if ((nbStruct = _segmentPub->watershedIntelligent(*_outputPub, *_segPub, 1e-8f, (float)_watershedRatio)) < 0) return (nbStruct);
		//fprintf(f, "watershedIntelligent\n");
	}

	//FILE* f = fopen("C:\\log.txt", "a+");
	//fprintf(f, "%d\n", nbStruct);
	//fclose(f);
	//end = clock();

	//FILE* fp = fopen("C:\\Watershed.txt", "a+");
	//fprintf(fp, "Time: %f\n", (end - begin) / (float)CLOCKS_PER_SEC );
	//fclose(fp);

	if ((nbPointsSign = _segmentPub->computeStatistics(*_outputPub, *_segPub, _aTrousPub->wDecomp->data + _waveletNo * _sizeX * _sizeY, _gaussFit)) < 0) return (nbPointsSign);

	//fprintf(f, "computeStatistics\n");

	// for debug
	//FILE* f = fopen("test_CPU.txt", "w");
	//InfoStruct2D** ptr;
	//ptr = (InfoStruct2D**)_segmentPub->getTabInfo();
	//for(int i =1; i<= 36; i++)
	//	fprintf(f, "%f, %f, %f, %d \n", ptr[i]->data.moment.x, ptr[i]->data.moment.y, ptr[i]->data.intensity, (int)ptr[i]->data.surf);
	//fclose(f);
	// fin debug

	// for debug
	//FILE* f = fopen("test_CPU.txt", "w");
	//for(int i =0; i< 128; i++)
	//	for(int j=0; j<128; j++)
	//		fprintf(f, "%d\n", _segPub->data[i + 128*j]);
	//fclose(f);
	// fin debug

	begin = clock();

	if ((_volMin > 0.) || (_intMin > 0.)) if ((nbStruct = _segmentPub->filterData(*_segPub, _structFiltParam, _structFiltParamUse)) < 0) return (ret);
	//fprintf(f, "filterData\n");

	InfoStruct2D** p;
	p = (InfoStruct2D**)_segmentPub->getTabInfo();

	for (int i = 1; i <= min(nbStruct, _nbPointsMaxPerPlane - 2); i++) {
		switch (_gaussFit) {
			case 0: {
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = p[i]->data.moment.x;
				_pointsList[currentPoint++] = p[i]->data.moment.y;
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = p[i]->data.intensity;
				_pointsList[currentPoint++] = p[i]->data.surf;

				_pointsList[currentPoint++] = 0.; // Z coordinate
				_pointsList[currentPoint++] = 0.; // Spectral or other
				_pointsList[currentPoint++] = 0.; // Other

				break;
			}

			case 5: {
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = p[i]->data.moment.x;
				_pointsList[currentPoint++] = p[i]->data.moment.y;
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = p[i]->data.intensity;
				_pointsList[currentPoint++] = p[i]->data.surf;

				_pointsList[currentPoint++] = 0.; // Z coordinate
				_pointsList[currentPoint++] = 0.; // Spectral or other
				_pointsList[currentPoint++] = 0.; // Other

				break;
			}

			case 1: // Sigma is constant
			{
				float y0 = p[i]->data.moment.x;
				float x0 = p[i]->data.moment.y;
				float I0 = 0.;
				float sigmaX = _sigmaX;
				float offset = 0.;
				float RSS = -1.;

				/*begin = clock();*/

				// Gaussian fit of 2D
				GaussianFit2D* fit = new GaussianFit2D();
				fit->Run(_image, _sizeY, _sizeX, x0, y0, I0, sigmaX, offset, RSS, /*30*/10, _size);

				//end = clock();
				//fp = fopen("GaussianFit.txt", "a+");
				//fprintf(fp, "Time: %f\n", (end - begin) / (float)CLOCKS_PER_SEC );
				//fclose(fp);

				//if ((int)RSS != -1)
				//{
				//	sum += RSS;
				//	pointsNumber++;
				//}

				_pointsList[currentPoint++] = sigmaX;
				_pointsList[currentPoint++] = sigmaX;
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = y0;
				_pointsList[currentPoint++] = x0;
				_pointsList[currentPoint++] = I0;
				_pointsList[currentPoint++] = offset;
				_pointsList[currentPoint++] = RSS;
				_pointsList[currentPoint++] = p[i]->data.intensity;
				_pointsList[currentPoint++] = p[i]->data.surf;

				_pointsList[currentPoint++] = 0.; // Z coordinate
				_pointsList[currentPoint++] = 0.; // Spectral or other
				_pointsList[currentPoint++] = 0.; // Other

				break;
			}

			case 2: {
				float y0 = p[i]->data.moment.x;
				float x0 = p[i]->data.moment.y;
				float I0 = 0.;
				float sigmaX = _sigmaX;
				float offset;
				float RSS = -1.;

				/*begin = clock();*/

				// Gaussian fit of 2D
				GaussianFit2D_Sigma* fit = new GaussianFit2D_Sigma();
				fit->Run(_image, _sizeY, _sizeX, x0, y0, I0, sigmaX, offset, RSS, 10/*30*/, _size);

				//end = clock();
				//fp = fopen("GaussianFit.txt", "a+");
				//fprintf(fp, "Time: %f\n", (end - begin) / (float)CLOCKS_PER_SEC );
				//fclose(fp);

				//if ((int)RSS != -1)
				//{
				//	sum += RSS;
				//	pointsNumber++;

				//	FILE* f = fopen("RSS.txt", "a+");
				//	fprintf(f, "%f\n", RSS);
				//	fclose(f);
				//}

				_pointsList[currentPoint++] = sigmaX;
				_pointsList[currentPoint++] = sigmaX;
				_pointsList[currentPoint++] = 0.;
				_pointsList[currentPoint++] = y0;
				_pointsList[currentPoint++] = x0;
				_pointsList[currentPoint++] = I0;
				_pointsList[currentPoint++] = offset;
				_pointsList[currentPoint++] = RSS;
				_pointsList[currentPoint++] = p[i]->data.intensity;
				_pointsList[currentPoint++] = p[i]->data.surf;

				_pointsList[currentPoint++] = 0.; // Z coordinate
				_pointsList[currentPoint++] = 0.; // Spectral or other
				_pointsList[currentPoint++] = 0.; // Other

				break;
			}

			case 3: // Theta is constant
			{
				float y0 = p[i]->data.moment.x;
				float x0 = p[i]->data.moment.y;
				float I0 = 0./*p[i]->data.intensity*/; // p[i]->data.max.val - p[i]->data.min.val
				float sigmaX = _sigmaX;
				float sigmaY = _sigmaY;
				float angle = _angle;
				float offset = 0.; // p[i]->data.min.val
				float RSS = -1.;

				// Gaussian fit of 3D

				GaussianFit3D* fit = new GaussianFit3D();
				fit->Run(_image, _sizeY, _sizeX, x0, y0, I0, sigmaX, sigmaY, angle, offset, RSS, /*20*//*50*//*20 * 6*/10/*7 * 20*/, _size);

				//if ((int)RSS != -1)
				//{
				//	sum += RSS;
				//	pointsNumber++;
				//}

				_pointsList[currentPoint++] = sigmaX;
				_pointsList[currentPoint++] = sigmaY;
				_pointsList[currentPoint++] = angle;
				_pointsList[currentPoint++] = y0;
				_pointsList[currentPoint++] = x0;
				_pointsList[currentPoint++] = I0;
				_pointsList[currentPoint++] = offset;
				_pointsList[currentPoint++] = RSS;
				_pointsList[currentPoint++] = p[i]->data.intensity;
				_pointsList[currentPoint++] = p[i]->data.surf;

				_pointsList[currentPoint++] = 0.; // Z coordinate
				_pointsList[currentPoint++] = 0.; // Spectral or other
				_pointsList[currentPoint++] = 0.; // Other

				break;
			}

			case 4: {
				float y0 = p[i]->data.moment.x;
				float x0 = p[i]->data.moment.y;
				float I0 = 0./*p[i]->data.intensity*/; // p[i]->data.max.val - p[i]->data.min.val
				float sigmaX = _sigmaX;
				float sigmaY = _sigmaY;
				float angle = _angle;
				float offset = 0.; // p[i]->data.min.val
				float RSS = -1.;

				// Gaussian fit of 3D

				GaussianFit3D_Theta* fit = new GaussianFit3D_Theta();
				fit->Run(_image, _sizeY, _sizeX, x0, y0, I0, sigmaX, sigmaY, angle, offset, RSS, 20 * 7, _size);

				//if ((int)RSS != -1)
				//{
				//	sum += RSS;
				//	pointsNumber++;

				//	FILE* f = fopen("RSS.txt", "a+");
				//	fprintf(f, "%f\n", RSS);
				//	fclose(f);
				//}

				_pointsList[currentPoint++] = sigmaX;
				_pointsList[currentPoint++] = sigmaY;
				_pointsList[currentPoint++] = angle;
				_pointsList[currentPoint++] = y0;
				_pointsList[currentPoint++] = x0;
				_pointsList[currentPoint++] = I0;
				_pointsList[currentPoint++] = offset;
				_pointsList[currentPoint++] = RSS;
				_pointsList[currentPoint++] = p[i]->data.intensity;
				_pointsList[currentPoint++] = p[i]->data.surf;

				_pointsList[currentPoint++] = 0.; // Z coordinate
				_pointsList[currentPoint++] = 0.; // Spectral or other
				_pointsList[currentPoint++] = 0.; // Other

				break;
			}

			case 6: // Theta is constant --> for 3D calibration
			{
				float y0 = p[i]->data.moment.x;
				float x0 = p[i]->data.moment.y;
				float I0 = 0./*p[i]->data.intensity*/; // p[i]->data.max.val - p[i]->data.min.val
				float sigmaX = _sigmaX;
				float sigmaY = _sigmaY;
				float angle = _angle;
				float offset = 0.; // p[i]->data.min.val
				float RSS = -1.;

				// Gaussian fit of 3D

				GaussianFit3D* fit = new GaussianFit3D();
				fit->Run(_image, _sizeY, _sizeX, x0, y0, I0, sigmaX, sigmaY, angle, offset, RSS, 7 * 20, _size);

				//if ((int)RSS != -1)
				//{
				//	sum += RSS;
				//	pointsNumber++;
				//}

				_pointsList[currentPoint++] = sigmaX;
				_pointsList[currentPoint++] = sigmaY;
				_pointsList[currentPoint++] = angle;
				_pointsList[currentPoint++] = y0;
				_pointsList[currentPoint++] = x0;
				_pointsList[currentPoint++] = I0;
				_pointsList[currentPoint++] = offset;
				_pointsList[currentPoint++] = RSS;
				_pointsList[currentPoint++] = p[i]->data.intensity;
				_pointsList[currentPoint++] = p[i]->data.surf;

				_pointsList[currentPoint++] = 0.; // Z coordinate
				_pointsList[currentPoint++] = 0.; // Spectral or other
				_pointsList[currentPoint++] = 0.; // Other

				break;
			}
		}
	}

	_pointsList[currentPoint++] = -1;
	_pointsList[currentPoint++] = -1;
	_pointsList[currentPoint++] = -1;
	_pointsList[currentPoint++] = -1;
	_pointsList[currentPoint++] = -1;
	_pointsList[currentPoint++] = -1;
	_pointsList[currentPoint++] = -1;
	_pointsList[currentPoint++] = -1;
	_pointsList[currentPoint++] = -1;
	_pointsList[currentPoint++] = -1;

	_pointsList[currentPoint++] = -1; // Z coordinate
	_pointsList[currentPoint++] = -1; // Spectral or other
	_pointsList[currentPoint++] = -1; // Other

	//fprintf(f, "End\n");
	//fclose(f);

	return (currentPoint);
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
