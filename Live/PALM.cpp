// PALM.cpp : définit le point d'entrée pour l'application DLL.
//

//#include <GL\glut.h>
#include <cuda_gl_interop.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <helper_cuda.h>
#include <helper_cuda_gl.h>
#include <helper_functions.h>
#include <math_constants.h>

//#include "PALM.h"
#include "GaussianFit2D.h"
#include "GaussianFit2D_Sigma.h"
#include "GaussianFit3D.h"
#include "GaussianFit3D_Theta.h"
#include "PALM_Class.h"
#include "palm/parameters.h"

#ifdef _MANAGED
#pragma managed(push, off)
#endif

#define min(a,b) (((a) < (b)) ? (a) : (b))


#ifdef _MANAGED
#pragma managed(pop)
#endif


CPALM* c;

double sum;
unsigned int pointsNumber;

// for debug
//unsigned int* tab;
//unsigned int numImage;
/**********************/


CPALM::CPALM(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short roi)
{
	long ret;
	nbCalls = 0;
	currentPoint = 0;
	_image = image;
	_pointsList = pointsList;
	_pointsListSize = pointsListSize;
	_nbSegParams = 13;
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

	for (int i = 0; i < 13; i++) { _structFiltParam[i] = 0; }
	for (int i = 0; i < 13; i++) { _structFiltParamUse[i] = 0; }
	for (int i = 0; i < 10; i++) { _planesToUse[i] = 0; }
	for (int i = 0; i < 9; i++) { _operators[i] = 0; }
	for (int i = 0; i < 9; i++) { _threshold[i] = 0; }


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
	if ((ret = _input->alloc(_image, _sizeX, _sizeY)) < 0) { return; }
	if ((ret = _aTrousPub->wDecomp->alloc(_sizeX, _sizeY, waveletNo + 2)) < 0) { return; }
	if ((ret = _segPub->realloc(_sizeX, _sizeY)) < 0) { return; }
	if ((ret = _outputPub->alloc(_sizeX, _sizeY)) < 0) { return; }


	_planesToUseB = new BlArray(9);

	_nbPlane = 0;
	for (int i = 0; i < 9; i++) {
		(*_planesToUseB)(i) = (_planesToUse[i] != 0);
		if ((*_planesToUseB)(i)) { _nbPlane = i + 2; }
	}
	(*_planesToUseB)(8) = false;
}

unsigned long twoToThe_(unsigned long val)
{
	if (val == 0) { return (1); }
	return (((unsigned long)1) << val);
}

int aTrousDecomposition2D(unsigned short* d_inputData, float* d_outputData, unsigned short sizeDataX, unsigned short sizeDataY, unsigned short sizeStk, /*BlArray &planesToUse,*/ unsigned long nbPlane, bool realProj, bool planC)
{
	unsigned short step;
	unsigned short currentPlane;

	//
	float* d_buffer0;
	checkCudaErrors(cudaMalloc((void**)&d_buffer0, sizeStk * sizeDataX * sizeDataY * sizeof(float)));

	cudaProjectionV02DGPU(d_inputData, d_buffer0, sizeStk, sizeDataX, sizeDataY, realProj);

	float* d_buffer1;
	checkCudaErrors(cudaMalloc((void**)&d_buffer1, sizeStk * sizeDataX * sizeDataY * sizeof(float)));

	for (currentPlane = 0; currentPlane < nbPlane - 1; currentPlane++) {
		step = twoToThe_(currentPlane);

		cudaTrous2DGPU(d_buffer0, d_buffer1, sizeStk, sizeDataX, sizeDataY, step);

		checkCudaErrors(cudaMemcpy(d_outputData + currentPlane * sizeStk * sizeDataX * sizeDataY, d_buffer0, sizeStk * sizeDataX * sizeDataY * sizeof(float), cudaMemcpyDeviceToDevice));
		checkCudaErrors(cudaMemcpy(d_buffer0, d_buffer1, sizeStk * sizeDataX * sizeDataY * sizeof(float), cudaMemcpyDeviceToDevice));
	}

	if (planC == true) {
		// Copier le dernier plan lissé Cn
		checkCudaErrors(cudaMemcpy(d_outputData + currentPlane * sizeStk * sizeDataX * sizeDataY, d_buffer1, sizeStk * sizeDataX * sizeDataY * sizeof(float), cudaMemcpyDeviceToDevice));

		float* data = (float*)malloc(sizeDataX * sizeDataY * sizeof(float));
		checkCudaErrors(cudaMemcpy(data, d_buffer1, sizeDataX * sizeDataY * sizeof(float), cudaMemcpyDeviceToHost));
		free(data);
	}

	// cleanup memory
	checkCudaErrors(cudaFree(d_buffer0));
	checkCudaErrors(cudaFree(d_buffer1));

	return 1;
}


long CPALM::CPALMProcessing()
{
	int ret;
	nbCalls++;
	long nbStruct;
	int nbPointsSign;

	float begin = clock();
	// copie des images : host --> device
	unsigned short* d_images;
	checkCudaErrors(cudaMalloc((void**)&d_images, _sizeX * _sizeY * sizeof(unsigned short)));
	checkCudaErrors(cudaMemcpy(d_images, _input->data, _sizeX * _sizeY * sizeof(unsigned short), cudaMemcpyHostToDevice));

	// memory allocation for : d_output
	float* d_data;
	checkCudaErrors(cudaMalloc((void**)&d_data, (_nbPlane) * _sizeX * _sizeY * sizeof(float)));

	// Decompostition : wavelet algorithm
	aTrousDecomposition2D(d_images, d_data, _sizeY, _sizeX, 1, /**_planesToUseB,*/ _nbPlane, false/*true*/, false);

	checkCudaErrors(cudaMemcpy(_aTrousPub->wDecomp->data + _waveletNo * _sizeX * _sizeY, d_data + _waveletNo * _sizeX * _sizeY, _sizeX * _sizeY * sizeof(float), cudaMemcpyDeviceToHost));

	// cleanup memory
	checkCudaErrors(cudaFree(d_images));
	checkCudaErrors(cudaFree(d_data));

	float end = clock();

	begin = clock();
	for (long j = 0; j < long(_sizeX * _sizeY); j++) {
		long offs = _waveletNo * _sizeX * _sizeY;

		if (_aTrousPub->wDecomp->data[offs + j] > _thresholdVal) {
			_outputPub->data[j] = _aTrousPub->wDecomp->data[offs + j];
			_segPub->data[j] = 1;
		}
		else {
			_segPub->data[j] = 0;
			_outputPub->data[j] = 0;
		}
	}
	end = clock();
	begin = clock();

	if (_watershedRatio == 1.) { if ((nbStruct = _segmentPub->segmentation(*_outputPub, *_segPub, 1e-8f, true, false)) < 0) { return (nbStruct); } }
	else { if ((nbStruct = _segmentPub->watershedIntelligent(*_outputPub, *_segPub, 1e-8f, (float)_watershedRatio)) < 0) { return (nbStruct); } }

	end = clock();;
	begin = clock();


	if ((nbPointsSign = _segmentPub->computeStatistics(*_outputPub, *_segPub, _aTrousPub->wDecomp->data + _waveletNo * _sizeX * _sizeY, _gaussFit)) < 0) { return (nbPointsSign); }

	begin = clock();

	if ((_volMin > 0.) || (_intMin > 0.)) { if ((nbStruct = _segmentPub->filterData(*_segPub, _structFiltParam, _structFiltParamUse)) < 0) { return (ret); } }

	end = clock();

	InfoStruct2D** p;
	p = (InfoStruct2D**)_segmentPub->getTabInfo();

	// le nombre de points détectés ne doit pas dépasser _nbPointsMaxPerPlane --> problème d'espace mémoire
	for (int i = 1; i <= min(nbStruct, _nbPointsMaxPerPlane - 2); i++) {
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

	_pointsList[currentPoint++] = -1;
	_pointsList[currentPoint++] = -1;
	_pointsList[currentPoint++] = -1;

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

_declspec(dllexport) void _stdcall __OpenPALMProcessing(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size)
{
	sum = 0.;
	pointsNumber = 0;

	c = new CPALM(image, pointsList, pointsListSize, sizeX, sizeY, waveletNo, thresholdVal, watershedRatio, volMin, intMin, gaussFit, sigmaX, sigmaY, angle, size);
}

_declspec(dllexport) long _stdcall __PALMProcessing() { return c->CPALMProcessing(); }

_declspec(dllexport) void _stdcall __closePALMProcessing() { delete c; }


unsigned short* stack;

_declspec(dllexport) void _stdcall AllocateMemory(unsigned short frameCount, unsigned short width, unsigned short height)
{
	//stack = (unsigned short*)malloc(frameCount * width * height * sizeof(unsigned short));
	stack = new unsigned short[frameCount * width * height];
}

_declspec(dllexport) void _stdcall CopyFrame(unsigned short numImage, unsigned short* image, unsigned short frameCount, unsigned short width, unsigned short height) { memcpy(stack + numImage * width * height, image, width * height * sizeof(unsigned short)); }

_declspec(dllexport) void _stdcall InitializeMemory(unsigned short* frames) { stack = frames; }


_declspec(dllexport) void _stdcall LaunchGaussianFit(double* pointsList, unsigned long pointsListSize, unsigned short sizeStk, unsigned short sizeDataX, unsigned short sizeDataY, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size, unsigned short iterNum, unsigned short fromcplusplus, unsigned short computeZ, double* SigmaXCoefs, double* SigmaYCoefs, double ThreeDThickness/*µm*/, double ThreeDSampling/*µm*/, double pixelSize /*µm*/)
{
	if ((gaussFit >= 1) && (gaussFit <= 4)) {
		// Checking the available memory
		cudaDeviceProp deviceProp;
		cudaGetDeviceProperties(&deviceProp, 0);
		int kernelExecTimeoutEnabled = deviceProp.kernelExecTimeoutEnabled;
		int deviceVersion = deviceProp.major;
		unsigned long long totalGlobalMem = (unsigned long long)deviceProp.totalGlobalMem;

		int newSizeStk, iter, r;

		newSizeStk = (totalGlobalMem * 0.8) / (sizeDataX * sizeDataY * (/*sizeof(float) + */2 * sizeof(unsigned short)));
		iter = 0;
		r = 0;

		// A vérifier
		if (kernelExecTimeoutEnabled == 1) { newSizeStk = 80; }

		//newSizeStk = 80;

		{
			iter = sizeStk / newSizeStk;
			r = sizeStk % newSizeStk;
		}

		// Fill nbStruct
		unsigned long NbSegParams = 13;

		unsigned int* cumNbStruct = (unsigned int*)malloc((sizeStk + 1) * sizeof(unsigned int));
		unsigned short* nbStruct = (unsigned short*)malloc(sizeStk * sizeof(unsigned short));

		unsigned short numFrame = 0;
		nbStruct[numFrame] = 0;
		cumNbStruct[numFrame] = 0;

		unsigned int nbPoints = 0;

		for (int i = 0; i < (pointsListSize / NbSegParams); i++) {
			if ((int)pointsList[NbSegParams * i] == -1) {
				nbPoints += nbStruct[numFrame];

				numFrame++;
				if (numFrame == sizeStk) { break; }
				nbStruct[numFrame] = 0;
			}
			else { nbStruct[numFrame]++; }
		}

		// Checking the available maximum sizes of each dimension of a grid (GPU launching)
		if ((int)ceil((float)(nbPoints) / (float)FIT_BLOCK_SIZE) > deviceProp.maxGridSize[0]) {
			newSizeStk = (int)(ceil((float)(deviceProp.maxGridSize[0] * FIT_BLOCK_SIZE) / (float)(nbPoints / sizeStk)));
			newSizeStk -= 10;
			iter = sizeStk / newSizeStk;
			r = sizeStk % newSizeStk;
		}

		float* list = (float*)malloc(nbPoints * NbSegParams * sizeof(float));

		int ind_src;
		int ind_dst;

		if (fromcplusplus == 1) {
			for (int plane = 0; plane < sizeStk; plane++) {
				for (int i = 0; i < nbStruct[plane]; i++) {
					ind_src = cumNbStruct[plane] + plane + i;
					ind_dst = cumNbStruct[plane] + i;

					list[ind_dst * NbSegParams + 0] = (float)pointsList[NbSegParams * ind_src + 3];
					list[ind_dst * NbSegParams + 1] = (float)pointsList[NbSegParams * ind_src + 4];
					list[ind_dst * NbSegParams + 2] = (float)pointsList[NbSegParams * ind_src + 8];
					list[ind_dst * NbSegParams + 3] = (float)pointsList[NbSegParams * ind_src + 9];
					list[ind_dst * NbSegParams + 4] = 0.;
					list[ind_dst * NbSegParams + 5] = 0.;
					list[ind_dst * NbSegParams + 6] = 0.;
					list[ind_dst * NbSegParams + 7] = 0.;
					list[ind_dst * NbSegParams + 8] = 0.;
					list[ind_dst * NbSegParams + 9] = (float)i;
					list[ind_dst * NbSegParams + 10] = (float)(plane % newSizeStk);
					list[ind_dst * NbSegParams + 11] = 0.;
					list[ind_dst * NbSegParams + 12] = 0.;
				}

				if ((plane + 1) <= sizeStk) { cumNbStruct[plane + 1] = cumNbStruct[plane] + nbStruct[plane]; }
			}
		}
		else {
			for (int plane = 0; plane < sizeStk; plane++) {
				for (int i = 0; i < nbStruct[plane]; i++) {
					ind_src = cumNbStruct[plane] + plane + i;
					ind_dst = cumNbStruct[plane] + i;

					list[ind_dst * NbSegParams + 0] = (float)pointsList[NbSegParams * ind_src + 4];
					list[ind_dst * NbSegParams + 1] = (float)pointsList[NbSegParams * ind_src + 3];
					list[ind_dst * NbSegParams + 2] = (float)pointsList[NbSegParams * ind_src + 8];
					list[ind_dst * NbSegParams + 3] = (float)pointsList[NbSegParams * ind_src + 9];
					list[ind_dst * NbSegParams + 4] = 0.;
					list[ind_dst * NbSegParams + 5] = 0.;
					list[ind_dst * NbSegParams + 6] = 0.;
					list[ind_dst * NbSegParams + 7] = 0.;
					list[ind_dst * NbSegParams + 8] = 0.;
					list[ind_dst * NbSegParams + 9] = (float)i;
					list[ind_dst * NbSegParams + 10] = (float)(plane % newSizeStk);
					list[ind_dst * NbSegParams + 11] = 0.;
					list[ind_dst * NbSegParams + 12] = 0.;
				}

				if ((plane + 1) <= sizeStk) { cumNbStruct[plane + 1] = cumNbStruct[plane] + nbStruct[plane]; }
			}
		}

		float* d_SigmaXCoefs;
		float* d_SigmaYCoefs;

		// 5: number of parameters of calibration equation
		float _SigmaXCoefs[5];
		float _SigmaYCoefs[5];

		for (int i = 0; i < 5; i++) {
			_SigmaXCoefs[i] = (float)SigmaXCoefs[i];
			_SigmaYCoefs[i] = (float)SigmaYCoefs[i];
		}

		checkCudaErrors(cudaMalloc((void**)&d_SigmaXCoefs, 5 * sizeof(float)));
		checkCudaErrors(cudaMalloc((void**)&d_SigmaYCoefs, 5 * sizeof(float)));

		checkCudaErrors(cudaMemcpy(d_SigmaXCoefs, _SigmaXCoefs, 5 * sizeof(float), cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_SigmaYCoefs, _SigmaYCoefs, 5 * sizeof(float), cudaMemcpyHostToDevice));

		// Loop for computing
		for (int i = 1; i <= iter; i++) {
			// For test
			unsigned short *d_images, *d_nbStruct;
			float* d_list;
			checkCudaErrors(cudaMalloc((void**)&d_images, newSizeStk * sizeDataX * sizeDataY * sizeof(unsigned short)));
			checkCudaErrors(cudaMalloc((void**)&d_nbStruct, newSizeStk * sizeof(unsigned short)));
			checkCudaErrors(cudaMalloc((void**)&d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]) * NbSegParams * sizeof(float)));

			checkCudaErrors(cudaMemcpy(d_images, stack + (i - 1) * newSizeStk * sizeDataX * sizeDataY, newSizeStk * sizeDataX * sizeDataY * sizeof(unsigned short), cudaMemcpyHostToDevice));
			checkCudaErrors(cudaMemcpy(d_nbStruct, nbStruct + (i - 1) * newSizeStk, newSizeStk * sizeof(unsigned short), cudaMemcpyHostToDevice));
			checkCudaErrors(cudaMemcpy(d_list, list + cumNbStruct[(i - 1) * newSizeStk] * NbSegParams, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]) * NbSegParams * sizeof(float), cudaMemcpyHostToDevice));

			// cleanup memory

			if (gaussFit == 1) // 2D with constant Sigma
			{
				if (size < 7) { cuda2DFittingGPU_5x5(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*30*/, size, sigmaX); }

				if (size >= 7) { cuda2DFittingGPU_7x7(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*30*/, size, sigmaX); }
			}

			if (gaussFit == 2) // 2D with Sigma
			{
				if (size < 7) { cuda2DFittingGPU_Sigma_5x5(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*30*/, size, sigmaX); }

				if (size >= 7) { cuda2DFittingGPU_Sigma_7x7(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*30*/, size, sigmaX); }
			}

			if (gaussFit == 3) // 3D with constant Theta
			{
				if ((size <= 7) || (size == 8)) { cuda3DFittingGPU_7x7(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/); }

				if ((size == 9) || (size == 10)) { cuda3DFittingGPU_9x9(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/); }

				if ((size == 11) || (size == 12)) { cuda3DFittingGPU_11x11(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/); }

				if (size >= 13) { cuda3DFittingGPU_13x13(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/); }
			}
			if (gaussFit == 4) // 3D with Theta
			{
				if (size <= 9) { cuda3DFittingGPU_Theta_9x9(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum /*20*7*/, size, sigmaX, sigmaY, angle); }
				if (size > 9) { cuda3DFittingGPU_Theta_13x13(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum /*20*7*/, size, sigmaX, sigmaY, angle); }
			}

			cudaError err = cudaThreadSynchronize();
			checkCudaErrors(cudaMemcpy(list + cumNbStruct[(i - 1) * newSizeStk] * NbSegParams, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]) * NbSegParams * sizeof(float), cudaMemcpyDeviceToHost));

			// Cleanup memory
			checkCudaErrors(cudaFree(d_images));
			checkCudaErrors(cudaFree(d_list));
			checkCudaErrors(cudaFree(d_nbStruct));
		}

		if (r != 0) {
			unsigned short *d_images, *d_nbStruct;
			float* d_list;
			checkCudaErrors(cudaMalloc((void**)&d_images, (sizeStk % newSizeStk) * sizeDataX * sizeDataY * sizeof(unsigned short)));
			checkCudaErrors(cudaMalloc((void**)&d_nbStruct, (sizeStk % newSizeStk) * sizeof(unsigned short)));
			checkCudaErrors(cudaMalloc((void**)&d_list, (nbPoints - cumNbStruct[iter * newSizeStk]) * NbSegParams * sizeof(float)));

			checkCudaErrors(cudaMemcpy(d_images, stack + iter * newSizeStk * sizeDataX * sizeDataY, (sizeStk % newSizeStk) * sizeDataX * sizeDataY * sizeof(unsigned short), cudaMemcpyHostToDevice));
			checkCudaErrors(cudaMemcpy(d_nbStruct, nbStruct + iter * newSizeStk, (sizeStk % newSizeStk) * sizeof(unsigned short), cudaMemcpyHostToDevice));
			checkCudaErrors(cudaMemcpy(d_list, list + cumNbStruct[iter * newSizeStk] * NbSegParams, (nbPoints - cumNbStruct[iter * newSizeStk]) * NbSegParams * sizeof(float), cudaMemcpyHostToDevice));

			//convertUSrtToFloat(d_images, d_data, (sizeStk % newSizeStk), sizeDataX, sizeDataY);

			// cleanup memory
			//checkCudaErrors( cudaFree(d_images) );

			if (gaussFit == 1) // 2D with constant Sigma
			{
				if (size < 7) { cuda2DFittingGPU_5x5(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*30*/, size, sigmaX); }

				if (size >= 7) { cuda2DFittingGPU_7x7(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*30*/, size, sigmaX); }
			}

			if (gaussFit == 2) // 2D with Sigma
			{
				if (size < 7) { cuda2DFittingGPU_Sigma_5x5(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*30*/, size, sigmaX); }

				if (size >= 7) { cuda2DFittingGPU_Sigma_7x7(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*30*/, size, sigmaX); }
			}

			if (gaussFit == 3) // 3D with constant Theta
			{
				if ((size <= 7) || (size == 8)) { cuda3DFittingGPU_7x7(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/); }
				if ((size == 9) || (size == 10)) { cuda3DFittingGPU_9x9(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/); }
				if ((size == 11) || (size == 12)) { cuda3DFittingGPU_11x11(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/); }
				if (size >= 13) { cuda3DFittingGPU_13x13(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/); }
			}
			if (gaussFit == 4) // 3D with Theta
			{
				if (size <= 9) { cuda3DFittingGPU_Theta_9x9(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*20*7*/, size, sigmaX, sigmaY, angle); }
				if (size > 9) { cuda3DFittingGPU_Theta_13x13(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*20*7*/, size, sigmaX, sigmaY, angle); }
			}

			cudaError err = cudaThreadSynchronize();
			checkCudaErrors(cudaMemcpy(list + cumNbStruct[iter * newSizeStk] * NbSegParams, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]) * NbSegParams * sizeof(float), cudaMemcpyDeviceToHost));

			// Cleanup memory
			checkCudaErrors(cudaFree(d_images));
			checkCudaErrors(cudaFree(d_list));
			checkCudaErrors(cudaFree(d_nbStruct));
		};

		// Cleanup memory
		checkCudaErrors(cudaFree(d_SigmaXCoefs));
		checkCudaErrors(cudaFree(d_SigmaYCoefs));

		int ind_list = 0;
		int ind_pointsList = 0;

		if (fromcplusplus != 1) {
			for (int i = 0; i < sizeStk /*host*/; i++) {
				for (int s = 0; s < nbStruct[i]; s++) {
					pointsList[ind_pointsList++] = (double)list[ind_list++];   // sigmaX
					pointsList[ind_pointsList++] = (double)list[ind_list++];   // sigmaY
					pointsList[ind_pointsList++] = (double)list[ind_list++];   // Theta
					pointsList[ind_pointsList++] = (double)list[ind_list + 1]; // X
					pointsList[ind_pointsList++] = (double)list[ind_list];     // Y
					ind_list += 2;
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Intensity
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Offset
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Error
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // I
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // S

					pointsList[ind_pointsList++] = (double)list[ind_list++]; //0.; // Z coordinate
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Spectral or other
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Other
				}
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;

				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
			}
		}
		else {
			for (int i = 0; i < sizeStk /*host*/; i++) {
				//ind_list = (i * maxNbStruct/* + 1*/) * NbSegParams;
				//while((int)list[ind_list + 2] != 0)
				for (int s = 0; s < nbStruct[i]; s++) {
					//pointsList[ind_pointsList++] = list[ind_list++];  // x
					//pointsList[ind_pointsList++] = list[ind_list++];  // y

					pointsList[ind_pointsList++] = (double)list[ind_list++]; // sigmaX
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // sigmaY

					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Theta


					pointsList[ind_pointsList++] = (double)list[ind_list++]; // X
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Y
					//ind_list += 2;
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Intensity
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Offset
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Error
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // I
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // S

					pointsList[ind_pointsList++] = (double)list[ind_list++]; //0.; // Z coordinate
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Spectral or other
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Other
					//ind_list++; // Z coordinate

					/*ind_list += 4;*/
				}
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;

				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
				pointsList[ind_pointsList++] = -1;
			}
		}

		// cleanup memory
		free(cumNbStruct);
		free(nbStruct);
		free(list);

		if (fromcplusplus != 1) { delete[] stack; }
	}
}

_declspec(dllexport) void _stdcall GaussianVisuGPU(double* gaussianPointsList/*(sizeReg*sizeReg*nbPoints)*/, double* pointsList, unsigned long pointsListSize, unsigned short sizeStk, unsigned short width, unsigned short height, unsigned short sizeReg, unsigned short zoom, double conversionRatio /*from: ADU to photons*/, double pixelSize /* unit: nm */)
{
	unsigned int nbPoints = 0;
	unsigned long NbSegParams = 13;

	unsigned int* cumNbStruct = (unsigned int*)malloc((sizeStk + 1) * sizeof(unsigned int));
	unsigned short* nbStruct = (unsigned short*)malloc(sizeStk * sizeof(unsigned short));

	unsigned short numFrame = 0;
	nbStruct[numFrame] = 0;
	cumNbStruct[numFrame] = 0;

	for (int i = 0; i < (pointsListSize / NbSegParams); i++) {
		if ((int)pointsList[NbSegParams * i] == -1) {
			nbPoints += nbStruct[numFrame];

			numFrame++;
			if (numFrame == sizeStk) { break; }
			nbStruct[numFrame] = 0;
		}
		else { nbStruct[numFrame]++; }
	}

	double* list = (double*)malloc(nbPoints * NbSegParams * sizeof(double));

	int ind_src;
	int ind_dst;
	for (int plane = 0; plane < sizeStk; plane++) {
		for (int i = 0; i < nbStruct[plane]; i++) {
			ind_src = cumNbStruct[plane] + plane + i;
			ind_dst = cumNbStruct[plane] + i;

			list[ind_dst * NbSegParams + 0] = pointsList[NbSegParams * ind_src + 4];
			list[ind_dst * NbSegParams + 1] = pointsList[NbSegParams * ind_src + 3];
			list[ind_dst * NbSegParams + 2] = 0.;
			list[ind_dst * NbSegParams + 3] = 10.;
			list[ind_dst * NbSegParams + 4] = 0.;
			list[ind_dst * NbSegParams + 5] = 0.;
			list[ind_dst * NbSegParams + 6] = 0.;
			list[ind_dst * NbSegParams + 7] = 0.;
			list[ind_dst * NbSegParams + 8] = 0.;
			list[ind_dst * NbSegParams + 9] = pointsList[NbSegParams * ind_src + 0];  // sigamX
			list[ind_dst * NbSegParams + 10] = pointsList[NbSegParams * ind_src + 5]; // Gaussian intensity
		}

		if ((plane + 1) <= sizeStk) { cumNbStruct[plane + 1] = cumNbStruct[plane] + nbStruct[plane]; }
	}

	double *d_gaussianPointsList, *d_pointsList;

	checkCudaErrors(cudaMalloc((void**)&d_gaussianPointsList, sizeReg * sizeReg * nbPoints * sizeof(double)));
	checkCudaErrors(cudaMalloc((void**)&d_pointsList, nbPoints * NbSegParams * sizeof(double)));

	checkCudaErrors(cudaMemcpy(d_gaussianPointsList, gaussianPointsList, sizeReg * sizeReg * nbPoints * sizeof(double), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_pointsList, list, nbPoints * NbSegParams * sizeof(double), cudaMemcpyHostToDevice));

	GaussianVisu(d_gaussianPointsList /*(sizeReg*sizeReg*nbPoints)*/, d_pointsList, width, height, nbPoints, sizeReg, zoom, conversionRatio, pixelSize);

	checkCudaErrors(cudaMemcpy(gaussianPointsList, d_gaussianPointsList, sizeReg * sizeReg * nbPoints * sizeof(double), cudaMemcpyDeviceToHost));

	checkCudaErrors(cudaFree(d_gaussianPointsList));
	checkCudaErrors(cudaFree(d_pointsList));
}
