// PALM.cpp : définit le point d'entrée pour l'application DLL.
//

//#include <GL\glut.h>
#include "stdafx.h"

#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include <cufft.h>

#include <helper_cuda.h>
#include <helper_cuda_gl.h>

#include <helper_functions.h>
#include <math_constants.h>

//#include "PALM.h"
#include "PALM_Class.h"

#include "GaussianFit2D.h"
#include "GaussianFit2D_Sigma.h"
#include "GaussianFit3D.h"
#include "GaussianFit3D_Theta.h"

#include "parameters.h"

#ifdef _MANAGED
#pragma managed(push, off)
#endif

#define min(a,b) (((a) < (b)) ? (a) : (b))

BOOL APIENTRY DllMain(HMODULE hModule, DWORD ul_reason_for_call, LPVOID lpReserved)
{
	switch (ul_reason_for_call) {
		case DLL_PROCESS_ATTACH:
		case DLL_THREAD_ATTACH:
		case DLL_THREAD_DETACH:
		case DLL_PROCESS_DETACH: break;
	}
	return TRUE;
}

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


CPALM::CPALM(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size)
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

	/*********************/
	//tab = (unsigned int*)malloc(10000 * sizeof(unsigned int));
	//numImage = 0;
	/*********************/
}

unsigned long twoToThe_(unsigned long val)
{
	if (val == 0) return (1);
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
		//FILE *f = fopen("C:\\C1_GPU.txt", "w");
		//for(int i=0; i<sizeDataX*sizeDataY; i++)
		//	fprintf(f, "%f\n", data[i]);
		//fclose(f);
		free(data);
	}

	// cleanup memory
	checkCudaErrors(cudaFree(d_buffer0));
	checkCudaErrors(cudaFree(d_buffer1));

	return 1;
}


long CPALM::CPALMProcessing()
{
	//
	//FILE* f = fopen("C:\\debug.txt", "w");
	//fprintf(f, "After CPALMProcessing\n");
	//fclose(f);

	//FILE* f = fopen("debug.txt", "a+");
	//fprintf(f, "Begin\n");

	/*********************/
	//numImage++;
	//unsigned int _sum = 0;
	//   for(int i = 0; i < _sizeY; i++) 
	//{
	//	for(int j = 0; j < _sizeX; j++)
	//	{
	//		_sum += _input->data[i * _sizeX + j];
	//	}
	//}
	//tab[numImage-1] = _sum;
	/*********************/

	int ret;
	nbCalls++;
	long nbStruct;
	int nbPointsSign;

	float begin = clock();

	//if( (ret=_aTrousPub->aTrousDecomposition2D(*_input, *(_aTrousPub->wDecomp), *_planesToUseB, _nbPlane, false))<0) 
	//	return(ret);

	//goto suite;
	;
	/*********************************************************************/
	// Wavelet on GPU

	// copie des images : host --> device
	unsigned short* d_images;
	checkCudaErrors(cudaMalloc((void**)&d_images, _sizeX * _sizeY * sizeof(unsigned short)));
	//f = fopen("C:\\debug.txt", "a+");
	//fprintf(f, "Begin\n");
	//fclose(f);
	checkCudaErrors(cudaMemcpy(d_images, _input->data, _sizeX * _sizeY * sizeof(unsigned short), cudaMemcpyHostToDevice));

	// memory allocation for : d_output
	float* d_data;
	//checkCudaErrors( cudaMalloc( (void**) &d_data, (_nbPlane-1) * _sizeX * _sizeY * sizeof(float) ) );
	checkCudaErrors(cudaMalloc((void**)&d_data, (_nbPlane) * _sizeX * _sizeY * sizeof(float)));

	// Decompostition : wavelet algorithm
	aTrousDecomposition2D(d_images, d_data, _sizeY, _sizeX, 1, /**_planesToUseB,*/ _nbPlane, false/*true*/, false);

	checkCudaErrors(cudaMemcpy(_aTrousPub->wDecomp->data + _waveletNo * _sizeX * _sizeY, d_data + _waveletNo * _sizeX * _sizeY, _sizeX * _sizeY * sizeof(float), cudaMemcpyDeviceToHost));

	////////checkCudaErrors( cudaMemcpy( _aTrousPub->wDecomp->data, d_data, (_nbPlane-1)* _sizeX * _sizeY * sizeof(float), cudaMemcpyDeviceToHost) );

	// cleanup memory
	checkCudaErrors(cudaFree(d_images));
	checkCudaErrors(cudaFree(d_data));

	/*********************************************************************/
	//fprintf(f, "wavelet\n");

	float end = clock();
	//FILE* fp = fopen("C:\\WaveletTimeGPU.txt", "a+");
	//fprintf(fp, "Time: %f\n", (end - begin) / (float)CLOCKS_PER_SEC );
	//fclose(fp);

suite:
	//return 0;

	begin = clock();
	//FILE *f = fopen("C:\\image_GPU.txt", "a+");
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

		//if(_input->data[j]>(unsigned short)_thresholdVal){
		//	_outputPub->data[j]=(float)_input->data[j];
		//	_segPub->data[j]=1;
		//}else{
		//	_segPub->data[j]=0;
		//	_outputPub->data[j]=0;
		//}
	}
	//fclose(f);
	end = clock();

	/**********************/
	//float* d_wav = (float*)malloc(_sizeX*_sizeY * sizeof(float));
	//FILE* 
	/**********************/

	///*FILE* */fp = fopen("C:\\Throsholding.txt", "a+");
	//fprintf(fp, "Time: %f\n", (end - begin) / (float)CLOCKS_PER_SEC );
	//fclose(fp)	;

	begin = clock();

	if (_watershedRatio == 1.) {
		if ((nbStruct = _segmentPub->segmentation(*_outputPub, *_segPub, 1e-8f, true, false)) < 0) return (nbStruct);
		//fprintf(f, "segmentation\n");
	}
	else {
		if ((nbStruct = _segmentPub->watershedIntelligent(*_outputPub, *_segPub, 1e-8f, (float)_watershedRatio)) < 0) return (nbStruct);
		//fprintf(f, "watershedIntelligent\n");
	}

	end = clock();

	//fp = fopen("C:\\Watershed.txt", "a+");
	//fprintf(fp, "Time: %f\n", (end - begin) / (float)CLOCKS_PER_SEC );
	//fclose(fp)	;

	/***********/
	//_image = _segPub->data;
	//FILE* file = fopen("C:\\seg.txt", "w");

	//int cos = 1;
	//float max = -99999999;
	//float min = +99999999;
	//for(int i=0; i<_sizeX*_sizeY;i++)
	//{
	//	if (_aTrousPub->wDecomp->data[cos*_waveletNo*_sizeX*_sizeY+i]>max)
	//		max = _aTrousPub->wDecomp->data[cos*_waveletNo*_sizeX*_sizeY+i];

	//	if (_aTrousPub->wDecomp->data[cos*_waveletNo*_sizeX*_sizeY+i]<min)
	//		min = _aTrousPub->wDecomp->data[cos*_waveletNo*_sizeX*_sizeY+i];
	//}

	//fprintf(file, "max: %f\n", max);
	//fprintf(file, "min: %f\n", min);
	//for(int i=0; i<_sizeX*_sizeY;i++)
	//{
	//fprintf(file, "%d\n", _segPub->data[i]);
	//fprintf(file, "%d\n", (unsigned short)(((_aTrousPub->wDecomp->data[(cos*_waveletNo*_sizeX*_sizeY)+i] -min) * 65536.)/(max-min)));
	//fprintf(file, "%f\n", _aTrousPub->wDecomp->data[(_waveletNo*_sizeX*_sizeY)+i]);

	//_aTrousPub->wDecomp->data[(cos*_waveletNo*_sizeX*_sizeY)+i] = _aTrousPub->wDecomp->data[(cos*_waveletNo*_sizeX*_sizeY)+i] -min + 1;
	//}
	//fclose(file);
	/***********/
	;
	begin = clock();


	/*long offs=_waveletNo*_sizeX*_sizeY;*/
	if ((nbPointsSign = _segmentPub->computeStatistics(*_outputPub, *_segPub, _aTrousPub->wDecomp->data + _waveletNo * _sizeX * _sizeY, _gaussFit)) < 0) return (nbPointsSign);

	//end = clock();
	//fp = fopen("C:\\computeStatistics.txt", "a+");
	//fprintf(fp, "Time: %f\n", (end - begin) / (float)CLOCKS_PER_SEC );
	//fclose(fp);

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

	end = clock();
	//fp = fopen("C:\\filterData.txt", "a+");
	//fprintf(fp, "Time: %f\n", (end - begin) / (float)CLOCKS_PER_SEC );
	//fclose(fp);

	//FILE* f = fopen("C:\\Nbrstruct.txt", "a+");
	//fprintf(f, "%d\n", nbStruct);
	//fclose(f);

	///*FILE**/ f = fopen("C:\\Trace_C++.txt", "a+");
	InfoStruct2D** p;
	p = (InfoStruct2D**)_segmentPub->getTabInfo();
	//fprintf(f, "getTabInfo\n");
	//fclose(f);

	//return 0;

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

	//fprintf(f, "End\n");
	//fclose(f);

	//FILE* fp= fopen("C:\\test.txt", "a+");
	//fprintf(fp, "%d\n", currentPoint);
	//fclose(fp);

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
	//FILE* f = fopen("C:\\debug.txt", "w");
	//fprintf(f, "Before CPALM\n");
	//fclose(f);

	sum = 0.;
	pointsNumber = 0;

	//unsigned short* _image = (unsigned short*) malloc(sizeX * sizeY * sizeof(unsigned short));
	//for(int i=0; i< sizeX* sizeY; i++)
	//	_image[i] = (unsigned short) ((int*)image)[i];

	c = new CPALM(image, pointsList, pointsListSize, sizeX, sizeY, waveletNo, thresholdVal, watershedRatio, volMin, intMin, gaussFit, sigmaX, sigmaY, angle, size);

	//f = fopen("C:\\debug.txt", "w");
	//fprintf(f, "After CPALM\n");
	//fclose(f);
}

_declspec(dllexport) long _stdcall __PALMProcessing()
{
	//FILE* f = fopen("C:\\LivePALM.txt", "a+");
	//for(int i = 0; i < c->_sizeX; i++)
	//	for(int j = 0; j < c->_sizeY; j++)
	//		fprintf(f, " %d\n", c->_image[i + j * c->_sizeX]);
	//fclose(f);

	//
	//FILE* f = fopen("C:\\debug.txt", "w");
	//fprintf(f, "Before CPALMProcessing\n");
	//fclose(f);

	long nbr = c->CPALMProcessing();

	//FILE* f = fopen("C:\\Points_C++.txt", "a+");
	//for(int i = 0; i < ((nbr / 10) - 1); i++)
	//	fprintf(f, "X: %f, Y: %f, Surface: %d\n", c->_pointsList[i * 10 + 4], c->_pointsList[i * 10 + 3], (int)c->_pointsList[i * 10 + 9]);
	//fprintf(f, "********************************************************\n");
	//fclose(f);

	//FILE* f = fopen("C:\\debug.txt", "a+");
	//fprintf(f, "hhhhhhhhh\n");
	//fclose(f);

	return (nbr);
	return 0;
}

_declspec(dllexport) void _stdcall __closePALMProcessing()
{
	//printf("Sum: %f\n", sum);	
	//printf("PointsNumber: %d\n", pointsNumber);	
	//printf("RSS: %f\n", sum/pointsNumber);

	/*********************/
	//FILE* fp = fopen("c:\\trace.txt", "w");
	//   for(int i = 0; i < 10000; i++) 
	//{
	//	fprintf(fp, "%d\n", tab[i]);
	//}
	//fclose(fp);
	/*********************/
	delete c;
}


unsigned short* stack;

_declspec(dllexport) void _stdcall AllocateMemory(unsigned short frameCount, unsigned short width, unsigned short height)
{
	//stack = (unsigned short*)malloc(frameCount * width * height * sizeof(unsigned short));
	stack = new unsigned short[frameCount * width * height];

	//FILE* f = fopen("C:\\debug.txt", "a+");
	//fprintf(f, "stack: %d\n", stack);
	//fclose(f);
	//checkCudaErrors( cudaMalloc( (void**) &d_stack, frameCount * width * height * sizeof(unsigned short) ) );
}

_declspec(dllexport) void _stdcall CopyFrame(unsigned short numImage, unsigned short* image, unsigned short frameCount, unsigned short width, unsigned short height)
{
	memcpy(stack + numImage * width * height, image, width * height * sizeof(unsigned short));

	//FILE* f = fopen("C:\\debug.txt", "a+");
	//fprintf(f, "%d\n", numImage);
	//fclose(f);
	//checkCudaErrors( cudaMemcpy( d_stack + numImage * width * height, image, width * height * sizeof(unsigned short), cudaMemcpyHostToDevice) );
}

_declspec(dllexport) void _stdcall InitializeMemory(unsigned short* frames) { stack = frames; }


_declspec(dllexport) void _stdcall LaunchGaussianFit(double* pointsList, unsigned long pointsListSize, unsigned short sizeStk, unsigned short sizeDataX, unsigned short sizeDataY, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size, unsigned short iterNum, unsigned short fromcplusplus, unsigned short computeZ, double* SigmaXCoefs, double* SigmaYCoefs, double ThreeDThickness/*µm*/, double ThreeDSampling/*µm*/, double pixelSize /*µm*/)
{
	if ((gaussFit >= 1) && (gaussFit <= 4)) {
		// debug
		//FILE *fp = fopen("D:\\Image.txt", "w");
		//		for(int i = 0; i < 15 * sizeDataX * sizeDataY; i++)
		//{
		//	fprintf(fp, "%d,", stack[i]);
		//}
		//fclose(fp);

		//fp = fopen("D:\\pointsList.txt", "w");
		//for(int i = 0; i < pointsListSize; i++)
		//{
		//	fprintf(fp, "%f,", pointsList[i]);
		//}
		//fclose(fp);

		//fp = fopen("D:\\Params.txt", "w");
		//fprintf(fp, "sizeStk: %d, sizeDataX: %d, sizeDataY: %d", sizeStk, sizeDataX, sizeDataY);
		//fclose(fp);

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
		if (kernelExecTimeoutEnabled == 1) newSizeStk = 80;

		//newSizeStk = 80;

		//if (sizeStk > newSizeStk)
		{
			iter = sizeStk / newSizeStk;
			r = sizeStk % newSizeStk;
		}

		//FILE* f = fopen("debug.txt", "w");
		//fprintf(f, "LaunchGaussianFit\n");

		//fprintf(f, "iter: %d\n", iter);
		//fprintf(f, "r: %d\n", r);
		//fprintf(f, "newSizeStk: %d\n", newSizeStk);
		//fprintf(f, "sizeStk: %d", sizeStk);

		//fclose(f);

		//FILE* f = fopen("C:\\log.txt", "w");
		//fprintf(f, "iter: %d\n", iter);
		//fprintf(f, "r: %d\n", r);
		//fprintf(f, "newSizeStk: %d\n", newSizeStk);
		//fprintf(f, "sizeStk: %d\n", sizeStk);
		//fclose(f);

		//
		//printf("LaunchGaussianFit\n");

		//FILE* f = fopen("C:\\debug.txt", "w");
		//fprintf(f, "%d, %d, %d, %d, %f, %f, %f, %d\n", sizeStk, sizeDataX, sizeDataY, gaussFit, sigmaX, sigmaY, angle, size);
		//fclose(f);

		// Fill nbStruct
		unsigned long NbSegParams = 13;

		unsigned int* cumNbStruct = (unsigned int*)malloc((sizeStk + 1) * sizeof(unsigned int));
		unsigned short* nbStruct = (unsigned short*)malloc(sizeStk * sizeof(unsigned short));

		//f = fopen("C:\\debug.txt", "a+");
		//fprintf(f, "520\n");
		//fclose(f);

		unsigned short numFrame = 0;
		nbStruct[numFrame] = 0;
		cumNbStruct[numFrame] = 0;

		unsigned int nbPoints = 0;

		//f = fopen("C:\\debug.txt", "a+");
		//fprintf(f, "529\n");
		//fclose(f);

		for (int i = 0; i < (pointsListSize / NbSegParams); i++) {
			if ((int)pointsList[NbSegParams * i] == -1) {
				nbPoints += nbStruct[numFrame];

				numFrame++;
				if (numFrame == sizeStk) break;
				nbStruct[numFrame] = 0;

				//printf("numFrame: %d\n", numFrame);
			}
			else //if((int)pointsList[NbSegParams * i] != -1)
				nbStruct[numFrame]++;
		}

		// Checking the available maximum sizes of each dimension of a grid (GPU launching)
		if ((int)ceil((float)(nbPoints) / (float)FIT_BLOCK_SIZE) > deviceProp.maxGridSize[0]) {
			newSizeStk = (int)(ceil((float)(deviceProp.maxGridSize[0] * FIT_BLOCK_SIZE) / (float)(nbPoints / sizeStk)));
			newSizeStk -= 10;
			iter = sizeStk / newSizeStk;
			r = sizeStk % newSizeStk;
		}

		float* list = (float*)malloc(nbPoints * NbSegParams * sizeof(float));

		//for(int i = 0; i < nbPoints * NbSegParams; i++)
		//	list[i] = 0.;

		int ind_src;
		int ind_dst;

		if (fromcplusplus == 1)
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

				if ((plane + 1) <= sizeStk) cumNbStruct[plane + 1] = cumNbStruct[plane] + nbStruct[plane];
			}
		else
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

				if ((plane + 1) <= sizeStk) cumNbStruct[plane + 1] = cumNbStruct[plane] + nbStruct[plane];
			}

		// debug
		//FILE *fp = fopen("C:\\nbStruct.txt", "w");
		//for(int i = 0; i < sizeStk; i++)
		//{
		//	fprintf(fp, "%d: %d: %d\n", i, nbStruct[i], cumNbStruct[i]);
		//}
		//fclose(fp);

		///*FILE**/ f = fopen("C:\\debug.txt", "a+");
		//for(int i = 0; i < nbPoints; i++)
		//{
		//	fprintf(f, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n", list[NbSegParams * i + 0],  list[NbSegParams * i + 1], list[NbSegParams * i + 2], list[NbSegParams * i + 3], list[NbSegParams * i + 4], list[NbSegParams * i + 5], list[NbSegParams * i + 6], list[NbSegParams * i + 7], list[NbSegParams * i + 8], list[NbSegParams * i + 9], list[NbSegParams * i + 10]);
		//}
		//fclose(f);


		//FILE* f = fopen("C:\\MM\\Processing DLLs\\FittingPostProcessing.txt", "w+");
		//fprintf(f, "%d\n", gaussFit);

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
			//float* d_data;
			unsigned short *d_images, *d_nbStruct;
			float* d_list;
			//checkCudaErrors( cudaMalloc( (void**) &d_data, newSizeStk * sizeDataX * sizeDataY * sizeof(float) ) );
			checkCudaErrors(cudaMalloc((void**)&d_images, newSizeStk * sizeDataX * sizeDataY * sizeof(unsigned short)));
			checkCudaErrors(cudaMalloc((void**)&d_nbStruct, newSizeStk * sizeof(unsigned short)));
			checkCudaErrors(cudaMalloc((void**)&d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]) * NbSegParams * sizeof(float)));

			checkCudaErrors(cudaMemcpy(d_images, stack + (i - 1) * newSizeStk * sizeDataX * sizeDataY, newSizeStk * sizeDataX * sizeDataY * sizeof(unsigned short), cudaMemcpyHostToDevice));
			checkCudaErrors(cudaMemcpy(d_nbStruct, nbStruct + (i - 1) * newSizeStk, newSizeStk * sizeof(unsigned short), cudaMemcpyHostToDevice));
			checkCudaErrors(cudaMemcpy(d_list, list + cumNbStruct[(i - 1) * newSizeStk] * NbSegParams, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]) * NbSegParams * sizeof(float), cudaMemcpyHostToDevice));

			//convertUSrtToFloat(d_images, d_data, newSizeStk, sizeDataX, sizeDataY);	

			//FILE* fp = fopen("C:\\log.txt", "a+");
			//fprintf(fp, "i: %d\n", i);
			//fclose(fp);

			// cleanup memory
			//checkCudaErrors( cudaFree(d_images) );

			if (gaussFit == 1) // 2D with constant Sigma
			{
				if (size < 7) cuda2DFittingGPU_5x5(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*30*/, size, sigmaX);

				if (size >= 7) cuda2DFittingGPU_7x7(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*30*/, size, sigmaX);
			}

			if (gaussFit == 2) // 2D with Sigma
			{
				if (size < 7) cuda2DFittingGPU_Sigma_5x5(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*30*/, size, sigmaX);

				if (size >= 7) cuda2DFittingGPU_Sigma_7x7(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*30*/, size, sigmaX);
			}

			if (gaussFit == 3) // 3D with constant Theta
			{
				if ((size <= 7) || (size == 8)) cuda3DFittingGPU_7x7(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/);

				if ((size == 9) || (size == 10)) cuda3DFittingGPU_9x9(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/);

				if ((size == 11) || (size == 12)) cuda3DFittingGPU_11x11(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/);

				if (size >= 13) cuda3DFittingGPU_13x13(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/);
			}
			if (gaussFit == 4) // 3D with Theta
			{
				if (size <= 9) cuda3DFittingGPU_Theta_9x9(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum /*20*7*/, size, sigmaX, sigmaY, angle);
				if (size > 9) cuda3DFittingGPU_Theta_13x13(d_images, newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]), iterNum /*20*7*/, size, sigmaX, sigmaY, angle);
			}

			cudaError err = cudaThreadSynchronize();
			//printf("Error : %d \n", err);
			checkCudaErrors(cudaMemcpy(list + cumNbStruct[(i - 1) * newSizeStk] * NbSegParams, d_list, (cumNbStruct[i * newSizeStk] - cumNbStruct[(i - 1) * newSizeStk]) * NbSegParams * sizeof(float), cudaMemcpyDeviceToHost));

			// Cleanup memory
			checkCudaErrors(cudaFree(d_images));
			checkCudaErrors(cudaFree(d_list));
			checkCudaErrors(cudaFree(d_nbStruct));

			//fprintf(f, "%d\n", i*newSizeStk);

			//fp = fopen("C:\\log.txt", "a+");
			//fprintf(fp, "i: %d\n", i);
			//fclose(fp);
		}

		if (r != 0) {
			//newSizeStk = sizeStk % newSizeStk;

			//float* d_data;
			unsigned short *d_images, *d_nbStruct;
			float* d_list;
			//checkCudaErrors( cudaMalloc( (void**) &d_data, (sizeStk % newSizeStk) * sizeDataX * sizeDataY * sizeof(float) ) );
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
				if (size < 7) cuda2DFittingGPU_5x5(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*30*/, size, sigmaX);

				if (size >= 7) cuda2DFittingGPU_7x7(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*30*/, size, sigmaX);
			}

			if (gaussFit == 2) // 2D with Sigma
			{
				if (size < 7) cuda2DFittingGPU_Sigma_5x5(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*30*/, size, sigmaX);

				if (size >= 7) cuda2DFittingGPU_Sigma_7x7(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*30*/, size, sigmaX);
			}

			if (gaussFit == 3) // 3D with constant Theta
			{
				if ((size <= 7) || (size == 8)) cuda3DFittingGPU_7x7(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/);

				if ((size == 9) || (size == 10)) cuda3DFittingGPU_9x9(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/);

				if ((size == 11) || (size == 12)) cuda3DFittingGPU_11x11(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/);

				if (size >= 13) cuda3DFittingGPU_13x13(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*20*//*50*//*20*6*/, size, sigmaX, sigmaY, angle, computeZ, d_SigmaXCoefs, d_SigmaYCoefs, ThreeDThickness/*µm*/, ThreeDSampling/*µm*/, pixelSize/*µm*/);
			}
			if (gaussFit == 4) // 3D with Theta
			{
				if (size <= 9) cuda3DFittingGPU_Theta_9x9(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*20*7*/, size, sigmaX, sigmaY, angle);
				if (size > 9) cuda3DFittingGPU_Theta_13x13(d_images, sizeStk % newSizeStk, sizeDataY, sizeDataX, d_nbStruct, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]), iterNum/*20*7*/, size, sigmaX, sigmaY, angle);
			}

			cudaError err = cudaThreadSynchronize();
			//printf("Error : %d \n", err);
			checkCudaErrors(cudaMemcpy(list + cumNbStruct[iter * newSizeStk] * NbSegParams, d_list, (nbPoints - cumNbStruct[iter * newSizeStk]) * NbSegParams * sizeof(float), cudaMemcpyDeviceToHost));

			// Cleanup memory
			checkCudaErrors(cudaFree(d_images));
			checkCudaErrors(cudaFree(d_list));
			checkCudaErrors(cudaFree(d_nbStruct));

			//fprintf(f, "%d\n", sizeStk);
		};
		//fclose(f);

		// Cleanup memory
		checkCudaErrors(cudaFree(d_SigmaXCoefs));
		checkCudaErrors(cudaFree(d_SigmaYCoefs));

		// debug
		//FILE *f = fopen("C:\\debug.txt", "w");
		//for(int i = 0; i < sizeStk * MAX_NB_STRUCTURES; i++)
		//{
		//	fprintf(f, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n", list[10 * i + 0],  list[10 * i + 1], list[10 * i + 2], list[10 * i + 3], list[10 * i + 4], list[10 * i + 5], list[10 * i + 6], list[10 * i + 7], list[10 * i + 8], list[10 * i + 9]);
		//}
		//fclose(f);

		//FILE* f = fopen("C:\\debug.txt", "w");
		//fprintf(f , "nbPoints: %d\n", nbPoints);
		//for(int i = 0; i < nbPoints; i++)
		//{
		//	fprintf(f, "%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f \n", list[NbSegParams * i + 0],  list[NbSegParams * i + 1], list[NbSegParams * i + 2], list[NbSegParams * i + 3], list[NbSegParams * i + 4], list[NbSegParams * i + 5], list[NbSegParams * i + 6], list[NbSegParams * i + 7], list[NbSegParams * i + 8], list[NbSegParams * i + 9], list[NbSegParams * i + 10]);
		//}
		//fclose(f);

		int ind_list = 0;
		int ind_pointsList = 0;

		if (fromcplusplus != 1)
			for (int i = 0; i < sizeStk /*host*/; i++) {
				//ind_list = (i * maxNbStruct/* + 1*/) * NbSegParams;
				//while((int)list[ind_list + 2] != 0)
				for (int s = 0; s < nbStruct[i]; s++) {
					//pointsList[ind_pointsList++] = list[ind_list++];  // x
					//pointsList[ind_pointsList++] = list[ind_list++];  // y

					pointsList[ind_pointsList++] = (double)list[ind_list++]; // sigmaX
					pointsList[ind_pointsList++] = (double)list[ind_list++]; // sigmaY

					pointsList[ind_pointsList++] = (double)list[ind_list++]; // Theta


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
		else
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

		// cleanup memory
		free(cumNbStruct);
		free(nbStruct);
		free(list);
		//free(stack);

		if (fromcplusplus != 1) delete[] stack;

		/********************************************************/
		//fp = fopen("D:\\Sigma.txt", "w");
		//	for(int i = 0; i < pointsListSize/NbSegParams; i++)
		//	{
		//		fprintf(fp, "%f,%f\n", pointsList[i*NbSegParams], pointsList[i*NbSegParams+1]);
		//	}
		//	fclose(fp);
		/********************************************************/

		//FILE* f=fopen("C:\\end.txt", "a+");
		//fprintf(f, "End of Gaussian fitting\n");
		//fclose(f);	
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

	//FILE* f = fopen("c:\\GaussianVisuGPU.txt", "w");
	//fprintf(f, "%d\n", zoom);
	//fclose(f);

	//FILE* f = fopen("C:\\log.txt", "w");
	//fprintf(f, "Hello !!!!!!!!!!!!!!\n%d\n", pointsListSize);
	//fclose(f);

	//FILE* f= fopen("C:\\pointsList.txt", "a+");
	//fprintf(f, "%d\n", pointsListSize);
	//for(int i = 0; i < (pointsListSize/NbSegParams); i++)
	//{
	//	for(int j = 0; j<NbSegParams; j++)
	//	{
	//		fprintf(f, "%f\t", pointsList[NbSegParams * i + j]);
	//	}
	//	fprintf(f, "\n");
	//}
	//fclose(f);

	for (int i = 0; i < (pointsListSize / NbSegParams); i++) {
		if ((int)pointsList[NbSegParams * i] == -1) {
			nbPoints += nbStruct[numFrame];

			numFrame++;
			if (numFrame == sizeStk) break;
			nbStruct[numFrame] = 0;

			//printf("numFrame: %d\n", numFrame);
		}
		else //if((int)pointsList[NbSegParams * i] != -1)
			nbStruct[numFrame]++;

		//fprintf(f, "%d\n", i);
	}

	/*FILE* f = fopen("C:\\log.txt", "w");
	fprintf(f, "%d\n", nbPoints);
	fclose(f);*/

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

		if ((plane + 1) <= sizeStk) cumNbStruct[plane + 1] = cumNbStruct[plane] + nbStruct[plane];
	}

	//FILE* f= fopen("C:\\log.txt", "a+");
	//for (int i=0; i<nbPoints;i++)
	//	fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", list[i * NbSegParams + 0], list[i * NbSegParams + 1], list[i * NbSegParams + 2], list[i * NbSegParams + 3], list[i * NbSegParams + 4], list[i * NbSegParams + 5], list[i * NbSegParams + 6], list[i * NbSegParams + 7], list[i * NbSegParams + 8], list[i * NbSegParams + 9], list[i * NbSegParams + 10]);
	//fclose(f);

	double *d_gaussianPointsList, *d_pointsList;

	//for(int i = 0; i<=nbPoints * sizeReg * sizeReg; i++)
	//	gaussianPointsList[i]=321;

	//cudaDeviceProp deviceProp;
	//cudaGetDeviceProperties(&deviceProp, 0);
	//FILE* fp = fopen("c:\\GraphicCards.txt", "w");
	//fprintf(fp, "%s\n", deviceProp.name);
	//fclose(fp);

	checkCudaErrors(cudaMalloc((void**)&d_gaussianPointsList, sizeReg * sizeReg * nbPoints * sizeof(double)));
	checkCudaErrors(cudaMalloc((void**)&d_pointsList, nbPoints * NbSegParams * sizeof(double)));

	checkCudaErrors(cudaMemcpy(d_gaussianPointsList, gaussianPointsList, sizeReg * sizeReg * nbPoints * sizeof(double), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(d_pointsList, list, nbPoints * NbSegParams * sizeof(double), cudaMemcpyHostToDevice));

	//******************************************
	//double *list_sup = (double*)malloc(nbPoints * NbSegParams * sizeof(double));
	//checkCudaErrors( cudaMemcpy( list_sup, d_pointsList, nbPoints * NbSegParams * sizeof(double), cudaMemcpyDeviceToHost) );
	//	FILE* f= fopen("C:\\log.txt", "a+");
	//for (int i=0; i<nbPoints;i++)
	//	fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", list_sup[i * NbSegParams + 0], list_sup[i * NbSegParams + 1], list_sup[i * NbSegParams + 2], list_sup[i * NbSegParams + 3], list_sup[i * NbSegParams + 4], list_sup[i * NbSegParams + 5], list_sup[i * NbSegParams + 6], list_sup[i * NbSegParams + 7], list_sup[i * NbSegParams + 8], list_sup[i * NbSegParams + 9], list_sup[i * NbSegParams + 10]);
	//fclose(f);
	//******************************************

	GaussianVisu(d_gaussianPointsList /*(sizeReg*sizeReg*nbPoints)*/, d_pointsList, width, height, nbPoints, sizeReg, zoom, conversionRatio, pixelSize);

	checkCudaErrors(cudaMemcpy(gaussianPointsList, d_gaussianPointsList, sizeReg * sizeReg * nbPoints * sizeof(double), cudaMemcpyDeviceToHost));

	checkCudaErrors(cudaFree(d_gaussianPointsList));
	checkCudaErrors(cudaFree(d_pointsList));

	//FILE* f = fopen("C:\\GPU_log.txt", "a+");
	//for(int i = 0; i<=nbPoints * sizeReg * sizeReg; i++)
	//	fprintf(f, "%f\n", gaussianPointsList[i]);
	//fclose(f);

	//FILE* f= fopen("C:\\pointsList.txt", "w");
	//fprintf(f, "%d\n", pointsListSize);
	//for(int i = 0; i < (pointsListSize/NbSegParams); i++)
	//{
	//	for(int j = 0; j<NbSegParams; j++)
	//	{
	//		fprintf(f, "%f\t", pointsList[NbSegParams * i + j]);
	//	}
	//	fprintf(f, "\n");
	//}
	//fclose(f);
}
