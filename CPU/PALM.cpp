// PALM.cpp : définit le point d'entrée pour l'application DLL.
//


#include "PALM.h"
#include "palm/CPALM.h"

#include "fit/CalibrationFit.h"

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

_declspec(dllexport) void _stdcall _OpenPALMProcessing(unsigned short* image, double* pointsList, unsigned long pointsListSize, unsigned long sizeX, unsigned long sizeY, unsigned long waveletNo, double thresholdVal, double watershedRatio, double volMin, double intMin, unsigned short gaussFit, double sigmaX, double sigmaY, double angle, unsigned short size)
{
	/* DEBUG */
	//FILE* f = fopen("C:\\debug.txt", "w");
	//fprintf(f, "Hello !!!!!!!!!!!\n");	
	//fprintf(f, "pointsListSize: %d,\n sizeX: %d,\n sizeY: %d,\n waveletNo: %d,\n thresholdVal: %f,\n watershedRatio: %f,\n volMin: %f,\n intMin: %f,\n gaussFit: %d,\n sigmaX: %f,\n sigmaY: %f,\n angle: %f,\n size: %d \n", pointsListSize, sizeX, sizeY, waveletNo, thresholdVal, watershedRatio, volMin, intMin, gaussFit, sigmaX, sigmaY, angle, size);
	//		  //maxNbPoints,     sizeX,     sizeY,     waveletNo,     thresholdVal,     watershedRatio,     volMin,     intMin,     threeD,     fit,     sigmaX,     sigmaY,     angle,     size);
	//
	//if (image  == NULL)
	//	fprintf(f, "%s\n", image);
	//else
	//	fprintf(f, "image: %d\n", image);
	//
	//if (pointsList == NULL)
	//	fprintf(f, "%s\n", pointsList);
	//else
	//	fprintf(f, "pointsList: %d\n", pointsList);

	//fclose(f);

	//maxNbPoints = /*250000*/500000;
	//sizeX = 128;
	//sizeY = 128;
	//waveletNo = 1;
	//thresholdVal = 100.;
	//watershedRatio = 0.;
	//volMin = 4.;
	//intMin = 0.;
	//threeD = 0;
	//fit = 0;
	//sigmaX = 1.;
	//sigmaY = 1.;
	//angle = 0.;
	//size = 7;

	/*************************************************************/
	//glutInit(NULL, NULL); /* initialise la bibliothèque GLUT */
	//glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	//glutInitWindowSize(500, 500);
	//glutInitWindowPosition(100, 100);
	//glutMainLoop();
	/*************************************************************/

	//FILE* f = fopen("C:\\debug.txt", "w");
	//fprintf(f, "Before CPALM\n");
	//fclose(f);

	sum = 0.;
	pointsNumber = 0;

	c = new CPALM(image, pointsList, pointsListSize, sizeX, sizeY, waveletNo, thresholdVal, watershedRatio, volMin, intMin, gaussFit, sigmaX, sigmaY, angle, size);

	//f = fopen("C:\\debug.txt", "w");
	//fprintf(f, "After CPALM\n");
	//fclose(f);
}

_declspec(dllexport) long _stdcall _PALMProcessing()
{
	//FILE* f = fopen("C:\\CPUPALM.txt", "a+");
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

	//FILE* f = fopen("C:\\Nbre.txt", "a+");
	//fprintf(f, "%d\n", nbr);
	//fclose(f);

	return (nbr);
	return 0;
}

_declspec(dllexport) void _stdcall _closePALMProcessing()
{
	//printf("Sum: %f\n", sum);	
	//printf("PointsNumber: %d\n", pointsNumber);	
	//printf("RSS: %f\n", sum/pointsNumber);

	delete c;
}

_declspec(dllexport) double _stdcall RunCalibrationFit(double* sigma, double* Z, unsigned short pointsNbr, double* C)
{
	//FILE* f= fopen("C:\\TraceC++.txt", "w");
	//fprintf(f, "%f, %f, %f, %f, %f\n", C[0], C[1], C[2], C[3], C[4]);
	float* _sigma = (float*)malloc(pointsNbr * sizeof(float));
	float* _Z = (float*)malloc(pointsNbr * sizeof(float));

	for (int i = 0; i < pointsNbr; i++) {
		_sigma[i] = (float)sigma[i];
		_Z[i] = (float)Z[i];
		//fprintf(f, "sigma: %f, Z: %f\n", _sigma[i], _Z[i]);
	}

	//fclose(f);

	float C0 = (float)C[0];
	float C1 = (float)C[1];
	float C2 = (float)C[2];
	float C3 = (float)C[3];
	float C4 = (float)C[4];

	CalibrationFit* fit = new CalibrationFit();
	float RSS = fit->Run(_sigma, _Z, pointsNbr, 150, C0, C1, C2, C3, C4);

	C[0] = (double)C0;
	C[1] = (double)C1;
	C[2] = (double)C2;
	C[3] = (double)C3;
	C[4] = (double)C4;


	//FILE* f = fopen("C:\\Results.txt", "a+");
	//for(int i=0; i < pointsNbr; i++)
	{
		//fprintf(f, "%f %f\n", _Z[i], fit->sigmaFunc(_Z[i], C0, C1, C2, C3, C4));
		//fprintf(f, "%f %f\n", _Z[i], fit->sigmaFunc(_Z[i], -.120020702338513, .283238721128094, .248737266621778, 9.69480906847152E-02, .135361495198247));
	}

	//for(int i=0; i < pointsNbr; i++)
	//{
	//	fprintf(f, "%f %f\n", _Z[i], fit->sigmaFunc(_Z[i],  .120020702338513, .25011170282746, .170348876804454, 3.27843842401171E-02, .142517439595579));
	//}

	//fclose(f);

	// Cleanup memory
	free(_sigma);
	free(_Z);

	return RSS;
}

float GaussianFunction_2D(float x, float y, float I0, float sigma) { return I0 * expf((-powf(x, 2) - powf(y, 2)) / (2 * powf(sigma, 2))); }

_declspec(dllexport) void _stdcall GaussianVisuCPU(double* gaussianPointsList/*(sizeReg*sizeReg*nbPoints)*/, double* pointsList, unsigned long pointsListSize, unsigned short sizeStk, unsigned short width, unsigned short height, unsigned short sizeReg, unsigned short zoom, double conversionRatio, double pixelSize)
{
	unsigned int index = 0;
	float x0, y0;
	unsigned long NbSegParams = 13;

	//FILE* f = fopen("C:\\log.txt", "w");
	//fprintf(f, "%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f\n", pointsListSize, sizeStk, width, height, sizeReg, zoom, conversionRatio, pixelSize);
	//for(int i = 0; i < (pointsListSize/NbSegParams); i++)
	//{
	//	fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", pointsList[NbSegParams * i + 0], pointsList[NbSegParams * i + 1],pointsList[NbSegParams * i + 2],pointsList[NbSegParams * i + 3],pointsList[NbSegParams * i + 4],pointsList[NbSegParams * i + 5],pointsList[NbSegParams * i + 6],pointsList[NbSegParams * i + 7],pointsList[NbSegParams * i + 8],pointsList[NbSegParams * i + 9], pointsList[NbSegParams * i + 10]);
	//}
	//fclose(f);

	//FILE* f= fopen("C:\\log.txt", "a+");
	for (int i = 0; i < (pointsListSize / NbSegParams); i++) {
		//fprintf(f, "%d\n", (pointsListSize/NbSegParams));	
		if ((int)pointsList[NbSegParams * i + 9] > 0) {
			if ((int)pointsList[NbSegParams * i + 3] != -1) {
				int sizeRegH = sizeReg / 2;

				x0 = (float)pointsList[NbSegParams * i + 4];
				y0 = (float)pointsList[NbSegParams * i + 3];

				x0 = x0 * (float)zoom + 0.5 * ((float)zoom - 1.0);
				y0 = y0 * (float)zoom + 0.5 * ((float)zoom - 1.0);

				//x0 = x0 * pixelSize;
				//y0 = y0 * pixelSize;

				float photonsNbr = 2 * 3.14 * (float)pointsList[NbSegParams * i + 5] * (float)pointsList[NbSegParams * i + 0] * (float)pointsList[NbSegParams * i + 0];
				photonsNbr = photonsNbr * conversionRatio;
				float sigma = ((float)pointsList[NbSegParams * i + 0] * pixelSize) / sqrt(photonsNbr); // en nm
				sigma = sigma / (pixelSize / zoom);                                                    // in nm / (pixelSize / zoom); // en pixel SR space
				float A = 10000 / (2 * 3.14 * sigma * sigma);                                          //photonsNbr / (2 * 3.14 * sigma * sigma);

				/*if ((int)pointsList[index * NbSegParams + 3] != 0)*/
				if (((unsigned int)y0 >= sizeRegH) && ((unsigned int)y0 <= (height * zoom) - 1 - sizeRegH) && ((unsigned int)x0 >= sizeRegH) && ((unsigned int)x0 <= (width * zoom) - 1 - sizeRegH))
					for (int jj = -sizeRegH; jj <= sizeRegH; jj++)
						for (int ii = -sizeRegH; ii <= sizeRegH; ii++) gaussianPointsList[index * sizeReg * sizeReg + (ii + sizeRegH) + (jj + sizeRegH) * sizeReg] = GaussianFunction_2D(x0 - (int)(x0 - ii), y0 - (int)(y0 - jj), A, sigma); // 1000.0f, /*I0*/, 1.0f/*sigma*/);			
			}
			index++;
		}
	}
	//fprintf(f, "End");
	//fclose(f);

	/*FILE**/ /*f = fopen("C:\\CPU_log.txt", "a+");
	for(int i = 0; i<=index * sizeReg * sizeReg; i++)
		fprintf(f, "%f\n", gaussianPointsList[i]);
	fclose(f);*/

	//f = fopen("C:\\log.txt", "a+");
	//fprintf(f, "%d\n", nbPoints);
	//fclose(f);

	//FILE* f= fopen("C:\\log.txt", "w");
	//for (int i=0; i<nbPoints;i++)
	//	fprintf(f, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", list[i * NbSegParams + 0], list[i * NbSegParams + 1], list[i * NbSegParams + 2], list[i * NbSegParams + 3], list[i * NbSegParams + 4], list[i * NbSegParams + 5], list[i * NbSegParams + 6], list[i * NbSegParams + 7], list[i * NbSegParams + 8], list[i * NbSegParams + 9], list[i * NbSegParams + 10]);
	//fclose(f);
}
