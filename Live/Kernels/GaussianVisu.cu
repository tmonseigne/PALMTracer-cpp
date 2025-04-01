#pragma once

#include <helper_cuda.h>
#include <helper_timer.h>
#include <math_constants.h>
#include "../parameters.h"

#define NbSegParams 13

__device__ float GaussianFunction_2D(float x, float y, float I0, float sigma) { return I0 * expf((-powf(x, 2) - powf(y, 2)) / (2 * powf(sigma, 2))); }

__global__ void GaussianVisu_kernel(double* gaussianPointsList, double* pointsList, unsigned short width, unsigned short height, unsigned int nbPoints, unsigned short sizeReg, unsigned short zoom, double conversionRatio /*from: ADU to photons: photons/ADU */, double pixelSize)
{
	unsigned int index = threadIdx.x + blockIdx.x * blockDim.x;

	if ((index < nbPoints) && ((int)pointsList[NbSegParams * index + 1] != -1)) {
		float x0, y0;
		//int sizeRegH = (int)(sizeReg/2);
		int sizeRegH = sizeReg / 2;

		x0 = (float)pointsList[NbSegParams * index + 0];
		y0 = (float)pointsList[NbSegParams * index + 1];

		x0 = x0 * (float)zoom + 0.5 * ((float)zoom - 1.0);
		y0 = y0 * (float)zoom + 0.5 * ((float)zoom - 1.0);

		//x0 = x0 * pixelSize;
		//y0 = y0 * pixelSize;

		float photonsNbr = 2 * 3.14 * (float)pointsList[NbSegParams * index + 10] * (float)pointsList[NbSegParams * index + 9] * (float)pointsList[NbSegParams * index + 9];
		photonsNbr = photonsNbr * conversionRatio;
		float sigma = ((float)pointsList[NbSegParams * index + 9] * pixelSize) / sqrt(photonsNbr); // en nm
		sigma = sigma / (pixelSize / zoom);                                                        // en pixel SR space
		//int sizeRegH = (int)(3*sigma/2);
		float A = 10000 / (2 * 3.14 * sigma * sigma); //photonsNbr / (2 * 3.14 * sigma * sigma);10000 / (2 * 3.14 * sigma * sigma)

		/*if ((int)pointsList[index * NbSegParams + 3] != 0)*/
		if (((unsigned int)y0 >= sizeRegH) && ((unsigned int)y0 <= (height * zoom) - 1 - sizeRegH) && ((unsigned int)x0 >= sizeRegH) && ((unsigned int)x0 <= (width * zoom) - 1 - sizeRegH)) {
			for (int j = -sizeRegH; j <= sizeRegH; j++) {
				for (int i = -sizeRegH; i <= sizeRegH; i++) {
					gaussianPointsList[index * sizeReg * sizeReg + (i + sizeRegH) + (j + sizeRegH) * sizeReg] = GaussianFunction_2D(x0 - (int)(x0 - i), y0 - (int)(y0 - j), A, sigma); //1000.0f/*I0*/, 1.0f/*sigma*/);
				}
			}
		}
	}
}


extern "C"void GaussianVisu(double* gaussianPointsList /*(sizeReg*sizeReg*nbPoints)*/, double* pointsList, unsigned short width, unsigned short height, unsigned int nbPoints, unsigned short sizeReg, unsigned short zoom, double conversionRatio /*from: ADU to photons: photons/ADU*/, double pixelSize)
{
	dim3 grid(1, 1, 1);
	dim3 threads(1, 1, 1);

	// allocate memory
	grid.x = (int)ceil((float)(nbPoints) / (float)FIT_BLOCK_SIZE);
	grid.y = 1;
	grid.z = 1;
	threads.x = FIT_BLOCK_SIZE;
	threads.y = 1;
	threads.z = 1; // Par rapport aux contraintes de la mémoire partagée (taille limitée à 16 Ko pour chaque block)
	// '64' au lieu de '256'
	GaussianVisu_kernel<<< grid, threads >>>(gaussianPointsList, pointsList, width, height, nbPoints, sizeReg, zoom, conversionRatio /*from: ADU to photons*/, pixelSize);

	cudaError err = cudaThreadSynchronize();
}

