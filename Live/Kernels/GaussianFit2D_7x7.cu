#ifndef _FITTING_KERNEL_H_
#define _FITTING_KERNEL_H_

#include <math_constants.h>
#include <helper_cuda.h>
#include <helper_timer.h>
#include "../parameters.h"

extern "C"unsigned int nextPow2(unsigned int x)
{
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

extern "C"bool isPow2(unsigned int x);

// Levenberg fit --> simple precision (float)
#define dimFit 7
#define Npar 4 // I0, x0, y0, offset
#define Npnt dimFit*dimFit
#define ref (dimFit/2)
#define NbSegParams 13

//__constant__ int MaxIter = /*20 * Npar*//*20*/40;
__constant__ float epsilon_1 = 1e-6;
__constant__ float epsilon_2 = 1e-6;
__constant__ float epsilon_3 = 1e-5;
__constant__ float epsilon_4 = 1e-2;
__constant__ float lambda_0 = 1e-2;
__constant__ float lambda_UP_fac = 11;
__constant__ float lambda_DN_fac = 9;

__device__ float PSF_7x7_2D_(float x, float y, float I0, float x0, float y0, float sigma, float offset) { return I0 * expf((-powf(x - x0, 2) - powf(y - y0, 2)) / (2 * powf(sigma, 2))) + offset; }

__device__ float deriv1I0_7x7_2D_(float x, float y, float I0, float x0, float y0, float sigma, float cons) { return cons / I0; }

__device__ float deriv1x0_7x7_2D_(float x, float y, float I0, float x0, float y0, float sigma, float cons) { return cons * ((x - x0) / powf(sigma, 2)); }

__device__ float deriv1y0_7x7_2D_(float x, float y, float I0, float x0, float y0, float sigma, float cons) { return cons * ((y - y0) / powf(sigma, 2)); }

__device__ float lm_matx_7x7_2D_(/*unsigned int* tIN,*/ float* p/*IN*/, float* y_dat/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/, float sigma)
{
	// Computing Jacobian matrix (dydp) & y_hat, delta_y
	float delta_y[Npnt];
	float dydp[Npnt * Npar];
	float cons;
	unsigned int x, y;
	float Chi_sq = 0.;

	for (int elt = 0; elt < Npnt; elt++) {
		x = (elt % dimFit);
		y = (elt / dimFit);
		cons = p[0] * expf((-powf((float)x - p[1], 2) - powf((float)y - p[2], 2)) / (2 * powf(sigma, 2)));
		dydp[(0) + elt * Npar] = deriv1I0_7x7_2D_((float)x, (float)y, p[0], p[1], p[2], sigma, cons); // I0
		dydp[(1) + elt * Npar] = deriv1x0_7x7_2D_((float)x, (float)y, p[0], p[1], p[2], sigma, cons); // x0
		dydp[(2) + elt * Npar] = deriv1y0_7x7_2D_((float)x, (float)y, p[0], p[1], p[2], sigma, cons); // y0
		dydp[(3) + elt * Npar] = 1;                                                                   // offset

		y_hat[elt] = cons + p[3];
		delta_y[elt] = y_dat[elt] - y_hat[elt];

		// Computing Chi_sq
		Chi_sq += delta_y[elt] * delta_y[elt];
	}

	// Computing alpha
	for (int j = 0; j < Npar; j++) {
		for (int i = 0; i < Npar; i++) {
			alpha[Npar * j + i] = 0;
			for (int cursor = 0; cursor < Npnt; cursor++) alpha[Npar * j + i] += dydp[cursor * Npar + j] * dydp[cursor * Npar + i]; // dydp' * inter : dydp'[i*Npnt + cursor]
		}

		// Computing beta
		beta[j] = 0;
		for (int cursor = 0; cursor < Npnt; cursor++) beta[j] += dydp[cursor * Npar + j] * delta_y[cursor]; // dydp' * inter : dydp'[j*Npnt + cursor]
	}

	return Chi_sq;
}

__device__ void inverse_7x7_2D_(float* MM, float* Minv)
{
	//int shift_ = (int)(threadIdx.x/16/*32*/) * (Npar * Npar * 16);
	//int id = threadIdx.x % 16/*32*/;

	float M1[Npar * Npar];

	// Copy MM to M1
	for (int i = 0; i < Npar * Npar; i++) M1[i] = MM[i];
	//M1[i] = MM[shift_ + (id + i * 16)];

	// Initializatiion of identity matrix
	float Mi1[Npar * Npar];
	for (int i = 0; i < Npar * Npar; i++) Mi1[i] = 0.;
	Mi1[0 * Npar + 0] = 1;
	Mi1[1 * Npar + 1] = 1;
	Mi1[2 * Npar + 2] = 1;
	Mi1[3 * Npar + 3] = 1;

	// Copy Mi1 to Minv
	for (int i = 0; i < Npar * Npar; i++) Minv[i] = Mi1[i];

	for (int i = 0; i < Npar; i++) {
		for (int j = 0; j < Npar; j++) {
			M1[i * Npar + j] = MM[i * Npar + j] / MM[i * Npar + i];
			Mi1[i * Npar + j] = Minv[i * Npar + j] / MM[i * Npar + i];
			//M1[i*Npar + j] = MM[shift_ + (id + (i*Npar + j) * 16)] / MM[shift_ + (id + (i*Npar + i) * 16)];
			//Mi1[i*Npar + j] = Minv[i*Npar + j] / MM[shift_ + (id + (i*Npar + i) * 16)];
		}

		// Copy M1 to MM & Copy Mi1 to Minv
		for (int index = 0; index < Npar * Npar; index++) {
			MM[index] = M1[index];
			//MM[shift_ + (id + index * 16)] = M1[index];
			Minv[index] = Mi1[index];
		}

		for (int k = 0; k < Npar; k++) // mise zéros de la colonne
			if (k != i)
				for (int j = 0; j < Npar; j++) {
					M1[k * Npar + j] = MM[k * Npar + j] - MM[i * Npar + j] * MM[k * Npar + i];
					Mi1[k * Npar + j] = Minv[k * Npar + j] - Minv[i * Npar + j] * MM[k * Npar + i];
					//M1[k*Npar+j] = MM[shift_ + (id + (k*Npar+j) * 16)] - MM[shift_ + (id + (i*Npar+j) * 16)] * MM[shift_ + (id + (k*Npar+i) * 16)];
					//Mi1[k*Npar+j] = Minv[k*Npar+j] - Minv[i*Npar+j] * MM[shift_ + (id + (k*Npar+i) * 16)];
				}

		// Copy M1 to MM & Copy Mi1 to Minv
		for (int index = 0; index < Npar * Npar; index++) {
			MM[index] = M1[index];
			//MM[shift_ + (id + index * 16)] = M1[index];
			Minv[index] = Mi1[index];
		}
	}
}

__global__ void GPUFitting_7x7_2D(unsigned int nbPoints, float* pointsList, unsigned short* data, unsigned short* nbStruct, unsigned short WIDTH, unsigned short HEIGHT, unsigned int iterNum, unsigned short _size, float sigmaX)
{
	unsigned int index = threadIdx.x + blockIdx.x * blockDim.x;

	if (index < nbPoints) {
		unsigned short idImage = (unsigned short)pointsList[NbSegParams * index + 10];
		unsigned int imageSize = WIDTH * HEIGHT;
		//unsigned short _seg = (unsigned short)pointsList[NbSegParams * index + 9];

		if ((int)pointsList[index * NbSegParams + 3] != 0) {
			if ((pointsList[NbSegParams * index + 1] >= 3) && (pointsList[NbSegParams * index + 1] <= HEIGHT - 4) && (pointsList[NbSegParams * index] >= 3) && (pointsList[NbSegParams * index] <= WIDTH - 4)) {
				//int shift = (int)(threadIdx.x/16/*32*/) * (Npar * 16);
				//int shift_ = (int)(threadIdx.x/16/*32*/) * (Npar * Npar * 16);
				//int id = threadIdx.x % 16/*32*/;

				int x = /*(int)*/floor(pointsList[NbSegParams * index]);
				int y = /*(int)*/floor(pointsList[NbSegParams * index + 1]);

				float _x = (pointsList[NbSegParams * index] - x) + ref;
				float _y = (pointsList[NbSegParams * index + 1] - y) + ref;

				// for real data
				float p_min[Npar] = {0., _x - 1, _y - 1, 0.};
				float p[Npar] = {1000., _x, _y, 1000.};
				float p_max[Npar] = {65535., _x + 1, _y + 1, 65535.};

				float y_dat[Npnt];

				// Initialization of y_dat & t
				float max_intensity = -70000;
				float min_intensity = +70000;

				for (int jj = -dimFit / 2; jj <= dimFit / 2; jj++)
					for (int ii = -dimFit / 2; ii <= dimFit / 2; ii++) {
						y_dat[(ii + dimFit / 2) + (jj + dimFit / 2) * dimFit] = (float)data[idImage * imageSize + (x + ii) + (y + jj) * WIDTH];

						max_intensity = fmaxf(max_intensity, y_dat[(ii + dimFit / 2) + (jj + dimFit / 2) * dimFit]);
						min_intensity = fminf(min_intensity, y_dat[(ii + dimFit / 2) + (jj + dimFit / 2) * dimFit]);
					}

				p[0] = max_intensity - min_intensity;
				p[3] = min_intensity;

				float lambda;

				float alpha[Npar * Npar];
				float beta[Npar];
				float y_hat[Npnt];

				float X2_old;
				float X2 = lm_matx_7x7_2D_(p/*IN*/, y_dat/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/, sigmaX);

				float max__;
				float max_ = 0;
				for (int ii = 0; ii < Npar; ii++) max_ = fmaxf(max_, fabsf(beta[ii]));

				unsigned int stop = 0;
				if (max_ < epsilon_1) stop = 1;

				lambda = lambda_0;
				X2_old = X2;

				unsigned short iteration = 0;

				float inter1[Npar * Npar];
				float inter[Npar * Npar];
				float inter3;

				float p_try[Npar];
				float delta_y;
				float X2_try;
				float m;
				float rho;

				while ((stop == 0) && (iteration < iterNum)) // Main Loop
				{
					iteration++;
					float delta_p[Npar] = {0., 0., 0., 0.};
					//delta_p[shift + (id + 0 * 16)] = 0.;
					//delta_p[shift + (id + 1 * 16)] = 0.;
					//delta_p[shift + (id + 2 * 16)] = 0.;
					//delta_p[shift + (id + 3 * 16)] = 0.;
					//delta_p[shift + (id + 4 * 16)] = 0.;

					// Copy alpha to inter
					for (int ii = 0; ii < Npar * Npar; ii++) inter[ii] = alpha[ii];
					//inter[shift_ + (id + ii * 16)] = alpha[ii];

					inter[0 * Npar + 0] += lambda * alpha[0 * Npar + 0];
					inter[1 * Npar + 1] += lambda * alpha[1 * Npar + 1];
					inter[2 * Npar + 2] += lambda * alpha[2 * Npar + 2];
					inter[3 * Npar + 3] += lambda * alpha[3 * Npar + 3];

					//inter[shift_ + (id + (0*Npar + 0) * 16)] += lambda * alpha[0*Npar + 0];
					//inter[shift_ + (id + (1*Npar + 1) * 16)] += lambda * alpha[1*Npar + 1];
					//inter[shift_ + (id + (2*Npar + 2) * 16)] += lambda * alpha[2*Npar + 2];
					//inter[shift_ + (id + (3*Npar + 3) * 16)] += lambda * alpha[3*Npar + 3];
					//inter[shift_ + (id + (4*Npar + 4) * 16)] += lambda * alpha[4*Npar + 4];

					inverse_7x7_2D_(inter, inter1);

					m = 0;
					for (int ii = 0; ii < Npar; ii++) {
						for (int cursor = 0; cursor < Npar; cursor++) delta_p[ii] += inter1[ii * Npar + cursor] * beta[cursor];
						//for(int cursor=0; cursor < Npar; cursor++)
						//	delta_p[shift + (id + ii * 16)] += inter1[ii*Npar + cursor] * beta[cursor];

						// Computing p_try
						p_try[ii] = fminf(fmaxf(p_min[ii], p[ii] + delta_p[ii]), p_max[ii]);
						/*p_try[ii] = min(max(p_min[ii], p[ii] + delta_p[shift + (id + ii * 16)]), p_max[ii]);*/

						// Computing m
						inter3 = lambda * delta_p[ii] + beta[ii]; // lambda * delta_p + beta
						m += delta_p[ii] * inter3;
						//inter3 = lambda * delta_p[shift + (id + ii * 16)] + beta[ii];// lambda * delta_p + beta
						//m += delta_p[shift + (id + ii * 16)] * inter3;
					}
					m *= 2;

					X2_try = 0;
					for (int ii = 0; ii < Npnt; ii++) {
						delta_y = y_dat[ii] - PSF_7x7_2D_((float)(ii % dimFit), (float)(ii / dimFit), p_try[0], p_try[1], p_try[2], sigmaX, p_try[3]);
						X2_try += delta_y * delta_y;
					}

					rho = (X2 - X2_try) / m;

					if (rho > epsilon_4) // it IS significantly better
					{
						X2_old = X2;
						for (int ii = 0; ii < Npar; ii++) p[ii] = p_try[ii]; // accept p_try

						X2 = lm_matx_7x7_2D_(p/*IN*/, y_dat/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/, sigmaX);

						lambda = fmaxf(lambda / lambda_DN_fac, 1.e-7);
					}
					else // it IS NOT better
					{
						X2 = X2_old; // do not accept a_try
						lambda = fminf(lambda * lambda_UP_fac, 1.e7);
					}

					max_ = 0;
					max__ = 0;
					for (int ii = 0; ii < Npar; ii++) {
						max_ = fmaxf(max_, fabsf(delta_p[ii] / p[ii]));
						//max_ = max(max_, abs(delta_p[shift + (id + ii * 16)] / p[ii]));
						max__ = fmaxf(max__, fabsf(beta[ii]));
					}

					if (((max_ < epsilon_2) || (X2 / Npnt < epsilon_3) || (max__ < epsilon_1)) && (iteration > 2)) stop = 1;

					if (iteration == iterNum) stop = 1;
				}

				// Computing statistics
				float R_sq;
				{
					///*float covar[Npar*Npar];*//* replaced by inter */
					///*float dydp[Npar*Npnt];*//* replaced by inter1*/
					/*float sigma_p[Npar];*/
					////float corr[Npar*Npar];
					////float res;
					float y_dat_mean;
					float y_hat_mean;
					float t1, t2, t3;
					float diff_hat, diff_dat;

					// parameter covariance matrix
					//inverse(alpha, inter/*covar*/);

					// standard error of parameters
					//sigma_p[0] = sqrt(/*covar*/inter[0*Npar+0]);
					//sigma_p[1] = sqrt(/*covar*/inter[1*Npar+1]);
					//sigma_p[2] = sqrt(/*covar*/inter[2*Npar+2]);
					//sigma_p[3] = sqrt(/*covar*/inter[3*Npar+3]);
					//sigma_p[4] = sqrt(/*covar*/inter[4*Npar+4]);

					// parameter correlation matrix (A vérifier)
					/*float q[Npar*Npar];*/ /* replaced by inter1 */
					//for(int ii = 0; ii < Npar; ii++) 
					//	for(int jj = 0; jj < Npar; jj++)
					//		/*q*/inter1[ii + jj * Npar] = sigma_p[ii] * sigma_p[jj];

					/*float corr [Npar * Npar];*/ /* replaced by inter */
					//for(int ii = 0; ii < Npar*Npar; ii++)
					//	/*corr*/inter[ii] = /*covar*/inter[ii] / /*q*/inter1[ii]; 

					// computing 'y_dat' mean & 'y_hat' mean
					y_dat_mean = 0;
					y_hat_mean = 0;
					for (int ii = 0; ii < Npnt; ii++) {
						y_dat_mean += y_dat[ii];
						y_hat_mean += y_hat[ii];
					}
					y_dat_mean /= Npnt;
					y_hat_mean /= Npnt;

					// computing 't1' & 't2' & 't3'
					t1 = 0;
					t2 = 0;
					t3 = 0;
					for (int ii = 0; ii < Npnt; ii++) {
						diff_dat = y_dat[ii] - y_dat_mean;
						diff_hat = y_hat[ii] - y_hat_mean;
						t1 += (diff_dat) * (diff_hat);
						t2 += (diff_dat) * (diff_dat);
						t3 += (diff_hat) * (diff_hat);
					}
					t2 = sqrtf(t2);
					t3 = sqrtf(t3);

					// computing 'R_sq'
					R_sq = pow(t1 / (t2 * t3), 2);
				}

				pointsList[NbSegParams * index + 8] = pointsList[NbSegParams * index + 2];
				pointsList[NbSegParams * index + 9] = pointsList[NbSegParams * index + 3];
				pointsList[NbSegParams * index] = sigmaX;
				pointsList[NbSegParams * index + 1] = sigmaX;
				pointsList[NbSegParams * index + 2] = 0.; // Theta

				pointsList[NbSegParams * index + 3] = /*(double)*/(p[1] - ref) + (float)x;
				pointsList[NbSegParams * index + 4] = /*(double)*/(p[2] - ref) + (float)y;
				pointsList[NbSegParams * index + 5] = p[0]; // Intensity
				pointsList[NbSegParams * index + 6] = p[3]; // Offset

				if (isnan(R_sq) || (isinf(R_sq))) pointsList[NbSegParams * index + 7] = (float)0;
				else pointsList[NbSegParams * index + 7] = R_sq;

				pointsList[NbSegParams * index + 10] = 0.;
				pointsList[NbSegParams * index + 11] = 0.;
				pointsList[NbSegParams * index + 12] = 0.;
			}
			else {
				float inter = pointsList[NbSegParams * index + 3];
				pointsList[NbSegParams * index + 3] = pointsList[NbSegParams * index + 0];
				pointsList[NbSegParams * index + 4] = pointsList[NbSegParams * index + 1];
				pointsList[NbSegParams * index] = sigmaX;
				pointsList[NbSegParams * index + 1] = sigmaX;
				pointsList[NbSegParams * index + 8] = pointsList[NbSegParams * index + 2];
				pointsList[NbSegParams * index + 9] = inter/*pointsList[NbSegParams * index + 3]*/;
				pointsList[NbSegParams * index + 2] = 0.;

				pointsList[NbSegParams * index + 5] = 0.;
				pointsList[NbSegParams * index + 6] = 0.;
				pointsList[NbSegParams * index + 7] = -1.;

				pointsList[NbSegParams * index + 10] = 0.;
				pointsList[NbSegParams * index + 11] = 0.;
				pointsList[NbSegParams * index + 12] = 0.;
			}
		}
	}
}

// version provisoire initSegBy0
__global__ void GPUconvert(float* data, unsigned short WIDTH, unsigned short HEIGHT, float minimum, float maximum, float new_minimum, float new_maximum)
{
	unsigned int index;
	unsigned short idImage = (blockIdx.x / (int)ceilf((float)WIDTH / (float)WIDTH_PORTION)); /* numero de l'image */
	unsigned int imageSize = WIDTH * HEIGHT;
	unsigned short bx = blockIdx.x % (int)ceilf((float)WIDTH / (float)WIDTH_PORTION);

	for (unsigned short y = 0; y < HEIGHT_PORTION; y++) {
		index = (threadIdx.x + bx * WIDTH_PORTION) + (blockIdx.y * HEIGHT_PORTION + y) * WIDTH;
		if (((threadIdx.x + bx * WIDTH_PORTION) < WIDTH) && ((blockIdx.y * HEIGHT_PORTION + y) < HEIGHT)) { data[idImage * imageSize + index] = (new_maximum - new_minimum) * ((data[idImage * imageSize + index] - minimum) / (maximum - minimum)); }
	}
}


extern "C"void cuda2DFittingGPU_7x7(unsigned short* data, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int iterNum, unsigned short sizeFit, float sigmaX)
{
	//{
	//	double* list = (double*)malloc(MAX_NB_STRUCTURES * numImages * 10 * sizeof(double));
	//	double* d_list;
	//	checkCudaErrors( cudaMalloc( (void**) &d_list, MAX_NB_STRUCTURES * numImages * 10 * sizeof(double) ) );
	//	checkCudaErrors( cudaMemcpy( list, d_list, MAX_NB_STRUCTURES * numImages * 10 * sizeof(double), cudaMemcpyDeviceToHost) );

	//	FILE* f = fopen("t1989.txt", "w");
	//	for(int i=0; i < MAX_NB_STRUCTURES * numImages; i++)
	//		fprintf(f, "%f, %f, %f, %f\n", list[i*10], list[i*10+1], list[i*10+2], list[i*10+3], list[i*10+4]);
	//	fclose(f);

	//	// cleanup memory
	//	free(list);
	//}

	printf("Optimized version\n");


	dim3 grid(1, 1, 1);
	dim3 threads(1, 1, 1);

	//{
	//	float* list = (float*)malloc(width * height * sizeof(float));
	//	checkCudaErrors( cudaMemcpy( list, data, width * height * sizeof(float), cudaMemcpyDeviceToHost) );

	//	FILE* f = fopen("t.txt", "w");
	//	for(int i=0; i < width * height; i++)
	//		fprintf(f, "(%d, %d) : %f\n", i%width, i/width, list[i]);
	//	fclose(f);

	//	// cleanup memory
	//	free(list);
	//}

	// Convert all pixels to [0, 1000] interval
	//GPUconvert<<<grid, threads>>>(data, width, height, gpu_min, gpu_max, 0, 1000);

	//{
	//	float* list = (float*)malloc(width * height * sizeof(float));
	//	checkCudaErrors( cudaMemcpy( list, data, width * height * sizeof(float), cudaMemcpyDeviceToHost) );

	//	FILE* f = fopen("t333.txt", "w");
	//	for(int i=0; i < width * height; i++)
	//		fprintf(f, "(%d, %d): %f\n", i%width, i/width, list[i]);
	//	fclose(f);

	//	// cleanup memory
	//	free(list);
	//}

	// allocate memory
	grid.x = (int)ceil((float)(nbPoints) / (float)FIT_BLOCK_SIZE);
	grid.y = 1;
	grid.z = 1;
	threads.x = FIT_BLOCK_SIZE;
	threads.y = 1;
	threads.z = 1; // Par rapport aux contraintes de la mémoire partagée (taille limitée à 16 Ko pour chaque block)
	// '64' au lieu de '256'
	GPUFitting_7x7_2D<<< grid, threads >>>(nbPoints, pointsList, data, nbStruct, width, height, iterNum, sizeFit, sigmaX);

	cudaError err = cudaThreadSynchronize();
	/*printf("Error : %d \n", err);*/
}

#endif // #ifndef _FITTING_KERNEL_H_

//{
//	float* list = (float*)malloc(MAX_NB_STRUCTURES * numImages * Npar * sizeof(float));
//	checkCudaErrors( cudaMemcpy( list, p, MAX_NB_STRUCTURES * numImages * Npar * sizeof(float), cudaMemcpyDeviceToHost) );

//	FILE* f = fopen("t1989.txt", "w");
//	for(int i=0; i < MAX_NB_STRUCTURES * numImages; i++)
//		fprintf(f, "%f, %f, %f, %f\n", list[i*Npar], list[i*Npar+1], list[i*Npar+2], list[i*Npar+3], list[i*Npar+4]);
//	fclose(f);

//	// cleanup memory
//	free(list);
//}
