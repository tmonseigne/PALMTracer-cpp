#pragma once

#include <helper_cuda.h>
#include <helper_timer.h>
#include <math_constants.h>
#include "../parameters.h"

extern "C"unsigned int nextPow2(unsigned int x);

extern "C"bool isPow2(unsigned int x);


#define opt
#ifdef opt
// Levenberg fit --> simple precision (float)
#define dimFit 7
#define N_PARAM 5 // I0, x0, y0, sigma, backGround
#define Npnt dimFit*dimFit
#define ref (dimFit/2)
#define NbSegParams 13

//__constant__ int MaxIter = /*20 * N_PARAM*//*20*/40;
__constant__ float epsilon_1 = 1e-6;
__constant__ float epsilon_2 = 1e-6;
__constant__ float epsilon_3 = 1e-5;
__constant__ float epsilon_4 = 1e-2;
__constant__ float lambda_0 = 1e-2;
__constant__ float lambda_UP_fac = 11;
__constant__ float lambda_DN_fac = 9;

__device__ float PSF_7x7_2D(float x, float y, float I0, float x0, float y0, float sigma, float backGround) { return I0 * expf((-powf(x - x0, 2) - powf(y - y0, 2)) / (2 * powf(sigma, 2))) + backGround; }

__device__ float deriv1I0_7x7_2D(float x, float y, float I0, float x0, float y0, float sigma, float cons) { return cons / I0; }

__device__ float deriv1x0_7x7_2D(float x, float y, float I0, float x0, float y0, float sigma, float cons) { return cons * ((x - x0) / powf(sigma, 2)); }

__device__ float deriv1y0_7x7_2D(float x, float y, float I0, float x0, float y0, float sigma, float cons) { return cons * ((y - y0) / powf(sigma, 2)); }

__device__ float deriv1sigma_7x7_2D(float x, float y, float I0, float x0, float y0, float sigma, float cons) { return (cons / powf(sigma, 3)) * (powf(x - x0, 2) + powf(y - y0, 2)); }

__device__ float lm_matx_7x7_2D(/*unsigned int* tIN,*/ float* p/*IN*/, float* y_dat/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/)
{
	// Computing Jacobian matrix (dydp) & y_hat, delta_y
	float delta_y[Npnt];
	float dydp[Npnt * N_PARAM];
	float cons;
	unsigned int x, y;
	float Chi_sq = 0.;

	for (int elt = 0; elt < Npnt; elt++) {
		x = (elt % dimFit);
		y = (elt / dimFit);
		cons = p[0] * expf((-powf((float)x - p[1], 2) - powf((float)y - p[2], 2)) / (2 * powf(p[3], 2)));
		dydp[(0) + elt * N_PARAM] = deriv1I0_7x7_2D((float)x, (float)y, p[0], p[1], p[2], p[3], cons);    // I0
		dydp[(1) + elt * N_PARAM] = deriv1x0_7x7_2D((float)x, (float)y, p[0], p[1], p[2], p[3], cons);    // x0
		dydp[(2) + elt * N_PARAM] = deriv1y0_7x7_2D((float)x, (float)y, p[0], p[1], p[2], p[3], cons);    // y0
		dydp[(3) + elt * N_PARAM] = deriv1sigma_7x7_2D((float)x, (float)y, p[0], p[1], p[2], p[3], cons); // sigma
		dydp[(4) + elt * N_PARAM] = 1;                                                                    // backGround

		y_hat[elt] = cons + p[4];
		delta_y[elt] = y_dat[elt] - y_hat[elt];

		// Computing Chi_sq
		Chi_sq += delta_y[elt] * delta_y[elt];
	}

	// Computing alpha
	for (int j = 0; j < N_PARAM; j++) {
		for (int i = 0; i < N_PARAM; i++) {
			alpha[N_PARAM * j + i] = 0;
			for (int cursor = 0; cursor < Npnt; cursor++) { alpha[N_PARAM * j + i] += dydp[cursor * N_PARAM + j] * dydp[cursor * N_PARAM + i]; } // dydp' * inter : dydp'[i*Npnt + cursor]
		}

		// Computing beta
		beta[j] = 0;
		for (int cursor = 0; cursor < Npnt; cursor++) { beta[j] += dydp[cursor * N_PARAM + j] * delta_y[cursor]; } // dydp' * inter : dydp'[j*Npnt + cursor]
	}

	return Chi_sq;
}

__device__ void inverse_7x7_2D(float* MM, float* Minv)
{
	//int shift_ = (int)(threadIdx.x/16/*32*/) * (N_PARAM * N_PARAM * 16);
	//int id = threadIdx.x % 16/*32*/;

	float M1[N_PARAM * N_PARAM];

	// Copy MM to M1
	for (int i = 0; i < N_PARAM * N_PARAM; i++) { M1[i] = MM[i]; }
	//M1[i] = MM[shift_ + (id + i * 16)];

	// Initializatiion of identity matrix
	float Mi1[N_PARAM * N_PARAM];
	for (int i = 0; i < N_PARAM * N_PARAM; i++) { Mi1[i] = 0.; }
	Mi1[0 * N_PARAM + 0] = 1;
	Mi1[1 * N_PARAM + 1] = 1;
	Mi1[2 * N_PARAM + 2] = 1;
	Mi1[3 * N_PARAM + 3] = 1;
	Mi1[4 * N_PARAM + 4] = 1;

	// Copy Mi1 to Minv
	for (int i = 0; i < N_PARAM * N_PARAM; i++) { Minv[i] = Mi1[i]; }

	for (int i = 0; i < N_PARAM; i++) {
		for (int j = 0; j < N_PARAM; j++) {
			M1[i * N_PARAM + j] = MM[i * N_PARAM + j] / MM[i * N_PARAM + i];
			Mi1[i * N_PARAM + j] = Minv[i * N_PARAM + j] / MM[i * N_PARAM + i];
			//M1[i*N_PARAM + j] = MM[shift_ + (id + (i*N_PARAM + j) * 16)] / MM[shift_ + (id + (i*N_PARAM + i) * 16)];
			//Mi1[i*N_PARAM + j] = Minv[i*N_PARAM + j] / MM[shift_ + (id + (i*N_PARAM + i) * 16)];
		}

		// Copy M1 to MM & Copy Mi1 to Minv
		for (int index = 0; index < N_PARAM * N_PARAM; index++) {
			MM[index] = M1[index];
			//MM[shift_ + (id + index * 16)] = M1[index];
			Minv[index] = Mi1[index];
		}

		for (int k = 0; k < N_PARAM; k++) { // mise zéros de la colonne
			if (k != i) {
				for (int j = 0; j < N_PARAM; j++) {
					M1[k * N_PARAM + j] = MM[k * N_PARAM + j] - MM[i * N_PARAM + j] * MM[k * N_PARAM + i];
					Mi1[k * N_PARAM + j] = Minv[k * N_PARAM + j] - Minv[i * N_PARAM + j] * MM[k * N_PARAM + i];
					//M1[k*N_PARAM+j] = MM[shift_ + (id + (k*N_PARAM+j) * 16)] - MM[shift_ + (id + (i*N_PARAM+j) * 16)] * MM[shift_ + (id + (k*N_PARAM+i) * 16)];
					//Mi1[k*N_PARAM+j] = Minv[k*N_PARAM+j] - Minv[i*N_PARAM+j] * MM[shift_ + (id + (k*N_PARAM+i) * 16)];
				}
			}
		}

		// Copy M1 to MM & Copy Mi1 to Minv
		for (int index = 0; index < N_PARAM * N_PARAM; index++) {
			MM[index] = M1[index];
			//MM[shift_ + (id + index * 16)] = M1[index];
			Minv[index] = Mi1[index];
		}
	}
}

__global__ void GPUFitting_Sigma_7x7_2D(unsigned int nbPoints, float* pointsList, unsigned short* data, unsigned short* nbStruct, unsigned short WIDTH, unsigned short HEIGHT, unsigned int iterNum, unsigned short _size, float sigmaX)
{
	unsigned int index = threadIdx.x + blockIdx.x * blockDim.x;

	if (index < nbPoints) {
		unsigned short idImage = (unsigned short)pointsList[NbSegParams * index + 10];
		unsigned int imageSize = WIDTH * HEIGHT;
		//unsigned short _seg = (unsigned short)pointsList[NbSegParams * index + 9];

		if ((int)pointsList[index * NbSegParams + 3] != 0) {
			if ((pointsList[NbSegParams * index + 1] >= 3) && (pointsList[NbSegParams * index + 1] <= HEIGHT - 4) && (pointsList[NbSegParams * index] >= 3) && (pointsList[NbSegParams * index] <= WIDTH - 4)) {
				//int shift = (int)(threadIdx.x/16/*32*/) * (N_PARAM * 16);
				//int shift_ = (int)(threadIdx.x/16/*32*/) * (N_PARAM * N_PARAM * 16);
				//int id = threadIdx.x % 16/*32*/;

				int x = /*(int)*/floor(pointsList[NbSegParams * index]);
				int y = /*(int)*/floor(pointsList[NbSegParams * index + 1]);

				float _x = (pointsList[NbSegParams * index] - x) + ref;
				float _y = (pointsList[NbSegParams * index + 1] - y) + ref;

				// for real data
				float p_min[N_PARAM] = {0., _x - 1, _y - 1, 0.1, 0.};
				float p[N_PARAM] = {1000., _x, _y, sigmaX, 1000.};
				float p_max[N_PARAM] = {65535., _x + 1, _y + 1, 20., 65535.};

				float y_dat[Npnt];

				// Initialization of y_dat & t
				float max_intensity = -70000;
				float min_intensity = +70000;

				for (int jj = -dimFit / 2; jj <= dimFit / 2; jj++) {
					for (int ii = -dimFit / 2; ii <= dimFit / 2; ii++) {
						y_dat[(ii + dimFit / 2) + (jj + dimFit / 2) * dimFit] = (float)data[idImage * imageSize + (x + ii) + (y + jj) * WIDTH];


						max_intensity = fmaxf(max_intensity, y_dat[(ii + dimFit / 2) + (jj + dimFit / 2) * dimFit]);
						min_intensity = fminf(min_intensity, y_dat[(ii + dimFit / 2) + (jj + dimFit / 2) * dimFit]);
					}
				}

				p[0] = max_intensity - min_intensity;
				p[4] = min_intensity;

				float lambda;

				float alpha[N_PARAM * N_PARAM];
				float beta[N_PARAM];
				float y_hat[Npnt];

				float X2_old;
				float X2 = lm_matx_7x7_2D(p/*IN*/, y_dat/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/);

				float max__;
				float max_ = 0;
				for (int ii = 0; ii < N_PARAM; ii++) { max_ = fmaxf(max_, fabsf(beta[ii])); }

				unsigned int stop = 0;
				if (max_ < epsilon_1) { stop = 1; }

				lambda = lambda_0;
				X2_old = X2;

				unsigned short iteration = 0;

				float inter1[N_PARAM * N_PARAM];
				float inter[N_PARAM * N_PARAM];
				float inter3;

				float p_try[N_PARAM];
				float delta_y;
				float X2_try;
				float m;
				float rho;

				while ((stop == 0) && (iteration < iterNum)) // Main Loop
				{
					iteration++;
					float delta_p[N_PARAM] = {0., 0., 0., 0., 0.};

					// Copy alpha to inter
					for (int ii = 0; ii < N_PARAM * N_PARAM; ii++) { inter[ii] = alpha[ii]; }
					//inter[shift_ + (id + ii * 16)] = alpha[ii];

					inter[0 * N_PARAM + 0] += lambda * alpha[0 * N_PARAM + 0];
					inter[1 * N_PARAM + 1] += lambda * alpha[1 * N_PARAM + 1];
					inter[2 * N_PARAM + 2] += lambda * alpha[2 * N_PARAM + 2];
					inter[3 * N_PARAM + 3] += lambda * alpha[3 * N_PARAM + 3];
					inter[4 * N_PARAM + 4] += lambda * alpha[4 * N_PARAM + 4];

					inverse_7x7_2D(inter, inter1);

					m = 0;
					for (int ii = 0; ii < N_PARAM; ii++) {
						for (int cursor = 0; cursor < N_PARAM; cursor++) { delta_p[ii] += inter1[ii * N_PARAM + cursor] * beta[cursor]; }

						// Computing p_try
						p_try[ii] = fminf(fmaxf(p_min[ii], p[ii] + delta_p[ii]), p_max[ii]);

						// Computing m
						inter3 = lambda * delta_p[ii] + beta[ii]; // lambda * delta_p + beta
						m += delta_p[ii] * inter3;
					}
					m *= 2;

					X2_try = 0;
					for (int ii = 0; ii < Npnt; ii++) {
						delta_y = y_dat[ii] - PSF_7x7_2D((float)(ii % dimFit), (float)(ii / dimFit), p_try[0], p_try[1], p_try[2], p_try[3], p_try[4]);
						X2_try += delta_y * delta_y;
					}

					rho = (X2 - X2_try) / m;

					if (rho > epsilon_4) // it IS significantly better
					{
						X2_old = X2;
						for (int ii = 0; ii < N_PARAM; ii++) {
							p[ii] = p_try[ii]; // accept p_try
						}

						X2 = lm_matx_7x7_2D(p/*IN*/, y_dat/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/);

						lambda = fmaxf(lambda / lambda_DN_fac, 1.e-7);
					}
					else // it IS NOT better
					{
						X2 = X2_old; // do not accept a_try
						lambda = fminf(lambda * lambda_UP_fac, 1.e7);
					}

					max_ = 0;
					max__ = 0;
					for (int ii = 0; ii < N_PARAM; ii++) {
						max_ = fmaxf(max_, fabsf(delta_p[ii] / p[ii]));
						max__ = fmaxf(max__, fabsf(beta[ii]));
					}

					if (((max_ < epsilon_2) || (X2 / Npnt < epsilon_3) || (max__ < epsilon_1)) && (iteration > 2)) { stop = 1; }

					if (iteration == iterNum) { stop = 1; }
				}

				// Computing statistics
				float R_sq;
				{
					float y_dat_mean;
					float y_hat_mean;
					float t1, t2, t3;
					float diff_hat, diff_dat;

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
				pointsList[NbSegParams * index] = p[3];
				pointsList[NbSegParams * index + 1] = p[3];
				pointsList[NbSegParams * index + 2] = 0.;

				pointsList[NbSegParams * index + 3] = /*(double)*/(p[1] - ref) + (float)x;
				pointsList[NbSegParams * index + 4] = /*(double)*/(p[2] - ref) + (float)y;
				pointsList[NbSegParams * index + 5] = p[0];
				pointsList[NbSegParams * index + 6] = p[4];

				if (isnan(R_sq) || (isinf(R_sq))) { pointsList[NbSegParams * index + 7] = (float)0; }
				else { pointsList[NbSegParams * index + 7] = R_sq; }

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

extern "C"void cuda2DFittingGPU_Sigma_7x7(unsigned short* data, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int iterNum, unsigned short sizeFit, float sigmaX)
{
	printf("Optimized version\n");
	dim3 grid(1, 1, 1);
	dim3 threads(WIDTH_PORTION, 1, 1);

	grid.x = (int)ceil((float)(nbPoints) / (float)FIT_BLOCK_SIZE);
	grid.y = 1;
	grid.z = 1;
	threads.x = FIT_BLOCK_SIZE;
	threads.y = 1;
	threads.z = 1; // Par rapport aux contraintes de la mémoire partagée (taille limitée à 16 Ko pour chaque block)
	// '64' au lieu de '256'
	GPUFitting_Sigma_7x7_2D<<< grid, threads >>>(nbPoints, pointsList, data, nbStruct, width, height, iterNum, sizeFit, sigmaX);
	cudaError err = cudaThreadSynchronize();
}

#endif

