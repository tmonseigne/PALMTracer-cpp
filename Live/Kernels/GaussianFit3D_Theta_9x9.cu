#pragma once

#include <helper_cuda.h>
#include <helper_timer.h>
#include <math_constants.h>
#include "../parameters.h"

// Levenberg fit --> simple precision (float)
#define dimFit 9
#define N_PARAM 7 // I0, x0, y0, sigmaX, sigmaY, angle, offset
#define Npnt dimFit*dimFit
#define ref (dimFit/2)
#define NbSegParams 13

__constant__ float epsilon_1 = 1e-6;
__constant__ float epsilon_2 = 1e-6;
__constant__ float epsilon_3 = 1e-5;
__constant__ float epsilon_4 = 1e-2;
__constant__ float lambda_0 = 1e-2;
__constant__ float lambda_UP_fac = 11;
__constant__ float lambda_DN_fac = 9;

__device__ float PSF_9x9_(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float COS, float SIN) { return I0 * expf(-powf(((x - x0) * COS - (y - y0) * SIN) / sigmaX, 2) - powf(((x - x0) * SIN + (y - y0) * COS) / sigmaY, 2)) + offset; }

__device__ float deriv1I0_9x9_(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float cons, float COS, float SIN) { return cons / I0; }

__device__ float deriv1x0_9x9_(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float cons, float COS, float SIN)
{
	float alpha = -((-2 * powf(COS, 2) * (x - x0) + 2 * (y - y0) * COS * SIN) / powf(sigmaX, 2));
	float beta = -((-2 * powf(SIN, 2) * (x - x0) - 2 * (y - y0) * COS * SIN) / powf(sigmaY, 2));
	return cons * (alpha + beta);
}

__device__ float deriv1y0_9x9_(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float cons, float COS, float SIN)
{
	float alpha = -((2 * powf(SIN, 2) * (y - y0) - 2 * (x - x0) * COS * SIN) / powf(sigmaX, 2));
	float beta = -((-2 * powf(COS, 2) * (y - y0) - 2 * (x - x0) * COS * SIN) / powf(sigmaY, 2));
	return cons * (alpha + beta);
}


__device__ float deriv1sigmaX_9x9_(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float cons, float COS, float SIN)
{
	float alpha = powf((x - x0) * COS, 2) - 2 * (x - x0) * (y - y0) * COS * SIN + powf((y - y0) * SIN, 2);
	return cons * ((2 * alpha) / powf(sigmaX, 3));
}

__device__ float deriv1sigmaY_9x9_(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float cons, float COS, float SIN)
{
	float alpha = powf((x - x0) * SIN, 2) + 2 * (x - x0) * (y - y0) * COS * SIN + powf((y - y0) * COS, 2);
	return cons * ((2 * alpha) / powf(sigmaY, 3));
}

__device__ float deriv1angle_9x9_(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float cons, float COS, float SIN)
{
	float a = powf(x - x0, 2);
	float b = powf(y - y0, 2);
	float sX2 = powf(sigmaX, 2);
	float sY2 = powf(sigmaY, 2);
	return cons * (2 * sinf(angle) * cosf(angle) * ((a - b) / sX2 + (b - a) / sY2) - powf(sinf(angle), 2) * ((2 * (x - x0) * (y - y0)) / sX2 + (-2 * (x - x0) * (y - y0)) / sY2) + powf(cosf(angle), 2) * ((2 * (x - x0) * (y - y0)) / sX2 + (-2 * (x - x0) * (y - y0)) / sY2));
}

__device__ float lm_matx_9x9_(/*unsigned int* tIN,*/ float* p/*IN*/, float* y_dat/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/)
{
	// Computing Jacobian matrix (dydp) & y_hat, delta_y
	float delta_y[Npnt];
	float dydp[Npnt * N_PARAM];
	float cons;
	unsigned int x, y;
	float Chi_sq = 0.;

	float COS = cosf(p[5]);
	float SIN = sinf(p[5]);

	for (int elt = 0; elt < Npnt; elt++) {
		x = (elt % dimFit);
		y = (elt / dimFit);
		cons = p[0] * expf(-powf((((float)x - p[1]) * COS - ((float)y - p[2]) * SIN) / p[3], 2) - powf((((float)x - p[1]) * SIN + ((float)y - p[2]) * COS) / p[4], 2));
		dydp[(0) + elt * N_PARAM] = deriv1I0_9x9_((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], p[6], cons, COS, SIN);     // I0
		dydp[(1) + elt * N_PARAM] = deriv1x0_9x9_((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], p[6], cons, COS, SIN);     // x0
		dydp[(2) + elt * N_PARAM] = deriv1y0_9x9_((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], p[6], cons, COS, SIN);     // y0
		dydp[(3) + elt * N_PARAM] = deriv1sigmaX_9x9_((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], p[6], cons, COS, SIN); // sigmaX
		dydp[(4) + elt * N_PARAM] = deriv1sigmaY_9x9_((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], p[6], cons, COS, SIN); // sigmaY
		dydp[(5) + elt * N_PARAM] = deriv1angle_9x9_((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], p[6], cons, COS, SIN);  // angle
		dydp[(6) + elt * N_PARAM] = 1;                                                                                               // offset

		y_hat[elt] = cons + p[6]; // PSF
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

__device__ void inverse_matrix_9x9_(float* MM, float* Minv)
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
	Mi1[5 * N_PARAM + 5] = 1;
	Mi1[6 * N_PARAM + 6] = 1;

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

		for (int k = 0; k < N_PARAM; k++) // mise z�ros de la colonne
		{
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

__global__ void GPU3DFitting_Theta_9x9_(unsigned int nbPoints, float* pointsList, unsigned short* data, unsigned short* nbStruct, unsigned short WIDTH, unsigned short HEIGHT, unsigned int iterNum, unsigned short _size, float sigmaX, float sigmaY, float angle)
{
	unsigned int index = threadIdx.x + blockIdx.x * blockDim.x;

	if (index < nbPoints) {
		unsigned short idImage = (unsigned short)pointsList[NbSegParams * index + 10];
		unsigned int imageSize = WIDTH * HEIGHT;

		if ((int)pointsList[index * NbSegParams + 3] != 0) {
			if ((pointsList[NbSegParams * index + 1] >= (int)ref) && (pointsList[NbSegParams * index + 1] <= HEIGHT - (int)ref - 1) && (pointsList[NbSegParams * index] >= (int)ref) && (pointsList[NbSegParams * index] <= WIDTH - (int)ref - 1)) {
				int x = (int)pointsList[NbSegParams * index];
				int y = (int)pointsList[NbSegParams * index + 1];

				float _x = (pointsList[NbSegParams * index] - x) + ref;
				float _y = (pointsList[NbSegParams * index + 1] - y) + ref;

				// for real data
				float p_min[N_PARAM] = {0., _x - 1, _y - 1, 0.1, 0.1, -150., 0.};
				float p[N_PARAM] = {2000., _x, _y, sqrtf(2.) * sigmaX, sqrtf(2.) * sigmaY, angle, 700.};
				float p_max[N_PARAM] = {65535., _x + 1, _y + 1, 20., 20., +150., 65535.};

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
				p[6] = min_intensity;

				float lambda;

				float alpha[N_PARAM * N_PARAM];
				float beta[N_PARAM];
				float y_hat[Npnt];

				float X2_old;
				float X2 = lm_matx_9x9_(p/*IN*/, y_dat/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/);

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

				float COS = cosf(p[5]);
				float SIN = sinf(p[5]);


				while ((stop == 0) && (iteration <= iterNum)) // Main Loop
				{
					iteration++;
					float delta_p[N_PARAM] = {0., 0., 0., 0., 0., 0., 0.};

					// Copy alpha to inter
					for (int ii = 0; ii < N_PARAM * N_PARAM; ii++) { inter[ii] = alpha[ii]; }

					inter[0 * N_PARAM + 0] += lambda * alpha[0 * N_PARAM + 0];
					inter[1 * N_PARAM + 1] += lambda * alpha[1 * N_PARAM + 1];
					inter[2 * N_PARAM + 2] += lambda * alpha[2 * N_PARAM + 2];
					inter[3 * N_PARAM + 3] += lambda * alpha[3 * N_PARAM + 3];
					inter[4 * N_PARAM + 4] += lambda * alpha[4 * N_PARAM + 4];
					inter[5 * N_PARAM + 5] += lambda * alpha[5 * N_PARAM + 5];
					inter[6 * N_PARAM + 6] += lambda * alpha[6 * N_PARAM + 6];

					inverse_matrix_9x9_(inter, inter1);

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
						delta_y = y_dat[ii] - PSF_9x9_((float)(ii % dimFit), (float)(ii / dimFit), p_try[0], p_try[1], p_try[2], p_try[3], p_try[4], p_try[5], p_try[6], COS, SIN);
						X2_try += delta_y * delta_y;
					}

					rho = (X2 - X2_try) / m;

					if (rho > epsilon_4) // it IS significantly better
					{
						X2_old = X2;
						for (int ii = 0; ii < N_PARAM; ii++) {
							p[ii] = p_try[ii]; // accept p_try
						}

						X2 = lm_matx_9x9_(p/*IN*/, y_dat/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/);

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
				float R_sq = 0.;
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
					R_sq = powf(t1 / (t2 * t3), 2);
				}

				pointsList[NbSegParams * index] = p[3] / sqrtf(2.);
				pointsList[NbSegParams * index + 1] = p[4] / sqrtf(2.);
				pointsList[NbSegParams * index + 8] = pointsList[NbSegParams * index + 2];
				pointsList[NbSegParams * index + 9] = pointsList[NbSegParams * index + 3];
				pointsList[NbSegParams * index + 2] = p[5];
				pointsList[NbSegParams * index + 3] = (p[1] - ref) + (float)x;
				pointsList[NbSegParams * index + 4] = (p[2] - ref) + (float)y;
				pointsList[NbSegParams * index + 5] = p[0];
				pointsList[NbSegParams * index + 6] = p[6];

				if (isnan(R_sq) || (isinf(R_sq))) { pointsList[NbSegParams * index + 7] = (float)0; }
				else { pointsList[NbSegParams * index + 7] = R_sq; }

				pointsList[NbSegParams * index + 11] = 0.;
				pointsList[NbSegParams * index + 12] = 0.;
			}
			else {
				float inter = pointsList[NbSegParams * index + 3];
				pointsList[NbSegParams * index + 3] = pointsList[NbSegParams * index + 0];
				pointsList[NbSegParams * index + 4] = pointsList[NbSegParams * index + 1];
				pointsList[NbSegParams * index] = sigmaX;
				pointsList[NbSegParams * index + 1] = sigmaY;
				pointsList[NbSegParams * index + 8] = pointsList[NbSegParams * index + 2];
				pointsList[NbSegParams * index + 9] = inter/*pointsList[NbSegParams * index + 3]*/;
				pointsList[NbSegParams * index + 2] = angle;

				pointsList[NbSegParams * index + 5] = 0.;
				pointsList[NbSegParams * index + 6] = 0.;
				pointsList[NbSegParams * index + 7] = -1.;

				pointsList[NbSegParams * index + 11] = 0.;
				pointsList[NbSegParams * index + 12] = 0.;
			}
		}
	}
}

extern "C"void cuda3DFittingGPU_Theta_9x9(unsigned short* data, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints, unsigned int iterNum, unsigned short sizeFit, float sigmaX, float sigmaY, float angle)
{
	printf("3D fitting\n");

	dim3 grid(1, 1, 1);
	dim3 threads(1, 1, 1);


	// allocate memory
	grid.x = (int)ceil((float)(nbPoints) / (float)FIT_BLOCK_SIZE);
	grid.y = 1;
	grid.z = 1;
	threads.x = FIT_BLOCK_SIZE;
	threads.y = 1;
	threads.z = 1; // Par rapport aux contraintes de la m�moire partag�e (taille limit�e � 16 Ko pour chaque block)
	// '64' au lieu de '256'
	GPU3DFitting_Theta_9x9_<<< grid, threads >>>(nbPoints, pointsList, data, nbStruct, width, height, iterNum, sizeFit, sigmaX, sigmaY, angle);
}


