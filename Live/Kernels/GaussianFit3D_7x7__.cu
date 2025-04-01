#pragma once

#include <math_constants.h>
#include <cutil_inline.h>
#include "../parameters.h"

// Levenberg fit --> simple precision (float)
#ifndef dimFit
	#define dimFit 7
#endif
#define Npar 6 // I0, x0, y0, sigmaX, sigmaY, angle, backGround
#define Npnt dimFit*dimFit
#define ref (dimFit/2)
#define NbSegParams 11

__constant__ float epsilon_1 = 1e-6;
__constant__ float epsilon_2 = 1e-6;
__constant__ float epsilon_3 = 1e-5;
__constant__ float epsilon_4 = 1e-2;
__constant__ float lambda_0  = 1e-2;
__constant__ float lambda_UP_fac = 11;
__constant__ float lambda_DN_fac = 9;

__device__ float PSF_7x7(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float backGround, float COS, float SIN)
{
	return I0 * expf( -powf(((x-x0)*COS-(y-y0)*SIN)/sigmaX, 2) - powf(((x-x0)*SIN+(y-y0)*COS)/sigmaY, 2) ) + backGround;
 }

__device__ float deriv1I0_7x7(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float backGround, float cons, float COS, float SIN)
{
	return cons/I0;
}

__device__ float deriv1x0_7x7(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float backGround, float cons, float COS, float SIN)
{
	float alpha = -( (-2 * powf(COS, 2) * (x - x0) +2 * (y - y0) * COS * SIN ) / powf(sigmaX, 2) );
	float beta = -( (-2 * powf(SIN, 2) * (x - x0) -2 * (y - y0) * COS * SIN ) / powf(sigmaY, 2) );
	return cons * (alpha + beta);
}

__device__ float deriv1y0_7x7(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float backGround, float cons, float COS, float SIN)
{
	float alpha = -( (2 * powf(SIN, 2) * (y - y0) -2 * (x - x0) * COS * SIN ) / powf(sigmaX, 2) );
	float beta = -( (-2 * powf(COS, 2) * (y - y0) -2 * (x - x0) * COS * SIN ) / powf(sigmaY, 2) );
	return cons * (alpha + beta);
}

__device__ float deriv1sigmaX_7x7(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float backGround, float cons, float COS, float SIN)
{
	float alpha = powf((x-x0)*COS, 2) - 2 * (x-x0) * (y-y0) * COS * SIN + powf((y-y0)*SIN, 2); 
	return cons * ( (2*alpha) / powf(sigmaX, 3) );
}

__device__ float deriv1sigmaY_7x7(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float backGround, float cons, float COS, float SIN)
{
	float alpha = powf((x-x0)*SIN, 2) + 2 * (x-x0) * (y-y0) * COS * SIN + powf((y-y0)*COS, 2); 
	return cons * ( (2*alpha) / powf(sigmaY, 3) );
}

__device__ float lm_matx_7x7(/*unsigned int* tIN,*/ float* p/*IN*/, float* y_dat/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/, float angle)
{
	// Computing Jacobian matrix (dydp) & y_hat, delta_y
	float delta_y[Npnt];
	float dydp[Npnt * Npar];
	float cons;
	unsigned int x, y;
	float Chi_sq = 0.;

	float COS = cosf(angle);
	float SIN = sinf(angle);

	for(int elt=0; elt < Npnt; elt++)
	{
		x = (elt%dimFit);
		y = (elt/dimFit);
		cons = p[0] * expf( -powf((((float)x - p[1])*COS-((float)y - p[2])*SIN)/p[3], 2) - powf((((float)x - p[1])*SIN+((float)y - p[2])*COS)/p[4], 2) );
		dydp[(0) + elt * Npar] = deriv1I0_7x7   ((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], cons, COS, SIN);  // I0
		dydp[(1) + elt * Npar] = deriv1x0_7x7   ((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], cons, COS, SIN);  // x0
		dydp[(2) + elt * Npar] = deriv1y0_7x7   ((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], cons, COS, SIN);  // y0
		dydp[(3) + elt * Npar] = deriv1sigmaX_7x7((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], cons, COS, SIN); // sigmaX
		dydp[(4) + elt * Npar] = deriv1sigmaY_7x7((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], cons, COS, SIN); // sigmaY
		dydp[(5) + elt * Npar] = 1;																			           // backGround

		y_hat[elt] = cons + p[5]; // PSF
		delta_y[elt] = y_dat[elt] - y_hat[elt];

		// Computing Chi_sq
		Chi_sq += delta_y[elt] * delta_y[elt];
	}

	// Computing alpha
	for(int j=0; j < Npar; j++)
	{
		for(int i=0; i < Npar; i++)
		{
			alpha[Npar * j + i] = 0;
			for(int cursor=0; cursor < Npnt; cursor++)
				alpha[Npar * j + i] += dydp[cursor*Npar + j] * dydp[cursor*Npar + i];  // dydp' * inter : dydp'[i*Npnt + cursor]
		}

		// Computing beta
		beta[j] = 0;
		for(int cursor=0; cursor < Npnt; cursor++)
			beta[j] += dydp[cursor*Npar + j] * delta_y[cursor];  // dydp' * inter : dydp'[j*Npnt + cursor]
	}

	return Chi_sq;
}

__device__ void inverse_7x7(float* MM, float* Minv)
{
	//int shift_ = (int)(threadIdx.x/16/*32*/) * (Npar * Npar * 16);
	//int id = threadIdx.x % 16/*32*/;

	float M1[Npar * Npar];

	// Copy MM to M1
	for(int i=0; i< Npar*Npar; i++)
		M1[i] = MM[i];
		//M1[i] = MM[shift_ + (id + i * 16)];

	// Initializatiion of identity matrix
	float Mi1[Npar * Npar];
	for(int i=0; i < Npar*Npar; i++)    Mi1[i] = 0.;
	Mi1[0*Npar + 0] = 1;
	Mi1[1*Npar + 1] = 1;
	Mi1[2*Npar + 2] = 1;
	Mi1[3*Npar + 3] = 1;
	Mi1[4*Npar + 4] = 1;
	Mi1[5*Npar + 5] = 1;

	// Copy Mi1 to Minv
	for(int i=0; i< Npar*Npar; i++)
		Minv[i] = Mi1[i];

	for(int i=0; i<Npar; i++)
	{
		for(int j=0; j<Npar; j++)
		{
			M1[i*Npar + j] = MM[i*Npar + j] / MM[i*Npar + i];
			Mi1[i*Npar + j] = Minv[i*Npar + j] / MM[i*Npar + i];
			//M1[i*Npar + j] = MM[shift_ + (id + (i*Npar + j) * 16)] / MM[shift_ + (id + (i*Npar + i) * 16)];
			//Mi1[i*Npar + j] = Minv[i*Npar + j] / MM[shift_ + (id + (i*Npar + i) * 16)];
		}

		// Copy M1 to MM & Copy Mi1 to Minv
		for(int index=0; index < Npar*Npar; index++)
		{
			MM[index] = M1[index];
			//MM[shift_ + (id + index * 16)] = M1[index];
			Minv[index] = Mi1[index];
		}

		for(int k=0; k<Npar; k++)// mise zéros de la colonne
			if (k != i)
				for(int j=0; j<Npar; j++)
				{
					M1[k*Npar+j] = MM[k*Npar+j] - MM[i*Npar+j] * MM[k*Npar+i];
					Mi1[k*Npar+j] = Minv[k*Npar+j] - Minv[i*Npar+j] * MM[k*Npar+i];
					//M1[k*Npar+j] = MM[shift_ + (id + (k*Npar+j) * 16)] - MM[shift_ + (id + (i*Npar+j) * 16)] * MM[shift_ + (id + (k*Npar+i) * 16)];
					//Mi1[k*Npar+j] = Minv[k*Npar+j] - Minv[i*Npar+j] * MM[shift_ + (id + (k*Npar+i) * 16)];
				}
		
		// Copy M1 to MM & Copy Mi1 to Minv
		for(int index=0; index < Npar*Npar; index++)
		{
			MM[index] = M1[index];
			//MM[shift_ + (id + index * 16)] = M1[index];
			Minv[index] = Mi1[index];
		}
	}
}

__global__ void GPU3DFitting_7x7(unsigned int nbPoints, float* pointsList, unsigned short* data, unsigned short* nbStruct, unsigned short WIDTH, unsigned short HEIGHT, unsigned int iterNum,
							 unsigned short _size, float sigmaX, float sigmaY, float angle, unsigned short computeZ, float* SigmaXCoefs, float* SigmaYCoefs, float ThreeDThickness/*um*/, float ThreeDSampling/*um*/)
{
	unsigned int index = threadIdx.x + blockIdx.x * blockDim.x;

	if (index < nbPoints)
	{
		unsigned short idImage = (unsigned short)pointsList[NbSegParams * index + 10];
		unsigned int imageSize = WIDTH * HEIGHT;
		//unsigned short _seg = (unsigned short)pointsList[NbSegParams * index + 9];

		if ((int)pointsList[index * NbSegParams + 3] != 0)
		{
			{
				int x = (int)pointsList[NbSegParams * index];
				int y = (int)pointsList[NbSegParams * index + 1];

				float _x = (pointsList[NbSegParams * index]     - x) + ref;
				float _y = (pointsList[NbSegParams * index + 1] - y) + ref;			

				// for real data
				float p_min[Npar] = {0.    , _x - 1, _y - 1, 0.1              , 0.1              , 0.    };
				float p[Npar]     = {2000. , _x    , _y    , sqrtf(2.) * sigmaX, sqrtf(2.) * sigmaY, 700.  };
				float p_max[Npar] = {65535., _x + 1, _y + 1, 20.              , 20.              , 65535.};

				float y_dat[Npnt];
				//float* y_dat = (float*)malloc(Npnt * sizeof(float));

				// Initialization of y_dat & t
				float max_intensity = -70000;
				float min_intensity = +70000;
				int shiftX = 0;
				int shiftY = 0;

				for(int jj = -dimFit/2; jj <= dimFit/2; jj++)
					for(int ii = -dimFit/2; ii <= dimFit/2; ii++)
					{
						shiftX = x + ii;
						shiftY = y + jj;

						if ((x + ii) < 0 )      shiftX = -ii;
						if ((x + ii) >= WIDTH)  shiftX = -ii;
						if ((y + jj) < 0 )      shiftY = +jj;
						if ((y + jj) >= HEIGHT) shiftY = -jj;

						y_dat[(ii + dimFit/2) + (jj + dimFit/2) * dimFit] = (float)data[idImage * imageSize + (shiftX) + (shiftY) * WIDTH];

						max_intensity = fmaxf(max_intensity, y_dat[(ii + dimFit/2) + (jj + dimFit/2) * dimFit]);
						min_intensity = fminf(min_intensity, y_dat[(ii + dimFit/2) + (jj + dimFit/2) * dimFit]);
					}

				p[0] = max_intensity - min_intensity;
				p[5] = min_intensity;

				float lambda;

				float alpha[Npar * Npar];
				float beta[Npar];
				float y_hat[Npnt];
				
				float X2_old;
				float X2 = lm_matx_7x7(p/*IN*/, y_dat/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/, angle);

				float max__;
				float max_ = 0;
				for(int ii=0; ii < Npar; ii++)
						max_ = fmaxf(max_, fabsf(beta[ii]));

				unsigned int stop = 0;
				if(max_ < epsilon_1) 
					stop = 1;

				lambda = lambda_0;
				X2_old = X2;

				unsigned short iteration = 0;

				float inter1[Npar * Npar];
				float inter[Npar * Npar];
				float inter3;

				float p_try  [Npar];
				float delta_y;
				float X2_try;
				float m;
				float rho;

				float COS = cosf(angle);
				float SIN = sinf(angle);


				while ((stop == 0) && (iteration <= iterNum)) // Main Loop
				{
					iteration++;
					float delta_p[Npar] = {0., 0., 0., 0., 0., 0.};

					// Copy alpha to inter
					for(int ii=0; ii < Npar*Npar; ii++)
						inter[ii] = alpha[ii];

					inter[0*Npar + 0] += lambda * alpha[0*Npar + 0];
					inter[1*Npar + 1] += lambda * alpha[1*Npar + 1];
					inter[2*Npar + 2] += lambda * alpha[2*Npar + 2];
					inter[3*Npar + 3] += lambda * alpha[3*Npar + 3];
					inter[4*Npar + 4] += lambda * alpha[4*Npar + 4];
					inter[5*Npar + 5] += lambda * alpha[5*Npar + 5];


					inverse_7x7(inter, inter1);

					m = 0;
					for(int ii=0; ii < Npar; ii++)
					{
						for(int cursor=0; cursor < Npar; cursor++)
							delta_p[ii] += inter1[ii*Npar + cursor] * beta[cursor];

						// Computing p_try
						p_try[ii] = fminf(fmaxf(p_min[ii], p[ii] + delta_p[ii]), p_max[ii]);

						// Computing m
						inter3 = lambda * delta_p[ii] + beta[ii];// lambda * delta_p + beta
						m += delta_p[ii] * inter3;
					}
					m *= 2;

					X2_try = 0;
					for(int ii=0; ii < Npnt; ii++)
					{
						delta_y = y_dat[ii] - PSF_7x7((float)(ii%dimFit), (float)(ii/dimFit), p_try[0], p_try[1], p_try[2], p_try[3], p_try[4], p_try[5], COS, SIN);
						X2_try += delta_y * delta_y;
					}				

					rho = (X2 - X2_try) / m;

					if(rho > epsilon_4) // it IS significantly better
					{
						X2_old = X2;
						for(int ii=0; ii < Npar; ii++)
							p[ii] = p_try[ii]; // accept p_try

						X2 = lm_matx_7x7(p/*IN*/, y_dat/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/, angle);

						lambda = fmaxf(lambda/lambda_DN_fac, 1.e-7);
					}
					else  // it IS NOT better
					{
						X2 = X2_old;  // do not accept a_try
						lambda = fminf(lambda*lambda_UP_fac, 1.e7);
					}

					max_ = 0;
					max__ = 0;
					for(int ii=0; ii < Npar; ii++)
					{
						max_ = fmaxf(max_, fabsf(delta_p[ii] / p[ii]));
						max__ = fmaxf(max__, fabsf(beta[ii]));
					}

					if( ( (max_ < epsilon_2) || (X2/Npnt < epsilon_3) || (max__ < epsilon_1) ) && (iteration > 2) ) 
						stop = 1;

					if(iteration == iterNum)
						stop = 1;
				}

				// Computing statistics
				float R_sq = 0.;
				float R = 0;
				{
					float y_dat_mean;
					float y_hat_mean;
					float t1 ,t2, t3;
					float diff_hat, diff_dat;

					// computing 'y_dat' mean & 'y_hat' mean
					y_dat_mean = 0;
					y_hat_mean = 0;
					for(int ii=0; ii < Npnt; ii++)
					{
						y_dat_mean += y_dat[ii];
						y_hat_mean += y_hat[ii];
					}
					y_dat_mean /= Npnt;
					y_hat_mean /= Npnt;

					// computing 't1' & 't2' & 't3'
					t1 = 0;
					t2 = 0;	
					t3 = 0;
					for(int ii=0; ii < Npnt; ii++)
					{
						diff_dat = y_dat[ii] - y_dat_mean;
						diff_hat = y_hat[ii] - y_hat_mean;
						t1 += (diff_dat) * (diff_hat);
						t2 += (diff_dat) * (diff_dat);
						t3 += (diff_hat) * (diff_hat);
					}
					t2 = sqrtf(t2);
					t3 = sqrtf(t3);

					// computing 'R_sq'
					R_sq = (t1 / (t2 * t3)) * (t1 / (t2 * t3)); //powf(t1 / (t2 * t3), 2);
					R = t3 / t2;
				}

				// Z extraction
				if (computeZ == 1)
				{
					float Zval = 999999;
					float Zerror = 1000.0;
					
						float ZYfunc, ZXfunc;
						float a;
						float val;
						float delta;
						float val_pixel;
						float sX = p[3]/sqrtf(2.);
						float sY = p[4]/sqrtf(2.);

						for(delta = (int)-(ThreeDThickness * 500); delta <= (int)(ThreeDThickness * 500); delta += (int)(ThreeDSampling * 1000)) // Delta en nm pour valeurs entières
						{
							val = delta / 1000.0; // Retour en um pour calibration

							//INVERSION SIGMA_X et SIGMA_Y

							ZYfunc = SigmaXCoefs[4] * sqrtf(1.0 + powf((val - SigmaXCoefs[0]) / SigmaXCoefs[1], 2) + SigmaXCoefs[2] * powf((val - SigmaXCoefs[0]) / SigmaXCoefs[1], 3) + SigmaXCoefs[3] * powf((val - SigmaXCoefs[0]) / SigmaXCoefs[1], 4));
							ZXfunc = SigmaYCoefs[4] * sqrtf(1.0 + powf((val - SigmaYCoefs[0]) / SigmaYCoefs[1], 2) + SigmaYCoefs[2] * powf((val - SigmaYCoefs[0]) / SigmaYCoefs[1], 3) + SigmaYCoefs[3] * powf((val - SigmaYCoefs[0]) / SigmaYCoefs[1], 4));

							a = (powf(sX * 0.16 - ZXfunc, 2) + powf(sY * 0.16 - ZYfunc, 2)) / (powf(ZXfunc, 2) + powf(ZYfunc, 2) );

							if (a < Zerror)
							{
								Zerror = a;
								Zval = val;
							}
						}

					if (Zerror < 0.4)			
						pointsList[NbSegParams * index + 10] = Zval;
					else
						pointsList[NbSegParams * index + 10] = 999999;
				}

	                                  		
				pointsList[NbSegParams * index    ] = p[3]/sqrtf(2.)/*iteration*/;
				pointsList[NbSegParams * index + 1] = p[4]/sqrtf(2.);
				pointsList[NbSegParams * index + 8] = /*R*/pointsList[NbSegParams * index + 2];
				pointsList[NbSegParams * index + 9] = pointsList[NbSegParams * index + 3];
				pointsList[NbSegParams * index + 2] = angle;
				
					
				//pointsList[NbSegParams * index + 2] = Zerror;
				
				pointsList[NbSegParams * index + 3] = (p[1] - ref) + (float)x;
				pointsList[NbSegParams * index + 4] = (p[2] - ref) + (float)y;
				pointsList[NbSegParams * index + 5] = p[0];
				pointsList[NbSegParams * index + 6] = p[5];
				
				if (isnan(R_sq) || (isinf(R_sq)))
					pointsList[NbSegParams * index + 7] = (float)0;
				else
					pointsList[NbSegParams * index + 7] = R_sq;

			}
		}
	}
}

extern "C"
void cuda3DFittingGPU_7x7(unsigned short* data, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, float* pointsList, unsigned int nbPoints,
							 unsigned int iterNum, unsigned short sizeFit, float sigmaX, float sigmaY, float angle, unsigned short computeZ, float* SigmaXCoefs, float* SigmaYCoefs, float ThreeDThickness/*um*/, float ThreeDSampling/*um*/)
{
	printf("3D fitting\n");

	dim3  grid( 1, 1, 1);
    dim3  threads(1, 1, 1);


	// allocate memory
	grid.x = (int)ceil( (float)(nbPoints)/(float)FIT_BLOCK_SIZE);
	grid.y = 1;
	grid.z = 1;
    threads.x = FIT_BLOCK_SIZE;
	threads.y = 1;
	threads.z = 1;  // Par rapport aux contraintes de la mémoire partagée (taille limitée à 16 Ko pour chaque block)

								// '64' au lieu de '256'
	GPU3DFitting_7x7<<< grid, threads >>>(nbPoints, pointsList, data, nbStruct, width, height, iterNum, sizeFit, sigmaX, sigmaY, angle, computeZ, SigmaXCoefs, SigmaYCoefs, ThreeDThickness/*um*/, ThreeDSampling/*um*/);
}

//#endif // imFit




