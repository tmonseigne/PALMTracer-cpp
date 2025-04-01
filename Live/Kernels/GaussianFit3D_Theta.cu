#pragma once

#include <math_constants.h>
#include <cutil_inline.h>
#include "parameters.h"

// Levenberg fit --> simple precision (float)
#define dimFit 13
#define Npar 7 // I0, x0, y0, sigmaX, sigmaY, angle, offset
#define Npnt dimFit*dimFit
#define ref (dimFit/2)

__constant__ double epsilon_1 = 1e-6;
__constant__ double epsilon_2 = 1e-6;
__constant__ double epsilon_3 = 1e-5;
__constant__ double epsilon_4 = 1e-2;
__constant__ double lambda_0  = 1e-2;
__constant__ double lambda_UP_fac = 11;
__constant__ double lambda_DN_fac = 9;

__device__ float PSF(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float COS, float SIN)
{
	return I0 * exp( -pow(((x-x0)*COS-(y-y0)*SIN)/sigmaX, 2) - pow(((x-x0)*SIN+(y-y0)*COS)/sigmaY, 2) ) + offset;
 }

__device__ float deriv1I0(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float cons, float COS, float SIN)
{
	return cons/I0;
}

__device__ float deriv1x0(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float cons, float COS, float SIN)
{
	float alpha = -( (-2 * pow(COS, 2) * (x - x0) +2 * (y - y0) * COS * SIN ) / pow(sigmaX, 2) );
	float beta = -( (-2 * pow(SIN, 2) * (x - x0) -2 * (y - y0) * COS * SIN ) / pow(sigmaY, 2) );
	return cons * (alpha + beta);
}

__device__ float deriv1y0(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float cons, float COS, float SIN)
{
	float alpha = -( (2 * pow(SIN, 2) * (y - y0) -2 * (x - x0) * COS * SIN ) / pow(sigmaX, 2) );
	float beta = -( (-2 * pow(COS, 2) * (y - y0) -2 * (x - x0) * COS * SIN ) / pow(sigmaY, 2) );
	return cons * (alpha + beta);
}


__device__ float deriv1sigmaX(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float cons, float COS, float SIN)
{
	float alpha = pow((x-x0)*COS, 2) - 2 * (x-x0) * (y-y0) * COS * SIN + pow((y-y0)*SIN, 2); 
	return cons * ( (2*alpha) / pow(sigmaX, 3) );
}

__device__ float deriv1sigmaY(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float cons, float COS, float SIN)
{
	float alpha = pow((x-x0)*SIN, 2) + 2 * (x-x0) * (y-y0) * COS * SIN + pow((y-y0)*COS, 2); 
	return cons * ( (2*alpha) / pow(sigmaY, 3) );
}

__device__ float deriv1angle(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float angle, float offset, float cons, float COS, float SIN)
{
	float a = pow(x-x0, 2);
	float b = pow(y-y0, 2);
	float sX2 = pow(sigmaX, 2);
	float sY2 = pow(sigmaY, 2);
	return  cons * (2*sin(angle)*cos(angle)*((a-b)/sX2 + (b-a)/sY2 ) - pow(sin(angle), 2)*((2*(x-x0)*(y-y0))/sX2 + (-2*(x-x0)*(y-y0))/sY2) + pow(cos(angle), 2)*((2*(x-x0)*(y-y0))/sX2 + (-2*(x-x0)*(y-y0))/sY2) );
}

__device__ float lm_matx_(/*unsigned int* tIN,*/ float* p/*IN*/, float* y_dat/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/)
{
	// Computing Jacobian matrix (dydp) & y_hat, delta_y
	float delta_y[Npnt];
	float dydp[Npnt * Npar];
	float cons;
	unsigned int x, y;
	float Chi_sq = 0.;

	float COS = cos(p[5]);
	float SIN = sin(p[5]);

	for(int elt=0; elt < Npnt; elt++)
	{
		x = (elt%dimFit);
		y = (elt/dimFit);
		cons = p[0] * exp( -pow((((float)x - p[1])*COS-((float)y - p[2])*SIN)/p[3], 2) - pow((((float)x - p[1])*SIN+((float)y - p[2])*COS)/p[4], 2) );
		dydp[(0) + elt * Npar] = deriv1I0    ((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], p[6], cons, COS, SIN); // I0
		dydp[(1) + elt * Npar] = deriv1x0    ((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], p[6], cons, COS, SIN); // x0
		dydp[(2) + elt * Npar] = deriv1y0    ((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], p[6], cons, COS, SIN); // y0
		dydp[(3) + elt * Npar] = deriv1sigmaX((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], p[6], cons, COS, SIN); // sigmaX
		dydp[(4) + elt * Npar] = deriv1sigmaY((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], p[6], cons, COS, SIN); // sigmaY
		dydp[(5) + elt * Npar] = deriv1angle ((float)x, (float)y, p[0], p[1], p[2], p[3], p[4], p[5], p[6], cons, COS, SIN); // angle
		dydp[(6) + elt * Npar] = 1;																			                 // offset

		y_hat[elt] = cons + p[6]; // PSF
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

__device__ void inverse_matrix(float* MM, float* Minv)
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
	Mi1[6*Npar + 6] = 1;

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

		for(int k=0; k<Npar; k++)// mise z�ros de la colonne
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

__global__ void GPU3DFitting_Theta(unsigned int size, unsigned int MaxNbStruct, double* pointsList, float* data, unsigned short* nbStruct, unsigned short WIDTH, unsigned short HEIGHT, unsigned int iterNum,
							 unsigned short _size, float sigmaX, float sigmaY, float angle)
{
	unsigned int idThread = threadIdx.x + blockIdx.x * blockDim.x;
	unsigned short idImage = idThread / MaxNbStruct;
	unsigned int imageSize = WIDTH * HEIGHT;
	unsigned short _seg = idThread % MaxNbStruct;
	unsigned int index = idImage * MaxNbStruct + _seg;

	if ( ( idThread < size ) /*&& (_seg >= 1)*/ && (_seg <= nbStruct[idImage]) && ((int)pointsList[index * 10 + 3] != 0) )
	{
		if( (pointsList[10 * index + 1] >= (int)ref) && (pointsList[10 * index + 1] <= HEIGHT-(int)ref-1) &&  (pointsList[10 * index] >= (int)ref) && (pointsList[10 * index] <= WIDTH-(int)ref-1))
		{
			int x = (int)pointsList[10 * index];
			int y = (int)pointsList[10 * index + 1];

			float _x = (pointsList[10 * index]     - x) + ref;
			float _y = (pointsList[10 * index + 1] - y) + ref;			

			// for real data
			float p_min[Npar] = {0.    , _x - 1, _y - 1, 0.1              , 0.1              , -150., 0.    };
			float p[Npar]     = {2000. , _x    , _y    , sqrt(2.) * sigmaX, sqrt(2.) * sigmaY, angle, 700.  };
			float p_max[Npar] = {65535., _x + 1, _y + 1, 20.              , 20.              , +150., 65535.};

			float y_dat[Npnt];

			// Initialization of y_dat & t
			for(int jj = -dimFit/2; jj <= dimFit/2; jj++)
				for(int ii = -dimFit/2; ii <= dimFit/2; ii++)
					y_dat[(ii + dimFit/2) + (jj + dimFit/2) * dimFit] = data[idImage * imageSize + (x + ii) + (y + jj) * WIDTH];

			float lambda;

			float alpha[Npar * Npar];
			float beta[Npar];
			float y_hat[Npnt];
			
			float X2_old;
			float X2 = lm_matx_(p/*IN*/, y_dat/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/);

			float max__;
			float max_ = 0;
			for(int ii=0; ii < Npar; ii++)
					max_ = max(max_, abs(beta[ii]));

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

			float COS = cos(p[5]);
			float SIN = sin(p[5]);


			while ((stop == 0) && (iteration <= iterNum)) // Main Loop
			{
				iteration++;
				float delta_p[Npar] = {0., 0., 0., 0., 0., 0., 0.};

				// Copy alpha to inter
				for(int ii=0; ii < Npar*Npar; ii++)
					inter[ii] = alpha[ii];

				inter[0*Npar + 0] += lambda * alpha[0*Npar + 0];
				inter[1*Npar + 1] += lambda * alpha[1*Npar + 1];
				inter[2*Npar + 2] += lambda * alpha[2*Npar + 2];
				inter[3*Npar + 3] += lambda * alpha[3*Npar + 3];
				inter[4*Npar + 4] += lambda * alpha[4*Npar + 4];
				inter[5*Npar + 5] += lambda * alpha[5*Npar + 5];
				inter[6*Npar + 6] += lambda * alpha[6*Npar + 6];

				inverse_matrix(inter, inter1);

				m = 0;
				for(int ii=0; ii < Npar; ii++)
				{
					for(int cursor=0; cursor < Npar; cursor++)
						delta_p[ii] += inter1[ii*Npar + cursor] * beta[cursor];

					// Computing p_try
					p_try[ii] = min(max(p_min[ii], p[ii] + delta_p[ii]), p_max[ii]);

					// Computing m
					inter3 = lambda * delta_p[ii] + beta[ii];// lambda * delta_p + beta
					m += delta_p[ii] * inter3;
				}
				m *= 2;

				X2_try = 0;
				for(int ii=0; ii < Npnt; ii++)
				{
					delta_y = y_dat[ii] - PSF((float)(ii%dimFit), (float)(ii/dimFit), p_try[0], p_try[1], p_try[2], p_try[3], p_try[4], p_try[5], p_try[6], COS, SIN);
					X2_try += delta_y * delta_y;
				}				

				rho = (X2 - X2_try) / m;

				if(rho > epsilon_4) // it IS significantly better
				{
					X2_old = X2;
					for(int ii=0; ii < Npar; ii++)
						p[ii] = p_try[ii]; // accept p_try

					X2 = lm_matx_(p/*IN*/, y_dat/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/);

					lambda = max(lambda/lambda_DN_fac, 1.e-7);
				}
				else  // it IS NOT better
				{
					X2 = X2_old;  // do not accept a_try
					lambda = min(lambda*lambda_UP_fac, 1.e7);
				}

				max_ = 0;
				max__ = 0;
				for(int ii=0; ii < Npar; ii++)
				{
					max_ = max(max_, abs(delta_p[ii] / p[ii]));
					max__ = max(max__, abs(beta[ii]));
				}

				if( ( (max_ < epsilon_2) || (X2/Npnt < epsilon_3) || (max__ < epsilon_1) ) && (iteration > 2) ) 
					stop = 1;

				if(iteration == iterNum)
					stop = 1;
			}

			// Computing statistics
			float R_sq = 0.;
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
				t2 = sqrt(t2);
				t3 = sqrt(t3);

				// computing 'R_sq'
				R_sq = pow(t1 / (t2 * t3), 2);
			}

			pointsList[10 * index    ] = p[3]/sqrt(2.);
			pointsList[10 * index + 1] = p[4]/sqrt(2.);
			pointsList[10 * index + 8] = pointsList[10 * index + 2];
			pointsList[10 * index + 9] = pointsList[10 * index + 3];			
			pointsList[10 * index + 2] = p[5];
			pointsList[10 * index + 3] = (p[1] - ref) + (double)x;
			pointsList[10 * index + 4] = (p[2] - ref) + (double)y;
			pointsList[10 * index + 5] = p[0];
			pointsList[10 * index + 6] = p[6];
			pointsList[10 * index + 7] = R_sq;

		}
		else
		{
			double inter = pointsList[10 * index + 3];
			pointsList[10 * index + 3] = pointsList[10 * index + 0];
			pointsList[10 * index + 4] = pointsList[10 * index + 1];			
			pointsList[10 * index    ] = sigmaX;
			pointsList[10 * index + 1] = sigmaY;
			pointsList[10 * index + 8] = pointsList[10 * index + 2];
			pointsList[10 * index + 9] = inter/*pointsList[10 * index + 3]*/;			
			pointsList[10 * index + 2] = angle;

			pointsList[10 * index + 5] = 0.;
			pointsList[10 * index + 6] = 0.;
			pointsList[10 * index + 7] = -1.;
		}
	}
}

extern "C"
void cuda3DFittingGPU_Theta(float* data, unsigned int numImages, unsigned short width, unsigned short height, unsigned short* nbStruct, double* pointsList, unsigned int MaxNbStruct,
							 unsigned int iterNum, unsigned short size, float sigmaX, float sigmaY, float angle)
{
	printf("3D fitting\n");

	dim3  grid( numImages * (unsigned int)ceil( (float)width/(float)WIDTH_PORTION ), (unsigned int)ceil( (float)height/(float)HEIGHT_PORTION ), 1);
    dim3  threads(WIDTH_PORTION, 1, 1);


	// allocate memory
	grid.x = (int)ceil( (float)(numImages * MaxNbStruct)/(float)FIT_BLOCK_SIZE);
	grid.y = 1;
	grid.z = 1;
    threads.x = FIT_BLOCK_SIZE;
	threads.y = 1;
	threads.z = 1;  // Par rapport aux contraintes de la m�moire partag�e (taille limit�e � 16 Ko pour chaque block)
											 // '64' au lieu de '256'
	GPU3DFitting_Theta<<< grid, threads/*, FIT_BLOCK_SIZE * Npar * Npar * sizeof(float) * FIT_BLOCK_SIZE*/>>>(numImages * MaxNbStruct, MaxNbStruct, pointsList, data, nbStruct, width, height, iterNum, size, sigmaX, sigmaY, angle);
}






