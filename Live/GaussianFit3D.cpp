#include "stdafx.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>

#include "GaussianFit3D.h"
// CPU version of gaussian fit 3D
//#define dimFit 13
#define Npar 6 // I0, x0, y0, sigmaX, sigmaY, offset
//#define Npnt dimFit*dimFit
//#define ANGLE 0.f
//#define ref (dimFit/2)


float GaussianFit3D::PSF(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset)
{
	return I0 * exp( -pow(((x-x0)*cos(ANGLE)-(y-y0)*sin(ANGLE))/sigmaX, 2) - pow(((x-x0)*sin(ANGLE)+(y-y0)*cos(ANGLE))/sigmaY, 2) ) + offset;
}

float GaussianFit3D::deriv1I0(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset)
{
	return exp( -pow(((x-x0)*cos(ANGLE)-(y-y0)*sin(ANGLE))/sigmaX, 2) - pow(((x-x0)*sin(ANGLE)+(y-y0)*cos(ANGLE))/sigmaY, 2) );
}

//float GaussianFit3D::deriv1x0(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset)
//{
//	float p = I0 * exp( -pow(((x-x0)*cos(ANGLE)-(y-y0)*sin(ANGLE))/sigmaX, 2) - pow(((x-x0)*sin(ANGLE)+(y-y0)*cos(ANGLE))/sigmaY, 2) );
//	float alpha = -( (2 * pow(cos(ANGLE), 2) * (x - x0) -2 * (y - y0) * cos(ANGLE) * sin(ANGLE) ) / pow(sigmaX, 2) );
//	float beta = -( (2 * pow(sin(ANGLE), 2) * (x - x0) +2 * (y - y0) * cos(ANGLE) * sin(ANGLE) ) / pow(sigmaY, 2) );
//	return p * (alpha + beta);
//}
//
//float GaussianFit3D::deriv1y0(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset)
//{
//	float p = I0 * exp( -pow(((x-x0)*cos(ANGLE)-(y-y0)*sin(ANGLE))/sigmaX, 2) - pow(((x-x0)*sin(ANGLE)+(y-y0)*cos(ANGLE))/sigmaY, 2) );
//	float alpha = -( (2 * pow(sin(ANGLE), 2) * (y - y0) -2 * (x - x0) * cos(ANGLE) * sin(ANGLE) ) / pow(sigmaX, 2) );
//	float beta = -( (2 * pow(cos(ANGLE), 2) * (y - y0) +2 * (x - x0) * cos(ANGLE) * sin(ANGLE) ) / pow(sigmaY, 2) );
//	return p * (alpha + beta);
//}

float GaussianFit3D::deriv1x0(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset)
{
	float p = I0 * exp( -pow(((x-x0)*cos(ANGLE)-(y-y0)*sin(ANGLE))/sigmaX, 2) - pow(((x-x0)*sin(ANGLE)+(y-y0)*cos(ANGLE))/sigmaY, 2) );
	float alpha = -( (-2 * pow(cos(ANGLE), 2) * (x - x0) +2 * (y - y0) * cos(ANGLE) * sin(ANGLE) ) / pow(sigmaX, 2) );
	float beta = -( (-2 * pow(sin(ANGLE), 2) * (x - x0) -2 * (y - y0) * cos(ANGLE) * sin(ANGLE) ) / pow(sigmaY, 2) );
	return p * (alpha + beta);
}

float GaussianFit3D::deriv1y0(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset)
{
	float p = I0 * exp( -pow(((x-x0)*cos(ANGLE)-(y-y0)*sin(ANGLE))/sigmaX, 2) - pow(((x-x0)*sin(ANGLE)+(y-y0)*cos(ANGLE))/sigmaY, 2) );
	float alpha = -( (2 * pow(sin(ANGLE), 2) * (y - y0) -2 * (x - x0) * cos(ANGLE) * sin(ANGLE) ) / pow(sigmaX, 2) );
	float beta = -( (-2 * pow(cos(ANGLE), 2) * (y - y0) -2 * (x - x0) * cos(ANGLE) * sin(ANGLE) ) / pow(sigmaY, 2) );
	return p * (alpha + beta);
}

float GaussianFit3D::deriv1sigmaX(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset)
{
	float p = I0 * exp( -pow(((x-x0)*cos(ANGLE)-(y-y0)*sin(ANGLE))/sigmaX, 2) - pow(((x-x0)*sin(ANGLE)+(y-y0)*cos(ANGLE))/sigmaY, 2) );
	float alpha = pow((x-x0)*cos(ANGLE), 2) - 2 * (x-x0) * (y-y0) * cos(ANGLE) * sin(ANGLE) + pow((y-y0)*sin(ANGLE), 2); 
	return p * ( (2*alpha) / pow(sigmaX, 3) );
}

float GaussianFit3D::deriv1sigmaY(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset)
{
	float p = I0 * exp( -pow(((x-x0)*cos(ANGLE)-(y-y0)*sin(ANGLE))/sigmaX, 2) - pow(((x-x0)*sin(ANGLE)+(y-y0)*cos(ANGLE))/sigmaY, 2) );
	float alpha = pow((x-x0)*sin(ANGLE), 2) + 2 * (x-x0) * (y-y0) * cos(ANGLE) * sin(ANGLE) + pow((y-y0)*cos(ANGLE), 2); 
	return p * ( (2*alpha) / pow(sigmaY, 3) );
}

float GaussianFit3D::lm_matx(unsigned int* t/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight_sq/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/)
{
	//float delta_y[Npnt];	
	float* delta_y = (float*)malloc(Npnt * sizeof(float));
	for(int elt = 0; elt < Npnt; elt++)
	{
		y_hat[elt] = PSF((float)t[elt*2], (float)t[elt*2+1], p[0], p[1], p[2], p[3], p[4], p[5]);
		delta_y[elt] = y_dat[elt] - y_hat[elt];
	}

	// Computing Jacobian matrix (dydp)
	//float dydp[Npnt * Npar];
	float* dydp = (float*)malloc(Npnt * Npar * sizeof(float));
	for(int elt=0; elt < Npnt; elt++)
	{
		dydp[(0) + elt * Npar] = deriv1I0   ((float)t[elt*2], (float)t[elt*2+1], p[0], p[1], p[2], p[3], p[4], p[5]);  // I0
		dydp[(1) + elt * Npar] = deriv1x0   ((float)t[elt*2], (float)t[elt*2+1], p[0], p[1], p[2], p[3], p[4], p[5]);  // x0
		dydp[(2) + elt * Npar] = deriv1y0   ((float)t[elt*2], (float)t[elt*2+1], p[0], p[1], p[2], p[3], p[4], p[5]);  // y0
		dydp[(3) + elt * Npar] = deriv1sigmaX((float)t[elt*2], (float)t[elt*2+1], p[0], p[1], p[2], p[3], p[4], p[5]); // sigmaX
		dydp[(4) + elt * Npar] = deriv1sigmaY((float)t[elt*2], (float)t[elt*2+1], p[0], p[1], p[2], p[3], p[4], p[5]); // sigmaY
		dydp[(5) + elt * Npar] = 1;																					   // offset
	}

	//float inter[Npar * Npnt];
	float* inter = (float*)malloc(Npar * Npnt * sizeof(float));
	for(int elt=0; elt < Npnt; elt++)
		for(int i=0; i < Npar; i++)
			inter[elt*Npar + i] = weight_sq[elt] * dydp[elt*Npar + i];

	// Computing alpha (A vérifier)
	for(int j=0; j < Npar; j++)
		for(int i=0; i < Npar; i++)
		{
			alpha[Npar * j + i] = 0;
			for(int cursor=0; cursor < Npnt; cursor++)
				alpha[Npar * j + i] += dydp[cursor*Npar + j] * inter[cursor*Npar + i];  // dydp' * inter : dydp'[i*Npnt + cursor]
		}

	// Computing beta
	//float inter1[Npnt];
	float* inter1 = (float*)malloc(Npnt * sizeof(float));
	for(int elt=0; elt < Npnt; elt++)
		inter1[elt] = weight_sq[elt] * delta_y[elt];

	for(int j=0; j < Npar; j++)
	{
		beta[j] = 0;
		for(int cursor=0; cursor < Npnt; cursor++)
			beta[j] += dydp[cursor*Npar + j] * inter1[cursor];  // dydp' * inter : dydp'[j*Npnt + cursor]
	}

	// Computing Chi_sq
	float Chi_sq = 0.;
	for(int elt=0; elt < Npnt; elt++)
		Chi_sq += delta_y[elt] * inter1[elt];

	// cleanup memory
	free(delta_y);
	free(dydp);
	free(inter);
	free(inter1);

	return Chi_sq;
}

void GaussianFit3D::inverse(float* MM, float* Minv)
{
	float M1[Npar * Npar];

	// Copy MM to M1
	for(int i=0; i< Npar*Npar; i++)
		M1[i] = MM[i];

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
		}

		// Copy M1 to MM
		for(int index=0; index < Npar*Npar; index++)
			MM[index] = M1[index];

		// Copy Mi1 to Minv
		for(int index=0; index < Npar*Npar; index++)
			Minv[index] = Mi1[index];

		for(int k=0; k<Npar; k++)// mise zéros de la colonne
			if (k != i)
				for(int j=0; j<Npar; j++)
				{
					M1[k*Npar+j] = MM[k*Npar+j] - MM[i*Npar+j]*MM[k*Npar+i];
					Mi1[k*Npar+j] = Minv[k*Npar+j] - Minv[i*Npar+j]*MM[k*Npar+i];
				}
		
		// Copy M1 to MM
		for(int index=0; index < Npar*Npar; index++)
			MM[index] = M1[index];

		// Copy Mi1 to Minv
		for(int index=0; index < Npar*Npar; index++)
			Minv[index] = Mi1[index];
	}
}


void GaussianFit3D::Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float &x0, float &y0, float &I0, float &sigmaX, float &sigmaY, float &angle, float &offset, float &RSS, int iterNum, unsigned short size)
{
	dimFit = size;
	Npnt = dimFit*dimFit;
	ref = (dimFit/2);
	ANGLE = angle;

	float X2;
	float time_initial;
	float time_final;
	float time = 0;

	if((x0 >= (int)ref) && (x0 <= sizeX-(int)ref-1) &&  (y0 >= (int)ref) && (y0 <= sizeY-(int)ref-1))
	{
		// for simulation
		//float p_min[5] = {0  , pointsList[4 * index]-1, pointsList[4 * index + 1]-1, 0.5, 300};
		//float p[5]     = {250, pointsList[4 * index]  , pointsList[4 * index + 1]  , 1. , 1000};
		//float p_max[5] = {350, pointsList[4 * index]+1, pointsList[4 * index + 1]+1, 1.5, 1700};

		// for real data
		float p_min[Npar] = {0.   , x0-1, y0-1, 0.1    , 0.1    , 0.   };
		float p[Npar]     = {2000., x0  , y0  , sqrt(2.) * sigmaX, sqrt(2.) * sigmaY, 700. };
		float p_max[Npar] = {65535., x0+1, y0+1, 20.   , 20.   , 65535.};

		//float y_dat[Npnt];
		//float weight[Npnt];
		//float weight_sq[Npnt];
		//unsigned int t[2 * Npnt];

		float* y_dat = (float*)malloc(Npnt * sizeof(float));
		float* weight = (float*)malloc(Npnt * sizeof(float));
		float* weight_sq = (float*)malloc(Npnt * sizeof(float));
		unsigned int* t = (unsigned int*)malloc(2 * Npnt * sizeof(unsigned int));

		// Initialization of y_dat & t (A vérifier)
		for(int jj = -dimFit/2; jj <= dimFit/2; jj++)
			for(int ii = -dimFit/2; ii <= dimFit/2; ii++)
			{
				y_dat[(ii + dimFit/2) + (jj + dimFit/2) * dimFit] = (float)image[((int)x0 + ii) + ((int)y0 + jj) * sizeX];
				t[((ii + dimFit/2) + (jj + dimFit/2) * dimFit) * 2]     = (int)x0 + ii;
				t[((ii + dimFit/2) + (jj + dimFit/2) * dimFit) * 2 + 1] = (int)y0 + jj;
			}

		for(int ii = 0; ii < Npnt; ii++) 
		{
			weight[ii] = 1/*sqrt((Npnt-Npar+1) / (y_dat[ii]*y_dat[ii]))*/;
			weight_sq[ii] = pow(weight[ii], 2);
		}

		float epsilon_1 = 1e-6;
		float epsilon_2 = 1e-6;
		float epsilon_3 = 1e-5;
		float epsilon_4 = 1e-2;
		float lambda_0  = 1e-2;
		float lambda_UP_fac = 11;
		float lambda_DN_fac = 9;

		float lambda;

		float alpha[Npar * Npar];
		float beta[Npar];
		/*float y_hat[Npnt];*/
		float* y_hat = (float*)malloc(Npnt * sizeof(float));
		
		//float X2;
		float X2_old;
		X2 = lm_matx(t/*IN*/, p/*IN*/, y_dat/*IN*/, weight_sq/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/);

		float max_ = 0;
		for(int ii=0; ii < Npar; ii++)
				max_ = max(max_, abs(beta[ii]));

		unsigned int stop = 0;
		if(max_ < epsilon_1) 
			stop = 1;

		lambda = lambda_0;
		X2_old = X2;

		float p_old[Npar] = {2*p[0], 2*p[1], 2*p[2], 2*p[3], 2*p[4], 2*p[5]};
		unsigned short iteration = 0;

		float inter [Npar * Npar];
		float inter1[Npar * Npar];
		/*float inter2[Npnt];*/
		float* inter2 = (float*)malloc(Npnt * sizeof(float));
		float inter3[Npar];

		float p_try  [Npar];
		/*float delta_y[Npnt];*/
		float* delta_y = (float*)malloc(Npnt * sizeof(float));
		float X2_try;
		float m;
		float rho;

		while ( (stop == 0) && (iteration <= iterNum) ) // Main Loop
		{
			iteration++;
			float delta_p[Npar] = {0, 0, 0, 0, 0, 0};

			// Copy alpha to inter
			for(int ii=0; ii < Npar*Npar; ii++)
				inter[ii] = alpha[ii];

			inter[0*Npar + 0] += lambda * alpha[0*Npar + 0];
			inter[1*Npar + 1] += lambda * alpha[1*Npar + 1];
			inter[2*Npar + 2] += lambda * alpha[2*Npar + 2];
			inter[3*Npar + 3] += lambda * alpha[3*Npar + 3];
			inter[4*Npar + 4] += lambda * alpha[4*Npar + 4];
			inter[5*Npar + 5] += lambda * alpha[5*Npar + 5];

			inverse(inter, inter1);

			for(int ii=0; ii < Npar; ii++)
			{
				for(int cursor=0; cursor < Npar; cursor++)
					delta_p[ii] += inter1[ii*Npar + cursor] * beta[cursor];
			}

			for(int ii=0; ii < Npar; ii++)
			{
				p_try[ii] = p[ii] + delta_p[ii];
				p_try[ii] = min(max(p_min[ii], p_try[ii]), p_max[ii]);
			}

			for(int ii=0; ii < Npnt; ii++)
				delta_y[ii] = y_dat[ii] - PSF((float)t[ii*2], (float)t[ii*2+1], p_try[0], p_try[1], p_try[2], p_try[3], p_try[4], p_try[5]);

			for(int ii=0; ii < Npnt; ii++) // delta_y .* weight_sq
				inter2[ii] = delta_y[ii] * weight_sq[ii];

			X2_try = 0;
			for(int ii=0; ii < Npnt; ii++)
				X2_try += delta_y[ii] * inter2[ii];

			for(int ii=0; ii < Npar; ii++) // lambda * delta_p + beta
				inter3[ii] = lambda * delta_p[ii] + beta[ii];

			m = 0;
			for(int ii=0; ii < Npar; ii++)
				m += delta_p[ii] * inter3[ii];
			m *= 2;

			rho = (X2 - X2_try) / m;

			if(rho > epsilon_4) // it IS significantly better
			{
				X2_old = X2;
				for(int ii=0; ii < Npar; ii++)
				{
					p_old[ii] = p[ii];
					p[ii] = p_try[ii]; // accept p_try
				}

				X2 = lm_matx(t/*IN*/, p/*IN*/, y_dat/*IN*/, weight_sq/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/);

				lambda = max(lambda/lambda_DN_fac,1.e-7);
			}
			else  // it IS NOT better
			{
				X2 = X2_old;  // do not accept a_try
				lambda = min(lambda*lambda_UP_fac,1.e7);
			}

			max_ = 0;
			for(int ii=0; ii < Npar; ii++)
					max_ = max(max_, abs(delta_p[ii] / p[ii]));

			if( (max_ < epsilon_2) && (iteration > 2) ) 
				stop = 1;

			if( (X2/Npnt < epsilon_3) && (iteration > 2 ) )
				stop = 1;

			max_ = 0;
			for(int ii=0; ii < Npar; ii++)
					max_ = max(max_, abs(beta[ii]));

			if( (max_ < epsilon_1) && (iteration > 2) ) 
				stop = 1;

			if(iteration == iterNum)
				stop = 1;
		}

		// Computing standard error
		;
		//for(int ii = 0; ii < Npnt; ii++) 
		//	weight_sq[ii] = (Npnt-Npar+1)/(delta_y[ii]*delta_y[ii]);

		//X2 = lm_matx_(t/*IN*/, p/*IN*/, y_dat/*IN*/, weight_sq/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/);

		//float covar[Npar * Npar];
		//float sigma_p[Npar];
		//float corr [Npar * Npar];

		//X2 /= 2;
		//
		//inverse_(alpha, covar);

		//// standard error of parameters
		//sigma_p[0] = sqrt(covar[0*Npar+0]);
		//sigma_p[1] = sqrt(covar[1*Npar+1]);
		//sigma_p[2] = sqrt(covar[2*Npar+2]);
		//sigma_p[3] = sqrt(covar[3*Npar+3]);
		//sigma_p[4] = sqrt(covar[4*Npar+4]);	
		//sigma_p[5] = sqrt(covar[5*Npar+5]);

		//float q[Npar*Npar];
		//for(int ii = 0; ii < Npar; ii++) 
		//	for(int jj = 0; jj < Npar; jj++)
		//		q[ii + jj * Npar] = sigma_p[ii] * sigma_p[jj];

		//for(int ii = 0; ii < Npar*Npar; ii++)
		//	corr[ii] = covar[ii] / q[ii]; 

		// computing "y_dat" mean & "y_hat" mean
		float y_dat_mean = 0.;
		float y_hat_mean = 0.;
		for(int ii=0; ii < Npnt; ii++)
		{
			y_dat_mean += y_dat[ii];
			y_hat_mean += y_hat[ii];
		}
		y_dat_mean /= Npnt;
		y_hat_mean /= Npnt;

		// computing 't1' & 't2' & 't3'
		float t1 = 0.;
		float t2 = 0.;	
		float t3 = 0.;
		float diff_hat, diff_dat;

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
		float R_sq = pow(t1 / (t2 * t3), 2);

		x0 = p[1];
		y0 = p[2];
		I0 = p[0];
		sigmaX = p[3]/sqrt(2.);
		sigmaY = p[4]/sqrt(2.);
		offset = p[5];
		RSS = R_sq/*iteration*/;	

		// cleanup memory
		free(y_dat);
		free(weight);
		free(weight_sq);
		free(t);
		free(y_hat);
		free(inter2);
		free(delta_y);
	}
}
