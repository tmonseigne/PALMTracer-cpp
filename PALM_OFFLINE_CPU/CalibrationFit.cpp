#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctime>
#include <stdlib.h>
#include <float.h>
#include "CalibrationFit.h"

#define Npar 5 // C0, C1, C2, C3, C4
int Npnt; 

float CalibrationFit::sigmaFunc(float Z, float C0, float C1, float C2, float C3, float C4)
{
	return C4 * sqrt(1 + pow((Z-C0)/C1, 2) + C2 * pow((Z-C0)/C1, 3) + C3 * pow((Z-C0)/C1, 4));
}

float CalibrationFit::deriv1C0(float Z, float C0, float C1, float C2, float C3, float C4)
{
	float cons = 2 * sqrt(1 + pow((Z-C0)/C1, 2) + C2 * pow((Z-C0)/C1, 3) + C3 * pow((Z-C0)/C1, 4));
	return -(C4/pow(C1, 2)) * ((2*(Z-C0) + 3*(C2/C1)*pow(Z-C0, 2) + 4*(C3/pow(C1, 2))*pow(Z-C0, 3))/(cons));
}

float CalibrationFit::deriv1C1(float Z, float C0, float C1, float C2, float C3, float C4)
{
	float cons = 2 * sqrt(1 + pow((Z-C0)/C1, 2) + C2 * pow((Z-C0)/C1, 3) + C3 * pow((Z-C0)/C1, 4));
	return -(C4/pow(C1, 3)) * ((2*pow(Z-C0, 2) + 3*(C2/C1)*pow(Z-C0, 3) + 4*(C3/pow(C1, 2))*pow(Z-C0, 4))/(cons));
}

float CalibrationFit::deriv1C2(float Z, float C0, float C1, float C2, float C3, float C4)
{
	float cons = 2 * sqrt(1 + pow((Z-C0)/C1, 2) + C2 * pow((Z-C0)/C1, 3) + C3 * pow((Z-C0)/C1, 4));
	return (C4/pow(C1, 3)) * (pow(Z-C0, 3)/(cons));
}

float CalibrationFit::deriv1C3(float Z, float C0, float C1, float C2, float C3, float C4)
{
	float cons = 2 * sqrt(1 + pow((Z-C0)/C1, 2) + C2 * pow((Z-C0)/C1, 3) + C3 * pow((Z-C0)/C1, 4));
	return (C4/pow(C1, 4)) * (pow(Z-C0, 4)/(cons));
}

float CalibrationFit::deriv1C4(float Z, float C0, float C1, float C2, float C3, float C4)
{
	return sqrt(1 + pow((Z-C0)/C1, 2) + C2 * pow((Z-C0)/C1, 3) + C3 * pow((Z-C0)/C1, 4));
}

float CalibrationFit::lm_matx(float* Z/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight_sq/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/)
{
	/*float delta_y[Npnt];*/
	float* delta_y = (float*)malloc(Npnt * sizeof(float));

	//FILE* f = fopen("C:\\FitTrace.txt", "w+");
	for(int elt = 0; elt < Npnt; elt++)
	{
		y_hat[elt] = sigmaFunc(Z[elt], p[0], p[1], p[2], p[3], p[4]);
		//fprintf(f, "%f\n", y_hat[elt]);
		delta_y[elt] = y_dat[elt] - y_hat[elt];
	}

	//fclose(f);

	// Computing Jacobian matrix (dydp)
	//float dydp[Npnt * Npar];
	float* dydp = (float*)malloc(Npnt * Npar * sizeof(float));
	for(int elt=0; elt < Npnt; elt++)
	{
		dydp[(0) + elt * Npar] = deriv1C0(Z[elt], p[0], p[1], p[2], p[3], p[4]); // C0
		dydp[(1) + elt * Npar] = deriv1C1(Z[elt], p[0], p[1], p[2], p[3], p[4]); // C1
		dydp[(2) + elt * Npar] = deriv1C2(Z[elt], p[0], p[1], p[2], p[3], p[4]); // C2
		dydp[(3) + elt * Npar] = deriv1C3(Z[elt], p[0], p[1], p[2], p[3], p[4]); // C3
		dydp[(4) + elt * Npar] = deriv1C4(Z[elt], p[0], p[1], p[2], p[3], p[4]); // C4
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

void CalibrationFit::inverse(float* MM, float* Minv)
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


float CalibrationFit::Run(float* sigma, float* Z, unsigned short pointsNbr, unsigned short iterNum, float &C0, float &C1, float &C2, float &C3, float &C4)
{
	Npnt = pointsNbr;

	float X2;

	float zMin = +999999;
	float zMax = -999999;
	for(int i=0; i <pointsNbr; i++)
	{
		if (zMin > Z[i])
			zMin = Z[i];

		if (zMax < Z[i])
			zMax = Z[i];
	}

	zMax += 2.0;
	zMin -= 2.0;

	float p_min[Npar] = {zMin, -30., -30., -30., -30.};
	float p[Npar]     = {C0, C1, C2, C3, C4};
	float p_max[Npar] = {zMax, 30., 30., 30., 30.};

	float* weight = (float*)malloc(Npnt * sizeof(float));
	float* weight_sq = (float*)malloc(Npnt * sizeof(float));

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
	//float y_hat[Npnt];
	float* y_hat = (float*)malloc(Npnt * sizeof(float));

	//float X2;
	float X2_old;
	X2 = lm_matx(Z/*IN*/, p/*IN*/, sigma/*IN*/, weight_sq/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/);

	float max_ = 0;
	for(int ii=0; ii < Npar; ii++)
			max_ = max(max_, abs(beta[ii]));

	unsigned int stop = 0;
	if(max_ < epsilon_1) 
		stop = 1;

	lambda = lambda_0;
	X2_old = X2;

	float p_old[Npar] = {2*p[0], 2*p[1], 2*p[2], 2*p[3], 2*p[4]};
	unsigned short iteration = 0;

	float inter [Npar * Npar];
	float inter1[Npar * Npar];
	//float inter2[Npnt];
	float* inter2 = (float*)malloc(Npnt * sizeof(float));
	float inter3[Npar];

	float p_try  [Npar];
	//float delta_y[Npnt];
	float* delta_y = (float*)malloc(Npnt * sizeof(float));
	float X2_try;
	float m;
	float rho;

	while ((stop == 0) && (iteration < iterNum)) // Main Loop
	{
		iteration++;
		float delta_p[Npar] = {0, 0, 0, 0, 0};

		// Copy alpha to inter
		for(int ii=0; ii < Npar*Npar; ii++)
			inter[ii] = alpha[ii];

		inter[0*Npar + 0] += lambda * alpha[0*Npar + 0];
		inter[1*Npar + 1] += lambda * alpha[1*Npar + 1];
		inter[2*Npar + 2] += lambda * alpha[2*Npar + 2];
		inter[3*Npar + 3] += lambda * alpha[3*Npar + 3];
		inter[4*Npar + 4] += lambda * alpha[4*Npar + 4];

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
			delta_y[ii] = sigma[ii] - sigmaFunc(Z[ii], p_try[0], p_try[1], p_try[2], p_try[3], p_try[4]);

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

			X2 = lm_matx(Z/*IN*/, p/*IN*/, sigma/*IN*/, weight_sq/*IN*/, alpha/*OUT*/, beta/*OUT*/, y_hat/*OUT*/);

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

	// computing "y_dat" mean & "y_hat" mean
	float y_dat_mean = 0.;
	float y_hat_mean = 0.;
	for(int ii=0; ii < Npnt; ii++)
	{
		y_dat_mean += sigma[ii];
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
		diff_dat = sigma[ii] - y_dat_mean;
		diff_hat = y_hat[ii] - y_hat_mean;
		t1 += (diff_dat) * (diff_hat);
		t2 += (diff_dat) * (diff_dat);
		t3 += (diff_hat) * (diff_hat);
	}
	t2 = sqrt(t2);
	t3 = sqrt(t3);

	// computing 'R_sq'
	float R_sq = pow(t1 / (t2 * t3), 2);

	C0 = p[0];
	C1 = p[1];
	C2 = p[2];
	C3 = p[3];
	C4 = p[4];

	float RSS;
	if (_isnan((double)R_sq))
		RSS = (float)0;
	else
		RSS = R_sq;
	

	//FILE* f = fopen("C:\\Fit.txt", "w");
	//fprintf(f, "%d\n", iteration);
	//fclose(f);

	// cleanup memory
	free(weight);
	free(weight_sq);
	free(y_hat);
	free(inter2);
	free(delta_y);

	return RSS;
}


