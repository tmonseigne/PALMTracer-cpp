
#include "GaussianFit.h"


class GaussianFit3D : GaussianFit
{
private:
	float PSF(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset);
	float deriv1I0(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset);
	float deriv1x0(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset);
	float deriv1y0(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset);
	float deriv1sigmaX(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset);
	float deriv1sigmaY(float x, float y, float I0, float x0, float y0, float sigmaX, float sigmaY, float offset);
	float lm_matx(unsigned int* t/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight_sq/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/);
	void inverse(float* MM, float* Minv);

public:
	unsigned short dimFit;
	unsigned short Npnt;
	float ref;
	float ANGLE;

	void Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float &x0, float &y0, float &I0, float &sigmaX, float &sigmaY, float &angle, float &offset, float &RSS, int iterNum, unsigned short size);
	GaussianFit3D()
	{
		dimFit = 0;
		Npnt = 0;
		ref = 0.;
		ANGLE = 0.;
	}
};