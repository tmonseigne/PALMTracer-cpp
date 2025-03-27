#include "GaussianFit.h"


class GaussianFit2D_Sigma : GaussianFit
{
private:
	float PSF(float x, float y, float I0, float x0, float y0, float sigma, float offset);
	float deriv1I0(float x, float y, float I0, float x0, float y0, float sigma);
	float deriv1x0(float x, float y, float I0, float x0, float y0, float sigma);
	float deriv1y0(float x, float y, float I0, float x0, float y0, float sigma);
	float deriv1sigma(float x, float y, float I0, float x0, float y0, float sigma);
	float lm_matx(unsigned int* t/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight_sq/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/) override;
	void inverse(float* MM, float* Minv) override;

public :
	unsigned short dimFit;
	unsigned short Npnt;
	float ref;

	GaussianFit2D_Sigma() : dimFit(0), Npnt(0), ref(0) {}

	void Run(unsigned short* image, unsigned long sizeX, unsigned long sizeY, float& x0, float& y0, float& I0, float& sigmaX, float& offset, float& RSS, int iterNum, unsigned short size);
};
