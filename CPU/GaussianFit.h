
#ifndef _GaussianFit_H_
#define _GaussianFit_H_

class GaussianFit
{
	virtual float lm_matx(unsigned int* t/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight_sq/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/) = 0;
	virtual void inverse(float* MM, float* Minv){}
};

#endif // #ifndef _GaussianFit_H_
