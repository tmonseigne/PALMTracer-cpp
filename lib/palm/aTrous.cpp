#include "aTrous.h"


ATrous::ATrous()
{
	temp1 = new FltArray();
	temp2 = new FltArray();
	inputV0 = new FltArray();
	VN = new FltArray();
	VNplus1 = new FltArray();
	wDecomp = new FltArray();
}


ATrous::~ATrous()
{
	delete temp1;
	delete temp2;
	delete inputV0;
	delete VN;
	delete VNplus1;
	delete wDecomp;
}

/*
Gives the noise degeneration for the used wavelet in the a trous algorithm.
For example 0.3548 is the standard deviation in the first for a input image 
full filled of gaussian noise with a std=1.
*/
double ATrous::stdRef(long dim, unsigned long realProj, unsigned long plan)
{
	static double _1DRealProj[11] = {0.368862, 0.243017, 0.172412, 0.121539, 0.085503, 0.060303, 0.043978, 0.030842, 0.020708, 0.015012, 0.011179};
	static double _1D[11] = {0.723077, 0.285706, 0.179198, 0.122700, 0.085704, 0.060338, 0.043984, 0.030842, 0.020721, 0.015015, 0.011192};
	static double _2DRealProj[4] = {0.3548, 0.1669, 0.0816, 0.0590};
	static double _2D[4] = {0.8868, 0.2046, 0.0862, 0.0598};
	static double _3DRealProj[4] = {0.2900, 0.0930, 0.0340, 0.0160};
	static double _3D[4] = {0.9500, 0.1200, 0.0360, 0.0180};
	switch (dim) {
		case 1: {
			if (realProj) { return _1DRealProj[plan]; }
			return _1D[plan];
		}
		case 2: {
			if (realProj) { return _2DRealProj[plan]; }
			return _2D[plan];
		}
		case 3: {
			if (realProj) { return _3DRealProj[plan]; }
			return _3D[plan];
		}
		default: return 0.;
	}
}
