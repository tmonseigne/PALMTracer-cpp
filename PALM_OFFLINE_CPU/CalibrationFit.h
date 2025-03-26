#ifndef _CalibrationFit_H_
#define _CalibrationFit_H_

class CalibrationFit
{
private:
	float deriv1C0(float Z, float C0, float C1, float C2, float C3, float C4);
	float deriv1C1(float Z, float C0, float C1, float C2, float C3, float C4);
	float deriv1C2(float Z, float C0, float C1, float C2, float C3, float C4);
	float deriv1C3(float Z, float C0, float C1, float C2, float C3, float C4);
	float deriv1C4(float Z, float C0, float C1, float C2, float C3, float C4);
	float lm_matx(float* Z/*IN*/, float* p/*IN*/, float* y_dat/*IN*/, float* weight_sq/*IN*/, float* alpha/*OUT*/, float* beta/*OUT*/, float* y_hat/*OUT*/);
	void inverse(float* MM, float* Minv);

public :
	float sigmaFunc(float Z, float C0, float C1, float C2, float C3, float C4);
	float Run(float* sigma, float* Z, unsigned short pointsNbr, unsigned short iterNum, float &C0, float &C1, float &C2, float &C3, float &C4);
};

#endif // #ifndef _CalibrationFit_H_