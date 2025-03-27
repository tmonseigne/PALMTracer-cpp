#include "ATrous.h"


#define	OFFSET(i,j)	offset((i),(j),sizeDataY)
#define	OFFSETTEST(i,j)	offsetTest((i),(j),sizeDataX,sizeDataY)


int ATrous::offsetTest(int i, int j, int maxSizeI, int maxSizeJ)
{
	//Same as offset but if i or j are out of bounds, it changes them to be in bounds.
	if (i < 0) { i = 0; }
	if (i >= maxSizeI) { i = maxSizeI - 1; }
	if (j < 0) { j = 0; }
	if (j >= maxSizeJ) { j = maxSizeJ - 1; }
	return (i * maxSizeJ + j);
}


int ATrous::offset(int i, int j, int maxSizeJ)
{
	//Compute the offset translate   tab[i][j] to tab[i*maxSizeJ+j]
	return (i * maxSizeJ + j);
}


int ATrous::testSizes2D(FltArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long& nbPlane)
{
	//Test inputData, outputData, planesToUse, nbPlane
	unsigned long i;

	if (!inputData.isThereData()) { return (W_INIT_ERR); }
	if (!outputData.isThereData()) { return (W_INIT_ERR); }
	if (!planesToUse.isThereData()) { return (W_INIT_ERR); }
	if ((inputData.dim() + 1) != outputData.dim()) { return (W_DIM_ERR); }
	if ((inputData.size(0) != outputData.size(0)) || (inputData.size(1) != outputData.size(1))) { return (W_SIZE_ERR); }
	if (planesToUse.dim() != 1) { return (W_DIM_ERR); }
	for (i = outputData.size(2); i < (unsigned long)planesToUse.size(0); i++) { planesToUse(i - 1) = false; }


	long sizeWavelet;


	if ((inputData.size(0) < 9) || (inputData.size(1) < 9)) { return (W_SIZE_ERR); }
	if (nbPlane < 2) { return (W_SIZE_ERR); }
	sizeWavelet = twoToThe(nbPlane + 1);
	while ((inputData.size(0) <= sizeWavelet) || (inputData.size(1) <= sizeWavelet)) {
		nbPlane--;
		sizeWavelet = twoToThe(nbPlane + 1);
	}


	return (RET_OK);
}


int ATrous::testSizes2D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long& nbPlane)
{
	//Test inputData, outputData, planesToUse, nbPlane
	unsigned long i;

	if (!inputData.isThereData()) { return (W_INIT_ERR); }
	if (!outputData.isThereData()) { return (W_INIT_ERR); }
	if ((inputData.dim() + 1) != outputData.dim()) { return (W_DIM_ERR); }
	if ((inputData.size(0) != outputData.size(0)) || (inputData.size(1) != outputData.size(1))) { return (W_SIZE_ERR); }
	if (planesToUse.dim() != 1) { return (W_DIM_ERR); }
	for (i = outputData.size(2); i < (unsigned long)planesToUse.size(0); i++) { planesToUse(i - 1) = false; }


	long sizeWavelet;

	if ((inputData.size(0) < 9) || (inputData.size(1) < 9)) { return (W_SIZE_ERR); }
	if (nbPlane < 2) { return (W_SIZE_ERR); }
	sizeWavelet = twoToThe(nbPlane + 1);
	while ((inputData.size(0) <= sizeWavelet) || (inputData.size(1) <= sizeWavelet)) {
		nbPlane--;
		sizeWavelet = twoToThe(nbPlane + 1);
	}


	return (nbPlane);
}


void ATrous::projectionV02D(unsigned short* srcData, float* destData, unsigned long realProj, unsigned long sizeDataX, unsigned long sizeDataY)
{
	//Makes a projection in V0 space to be in the apropriate condition of the multiresolution, 
	//realProj make be set to 1 if realProj=0 it just copy the data from srcData to destData

	int i, j, cur_i, cur_j, offs;
	double tab[3] = {0.1666666666, 0.6666666666666, 0.16666666666};
	if (realProj == 1) {
		i = 0;
		for (j = 0; j < (long)sizeDataY; j++) {
			offs = OFFSET(i, j);
			destData[offs] = 0.0;
			for (cur_i = -1; cur_i < 2; cur_i++) { for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i+cur_i, j+cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); } }
		}


		for (i = 1; i < (long)sizeDataX - 1; i++) {
			j = 0;
			offs = OFFSET(i, j);
			destData[offs] = 0.0;
			for (cur_i = -1; cur_i < 2; cur_i++) { for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i+cur_i, j+cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); } }


			for (j = 1; j < (long)sizeDataY - 1; j++) {
				offs = OFFSET(i, j);
				destData[offs] = 0.0;
				for (cur_i = -1; cur_i < 2; cur_i++) { for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i+cur_i, j+cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); } }
			}

			j = sizeDataY - 1;
			offs = OFFSET(i, j);
			destData[offs] = 0.0;
			for (cur_i = -1; cur_i < 2; cur_i++) { for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i+cur_i, j+cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); } }
		}

		i = sizeDataX - 1;
		for (j = 0; j < (long)sizeDataY; j++) {
			offs = OFFSET(i, j);
			destData[offs] = 0.0;
			for (cur_i = -1; cur_i < 2; cur_i++) { for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i+cur_i, j+cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); } }
		}
	}
	else { for (i = 0; i < (long)(sizeDataX * sizeDataY); i++) { destData[i] = (float)srcData[i]; } }
}

void ATrous::projectionV02D(float* srcData, float* destData, unsigned long realProj, unsigned long sizeDataX, unsigned long sizeDataY)
{
	//Makes a projection in V0 space to be in the apropriate condition of the multiresolution,
	//realProj make be set to 1 if realProj=0 it just copy the data from srcData to destData

	int i, j, cur_i, cur_j, offs;
	double tab[3] = {0.1666666666, 0.6666666666666, 0.16666666666};
	if (realProj == 1) {
		i = 0;
		for (j = 0; j < (long)(sizeDataY); j++) {
			offs = OFFSET(i, j);
			destData[offs] = 0.0;
			for (cur_i = -1; cur_i < 2; cur_i++) { for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i+cur_i, j+cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); } }
		}


		for (i = 1; i < (long)sizeDataX - 1; i++) {
			j = 0;
			offs = OFFSET(i, j);
			destData[offs] = 0.0;
			for (cur_i = -1; cur_i < 2; cur_i++) { for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i+cur_i, j+cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); } }


			for (j = 1; j < (long)sizeDataY - 1; j++) {
				offs = OFFSET(i, j);
				destData[offs] = 0.0;
				for (cur_i = -1; cur_i < 2; cur_i++) { for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i+cur_i, j+cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); } }
			}

			j = sizeDataY - 1;
			offs = OFFSET(i, j);
			destData[offs] = 0.0;
			for (cur_i = -1; cur_i < 2; cur_i++) { for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i+cur_i, j+cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); } }
		}

		i = sizeDataX - 1;
		for (j = 0; j < (long)sizeDataY; j++) {
			offs = OFFSET(i, j);
			destData[offs] = 0.0;
			for (cur_i = -1; cur_i < 2; cur_i++) { for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i+cur_i, j+cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); } }
		}
	}
	else { for (i = 0; i < (long)(sizeDataX * sizeDataY); i++) { destData[i] = (float)srcData[i]; } }
}


int ATrous::b3SplineFiltering2D(float* srcData, float* destData, unsigned long sizeDataX, unsigned long sizeDataY, unsigned long currentPlane)
{
	//Makes the projection in V(x) with hte knowlage of V(x-1)
	int i, j, step, offsetij;
	float coeffH0 = 3. / 8.;
	float coeffH1 = 1. / 4.;
	float coeffH2 = 1. / 16.;

	float* buff;
	int ret;

	if ((ret = temp2->realloc(sizeDataX, sizeDataY)) < 0) { return (ret); }
	buff = temp2->data;


	//type_border Type = Bord;

	step = twoToThe(currentPlane);


	//convolution in Y axe
	for (i = 0; i < (long)sizeDataX; i++) {
		for (j = 0; j < 2 * step; j++) {
			offsetij = OFFSET(i, j);
			buff[offsetij] = coeffH0 * srcData[offsetij] + coeffH1 * (srcData[OFFSETTEST(i, j-step)] + srcData[OFFSETTEST(i, j+step)]) + coeffH2 * (srcData[OFFSETTEST(i, j-2*step)] + srcData[OFFSETTEST(i, j+2*step)]);
		}

		for (j = 2 * step; j < (long)sizeDataY - 2 * step; j++) {
			offsetij = OFFSET(i, j);
			buff[offsetij] = coeffH0 * srcData[offsetij] + coeffH1 * (srcData[OFFSET(i, j-step)] + srcData[OFFSET(i, j+step)]) + coeffH2 * (srcData[OFFSET(i, j-2*step)] + srcData[OFFSET(i, j+2*step)]);
		}

		for (j = sizeDataY - 2 * step; j < (long)sizeDataY; j++) {
			offsetij = OFFSET(i, j);
			buff[offsetij] = coeffH0 * srcData[offsetij] + coeffH1 * (srcData[OFFSETTEST(i, j-step)] + srcData[OFFSETTEST(i, j+step)]) + coeffH2 * (srcData[OFFSETTEST(i, j-2*step)] + srcData[OFFSETTEST(i, j+2*step)]);
		}
	}

	//convolution in X axe
	for (i = 0; i < 2 * step; i++) {
		for (j = 0; j < (long)sizeDataY; j++) {
			offsetij = OFFSET(i, j);
			destData[offsetij] = coeffH0 * buff[offsetij] + coeffH1 * (buff[OFFSETTEST(i-step, j)] + buff[OFFSETTEST(i+step, j)]) + coeffH2 * (buff[OFFSETTEST(i-2*step, j)] + buff[OFFSETTEST(i+2*step, j)]);
		}
	}

	for (i = 2 * step; i < (long)sizeDataX - 2 * step; i++) {
		for (j = 0; j < (long)sizeDataY; j++) {
			offsetij = OFFSET(i, j);
			destData[offsetij] = coeffH0 * buff[offsetij] + coeffH1 * (buff[OFFSET(i-step, j)] + buff[OFFSET(i+step, j)]) + coeffH2 * (buff[OFFSET(i-2*step, j)] + buff[OFFSET(i+2*step, j)]);
		}
	}

	for (i = sizeDataX - 2 * step; i < (long)sizeDataX; i++) {
		for (j = 0; j < (long)sizeDataY; j++) {
			offsetij = OFFSET(i, j);
			destData[offsetij] = coeffH0 * buff[offsetij] + coeffH1 * (buff[OFFSETTEST(i-step, j)] + buff[OFFSETTEST(i+step, j)]) + coeffH2 * (buff[OFFSETTEST(i-2*step, j)] + buff[OFFSETTEST(i+2*step, j)]);
		}
	}

	return (RET_OK);
}


int ATrous::aTrousDecomposition2D(FltArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long nbPlane, bool realProj)
{
	/**	
		WAVELET TRANSFORM WITH "A TROUS" ALGORITHM
	
		Makes the transformation of the 2D signal "inputData".
		"inputData" must be a FltArray array of 2 dimentions ( 2D image).
	
	
		The algorithm makes a decompostion of the input signal in "outputData.size(2)" wavelet planes ("outputData").
		"outputData" is a FltArray of 3 dimentions. "outputData(0,0,4)" is the first point of the fifth wavelet plane. 
		The last plan is not a wavelet plane but the soft plane. So "outputData.size(2)" must be >1. If "outputData.size(2)"=2, 
		it computes the first wavelet plane and the soft plane. If the images are too short the algorithm won't be age too compute 
		all the planes. So "planesToUse" is a boolean array to indicate if the plane is computed. "planesToUse" must 1D array of 
		"outputData.size(2)" points.
			!!! The "outputData" array and "planesToUse" array must be allocated before calling this function !!!
	
		The "a trous" algorithm is a multiresolution analysis. The initialize the process, according to the theorem of multiresolution
		ones have to project the signal to a Riesz basis of part of L2(R2) space (V0 basis). But you can also considere that you signal 
		is smooth enough to be already in this basis. If "realProj" is true the projection is done.
	**/

	//float **waveletPlanes;
	unsigned long i;
	int ret;
	unsigned long currentPlane;
	float** waveletPlanes;
	unsigned long sizeDataX;
	unsigned long sizeDataY;


	sizeDataX = inputData.size(0);
	sizeDataY = inputData.size(1);


	//test if nbPlane is not too big !!!
	if ((ret = testSizes2D(inputData, outputData, planesToUse, nbPlane)) < 0) { return (ret); }


	//allocation of a direct way to reach the differents planes
	if ((waveletPlanes = (float**)malloc(sizeof(float*) * (nbPlane))) == NULL) { return (W_INIT_ERR); }
	for (currentPlane = 0; currentPlane < nbPlane; currentPlane++) { waveletPlanes[currentPlane] = outputData.data + sizeDataX * sizeDataY * currentPlane; }

	//do the initialisation of multiresolution algorthm (initialisation of the first wavelet plane)
	//signal is projected in V0 space
	projectionV02D(inputData.data, outputData.data, (int)realProj, sizeDataX, sizeDataY);

	//############  sigma of input data ##########
	//Projection of the input signal in each wavelet space.
	//The signal waveletPlanes[0] is the projection of the input signal in the W1 space,
	//the signal waveletPlanes[x] is the projection of the input signal in the W(x+1) space...
	for (currentPlane = 0; currentPlane < (unsigned long)(nbPlane - 1); currentPlane++) {
		//Firstly it computes in waveletPlanes[x+1] the projection of V(x) with the knowlege of the projection in V(x-1)
		if ((ret = b3SplineFiltering2D(waveletPlanes[currentPlane], waveletPlanes[currentPlane + 1], sizeDataX, sizeDataY, currentPlane)) < 0) { return (ret); }
		//Secondly it computes in waveletPlanes[x] the projection of W(x+1)
		//As proj( W(x+1) ) = proj( V(x) ) - proj( V(x+1) ) ;	proj(Vx) is in waveletPlanes[x+1] 
		//and proj(Vx+1) is in waveletPlanes[x+2]
		//	waveletPlanes[x+1][*]=waveletPlanes[x+1][*]-waveletPlanes[x+2][*]
		for (i = 0; i < (unsigned long)(sizeDataX * sizeDataY); i++) { (waveletPlanes[currentPlane])[i] -= (waveletPlanes[currentPlane + 1])[i]; }
	}
	free(waveletPlanes);
	return (nbPlane);
}

int ATrous::aTrousDecomposition2D(USrtArray& inputData, FltArray& outputData, BlArray& planesToUse, unsigned long nbPlane, bool realProj)
{
	/**	
		WAVELET TRANSFORM WITH "A TROUS" ALGORITHM
	
		Makes the transformation of the 2D signal "inputData".
		"inputData" must be a USrtArray array of 2 dimentions ( 2D image).
	
	
		The algorithm makes a decompostion of the input signal in "outputData.size(2)" wavelet planes ("outputData").
		"outputData" is a FltArray of 3 dimentions. "outputData(0,0,4)" is the first point of the fifth wavelet plane. 
		The last plan is not a wavelet plane but the soft plane. So "outputData.size(2)" must be >1. If "outputData.size(2)"=2, 
		it computes the first wavelet plane and the soft plane. If the images are too short the algorithm won't be age too compute 
		all the planes. So "PlanesToUse" is a boolean array to indicate if the plane is computed. "planesToUse" must 1D array of 
		"outputData.size(2)" points.
			!!! The "outputData" array and "planesToUse" array must be allocated before calling this function !!!
	
		The "a trous" algorithm is a multiresolution analysis. The initialize the process, according to the theorem of multiresolution
		ones have to project the signal to a Riesz basis of part of L2(R2) space (V0 basis). But you can also considere that you signal 
		is smooth enough to be already in this basis. If "realProj" is true the projection is done.
	**/


	unsigned long i;
	int ret;
	unsigned long currentPlane;
	float** waveletPlanes;
	unsigned long sizeDataX;
	unsigned long sizeDataY;

	sizeDataX = inputData.size(0);
	sizeDataY = inputData.size(1);

	//test if nbPlane is not too big !!!
	if ((ret = testSizes2D(inputData, outputData, planesToUse, nbPlane)) < 0) { return (ret); }


	//allocation of a direct way to reach the differents planes
	if ((waveletPlanes = (float**)malloc(sizeof(float*) * (nbPlane))) == NULL) { return (W_INIT_ERR); }
	for (currentPlane = 0; currentPlane < nbPlane; currentPlane++) { waveletPlanes[currentPlane] = outputData.data + sizeDataX * sizeDataY * currentPlane; }

	//do the initialisation of multiresolution algorthm (initialisation of the first wavelet plane)
	//signal is projected in V0 space
	projectionV02D(inputData.data, outputData.data, (int)realProj, sizeDataX, sizeDataY);

	//Projection of the input signal in each wavelet space.
	//The signal waveletPlanes[0] is the projection of the input signal in the W1 space,
	//the signal waveletPlanes[x] is the projection of the input signal in the W(x+1) space...
	for (currentPlane = 0; currentPlane < (unsigned long)(nbPlane - 1); currentPlane++) {
		//Firstly it computes in waveletPlanes[x+1] the projection of V(x) with the knowlege of the projection in V(x-1)
		if ((ret = b3SplineFiltering2D(waveletPlanes[currentPlane], waveletPlanes[currentPlane + 1], sizeDataX, sizeDataY, currentPlane)) < 0) { return (ret); }
		//Secondly it computes in waveletPlanes[x] the projection of W(x+1)
		//As proj( W(x+1) ) = proj( V(x) ) - proj( V(x+1) ) ;	proj(Vx) is in waveletPlanes[x+1] 
		//and proj(Vx+1) is in waveletPlanes[x+2]
		//	waveletPlanes[x+1][*]=waveletPlanes[x+1][*]-waveletPlanes[x+2][*]

		for (i = 0; i < (unsigned long)(sizeDataX * sizeDataY); i++) { (waveletPlanes[currentPlane])[i] -= (waveletPlanes[currentPlane + 1])[i]; }
	}

	free(waveletPlanes);
	return (nbPlane);
}


#undef	OFFSET
#undef	OFFSETTEST
