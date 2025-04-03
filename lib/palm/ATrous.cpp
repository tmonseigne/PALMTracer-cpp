#include "ATrous.h"

#include "tools/allDataLib.h"
#include "tools/NoError.h"


#define	OFFSET(i,j)	offset((i),(j),width)
#define	OFFSETTEST(i,j)	offsetTest((i),(j),height,width)


size_t twoToThe(const size_t val)
{
	if (val == 0) { return 1; }           // 2^0 = 1
	return static_cast<size_t>(1) << val; // 1 décalé à gauche de 'val' bits
}

int ATrous::offsetTest(int i, int j, const int maxSizeI, const int maxSizeJ)
{
	//Same as offset but if i or j are out of bounds, it changes them to be in bounds.
	if (i < 0) { i = 0; }
	if (i >= maxSizeI) { i = maxSizeI - 1; }
	if (j < 0) { j = 0; }
	if (j >= maxSizeJ) { j = maxSizeJ - 1; }
	return offset(i, j, maxSizeJ);
}


//Compute the offset translate   tab[i][j] to tab[i*maxSizeJ+j]
int ATrous::offset(const int i, const int j, const int maxSizeJ) { return i * maxSizeJ + j; }


int ATrous::testSizes2D(const size_t height, const size_t width, FltArray& out, bool* planesToUse, size_t& nbPlane)
{
	//Test in, out, planesToUse, nbPlane

	if (!out.isThereData()) { return W_INIT_ERR; }
	for (size_t i = out.size(2); i < 9; ++i) { planesToUse[i - 1] = false; }


	if (nbPlane < 2) { return W_SIZE_ERR; }
	long sizeWavelet = twoToThe(nbPlane + 1);
	while (height <= sizeWavelet || width <= sizeWavelet) {
		nbPlane--;
		sizeWavelet = twoToThe(nbPlane + 1);
	}


	return nbPlane;
}


void ATrous::projectionV02D(const unsigned short* srcData, float* destData, const size_t realProj, const size_t height, const size_t width)
{
	//Makes a projection in V0 space to be in the apropriate condition of the multiresolution, 
	//realProj make be set to 1 if realProj=0 it just copy the data from srcData to destData

	int i, j, cur_i, cur_j, offs;
	const double tab[3] = {0.1666666666, 0.6666666666666, 0.16666666666};
	if (realProj == 1) {
		i = 0;
		for (j = 0; j < width; j++) {
			offs = OFFSET(i, j);
			destData[offs] = 0.0;
			for (cur_i = -1; cur_i < 2; cur_i++) {
				for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i + cur_i, j + cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); }
			}
		}


		for (i = 1; i < height - 1; i++) {
			j = 0;
			offs = OFFSET(i, j);
			destData[offs] = 0.0;
			for (cur_i = -1; cur_i < 2; cur_i++) {
				for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i + cur_i, j + cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); }
			}


			for (j = 1; j < width - 1; j++) {
				offs = OFFSET(i, j);
				destData[offs] = 0.0;
				for (cur_i = -1; cur_i < 2; cur_i++) {
					for (cur_j = -1; cur_j < 2; cur_j++) {
						destData[offs] += (float)(srcData[OFFSETTEST(i + cur_i, j + cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]);
					}
				}
			}

			j = width - 1;
			offs = OFFSET(i, j);
			destData[offs] = 0.0;
			for (cur_i = -1; cur_i < 2; cur_i++) {
				for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i + cur_i, j + cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); }
			}
		}

		i = height - 1;
		for (j = 0; j < (long)width; j++) {
			offs = OFFSET(i, j);
			destData[offs] = 0.0;
			for (cur_i = -1; cur_i < 2; cur_i++) {
				for (cur_j = -1; cur_j < 2; cur_j++) { destData[offs] += (float)(srcData[OFFSETTEST(i + cur_i, j + cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]); }
			}
		}
	}
	else { for (i = 0; i < (height * width); i++) { destData[i] = (float)srcData[i]; } }
}

int ATrous::b3SplineFiltering2D(const float* srcData, float* destData, const size_t height, const size_t width, const size_t currentPlane)
{
	//Makes the projection in V(x) with hte knowlage of V(x-1)
	int i, j, offsetij;
	const float coeffH0 = 3. / 8.;
	const float coeffH1 = 1. / 4.;
	const float coeffH2 = 1. / 16.;

	const int step = twoToThe(currentPlane);


	//convolution in Y axe
	for (i = 0; i < (long)height; i++) {
		for (j = 0; j < 2 * step; j++) {
			offsetij = OFFSET(i, j);
			tmp_[offsetij] = coeffH0 * srcData[offsetij] + coeffH1 * (srcData[OFFSETTEST(i, j - step)] + srcData[OFFSETTEST(i, j + step)]) + coeffH2 * (
				                 srcData[OFFSETTEST(i, j - 2 * step)] + srcData[OFFSETTEST(i, j + 2 * step)]);
		}

		for (j = 2 * step; j < (long)width - 2 * step; j++) {
			offsetij = OFFSET(i, j);
			tmp_[offsetij] = coeffH0 * srcData[offsetij] + coeffH1 * (srcData[OFFSET(i, j - step)] + srcData[OFFSET(i, j + step)]) + coeffH2 * (
				                 srcData[OFFSET(i, j - 2 * step)] + srcData[OFFSET(i, j + 2 * step)]);
		}

		for (j = width - 2 * step; j < (long)width; j++) {
			offsetij = OFFSET(i, j);
			tmp_[offsetij] = coeffH0 * srcData[offsetij] + coeffH1 * (srcData[OFFSETTEST(i, j - step)] + srcData[OFFSETTEST(i, j + step)]) + coeffH2 * (
				                 srcData[OFFSETTEST(i, j - 2 * step)] + srcData[OFFSETTEST(i, j + 2 * step)]);
		}
	}

	//convolution in X axe
	for (i = 0; i < 2 * step; i++) {
		for (j = 0; j < (long)width; j++) {
			offsetij = OFFSET(i, j);
			destData[offsetij] = coeffH0 * tmp_[offsetij] + coeffH1 * (tmp_[OFFSETTEST(i - step, j)] + tmp_[OFFSETTEST(i + step, j)]) + coeffH2 * (
				                     tmp_[OFFSETTEST(i - 2 * step, j)] + tmp_[OFFSETTEST(i + 2 * step, j)]);
		}
	}

	for (i = 2 * step; i < height - 2 * step; i++) {
		for (j = 0; j < width; j++) {
			offsetij = OFFSET(i, j);
			destData[offsetij] = coeffH0 * tmp_[offsetij] + coeffH1 * (tmp_[OFFSET(i - step, j)] + tmp_[OFFSET(i + step, j)]) + coeffH2 * (
				                     tmp_[OFFSET(i - 2 * step, j)] + tmp_[OFFSET(i + 2 * step, j)]);
		}
	}

	for (i = height - 2 * step; i < height; i++) {
		for (j = 0; j < width; j++) {
			offsetij = OFFSET(i, j);
			destData[offsetij] = coeffH0 * tmp_[offsetij] + coeffH1 * (tmp_[OFFSETTEST(i - step, j)] + tmp_[OFFSETTEST(i + step, j)]) + coeffH2 * (
				                     tmp_[OFFSETTEST(i - 2 * step, j)] + tmp_[OFFSETTEST(i + 2 * step, j)]);
		}
	}

	return RET_OK;
}

int ATrous::ATrousDecomposition2D(unsigned short* image, const size_t height, const size_t width, FltArray& outputData, bool* planesToUse, size_t nbPlane, const bool realProj)
{
	/**
		WAVELET TRANSFORM WITH "A TROUS" ALGORITHM

		Makes the transformation of the 2D signal "image".
		"image" must be a USrtArray array of 2 dimentions ( 2D image).


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


	int ret;
	size_t currentPlane;
	float** waveletPlanes;

	tmp_.resize(height * width);

	//test if nbPlane is not too big !!!
	if ((ret = testSizes2D(height, width, outputData, planesToUse, nbPlane)) < 0) { return ret; }


	//allocation of a direct way to reach the differents planes
	if ((waveletPlanes = (float**)malloc(sizeof(float*) * nbPlane)) == NULL) { return W_INIT_ERR; }
	for (currentPlane = 0; currentPlane < nbPlane; currentPlane++) { waveletPlanes[currentPlane] = outputData.data + height * width * currentPlane; }

	//do the initialisation of multiresolution algorthm (initialisation of the first wavelet plane)
	//signal is projected in V0 space
	projectionV02D(image, outputData.data, (int)realProj, height, width);

	//Projection of the input signal in each wavelet space.
	//The signal waveletPlanes[0] is the projection of the input signal in the W1 space,
	//the signal waveletPlanes[x] is the projection of the input signal in the W(x+1) space...
	for (currentPlane = 0; currentPlane < (nbPlane - 1); currentPlane++) {
		//Firstly it computes in waveletPlanes[x+1] the projection of V(x) with the knowlege of the projection in V(x-1)
		if ((ret = b3SplineFiltering2D(waveletPlanes[currentPlane], waveletPlanes[currentPlane + 1], height, width, currentPlane)) < 0) { return ret; }
		//Secondly it computes in waveletPlanes[x] the projection of W(x+1)
		//As proj( W(x+1) ) = proj( V(x) ) - proj( V(x+1) ) ;	proj(Vx) is in waveletPlanes[x+1] 
		//and proj(Vx+1) is in waveletPlanes[x+2]
		//	waveletPlanes[x+1][*]=waveletPlanes[x+1][*]-waveletPlanes[x+2][*]

		for (size_t i = 0; i < (height * width); i++) { waveletPlanes[currentPlane][i] -= waveletPlanes[currentPlane + 1][i]; }
	}

	free(waveletPlanes);
	return nbPlane;
}


#undef	OFFSET
#undef	OFFSETTEST
