#include "ATrous.h"

#include <algorithm>

#include "tools/NoError.h"


size_t twoToThe(const size_t val)
{
	if (val == 0) { return 1; }           // 2^0 = 1
	return static_cast<size_t>(1) << val; // 1 décalé à gauche de 'val' bits
}

size_t ATrous::offsetTest(int i, int j) const
{
	i = std::clamp(i, 0, static_cast<int>(height_) - 1);
	j = std::clamp(j, 0, static_cast<int>(width_) - 1);
	return offset(static_cast<size_t>(i), static_cast<size_t>(j));
}


//Compute the offset translate   tab[i][j] to tab[i*maxSizeJ+j]
size_t ATrous::offset(const size_t i, const size_t j) const { return i * width_ + j; }


int ATrous::testSizes2D(bool* planesToUse, size_t& nbPlane) const
{
	//Test in, out, planesToUse, nbPlane

	for (size_t i = wDecomp_.size(); i < 9; ++i) { planesToUse[i - 1] = false; }


	if (nbPlane < 2) { return W_SIZE_ERR; }
	long sizeWavelet = twoToThe(nbPlane + 1);
	while (height_ <= sizeWavelet || width_ <= sizeWavelet) {
		nbPlane--;
		sizeWavelet = twoToThe(nbPlane + 1);
	}

	return nbPlane;
}


void ATrous::projectionV02D(const unsigned short* image, const size_t realProj)
{
	//Makes a projection in V0 space to be in the apropriate condition of the multiresolution, 
	//realProj make be set to 1 if realProj=0 it just copy the data from image to destData

	constexpr double tab[3] = {0.1666666666, 0.6666666666666, 0.16666666666};
	if (realProj == 1) {
		int i = 0, j, cur_i, cur_j;
		size_t idx;
		for (j = 0; j < width_; ++j) {
			idx = offset(i, j);
			wDecomp_[0][idx] = 0.0;
			for (cur_i = -1; cur_i < 2; ++cur_i) {
				for (cur_j = -1; cur_j < 2; ++cur_j) { wDecomp_[0][idx] += image[offsetTest(i + cur_i, j + cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]; }
			}
		}


		for (i = 1; i < height_ - 1; ++i) {
			j = 0;
			idx = offset(i, j);
			wDecomp_[0][idx] = 0.0;
			for (cur_i = -1; cur_i < 2; ++cur_i) {
				for (cur_j = -1; cur_j < 2; ++cur_j) { wDecomp_[0][idx] += image[offsetTest(i + cur_i, j + cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]; }
			}


			for (j = 1; j < width_ - 1; ++j) {
				idx = offset(i, j);
				wDecomp_[0][idx] = 0.0;
				for (cur_i = -1; cur_i < 2; ++cur_i) {
					for (cur_j = -1; cur_j < 2; ++cur_j) { wDecomp_[0][idx] += image[offsetTest(i + cur_i, j + cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]; }
				}
			}

			j = width_ - 1;
			idx = offset(i, j);
			wDecomp_[0][idx] = 0.0;
			for (cur_i = -1; cur_i < 2; ++cur_i) {
				for (cur_j = -1; cur_j < 2; ++cur_j) { wDecomp_[0][idx] += image[offsetTest(i + cur_i, j + cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]; }
			}
		}

		i = height_ - 1;
		for (j = 0; j < width_; ++j) {
			idx = offset(i, j);
			wDecomp_[0][idx] = 0.0;
			for (cur_i = -1; cur_i < 2; ++cur_i) {
				for (cur_j = -1; cur_j < 2; ++cur_j) { wDecomp_[0][idx] += image[offsetTest(i + cur_i, j + cur_j)] * tab[cur_i + 1] * tab[cur_j + 1]; }
			}
		}
	}
	else { for (size_t i = 0; i < size_; ++i) { wDecomp_[0][i] = image[i]; } }
}

void ATrous::b3SplineFiltering2D(const std::vector<double>& in, std::vector<double>& out, const size_t plane)
{
	//Makes the projection in V(x) with hte knowlage of V(x-1)
	int i, j;
	size_t idx;
	constexpr double CoeffH0 = 3. / 8.;
	constexpr double CoeffH1 = 1. / 4.;
	constexpr double CoeffH2 = 1. / 16.;

	const size_t step = twoToThe(plane);


	//convolution in Y axe
	for (i = 0; i < height_; ++i) {
		for (j = 0; j < 2 * step; ++j) {
			idx = offset(i, j);
			tmp_[idx] = CoeffH0 * in[idx] + CoeffH1 * (in[offsetTest(i, j - step)] + in[offsetTest(i, j + step)])
			            + CoeffH2 * (in[offsetTest(i, j - 2 * step)] + in[offsetTest(i, j + 2 * step)]);
		}

		for (j = 2 * step; j < width_ - 2 * step; ++j) {
			idx = offset(i, j);
			tmp_[idx] = CoeffH0 * in[idx] + CoeffH1 * (in[offset(i, j - step)] + in[offset(i, j + step)])
			            + CoeffH2 * (in[offset(i, j - 2 * step)] + in[offset(i, j + 2 * step)]);
		}

		for (j = width_ - 2 * step; j < width_; ++j) {
			idx = offset(i, j);
			tmp_[idx] = CoeffH0 * in[idx] + CoeffH1 * (in[offsetTest(i, j - step)] + in[offsetTest(i, j + step)])
			            + CoeffH2 * (in[offsetTest(i, j - 2 * step)] + in[offsetTest(i, j + 2 * step)]);
		}
	}

	//convolution in X axe
	for (i = 0; i < 2 * step; ++i) {
		for (j = 0; j < width_; ++j) {
			idx = offset(i, j);
			out[idx] = CoeffH0 * tmp_[idx] + CoeffH1 * (tmp_[offsetTest(i - step, j)] + tmp_[offsetTest(i + step, j)])
			           + CoeffH2 * (tmp_[offsetTest(i - 2 * step, j)] + tmp_[offsetTest(i + 2 * step, j)]);
		}
	}

	for (i = 2 * step; i < height_ - 2 * step; ++i) {
		for (j = 0; j < width_; ++j) {
			idx = offset(i, j);
			out[idx] = CoeffH0 * tmp_[idx] + CoeffH1 * (tmp_[offset(i - step, j)] + tmp_[offset(i + step, j)])
			           + CoeffH2 * (tmp_[offset(i - 2 * step, j)] + tmp_[offset(i + 2 * step, j)]);
		}
	}

	for (i = height_ - 2 * step; i < height_; ++i) {
		for (j = 0; j < width_; ++j) {
			idx = offset(i, j);
			out[idx] = CoeffH0 * tmp_[idx] + CoeffH1 * (tmp_[offsetTest(i - step, j)] + tmp_[offsetTest(i + step, j)])
			           + CoeffH2 * (tmp_[offsetTest(i - 2 * step, j)] + tmp_[offsetTest(i + 2 * step, j)]);
		}
	}
}

int ATrous::ATrousDecomposition2D(const unsigned short* image, const size_t height, const size_t width, bool* planesToUse, size_t nbPlane, const bool realProj)
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

	height_ = height;
	width_ = width;
	size_ = height_ * width_;
	int ret;
	tmp_.resize(size_);

	//test if nbPlane is not too big !!!
	if ((ret = testSizes2D(planesToUse, nbPlane)) < 0) { return ret; }


	//allocation of a direct way to reach the differents planes

	//do the initialisation of multiresolution algorthm (initialisation of the first wavelet plane)
	//signal is projected in V0 space
	projectionV02D(image, realProj);

	//Projection of the input signal in each wavelet space.
	//the signal waveletPlanes[x] is the projection of the input signal in the W(x+1) space
	for (size_t i = 0; i < nbPlane - 1; ++i) {
		//Firstly it computes in waveletPlanes[x+1] the projection of V(x) with the knowlege of the projection in V(x-1)
		b3SplineFiltering2D(wDecomp_[i], wDecomp_[i + 1], i);
		//Secondly it computes in waveletPlanes[x] the projection of W(x+1)
		//As proj( W(x+1) ) = proj( V(x) ) - proj( V(x+1) ) ;	proj(Vx) is in waveletPlanes[x+1] 
		//and proj(Vx+1) is in waveletPlanes[x+2]
		//	waveletPlanes[x+1][*]=waveletPlanes[x+1][*]-waveletPlanes[x+2][*]
		for (size_t j = 0; j < size_; ++j) { wDecomp_[i][j] -= wDecomp_[i + 1][j]; }
	}

	return static_cast<int>(nbPlane);
}
