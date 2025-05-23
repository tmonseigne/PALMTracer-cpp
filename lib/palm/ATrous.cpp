#include "ATrous.hpp"

#include <algorithm>
#include <array>

#include "Constants.hpp"

//-------------------------------------------------------------------------------------------------
void ATrous::ATrousDecomposition2D(const unsigned short* image, const size_t height, const size_t width)
{
	/**
		WAVELET TRANSFORM WITH "A TROUS" ALGORITHM

		Makes the transformation of the 2D signal "image".
		"image" must be a unsigned short array of 2 dimentions ( 2D image).


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
	width_  = width;
	size_   = height_ * width_;
	tmp_.resize(size_);
	for (auto& wavelet : wDecomp_) { wavelet.resize(size_); }

	//do the initialisation of multiresolution algorthm (initialisation of the first wavelet plane) signal is projected in V0 space
	std::transform(image, image + size_, wDecomp_[0].begin(), [](const unsigned short val) { return static_cast<double>(val); });

	//Projection of the input signal in each wavelet space. the signal waveletPlanes[x] is the projection of the input signal in the W(x+1) space
	for (size_t i = 0; i < N_PLANES - 1; ++i) {
		//Firstly it computes in waveletPlanes[x+1] the projection of V(x) with the knowlege of the projection in V(x-1)
		b3SplineFiltering2D(wDecomp_[i], wDecomp_[i + 1], i);
		//Secondly it computes in waveletPlanes[x] the projection of W(x+1)
		//As proj(W(x+1)) = proj(W(x)) - proj(W(x+1));	proj(Vx) is in waveletPlanes[x+1] and proj(Vx+1) is in waveletPlanes[x+2]
		for (size_t j = 0; j < size_; ++j) { wDecomp_[i][j] -= wDecomp_[i + 1][j]; }
	}
}

//-------------------------------------------------------------------------------------------------
size_t ATrous::offsetTest(const int i, const int j) const
{
	const auto Clamped_I = static_cast<size_t>(std::clamp(i, 0, static_cast<int>(height_) - 1));
	const auto Clamped_J = static_cast<size_t>(std::clamp(j, 0, static_cast<int>(width_) - 1));
	return offset(Clamped_I, Clamped_J);
}

//-------------------------------------------------------------------------------------------------
void ATrous::b3SplineFiltering2D(const std::vector<double>& in, std::vector<double>& out, const size_t plane)
{
	//Makes the projection in V(x) with hte knowlage of V(x-1)
	constexpr double CoeffH0 = 3.0 / 8.0, CoeffH1 = 1.0 / 4.0, CoeffH2 = 1.0 / 16.0;

	const size_t Step       = twoToThe(plane);
	const size_t DoubleStep = 2 * Step;
	const int IStep         = static_cast<int>(Step);
	const int IDoubleStep   = static_cast<int>(DoubleStep);

	auto get = [&](const std::vector<double>& data, const int y, const int x) -> double { return data[offsetTest(y, x)]; };

	// Lambda pour simplifier l'accès au voisinage avec ou sans test de bord
	auto convolutionY = [&](const size_t idx, const int i, const int j) -> double {
		return CoeffH0 * in[idx] + CoeffH1 * (get(in, i, j - IStep) + get(in, i, j + IStep)) +
		       CoeffH2 * (get(in, i, j - IDoubleStep) + get(in, i, j + IDoubleStep));
	};

	// Lambda pour simplifier l'accès au voisinage avec ou sans test de bord
	auto convolutionX = [&](const size_t idx, const int i, const int j) -> double {
		return CoeffH0 * tmp_[idx] + CoeffH1 * (get(tmp_, i - IStep, j) + get(tmp_, i + IStep, j)) +
		       CoeffH2 * (get(tmp_, i - IDoubleStep, j) + get(tmp_, i + IDoubleStep, j));
	};

	// Convolution horizontale (Y-axis)
	#pragma omp parallel for
	for (int idx = 0; idx < static_cast<int>(size_); ++idx) { tmp_[idx] = convolutionY(idx, idx / static_cast<int>(width_), idx % static_cast<int>(width_)); }

	// Convolution verticale (X-axis)
	#pragma omp parallel for
	for (int idx = 0; idx < static_cast<int>(size_); ++idx) { out[idx] = convolutionX(idx, idx / static_cast<int>(width_), idx % static_cast<int>(width_)); }
}
