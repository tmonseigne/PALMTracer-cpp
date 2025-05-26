#include "palm/ATrous.hpp"

#include <algorithm>
#include <array>

#include "Constants.hpp"

extern "C"
void WaveletProjections(const unsigned short* image, size_t height, size_t width, size_t nbPlane, std::vector<std::vector<double>>& decomp);

size_t ATrous::offsetTest(const int i, const int j) const
{
	const auto Clamped_I = static_cast<size_t>(std::clamp(i, 0, static_cast<int>(height_) - 1));
	const auto Clamped_J = static_cast<size_t>(std::clamp(j, 0, static_cast<int>(width_) - 1));
	return offset(Clamped_I, Clamped_J);
}


int ATrous::testSizes2D(bool* planesToUse, size_t& nbPlane) const
{
	//Test in, out, planesToUse, nbPlane
	std::fill(planesToUse + wDecomp_.size(), planesToUse + 9, false);
	if (nbPlane < 2) { return W_SIZE_ERR; }
	while (height_ <= twoToThe(nbPlane + 1) || width_ <= twoToThe(nbPlane + 1)) {
		if (nbPlane-- <= 2) { return W_SIZE_ERR; } // s�curit� si boucle infinie
	}
	return static_cast<int>(nbPlane);
}

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

int ATrous::ATrousDecomposition2D(const unsigned short* image, const size_t height, const size_t width, bool* planesToUse, size_t nbPlane)
{
	height_ = height;
	width_  = width;
	size_   = height_ * width_;
	int ret;
	if ((ret = testSizes2D(planesToUse, nbPlane)) < 0) { return ret; }
	WaveletProjections(image, height, width, nbPlane, wDecomp_);
	return static_cast<int>(nbPlane);
}
