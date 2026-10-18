#include "ATrous.hpp"

#include <algorithm>
#include <array>

#include "Constants.hpp"

namespace palm_tracer::localization {
ATrous::ATrous(const uint16_t* image, const Vec2I& dim)
	: dim_(dim), size_(dim.GetSize())
{
	tmp_.resize(size_);
	for (auto& wavelet : wDecomp_) { wavelet.resize(size_); }
	decomposition2D(image);
}

//-------------------------------------------------------------------------------------------------
void ATrous::decomposition2D(const uint16_t* image)
{
	// Décomposition multirésolution à trous : le signal initial constitue l'approximation V0.
	// Chaque lissage B3 produit V(j+1), puis le détail W(j) = V(j) - V(j+1).
	// Les N_PLANES-1 premiers tableaux contiennent les détails ; le dernier conserve l'approximation résiduelle.
	// Aucun rééchantillonnage ni aucune projection initiale supplémentaire n'est effectué ici.
	std::transform(image, image + size_, wDecomp_[0].begin(), [](const uint16_t val) { return static_cast<double>(val); });

	// Construction successive des approximations et des détails, à résolution spatiale constante.
	for (size_t i = 0; i < N_PLANES - 1; ++i) {
		// Approximation suivante par filtrage B3 à l'échelle courante.
		b3SplineFiltering2D(wDecomp_[i], wDecomp_[i + 1], i);
		// Détail à cette échelle : différence entre les deux approximations successives.
		for (size_t j = 0; j < size_; ++j) { wDecomp_[i][j] -= wDecomp_[i + 1][j]; }
	}
}

//-------------------------------------------------------------------------------------------------
size_t ATrous::offsetTest(const int i, const int j) const
{
	const auto Clamped_I = static_cast<size_t>(std::clamp(i, 0, static_cast<int>(dim_.y) - 1));
	const auto Clamped_J = static_cast<size_t>(std::clamp(j, 0, static_cast<int>(dim_.x) - 1));
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
	for (int idx = 0; idx < static_cast<int>(size_); ++idx) { tmp_[idx] = convolutionY(idx, idx / static_cast<int>(dim_.x), idx % static_cast<int>(dim_.x)); }

	// Convolution verticale (X-axis)
	#pragma omp parallel for
	for (int idx = 0; idx < static_cast<int>(size_); ++idx) { out[idx] = convolutionX(idx, idx / static_cast<int>(dim_.x), idx % static_cast<int>(dim_.x)); }
}
} // namespace palm_tracer::localization
