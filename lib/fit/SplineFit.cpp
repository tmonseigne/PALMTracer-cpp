#include "SplineFit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <execution>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

#include "Utils.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
void SplineFit::parseParams()
{
	roi_         = static_cast<size_t>(fitParams_[0]);
	splineSizeX_ = fitParams_[1];
	splineSizeY_ = fitParams_[2];
	splineSizeZ_ = fitParams_[3];
}

//-------------------------------------------------------------------------------------------------
void SplineFit::fit(const SegmentationInfo& infos, std::array<double, N_LOCALIZATIONS>& point)
{
	double offset    = 0,
	       intensity = 0,
	       x         = infos.x_,
	       y         = infos.y_,
	       z         = infos.z_,
	       rss       = -1;

	// Vérification si la zone d'intérêt est entièrement dans l'image (évite les bords)
	const size_t Ux = static_cast<size_t>(x), Uy = static_cast<size_t>(y);
	if (roiHalf_ <= Ux && Ux <= width_ - roiHalf_ - 1 && roiHalf_ <= Uy && Uy <= height_ - roiHalf_ - 1) {
		// Création des variables
		double lambda = LAMBDA_0, minI, maxI;
		std::vector<double> alpha(nParams2_), alphaR(nParams2_), alphaRi(nParams2_),
		                    beta(nParams_), beta_l_deltaP(nParams_), p_try(nParams_), deltaP(nParams_),
		                    data(n_), hat(n_), weight(n_), deltaY(n_), deltaY_w(n_);
		std::vector<size_t> coords(2 * n_);

		// Initialization of data & t (A vérifier)
		preCompute(Ux, Uy, data, coords, minI, maxI);
		std::ranges::fill(weight, 1);

		// Initialisation du/des P :
		// - Calcul du centre de masse de la région pour x et y
		CenterOfMass(data, roi_, x, y);
	}
	fill(infos, point, offset, intensity, x, y, z);
}
}
