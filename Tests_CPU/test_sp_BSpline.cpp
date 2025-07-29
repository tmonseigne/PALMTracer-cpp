#include "gtest/gtest.h"

#include <algorithm>
#include <numeric>

#include "Structs.hpp"
#include "TiffManagement.hpp"
#include "palm/Palm.hpp"
#include "signal_processing/BSpline3D.hpp"
#include "signal_processing/Matrix.hpp"
#include "signal_processing/Spline.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
TEST(BSpline3DTests, DiracImpulse)
{
	constexpr Vec3I Dim{.x = 5, .y = 5, .z = 5};
	const size_t N = Dim.GetSize();

	// Volume : Dirac au centre
	std::vector volume(N, 0.0);
	const size_t CenterIdx = Dim.GetMid().GetIndex(Dim);
	volume[CenterIdx]      = 1.0;

	// Paramètres du filtre
	constexpr size_t Degree = 3;
	constexpr Vec3D Lambda{.x = 0.1, .y = 0.1, .z = 0.1}; // Régularisation modérée

	// Application du filtrage B-spline
	const BSpline3D Bspline(volume, Dim, Degree, Lambda);

	//std::cout << Bspline.dimPadded_ << "\n";
	//matrix::Print(Bspline.coeffsPadded_, Bspline.dimPadded_);
	//std::cout << Bspline.dim_ << "\n";
	//matrix::Print(Bspline.coeffsHDR_, Bspline.dim_);
	//std::cout << Bspline.dimFinal_ << "\n";
	//matrix::Print(Bspline.coeffsFinal_, Bspline.dimFinal_);
	//size_t idx = 0;
	//for (size_t i = 0; i < Bspline.dimFinal_.GetSize(); ++i) {
	//	for (size_t j = 0; j < 64; ++j, ++idx) { std::cout << Bspline.coeffsFinal_[idx] << ", "; }
	//	std::cout << "\n";
	//}
}
} // namespace palm_tracer
