#include <gtest/gtest.h>

#include <algorithm>
#include <numeric>
#include <type_traits>
#include <utility>

#include "Palm.hpp"
#include "Structs.hpp"
#include "TiffManagement.hpp"
#include "signal_processing/BSpline3D.hpp"
#include "signal_processing/Matrix.hpp"
#include "signal_processing/Spline.hpp"

namespace palm_tracer::signal_processing {
/// <summary>Vérifie le transfert effectif des allocations sans exécuter l'algorithme spline.</summary>
TEST(BSpline3DTests, MovesBuffers)
{
	static_assert(std::is_nothrow_move_constructible_v<BSpline3D>);
	static_assert(std::is_nothrow_move_assignable_v<BSpline3D>);
	BSpline3D source;
	source.dim_              = {.x = 2, .y = 3, .z = 4};
	source.coeffs_           = {1, 2, 3};
	source.coeffsPadded_     = {4, 5};
	source.coeffsHDR_        = {6, 7};
	source.coeffsFinal_      = {8, 9};
	source.filters_[0]       = {{0.25, 0.5}};
	const auto* const Coeffs = source.coeffs_.data();
	const auto* const Padded = source.coeffsPadded_.data();
	const auto* const HDR    = source.coeffsHDR_.data();
	const auto* const Final  = source.coeffsFinal_.data();
	const auto* const Filter = source.filters_[0].data();

	BSpline3D moved(std::move(source));
	EXPECT_EQ(moved.coeffs_.data(), Coeffs);
	EXPECT_EQ(moved.coeffsPadded_.data(), Padded);
	EXPECT_EQ(moved.coeffsHDR_.data(), HDR);
	EXPECT_EQ(moved.coeffsFinal_.data(), Final);
	EXPECT_EQ(moved.filters_[0].data(), Filter);
	EXPECT_EQ(moved.dim_.GetSize(), 24);

	BSpline3D destination;
	destination.coeffsFinal_ = {-1};
	destination              = std::move(moved);
	EXPECT_EQ(destination.coeffs_.data(), Coeffs);
	EXPECT_EQ(destination.coeffsPadded_.data(), Padded);
	EXPECT_EQ(destination.coeffsHDR_.data(), HDR);
	EXPECT_EQ(destination.coeffsFinal_.data(), Final);
	EXPECT_EQ(destination.filters_[0].data(), Filter);
	EXPECT_EQ(destination.coeffsFinal_, (std::vector<double>{8, 9}));

	// La copie reste disponible et ne partage pas les coefficients modifiables.
	BSpline3D copied(destination);
	copied.coeffsFinal_[0] = 10;
	EXPECT_DOUBLE_EQ(destination.coeffsFinal_[0], 8);
	source                 = copied;
	source.coeffsFinal_[0] = 11;
	EXPECT_DOUBLE_EQ(copied.coeffsFinal_[0], 10);
}

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
	//matrix::PrintVolume(Bspline.coeffsPadded_, Bspline.dimPadded_);
	//std::cout << Bspline.dim_ << "\n";
	//matrix::PrintVolume(Bspline.coeffsHDR_, Bspline.dim_);
	//std::cout << Bspline.dimFinal_ << "\n";
	//matrix::PrintVolume(Bspline.coeffsFinal_, Bspline.dimFinal_);
	//size_t idx = 0;
	//for (size_t i = 0; i < Bspline.dimFinal_.GetSize(); ++i) {
	//	for (size_t j = 0; j < 64; ++j, ++idx) { std::cout << Bspline.coeffsFinal_[idx] << ", "; }
	//	std::cout << "\n";
	//}
}
} // namespace palm_tracer::signal_processing
