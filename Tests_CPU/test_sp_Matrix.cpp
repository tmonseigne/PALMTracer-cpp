#include <gtest/gtest.h>

#include <array>
#include <limits>
#include <numeric>
#include <stdexcept>
#include "signal_processing/Matrix.hpp"

// Compatibilité avec les anciennes versions de GoogleTest.
#ifndef INSTANTIATE_TEST_SUITE_P
#define INSTANTIATE_TEST_SUITE_P INSTANTIATE_TEST_CASE_P
#endif

namespace palm_tracer::signal_processing {
constexpr double EPSILON = 1e-10;
const std::vector NAN_VECTOR(3, std::numeric_limits<double>::quiet_NaN());

/// <summary>Matrices inversibles qui exigent un pivot ou un équilibrage des lignes.</summary>
class GaussJordanPivotTests : public testing::TestWithParam<std::vector<double>> { };

/// <summary>Matrices dont l'inversion doit échouer explicitement.</summary>
class GaussJordanFailureTests : public testing::TestWithParam<std::vector<double>> { };

/// <summary>Contrats de taille et de séparation des buffers de l'inversion.</summary>
class GaussJordanArgumentTests : public testing::TestWithParam<size_t> { };

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, Print)
{
	const std::vector<std::vector<double>> M1 = {{1, 2}, {3, 4}};
	const std::vector<double> M2              = {2, 1, 2, 3, 3, 4, 5, -1};
	matrix::PrintMatrix(M1);
	matrix::PrintMatrix(M2, 4, 2);
	matrix::PrintVolume(M2, Vec3I(2, 2, 2));
	matrix::PrintVector(M2);
	EXPECT_TRUE(true);
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, SendIrregular)
{
	const std::vector<std::vector<double>> Data = {{1, 2}, {3, 4, 5}, {}};
	const std::vector<double> Ref               = {2, 1, 2, 3, 3, 4, 5, -1};

	std::vector res(Ref.size(), 0.0);
	const size_t N = matrix::SendIrregularArray(Data, res.data());

	ASSERT_EQ(N, Ref.size()) << "Nombre d'éléments écrits inattendu";
	for (size_t i = 0; i < Ref.size(); ++i) { EXPECT_DOUBLE_EQ(res[i], Ref[i]); }

	std::cout << " Res : [";
	for (size_t i = 0; i < N; ++i) { std::cout << res[i] << ", "; }
	std::cout << "]\n";
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, HasNaN)
{
	const auto Data = std::vector({std::numeric_limits<double>::quiet_NaN(), 2.0, 3.0});
	EXPECT_TRUE(matrix::HasNaN(Data));
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, HasInvalid)
{
	const auto Data = std::vector({std::numeric_limits<double>::infinity(), 2.0, 3.0});
	EXPECT_TRUE(matrix::HasInvalid(Data));
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, IsInvalid)
{
	const auto Data = std::vector({std::numeric_limits<double>::quiet_NaN(), 2.0, 3.0});
	EXPECT_FALSE(matrix::IsInvalid(Data));
	EXPECT_TRUE(matrix::IsInvalid(NAN_VECTOR));
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, GetValidIndex)
{
	const auto Data = std::vector({std::numeric_limits<double>::quiet_NaN(), 2.0, 3.0});
	auto valid      = matrix::GetValidIndex(Data);
	ASSERT_EQ(valid.size(), 3);
	EXPECT_FALSE(valid[0]);
	EXPECT_TRUE(valid[1]);
	EXPECT_TRUE(valid[2]);
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, CountValid)
{
	const auto Data = std::vector({std::numeric_limits<double>::quiet_NaN(), 2.0, 3.0, std::numeric_limits<double>::infinity()});
	EXPECT_EQ(matrix::CountValid(Data), 2);
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, MaxIgnoresNaN)
{
	const auto Data = std::vector({std::numeric_limits<double>::quiet_NaN(), 5.0, 3.0});
	EXPECT_DOUBLE_EQ(matrix::Max(Data), 5.0);
	const auto Data2 = std::vector({std::numeric_limits<double>::infinity(), 5.0, 3.0});
	EXPECT_DOUBLE_EQ(matrix::Max(Data2), 5.0);
	const auto Data3 = std::vector({std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()});
	EXPECT_DOUBLE_EQ(matrix::Max(Data3), std::numeric_limits<double>::infinity());
	EXPECT_TRUE(std::isnan(matrix::Max(NAN_VECTOR)));
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, SumValid)
{
	const auto Data = std::vector({std::numeric_limits<double>::quiet_NaN(), 2.0, 3.0});
	EXPECT_DOUBLE_EQ(matrix::SumValid(Data), 5.0);
	EXPECT_DOUBLE_EQ(matrix::SumValid(NAN_VECTOR), 0.0);
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, MeanValid)
{
	const auto Data = std::vector({std::numeric_limits<double>::quiet_NaN(), 4.0});
	EXPECT_DOUBLE_EQ(matrix::MeanValid(Data), 4.0);
	EXPECT_TRUE(std::isnan(matrix::MeanValid(NAN_VECTOR)));
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, MeanWithMask)
{
	const std::vector<double> Data{1, 2, 3, 4};
	const std::vector Mask{false, true, true, false};
	EXPECT_DOUBLE_EQ(matrix::MeanWithMask(Data, Mask), 2.5);
	EXPECT_TRUE(std::isnan(matrix::MeanWithMask(NAN_VECTOR, Mask)));
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, Median)
{
	EXPECT_DOUBLE_EQ(matrix::Median({ 1.0, 3.0, 2.0 }), 2.0);
	EXPECT_DOUBLE_EQ(matrix::Median({ 1.0, 2.0, 3.0, 4.0 }), 2.5);
	EXPECT_TRUE(std::isnan(matrix::Median(NAN_VECTOR)));
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, Quantile)
{
	EXPECT_DOUBLE_EQ(matrix::Quantile({ 1, 2, 3, 4 }, 0.75), 3.25);
	EXPECT_DOUBLE_EQ(matrix::Quantile({ 1, 2, 3, 4, 5 }, 0.5), 3.0);
	EXPECT_TRUE(std::isnan(matrix::Quantile(NAN_VECTOR, 0.5)));
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, RmsDifference)
{
	const std::vector<double> A{1, 2, 3};
	const std::vector<double> B{1, 4, 9};
	const double Expected = std::sqrt((0 + 4 + 36) / 3.0);
	EXPECT_DOUBLE_EQ(matrix::RmsDifference(A, B), Expected);
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, NormalisationWithShift)
{
	std::vector<double> vec(27);
	std::iota(vec.begin(), vec.end(), 0);
	matrix::Normalisation(vec, true);
	EXPECT_DOUBLE_EQ(vec.front(), 0.0);
	EXPECT_NEAR(std::accumulate(vec.begin(), vec.end(), 0.0), 1.0, EPSILON);
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, CropVolume)
{
	Vec3I size{.x = 3, .y = 3, .z = 3};
	std::vector<double> vol(27);
	std::iota(vol.begin(), vol.end(), 0);
	matrix::Crop(vol, size);
	EXPECT_EQ(size.x, 1);
	EXPECT_EQ(size.y, 1);
	EXPECT_EQ(size.z, 1);
	EXPECT_EQ(vol.size(), 1);
	EXPECT_EQ(vol[0], 13);
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, ExtractSubVolume)
{
	constexpr Vec3I Size{.x = 3, .y = 3, .z = 1};
	std::vector<double> data(9);
	std::iota(data.begin(), data.end(), 0);
	const auto Sub = matrix::ExtractSubVolume(data, Size, {1, 2}, {0, 1}, {0});
	EXPECT_EQ(Sub, (std::vector<double>{1, 2, 4, 5}));
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, CenterOfMass)
{
	std::vector<Pixel> pixels(9);
	for (size_t i = 0; i < 9; ++i) { pixels[i].v = static_cast<double>(i + 1); }
	double x, y;
	matrix::CenterOfMass(pixels, 3, x, y);
	EXPECT_NEAR(x, 1.14, 0.01);
	EXPECT_NEAR(y, 1.4, 0.01);
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, CholeskySuccess)
{
	const std::vector<double> A = {
		4, 12, -16,
		12, 37, -43,
		-16, -43, 98
	};
	std::vector<double> l(9), u(9);
	ASSERT_TRUE(matrix::Cholesky(A, 3, l, u));
	EXPECT_NEAR(l[0], 2.0, EPSILON);
	EXPECT_NEAR(u[8], 3.0, EPSILON);
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, CholeskyFail)
{
	const std::vector<double> A = {
		0, 1, 2,
		1, 0, 3,
		2, 3, 0
	};
	std::vector<double> l(9), u(9);
	EXPECT_FALSE(matrix::Cholesky(A, 3, l, u));
}

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, InverseGaussJordan)
{
	std::vector<double> a = {4, 7, 2, 6};
	std::vector<double> inv(4);
	matrix::InverseGaussJordan(a, inv, 2);
	EXPECT_NEAR(inv[0], 0.6, EPSILON);
	EXPECT_NEAR(inv[1], -0.7, EPSILON);
	EXPECT_NEAR(inv[2], -0.2, EPSILON);
	EXPECT_NEAR(inv[3], 0.4, EPSILON);
}

/// <summary>Vérifie A × inverse(A), y compris après permutation de plusieurs lignes.</summary>
TEST_P(GaussJordanPivotTests, RecoversIdentity)
{
	const auto& Original = GetParam();
	const size_t N       = Original.size() == 9 ? 3 : 2;
	auto matrix          = Original;
	std::vector<double> inverse(Original.size());
	ASSERT_NO_THROW(matrix::InverseGaussJordan(matrix, inverse, N));
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			double product = 0.0;
			for (size_t k = 0; k < N; ++k) { product += Original[(i * N) + k] * inverse[(k * N) + j]; }
			EXPECT_NEAR(product, i == j ? 1.0 : 0.0, EPSILON);
		}
	}
}

INSTANTIATE_TEST_SUITE_P(Matrices, GaussJordanPivotTests,
                         testing::Values(std::vector<double>{0, 1, 1, 0}, std::vector<double>{1e-20, 1, 1, 1},
	                         std::vector<double>{1e-100, 0, 0, 1e100}, std::vector<double>{0, 0, 2, 3, 0, 0, 0, 4, 0}),
                         ([](const testing::TestParamInfo<std::vector<double>>& info) {
	                         constexpr std::array Names{"ZeroPivot", "TinyPivot", "DifferentRowScales", "MultipleSwaps"};
	                         return Names[info.index];
	                         }));

/// <summary>L'échec efface l'inverse partiel au lieu de laisser des NaN, Inf ou coefficients trompeurs.</summary>
TEST_P(GaussJordanFailureTests, RejectsInvalidInverse)
{
	auto matrix = GetParam();
	std::vector inverse(4, 42.0);
	EXPECT_THROW(matrix::InverseGaussJordan(matrix, inverse, 2), std::domain_error);
	for (const double Value : inverse) { EXPECT_DOUBLE_EQ(Value, 0.0); }
}

INSTANTIATE_TEST_SUITE_P(Matrices, GaussJordanFailureTests,
                         testing::Values(std::vector<double>{1, 2, 2, 4}, std::vector<double>{1, 0, 0, 0},
	                         std::vector<double>{1, 1, 1, std::nextafter(1.0, 2.0)},
	                         std::vector<double>{NAN_VECTOR[0], 0, 0, 1},
	                         std::vector<double>{std::numeric_limits<double>::infinity(), 0, 0, 1},
	                         std::vector<double>{std::numeric_limits<double>::denorm_min(), 0, 0, 1}),
                         ([](const testing::TestParamInfo<std::vector<double>>& info) {
	                         constexpr std::array Names{
	                         "DependentRows", "ZeroRow", "NumericallySingular", "NonFiniteNaN", "NonFiniteInfinity", "InverseOverflow"
	                         };
	                         return Names[info.index];
	                         }));

/// <summary>Les buffers et la dimension sont validés avant tout accès ou mutation.</summary>
TEST_P(GaussJordanArgumentTests, RejectsInvalidDimensions)
{
	std::vector<double> matrix{1, 0, 0, 1}, inverse(4);
	const size_t N = GetParam() == 0 ? 0 : (GetParam() == 1 ? 3 : 2);
	if (GetParam() == 2) { inverse.resize(3); }
	auto& output = GetParam() == 3 ? matrix : inverse;
	EXPECT_THROW(matrix::InverseGaussJordan(matrix, output, N), std::invalid_argument);
}

INSTANTIATE_TEST_SUITE_P(Arguments, GaussJordanArgumentTests, testing::Values(size_t{0}, size_t{1}, size_t{2}, size_t{3}),
                         ([](const testing::TestParamInfo<size_t>& info) {
	                         constexpr std::array Names{"ZeroDimension", "InputSizeMismatch", "OutputSizeMismatch", "AliasedBuffers"};
	                         return Names[info.param];
	                         }));

//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, InverseLU)
{
	std::vector<double> a = {4, 7, 2, 6};
	std::vector<double> inv(4);
	matrix::InverseLU(a, inv, 2);
	EXPECT_NEAR(inv[0], 0.6, EPSILON);
	EXPECT_NEAR(inv[1], -0.7, EPSILON);
	EXPECT_NEAR(inv[2], -0.2, EPSILON);
	EXPECT_NEAR(inv[3], 0.4, EPSILON);
}
} // namespace palm_tracer::signal_processing
