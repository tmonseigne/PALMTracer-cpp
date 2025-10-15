#include "gtest/gtest.h"

#include <numeric>
#include "signal_processing/Matrix.hpp"


namespace palm_tracer {
constexpr double EPSILON = 1e-10;
const std::vector NAN_VECTOR(3, std::numeric_limits<double>::quiet_NaN());


//-------------------------------------------------------------------------------------------------
TEST(MatrixTest, SendIrregular)
{
	const std::vector<std::vector<double>> Data = {{1, 2}, {3, 4, 5}, {}};
	const std::vector<double> Ref               = {3, 1, 1, 2, 4, 2, 3, 4, 5, -1};

	std::vector<double> res(Ref.size(), 0.0);
	const size_t N = matrix::SendIrregularArray(Data, res.data());

	ASSERT_EQ(N, Ref.size()) << "Nombre d'elements ecrits inattendu";
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
} // namespace palm_tracer
