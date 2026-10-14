#include "TestUtils.hpp"
#include "localization/Segmentation2D.hpp"

#include <array>
#include <cmath>
#include <numeric>

namespace palm_tracer::localization {
namespace {
/// <summary>Injecte un potentiel connu pour isoler la segmentation de la decomposition en ondelettes.</summary>
Segmentation2D Segment(const Vec2I& dim, const std::vector<double>& potential, const double threshold = 0, const double ratio = 0)
{
	const std::vector<uint16_t> Image(dim.GetSize(), 0);
	ATrous decomposition(Image.data(), dim);
	decomposition.wDecomp_[PLANE_ID] = potential;
	return Segmentation2D(decomposition, threshold, ratio);
}

/// <summary>Scenarios independants de mesure des regions.</summary>
class SegmentationStatisticsTests : public testing::TestWithParam<bool> { };

/// <summary>Seuils incluant les valeurs exactement egales au potentiel.</summary>
class SegmentationThresholdTests : public testing::TestWithParam<double> { };

/// <summary>Nombre de pixels conserves autour de la surface minimale.</summary>
class SegmentationSurfaceTests : public testing::TestWithParam<size_t> { };

/// <summary>Fusion facultative de deux bassins relies par un col.</summary>
class SegmentationWatershedTests : public testing::TestWithParam<bool> { };
}

#pragma region ComputeStatistics
/// <summary>Verifie surface, intensite, centroide pondere et ecart-type de population.</summary>
TEST_P(SegmentationStatisticsTests, RegionMeasurements)
{
	constexpr Vec2I Dim{.x = 7, .y = 5};
	std::vector<double> potential(Dim.GetSize(), 0);
	double sum = 0, weightedX = 0, weightedY = 0;
	std::vector<double> values;
	for (size_t y = 1; y <= 3; ++y) {
		for (size_t x = 2; x <= 4; ++x) {
			const double Value       = GetParam() ? 12.0 : static_cast<double>(1 + (y - 1) * 3 + x - 2);
			potential[y * Dim.x + x] = Value;
			values.push_back(Value);
			sum += Value;
			weightedX += Value * static_cast<double>(x);
			weightedY += Value * static_cast<double>(y);
		}
	}
	const double Mean = sum / static_cast<double>(values.size());
	double variance   = 0;
	for (const double Value : values) { variance += (Value - Mean) * (Value - Mean); }
	const auto Infos = Segment(Dim, potential).Infos();
	ASSERT_EQ(Infos.size(), 1);
	EXPECT_DOUBLE_EQ(Infos[0].surf_, 9);
	EXPECT_DOUBLE_EQ(Infos[0].intensity_, sum);
	EXPECT_DOUBLE_EQ(Infos[0].mean_, Mean);
	EXPECT_NEAR(Infos[0].x_, weightedX / sum, 1e-12);
	EXPECT_NEAR(Infos[0].y_, weightedY / sum, 1e-12);
	EXPECT_NEAR(Infos[0].stdev_, std::sqrt(variance / 9), 1e-12);
}

INSTANTIATE_TEST_SUITE_P(Fields, SegmentationStatisticsTests, testing::Bool(),
                         [](const testing::TestParamInfo<bool>& info) { return info.param ? "ConstantPlateau" : "IntensityRamp"; });
#pragma endregion ComputeStatistics

#pragma region FilterData
/// <summary>Le seuil est strict : un plateau egal au seuil ne produit aucune region.</summary>
TEST_P(SegmentationThresholdTests, StrictThreshold)
{
	constexpr Vec2I Dim{.x = 7, .y = 5};
	std::vector<double> potential(Dim.GetSize(), 0);
	for (size_t y = 1; y <= 3; ++y) { for (size_t x = 2; x <= 4; ++x) { potential[y * Dim.x + x] = 10; } }
	EXPECT_EQ(Segment(Dim, potential, GetParam()).Infos().size(), GetParam() < 10 ? size_t{1} : size_t{0});
}

INSTANTIATE_TEST_SUITE_P(Thresholds, SegmentationThresholdTests, testing::Values(9.0, 10.0, 11.0),
                         [](const testing::TestParamInfo<double>& info) { return info.param < 10 ? "BelowPlateau" : info.param == 10 ? "EqualPlateau" :
                         "AbovePlateau"; });

/// <summary>Le filtre rejette les surfaces inferieures ou egales a VOLUME_MIN.</summary>
TEST_P(SegmentationSurfaceTests, MinimumSurface)
{
	constexpr Vec2I Dim{.x = 9, .y = 5};
	std::vector<double> potential(Dim.GetSize(), 0);
	for (size_t x = 1; x <= GetParam(); ++x) { potential[2 * Dim.x + x] = 10; }
	const auto Infos = Segment(Dim, potential).Infos();
	ASSERT_EQ(Infos.size(), GetParam() > VOLUME_MIN ? size_t{1} : size_t{0});
	if (!Infos.empty()) { EXPECT_DOUBLE_EQ(Infos[0].surf_, static_cast<double>(GetParam())); }
}

INSTANTIATE_TEST_SUITE_P(Surfaces, SegmentationSurfaceTests,
                         testing::Values(VOLUME_MIN - 1, VOLUME_MIN, VOLUME_MIN + 1),
                         [](const testing::TestParamInfo<size_t>& info) { return info.param < VOLUME_MIN ? "BelowMinimum" : info.param == VOLUME_MIN ?
                         "AtMinimum" : "AboveMinimum"; });
#pragma endregion FilterData

#pragma region WatershedIntelligent
/// <summary>Deux pics relies restent separes ou fusionnent selon le ratio, sans perdre de pixels.</summary>
TEST_P(SegmentationWatershedTests, ConnectedPeaks)
{
	constexpr Vec2I Dim{.x = 13, .y = 7};
	std::vector<double> potential(Dim.GetSize(), 0);
	// Le pic le plus fort est le second label : le premier bassin peut fusionner vers ce voisin.
	for (size_t y = 2; y <= 4; ++y) {
		for (size_t x = 2; x <= 4; ++x) { potential[y * Dim.x + x] = 60; }
		for (size_t x = 8; x <= 10; ++x) { potential[y * Dim.x + x] = 100; }
	}
	// Un col non plat evite de creer un troisieme germe au milieu du pont.
	potential[3 * Dim.x + 5] = 30;
	potential[3 * Dim.x + 6] = 20;
	potential[3 * Dim.x + 7] = 40;
	const auto Infos         = Segment(Dim, potential, 0, GetParam() ? 10 : 0).Infos();
	ASSERT_EQ(Infos.size(), GetParam() ? size_t{1} : size_t{2});
	double surface = 0, intensity = 0;
	for (const auto& Info : Infos) {
		surface += Info.surf_;
		intensity += Info.intensity_;
	}
	EXPECT_DOUBLE_EQ(surface, 21);
	EXPECT_DOUBLE_EQ(intensity, std::accumulate(potential.begin(), potential.end(), 0.0));
}

INSTANTIATE_TEST_SUITE_P(Modes, SegmentationWatershedTests, testing::Bool(),
                         [](const testing::TestParamInfo<bool>& info) { return info.param ? "Merged" : "Separated"; });
#pragma endregion WatershedIntelligent
} // namespace palm_tracer::localization
