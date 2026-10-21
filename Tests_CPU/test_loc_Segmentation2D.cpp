#include "TestUtils.hpp"
#include "localization/Segmentation2D.hpp"

#include <array>
#include <cmath>
#include <numeric>
#include <tuple>
#include <type_traits>

namespace palm_tracer::localization {
namespace {
/// <summary>Injecte un potentiel connu pour isoler la segmentation de la décomposition en ondelettes.</summary>
Segmentation2D Segment(const Vec2I& dim, const std::vector<double>& potential, const double threshold = 0, const double ratio = 0)
{
	const std::vector<uint16_t> Image(dim.GetSize(), 0);
	ATrous decomposition(Image.data(), dim);
	decomposition.wDecomp_[PLANE_ID] = potential;
	return Segmentation2D(decomposition, threshold, ratio);
}

/// <summary>Scénarios indépendants de mesure des régions.</summary>
class SegmentationStatisticsTests : public testing::TestWithParam<bool> { };

/// <summary>Seuils incluant les valeurs exactement égales au potentiel.</summary>
class SegmentationThresholdTests : public testing::TestWithParam<double> { };

/// <summary>Nombre de pixels conservés autour de la surface minimale.</summary>
class SegmentationSurfaceTests : public testing::TestWithParam<size_t> { };

/// <summary>Fusion facultative de deux bassins reliés par un col.</summary>
class SegmentationWatershedTests : public testing::TestWithParam<std::tuple<bool, bool>> { };
} // namespace

#pragma region Process
/// <summary>Les transitions dense/vide/dense et les nouveaux seuils et ratios remplacent tout l'état précédent.</summary>
TEST_P(SegmentationWatershedTests, ReusesWorkspaceAcrossFrames)
{
	const auto [merge, firstPeakHigher] = GetParam();
	constexpr Vec2I Dim{.x = 13, .y = 7};
	const std::vector<uint16_t> Image(Dim.GetSize(), 0);
	ATrous decomposition(Image.data(), Dim);
	Segmentation2D segmentation(Dim);
	for (size_t frame = 0; frame < 5; ++frame) {
		SCOPED_TRACE(frame);
		auto& potential = decomposition.wDecomp_[PLANE_ID];
		std::ranges::fill(potential, 0);
		if (frame != 1) {
			for (size_t y = 2; y <= 4; ++y) {
				for (size_t x = 2; x <= 4; ++x) { potential[(y * Dim.x) + x] = firstPeakHigher ? 100 : 60; }
				for (size_t x = 8; x <= 10; ++x) { potential[(y * Dim.x) + x] = firstPeakHigher ? 60 : 100; }
			}
			potential[(3 * Dim.x) + 5] = firstPeakHigher ? 40 : 30;
			potential[(3 * Dim.x) + 6] = 20;
			potential[(3 * Dim.x) + 7] = firstPeakHigher ? 30 : 40;
		}
		const double Threshold = frame == 3 ? 100 : 0;
		const double Ratio = (frame == 4 ? !merge : merge) ? 10 : 0;
		const Segmentation2D Expected(decomposition, Threshold, Ratio);
		segmentation.Process(decomposition, Threshold, Ratio);
		std::ranges::fill(potential, -999);
		ASSERT_EQ(segmentation.Infos().size(), Expected.Infos().size());
		for (size_t i = 0; i < Expected.Infos().size(); ++i) {
			const auto& Actual = segmentation.Infos()[i];
			const auto& Reference = Expected.Infos()[i];
			EXPECT_DOUBLE_EQ(Actual.x_, Reference.x_);
			EXPECT_DOUBLE_EQ(Actual.y_, Reference.y_);
			EXPECT_DOUBLE_EQ(Actual.z_, Reference.z_);
			EXPECT_DOUBLE_EQ(Actual.intensity_, Reference.intensity_);
			EXPECT_DOUBLE_EQ(Actual.surf_, Reference.surf_);
			EXPECT_DOUBLE_EQ(Actual.mean_, Reference.mean_);
			EXPECT_DOUBLE_EQ(Actual.stdev_, Reference.stdev_);
		}
	}
}

/// <summary>Une nouvelle décomposition redimensionne les cartes sans conserver les anciennes régions.</summary>
TEST(SegmentationTests, ReconfiguresDimensions)
{
	Segmentation2D segmentation(Vec2I{.x = 2, .y = 2});
	for (const Vec2I Dim : {Vec2I{.x = 13, .y = 7}, Vec2I{.x = 7, .y = 13}, Vec2I{.x = 3, .y = 3}}) {
		const std::vector<uint16_t> Image(Dim.GetSize(), 0);
		ATrous decomposition(Image.data(), Dim);
		std::ranges::fill(decomposition.wDecomp_[PLANE_ID], 10);
		segmentation.Process(decomposition, 0, 0);
		ASSERT_EQ(segmentation.Infos().size(), 1);
		EXPECT_DOUBLE_EQ(segmentation.Infos()[0].surf_, static_cast<double>(Dim.GetSize()));
		EXPECT_DOUBLE_EQ(segmentation.Infos()[0].intensity_, 10 * static_cast<double>(Dim.GetSize()));
		EXPECT_DOUBLE_EQ(segmentation.Infos()[0].x_, static_cast<double>(Dim.x - 1) / 2);
		EXPECT_DOUBLE_EQ(segmentation.Infos()[0].y_, static_cast<double>(Dim.y - 1) / 2);
	}
}
#pragma endregion Process

#pragma region SegmentationInfo
/// <summary>Une région commence vide et Reset restaure cet état après utilisation.</summary>
TEST(SegmentationInfoTests, DefaultStateAndReset)
{
	SegmentationInfo info;
	EXPECT_DOUBLE_EQ(info.x_, 0);
	EXPECT_DOUBLE_EQ(info.y_, 0);
	EXPECT_DOUBLE_EQ(info.z_, 0);
	EXPECT_DOUBLE_EQ(info.intensity_, 0);
	EXPECT_DOUBLE_EQ(info.surf_, 0);
	EXPECT_DOUBLE_EQ(info.mean_, 0);
	EXPECT_DOUBLE_EQ(info.stdev_, 0);
	info.x_ = info.y_ = info.z_ = info.intensity_ = info.surf_ = info.mean_ = info.stdev_ = 42;
	info.Reset();
	EXPECT_DOUBLE_EQ(info.x_, 0);
	EXPECT_DOUBLE_EQ(info.y_, 0);
	EXPECT_DOUBLE_EQ(info.z_, 0);
	EXPECT_DOUBLE_EQ(info.intensity_, 0);
	EXPECT_DOUBLE_EQ(info.surf_, 0);
	EXPECT_DOUBLE_EQ(info.mean_, 0);
	EXPECT_DOUBLE_EQ(info.stdev_, 0);
}
#pragma endregion SegmentationInfo

#pragma region ComputeStatistics
/// <summary>Les labels restent uniques au-dela de la capacité historique des entiers 16 bits.</summary>
TEST(SegmentationTests, ManyRegions)
{
	constexpr Vec2I Dim{.x = 1025, .y = 1025};
	std::vector<double> potential(Dim.GetSize(), 0);
	for (size_t y = 1; y < 1024; y += 4) {
		for (size_t x = 1; x < 1024; x += 4) {
			for (size_t dy = 0; dy < 3; ++dy) { for (size_t dx = 0; dx < 3; ++dx) { potential[((y + dy) * Dim.x) + x + dx] = 100; } }
		}
	}
	const auto Infos = Segment(Dim, potential).Infos();
	ASSERT_EQ(Infos.size(), 65536);
	for (const auto& Info : Infos) {
		EXPECT_DOUBLE_EQ(Info.surf_, 9);
		EXPECT_DOUBLE_EQ(Info.intensity_, 900);
	}
}

/// <summary>Vérifie surface, intensité, centroïde pondéré et écart-type de population.</summary>
TEST_P(SegmentationStatisticsTests, RegionMeasurements)
{
	constexpr Vec2I Dim{.x = 7, .y = 5};
	std::vector<double> potential(Dim.GetSize(), 0);
	double sum = 0, weightedX = 0, weightedY = 0;
	std::vector<double> values;
	for (size_t y = 1; y <= 3; ++y) {
		for (size_t x = 2; x <= 4; ++x) {
			const double Value         = GetParam() ? 12.0 : static_cast<double>(1 + ((y - 1) * 3) + x - 2);
			potential[(y * Dim.x) + x] = Value;
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
/// <summary>Le seuil est strict : un plateau égal au seuil ne produit aucune région.</summary>
TEST_P(SegmentationThresholdTests, StrictThreshold)
{
	constexpr Vec2I Dim{.x = 7, .y = 5};
	std::vector<double> potential(Dim.GetSize(), 0);
	for (size_t y = 1; y <= 3; ++y) { for (size_t x = 2; x <= 4; ++x) { potential[(y * Dim.x) + x] = 10; } }
	EXPECT_EQ(Segment(Dim, potential, GetParam()).Infos().size(), GetParam() < 10 ? size_t{1} : size_t{0});
}

INSTANTIATE_TEST_SUITE_P(Thresholds, SegmentationThresholdTests, testing::Values(9.0, 10.0, 11.0),
                         [](const testing::TestParamInfo<double>& info) { return info.param < 10 ? "BelowPlateau" : info.param == 10 ? "EqualPlateau" :
                         "AbovePlateau"; });

/// <summary>Le filtre rejette les surfaces inférieures ou égales à VOLUME_MIN.</summary>
TEST_P(SegmentationSurfaceTests, MinimumSurface)
{
	constexpr Vec2I Dim{.x = 9, .y = 5};
	std::vector<double> potential(Dim.GetSize(), 0);
	for (size_t x = 1; x <= GetParam(); ++x) { potential[(2 * Dim.x) + x] = 10; }
	const auto Infos = Segment(Dim, potential).Infos();
	ASSERT_EQ(Infos.size(), GetParam() > VOLUME_MIN ? size_t{1} : size_t{0});
	if (!Infos.empty()) { EXPECT_DOUBLE_EQ(Infos[0].surf_, static_cast<double>(GetParam())); }
}

INSTANTIATE_TEST_SUITE_P(Surfaces, SegmentationSurfaceTests,
                         testing::Values(VOLUME_MIN - 1, VOLUME_MIN, VOLUME_MIN + 1),
                         [](const testing::TestParamInfo<size_t>& info) { return info.param < VOLUME_MIN ? "BelowMinimum" : info.param == VOLUME_MIN ?
                         "AtMinimum" : "AboveMinimum"; });
#pragma endregion FilterData

#pragma region Infos
/// <summary>La référence désigne les mesures originales et une copie explicite conserve son indépendance.</summary>
TEST(SegmentationInfosTests, ReadOnlyReferenceAndIndependentSnapshot)
{
	constexpr Vec2I Dim{.x = 7, .y = 7};
	std::vector<double> potential(Dim.GetSize(), 0);
	for (size_t y = 2; y <= 4; ++y) { for (size_t x = 2; x <= 4; ++x) { potential[(y * Dim.x) + x] = 100; } }
	const auto Segmentation = Segment(Dim, potential);
	const auto& Infos       = Segmentation.Infos();
	static_assert(std::is_same_v<decltype(Segmentation.Infos()), const std::vector<SegmentationInfo>&>);
	ASSERT_EQ(Infos.size(), 1);
	EXPECT_EQ(Infos.data(), Segmentation.Infos().data());
	auto snapshot = Segmentation.Infos();
	ASSERT_EQ(snapshot.size(), Infos.size());
	EXPECT_NE(snapshot.data(), Infos.data());
	snapshot[0].intensity_ = -42;
	EXPECT_DOUBLE_EQ(Infos[0].intensity_, 900);
	snapshot.clear();
	EXPECT_EQ(Infos.size(), 1);
	EXPECT_DOUBLE_EQ(Segmentation.Infos()[0].intensity_, 900);
}

/// <summary>Une segmentation sans région fournit une référence vers un vecteur vide.</summary>
TEST(SegmentationInfosTests, EmptyReadOnlyReference)
{
	constexpr Vec2I Dim{.x = 7, .y = 7};
	const auto Segmentation = Segment(Dim, std::vector<double>(Dim.GetSize(), 0));
	EXPECT_TRUE(Segmentation.Infos().empty());
	EXPECT_TRUE(Segmentation.Infos().empty());
}
#pragma endregion Infos

#pragma region WatershedIntelligent
/// <summary>Deux pics reliés restent séparés ou fusionnent selon le ratio, sans perdre de pixels.</summary>
TEST_P(SegmentationWatershedTests, ConnectedPeaks)
{
	const auto [merge, firstPeakHigher] = GetParam();
	constexpr Vec2I Dim{.x = 13, .y = 7};
	std::vector<double> potential(Dim.GetSize(), 0);
	// Le bassin dominant doit absorber son voisin même s'il porte la première étiquette.
	for (size_t y = 2; y <= 4; ++y) {
		for (size_t x = 2; x <= 4; ++x) { potential[(y * Dim.x) + x] = firstPeakHigher ? 100 : 60; }
		for (size_t x = 8; x <= 10; ++x) { potential[(y * Dim.x) + x] = firstPeakHigher ? 60 : 100; }
	}
	// Un col non plat évite de créer un troisième germe au milieu du pont.
	potential[(3 * Dim.x) + 5] = firstPeakHigher ? 40 : 30;
	potential[(3 * Dim.x) + 6] = 20;
	potential[(3 * Dim.x) + 7] = firstPeakHigher ? 30 : 40;
	const auto Infos           = Segment(Dim, potential, 0, merge ? 10 : 0).Infos();
	ASSERT_EQ(Infos.size(), merge ? size_t{1} : size_t{2});
	double surface = 0, intensity = 0;
	for (const auto& Info : Infos) {
		surface += Info.surf_;
		intensity += Info.intensity_;
	}
	EXPECT_DOUBLE_EQ(surface, 21);
	EXPECT_DOUBLE_EQ(intensity, std::accumulate(potential.begin(), potential.end(), 0.0));
}

INSTANTIATE_TEST_SUITE_P(Modes, SegmentationWatershedTests, testing::Combine(testing::Bool(), testing::Bool()),
                         [](const testing::TestParamInfo<std::tuple<bool, bool>>& info) {
                         return std::string(std::get<0>(info.param) ? "Merged" : "Separated")
                         + (std::get<1>(info.param) ? "FirstPeakHigher" : "SecondPeakHigher");
                         });

/// <summary>Une région absorbée ne recrée pas un germe lors du réétiquetage des passages suivants.</summary>
TEST(SegmentationTests, AbsorbedPeakDoesNotReturn)
{
	constexpr Vec2I Dim{.x = 18, .y = 7};
	std::vector<double> potential(Dim.GetSize(), 0);
	for (size_t y = 2; y <= 4; ++y) {
		for (size_t x = 1; x <= 3; ++x) { potential[(y * Dim.x) + x] = 200; }
		for (size_t x = 7; x <= 9; ++x) { potential[(y * Dim.x) + x] = 100; }
		for (size_t x = 13; x <= 15; ++x) { potential[(y * Dim.x) + x] = 60; }
	}
	potential[(3 * Dim.x) + 10] = 40;
	potential[(3 * Dim.x) + 11] = 20;
	potential[(3 * Dim.x) + 12] = 30;
	// Le premier bassin est isolé : la fusion des deux suivants ne dépend pas du défaut de sentinelle 0.
	const auto Infos = Segment(Dim, potential, 0, 10).Infos();
	ASSERT_EQ(Infos.size(), 2);
	EXPECT_DOUBLE_EQ(Infos[0].surf_, 9);
	EXPECT_DOUBLE_EQ(Infos[0].intensity_, 1800);
	EXPECT_DOUBLE_EQ(Infos[1].surf_, 21);
	EXPECT_DOUBLE_EQ(Infos[1].intensity_, 1530);
}
#pragma endregion WatershedIntelligent
} // namespace palm_tracer::localization
