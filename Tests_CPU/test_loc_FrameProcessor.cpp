#include "TestUtils.hpp"
#include "TiffManagement.hpp"
#include "localization/FrameProcessor.hpp"

#include <algorithm>
#include <array>
#include <compare>
#include <limits>
#include <tuple>

namespace palm_tracer::localization {
namespace {
constexpr std::array Params{7.0, 1.0, 2.0, 0.0}; ///< Paramètres communs des ajustements.
using ProcessingParams = std::tuple<size_t, bool>;

/// <summary>Modes d'ajustement et variantes de fusion disposant de données de référence.</summary>
class FrameProcessorTests : public testing::TestWithParam<ProcessingParams> { };

/// <summary>Compare toutes les colonnes de résultats obtenus après réutilisation.</summary>
void ExpectLocalizationsEqual(const FrameLocalizations& actual, const FrameLocalizations& expected)
{
	ASSERT_EQ(actual.size(), expected.size());
	for (size_t row = 0; row < expected.size(); ++row) {
		SCOPED_TRACE("Localisation " + std::to_string(row));
		for (size_t column = 0; column < localization_columns::Count; ++column) {
			SCOPED_TRACE("Colonne " + std::to_string(column));
			EXPECT_DOUBLE_EQ(actual[row][column], expected[row][column]);
		}
	}
}
} // namespace

/// <summary>Réutilise le même processeur et le même résultat après changements d'image et de seuil.</summary>
TEST_P(FrameProcessorTests, ReusesWorkspaceAcrossImagesAndThresholds)
{
	const auto [Fit, Watershed] = GetParam();
	std::vector<uint16_t> pixels;
	Vec3I dim{};
	LoadTiff3D((InputPath / "stack.tif").string().c_str(), pixels, dim);
	ASSERT_GE(dim.z, 2);
	const Vec2I Dim{.x = dim.x, .y = dim.y};
	const auto Size = dim.GetPlaneSize();
	const std::vector<uint16_t> Black(Size, 0);
	const double Ratio    = Watershed ? 0.0 : 10.0;
	const size_t Capacity = Size;
	FrameProcessor processor(Dim, Capacity, Ratio, Fit, Params.data());
	FrameLocalizations points;
	processor.Process(pixels.data(), 103.6, points);
	ASSERT_FALSE(points.empty());
	const auto First = points;

	processor.Process(Black.data(), 103.6, points);
	EXPECT_TRUE(points.empty());

	processor.Process(pixels.data() + Size, 103.6, points);
	FrameProcessor fresh(Dim, Capacity, Ratio, Fit, Params.data());
	FrameLocalizations expected;
	fresh.Process(pixels.data() + Size, 103.6, expected);
	ASSERT_FALSE(expected.empty());
	ExpectLocalizationsEqual(points, expected);

	processor.Process(pixels.data(), std::numeric_limits<double>::max(), points);
	EXPECT_TRUE(points.empty());
	processor.Process(pixels.data(), 103.6, points);
	ExpectLocalizationsEqual(points, First);
}

INSTANTIATE_TEST_SUITE_P(, FrameProcessorTests,
                           testing::Combine(testing::Values(size_t{0}, size_t{1}, size_t{2}, size_t{3}, size_t{4}), testing::Bool()),
                           ([](const testing::TestParamInfo<ProcessingParams>& info) {
	                           constexpr std::array Names{"NoFit", "GaussianXY", "GaussianSigma", "GaussianSigmaXY", "GaussianTheta"};
	                           return std::string(Names[std::get<0>(info.param)]) + (std::get<1>(info.param) ? "WithWatershed" : "WithoutWatershed");
	                           }));

/// <summary>Capacités indépendantes du tampon de localisations d'une image.</summary>
class FrameCapacityTests : public testing::TestWithParam<size_t> { };

/// <summary>Limite la sortie, remplace les valeurs antérieures et conserve le tri Y/X.</summary>
TEST_P(FrameCapacityTests, LimitsAndSortsOutput)
{
	std::vector<uint16_t> pixels;
	Vec3I dim{};
	LoadTiff3D((InputPath / "stack.tif").string().c_str(), pixels, dim);
	const size_t Capacity = GetParam();
	FrameProcessor processor({.x = dim.x, .y = dim.y}, Capacity, 0.0, 0, Params.data());
	FrameLocalizations points(Capacity + 1);
	processor.Process(pixels.data(), 103.6, points);
	ASSERT_EQ(points.size(), Capacity);
	EXPECT_TRUE(std::ranges::is_sorted(points, [](const auto& a, const auto& b) {
		if (const auto Order = std::weak_order(a[localization_columns::Y], b[localization_columns::Y]); Order != 0) { return Order < 0; }
		return std::weak_order(a[localization_columns::X], b[localization_columns::X]) < 0;
		}));
}

INSTANTIATE_TEST_SUITE_P(, FrameCapacityTests, testing::Values(size_t{0}, size_t{1}, size_t{4}),
                           [](const testing::TestParamInfo<size_t>& info) { return "Capacity" + std::to_string(info.param); });
} // namespace palm_tracer::localization
