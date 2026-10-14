#include "TestUtils.hpp"
#include "localization/ATrous.hpp"

#include <array>
#include <numeric>

namespace palm_tracer::localization {
namespace {
/// <summary>Champs constants dont la décomposition doit être exactement nulle hors résidu.</summary>
class ATrousConstantTests : public testing::TestWithParam<uint16_t> { };

/// <summary>Images non carrées, minces et réduites pour vérifier aussi le traitement des bords.</summary>
struct ImageCase {
	const char* name; ///< Nom du scénario.
	Vec2I dim;        ///< Dimensions de l'image.
};

/// <summary>Conserve un nom lisible dans les diagnostics des tests paramétrés.</summary>
void PrintTo(const ImageCase& value, std::ostream* stream) { *stream << value.name; }

/// <summary>Scénarios indépendants de reconstruction et de convolution.</summary>
class ATrousImageTests : public testing::TestWithParam<ImageCase> { };
}

#pragma region Decomposition
/// <summary>La somme de tous les plans, résidu compris, reconstruit chaque pixel de l'image.</summary>
TEST_P(ATrousImageTests, ReconstructsInput)
{
	const Vec2I Dim = GetParam().dim;
	std::vector<uint16_t> input(Dim.GetSize());
	for (size_t i = 0; i < input.size(); ++i) { input[i] = static_cast<uint16_t>((i * 37 + 11) % 4096); }
	const ATrous Decomposition(input.data(), Dim);
	std::vector<double> reconstructed(input.size(), 0);
	for (const auto& Plane : Decomposition.wDecomp_) {
		ASSERT_EQ(Plane.size(), input.size());
		for (size_t i = 0; i < input.size(); ++i) { reconstructed[i] += Plane[i]; }
	}
	for (size_t i = 0; i < input.size(); ++i) { EXPECT_NEAR(reconstructed[i], input[i], 1e-10) << "Pixel " << i; }
}

/// <summary>Compare le premier plan à une convolution B3 directe en deux dimensions, avec prolongement constant aux bords.</summary>
TEST_P(ATrousImageTests, FirstLevelMatchesDirectConvolution)
{
	constexpr std::array<double, 5> Kernel{1.0 / 16, 4.0 / 16, 6.0 / 16, 4.0 / 16, 1.0 / 16};
	const Vec2I Dim = GetParam().dim;
	std::vector<uint16_t> input(Dim.GetSize());
	for (size_t i = 0; i < input.size(); ++i) { input[i] = static_cast<uint16_t>((i * 37 + 11) % 4096); }
	const ATrous Decomposition(input.data(), Dim);
	for (size_t y = 0; y < Dim.y; ++y) {
		for (size_t x = 0; x < Dim.x; ++x) {
			double smoothed = 0;
			for (int dy = -2; dy <= 2; ++dy) {
				for (int dx = -2; dx <= 2; ++dx) {
					const auto X = static_cast<size_t>(std::clamp(static_cast<int>(x) + dx, 0, static_cast<int>(Dim.x) - 1));
					const auto Y = static_cast<size_t>(std::clamp(static_cast<int>(y) + dy, 0, static_cast<int>(Dim.y) - 1));
					smoothed += Kernel[static_cast<size_t>(dx + 2)] * Kernel[static_cast<size_t>(dy + 2)] * input[Y * Dim.x + X];
				}
			}
			EXPECT_NEAR(Decomposition.wDecomp_[0][y * Dim.x + x], input[y * Dim.x + x] - smoothed, 1e-10);
		}
	}
}

INSTANTIATE_TEST_SUITE_P(Shapes, ATrousImageTests,
                         testing::Values(ImageCase{"SinglePixel", {.x = 1, .y = 1}}, ImageCase{"SingleRow", {.x = 9, .y = 1}},
	                         ImageCase{"SingleColumn", {.x = 1, .y = 9}}, ImageCase{"Rectangle", {.x = 13, .y = 7}}),
                         [](const testing::TestParamInfo<ImageCase>& info) { return info.param.name; });

/// <summary>Un champ constant produit uniquement un résidu constant, même pour les intensités limites uint16.</summary>
TEST_P(ATrousConstantTests, DetailsVanishAndResidualIsPreserved)
{
	constexpr Vec2I Dim{.x = 11, .y = 7};
	const std::vector<uint16_t> Input(Dim.GetSize(), GetParam());
	const ATrous Decomposition(Input.data(), Dim);
	for (size_t level = 0; level < N_PLANES; ++level) {
		const double Expected = level + 1 == N_PLANES ? GetParam() : 0;
		ASSERT_EQ(Decomposition.wDecomp_[level].size(), Input.size());
		for (const double Value : Decomposition.wDecomp_[level]) { EXPECT_DOUBLE_EQ(Value, Expected); }
	}
}

INSTANTIATE_TEST_SUITE_P(Intensities, ATrousConstantTests, testing::Values(uint16_t{0}, uint16_t{1}, uint16_t{1024}, uint16_t{65535}),
                         [](const testing::TestParamInfo<uint16_t>& info) { return "Intensity" + std::to_string(info.param); });
#pragma endregion Decomposition
} // namespace palm_tracer::localization
