#include <gtest/gtest-param-test.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <numeric>
#include <tuple>

#include "TestUtils.hpp"
#include "TiffManagement.hpp"
#include "signal_processing/Alignment.hpp"

namespace palm_tracer::signal_processing {
const std::filesystem::path PATH          = std::filesystem::path(__FILE__).parent_path();
const std::filesystem::path REF_FILE_PATH = PATH / "Input/Alignement_Reference.tiff";
const std::filesystem::path OUTPUT_PATH   = PATH / "Output";

//-------------------------------------------------------------------------------------------------
struct AlignTestParams {
	std::string outputPath;
	std::array<std::array<double, 10>, 2> factors;
	int upsampling;
};

/// <summary>Affiche le fichier du scénario au lieu du contenu binaire des coefficients.</summary>
void PrintTo(const AlignTestParams& value, std::ostream* stream) { *stream << std::filesystem::path(value.outputPath).stem().string(); }

//-------------------------------------------------------------------------------------------------
class AlignmentTests : public testing::TestWithParam<AlignTestParams> { };

/// <summary>Translations hors champ, indépendantes selon l'axe et le sens.</summary>
class AlignmentOutsideTests : public testing::TestWithParam<size_t> { };

//-------------------------------------------------------------------------------------------------
TEST(AlignmentTests, BlankCallWithReference)
{
	Alignment align;

	constexpr size_t Height     = 100, Width = 100, Planes = 10;
	constexpr size_t ImSize     = Height * Width;
	constexpr size_t Size       = ImSize * Planes;
	constexpr size_t Upsampling = 1;
	constexpr Vec3I Dim         = {.x = Width, .y = Height, .z = Planes};
	const std::vector<uint16_t> HighSNR(ImSize, 0), LowSNR(ImSize, 0), Input(Size, 0);
	std::vector<uint16_t> output(Size, 0);
	EXPECT_THROW(align.Align(HighSNR.data(), LowSNR.data(), Input.data(), output.data(), Dim, Upsampling), std::runtime_error);
}

//-------------------------------------------------------------------------------------------------
TEST(AlignmentTests, BlankCallWithFactor)
{
	Alignment align;

	constexpr size_t Height     = 100, Width = 100, Planes = 10;
	constexpr size_t ImSize     = Height * Width;
	constexpr size_t Size       = ImSize * Planes;
	constexpr size_t Upsampling = 1;
	constexpr Vec3I Dim         = {.x = Width, .y = Height, .z = Planes};
	const std::vector<uint16_t> Input(Size, 0);
	std::vector<uint16_t> output(Size, 0);

	constexpr std::array<std::array<double, 10>, 2> Factors = {};

	align.Align(Input.data(), output.data(), Dim, Factors, Upsampling);
	EXPECT_EQ(output, Input);
}

/// <summary>Vérifie l'identite et la réutilisation d'un objet sur une pile non carrée.</summary>
TEST(AlignmentTests, IdentityAndReuse)
{
	constexpr Vec3I Dim{.x = 7, .y = 5, .z = 3};
	std::array<std::array<double, 10>, 2> factors{};
	factors[0][7] = factors[1][8] = 1;
	std::vector<uint16_t> input(Dim.GetSize()), output(Dim.GetSize(), 0);
	std::iota(input.begin(), input.end(), uint16_t{1});
	Alignment align;
	align.Align(input.data(), output.data(), Dim, factors, 1);
	EXPECT_EQ(output, input);
	std::ranges::reverse(input);
	align.Align(input.data(), output.data(), Dim, factors, 1);
	EXPECT_EQ(output, input);
}

/// <summary>Le suréchantillonnage conserve aussi les derniers sous-pixels du bord droit et inférieur.</summary>
TEST(AlignmentTests, UpsampledIdentity)
{
	constexpr Vec3I Dim{.x = 7, .y = 5, .z = 2};
	std::array<std::array<double, 10>, 2> factors{};
	factors[0][7] = factors[1][8] = 1;
	std::vector<uint16_t> input(Dim.GetSize()), output(Dim.GetSize() * 4 + 2, 999);
	std::iota(input.begin(), input.end(), uint16_t{1});
	Alignment align;
	align.Align(input.data(), output.data() + 1, Dim, factors, 2);
	EXPECT_EQ(output.front(), 999);
	EXPECT_EQ(output.back(), 999);
	for (size_t z = 0; z < Dim.z; ++z) {
		for (size_t y = 0; y < Dim.y * 2; ++y) {
			for (size_t x = 0; x < Dim.x * 2; ++x) {
				EXPECT_EQ(output[1 + z * Dim.GetPlaneSize() * 4 + y * Dim.x * 2 + x], input[z * Dim.GetPlaneSize() + y / 2 * Dim.x + x / 2]);
			}
		}
	}
}

/// <summary>Une translation entièrement hors champ ne produit aucune contribution.</summary>
TEST_P(AlignmentOutsideTests, OutsideImage)
{
	constexpr Vec3I Dim{.x = 7, .y = 5, .z = 2};
	const std::vector<uint16_t> Input(Dim.GetSize(), 100);
	std::vector<uint16_t> output(Dim.GetSize() * 4 + 2, 999);
	Alignment align;
	constexpr std::array Shifts{
		Vec2D{.x = -20, .y = 0}, Vec2D{.x = 20, .y = 0},
		Vec2D{.x = 0, .y = -20}, Vec2D{.x = 0, .y = 20}
	};
	const auto& Shift = Shifts[GetParam()];
	std::array<std::array<double, 10>, 2> factors{};
	factors[0][7] = factors[1][8] = 1;
	factors[0][9] = Shift.x;
	factors[1][9] = Shift.y;
	align.Align(Input.data(), output.data() + 1, Dim, factors, 2);
	EXPECT_EQ(output.front(), 999);
	EXPECT_EQ(output.back(), 999);
	EXPECT_TRUE(std::all_of(output.begin() + 1, output.end() - 1, [](const auto value) { return value == 0; }));
}

INSTANTIATE_TEST_SUITE_P(Translations, AlignmentOutsideTests, testing::Values(size_t{0}, size_t{1}, size_t{2}, size_t{3}),
                         ([](const testing::TestParamInfo<size_t>& info) {
	                         constexpr std::array Names{"Left", "Right", "Above", "Below"};
	                         return Names[info.param];
	                         }));

/// <summary>Vérifie la translation entière et le remplissage des bords par des zéros.</summary>
TEST(AlignmentTests, IntegerTranslation)
{
	constexpr Vec3I Dim{.x = 6, .y = 5, .z = 2};
	std::array<std::array<double, 10>, 2> factors{};
	factors[0][7] = factors[1][8] = factors[0][9] = factors[1][9] = 1;
	std::vector<uint16_t> input(Dim.GetSize()), output(Dim.GetSize(), 999), expected(Dim.GetSize(), 0);
	std::iota(input.begin(), input.end(), uint16_t{1});
	for (size_t z = 0; z < Dim.z; ++z) {
		for (size_t y = 1; y < Dim.y; ++y) {
			for (size_t x = 1; x < Dim.x; ++x) { expected[z * Dim.GetPlaneSize() + y * Dim.x + x] = input[z * Dim.GetPlaneSize() + (y - 1) * Dim.x + x - 1]; }
		}
	}
	Alignment align;
	align.Align(input.data(), output.data(), Dim, factors, 1);
	EXPECT_EQ(output, expected);
}

/// <summary>Vérifie tous les pixels de la transposition sans mélanger les plans.</summary>
TEST(AlignmentTests, TransposePixels)
{
	constexpr Vec3I Dim{.x = 5, .y = 5, .z = 2};
	std::array<std::array<double, 10>, 2> factors{};
	factors[0][8] = factors[1][7] = 1;
	std::vector<uint16_t> input(Dim.GetSize()), output(Dim.GetSize(), 0), expected(Dim.GetSize(), 0);
	std::iota(input.begin(), input.end(), uint16_t{1});
	for (size_t z = 0; z < Dim.z; ++z) {
		for (size_t y = 0; y < Dim.y; ++y) {
			for (size_t x = 0; x < Dim.x; ++x) { expected[z * Dim.GetPlaneSize() + y * Dim.x + x] = input[z * Dim.GetPlaneSize() + x * Dim.x + y]; }
		}
	}
	Alignment align;
	align.Align(input.data(), output.data(), Dim, factors, 1);
	EXPECT_EQ(output, expected);
}

namespace {
/// <summary>Accès aux coordonnées polynomiales sans modifier l'API de production.</summary>
class AlignmentAccess final : public Alignment {
public:
	/// <summary>Évalue un seul coefficient du polynôme sur l'axe choisi.</summary>
	Vec2D Coordinate(const size_t coefficient, const bool vertical)
	{
		factors_                                = {};
		factors_[vertical ? 1 : 0][coefficient] = 2.5;
		Vec2D result{};
		computeCoord(2, 3, result.x, result.y);
		return result;
	}
};

/// <summary>Contrôle chaque monome et chacun des deux axes indépendamment.</summary>
class AlignmentPolynomialTests : public testing::TestWithParam<std::tuple<size_t, bool>> { };

/// <summary>Deformations non linéaires dont les recouvrements se calculent analytiquement.</summary>
class AlignmentWarpTests : public testing::TestWithParam<bool> { };
}

/// <summary>Vérifie l'ordre des dix coefficients cubiques, quadratiques et affines.</summary>
TEST_P(AlignmentPolynomialTests, CoordinateMonomials)
{
	const auto [coefficient, vertical] = GetParam();
	constexpr std::array<double, 10> Monomials{8, 27, 12, 18, 4, 9, 6, 2, 3, 1};
	AlignmentAccess align;
	const auto Result = align.Coordinate(coefficient, vertical);
	EXPECT_DOUBLE_EQ(vertical ? Result.y : Result.x, 2.5 * Monomials[coefficient]);
	EXPECT_DOUBLE_EQ(vertical ? Result.x : Result.y, 0);
}

INSTANTIATE_TEST_SUITE_P(Coefficients, AlignmentPolynomialTests,
                         testing::Combine(testing::Range(size_t{0}, size_t{10}), testing::Bool()),
                         ([](const testing::TestParamInfo<std::tuple<size_t, bool>>& info) {
	                         constexpr std::array Names{"XCubed", "YCubed", "XSquaredY", "XYSquared", "XSquared", "YSquared", "XY", "X", "Y", "Constant"};
	                         return std::string(std::get<1>(info.param) ? "Vertical" : "Horizontal") + Names[std::get<0>(info.param)];
	                         }));

/// <summary>Vérifie tous les poids d'une translation fractionnaire par une interpolation bilineaire indépendante.</summary>
TEST(AlignmentTests, FractionalTranslation)
{
	constexpr Vec3I Dim{.x = 6, .y = 5, .z = 2};
	std::array<std::array<double, 10>, 2> factors{};
	factors[0][7] = factors[1][8] = 1;
	factors[0][9] = 0.25;
	factors[1][9] = 0.5;
	std::vector<uint16_t> input(Dim.GetSize()), output(Dim.GetSize(), 999), expected(Dim.GetSize());
	for (size_t i = 0; i < input.size(); ++i) { input[i] = static_cast<uint16_t>(16 * (i + 1)); }
	for (size_t z = 0; z < Dim.z; ++z) {
		for (size_t y = 0; y < Dim.y; ++y) {
			for (size_t x = 0; x < Dim.x; ++x) {
				double value = 0;
				for (size_t dy = 0; dy < 2; ++dy) {
					for (size_t dx = 0; dx < 2; ++dx) {
						if (x < dx || y < dy) { continue; }
						value += (dx == 0 ? 0.75 : 0.25) * 0.5 * input[z * Dim.GetPlaneSize() + (y - dy) * Dim.x + x - dx];
					}
				}
				expected[z * Dim.GetPlaneSize() + y * Dim.x + x] = static_cast<uint16_t>(std::lround(value));
			}
		}
	}
	Alignment align;
	align.Align(input.data(), output.data(), Dim, factors, 1);
	EXPECT_EQ(output, expected);
}

/// <summary>Le pixel source devient un parallelogramme dont les aires de recouvrement sont connues.</summary>
TEST_P(AlignmentWarpTests, PolynomialPixelOverlap)
{
	constexpr Vec3I Dim{.x = 6, .y = 5, .z = 1};
	std::array<std::array<double, 10>, 2> factors{};
	factors[0][7] = factors[1][8] = 1;
	// X' = X + Y^2/4 ou X' = X + Y^3/8 ; Y' = Y.
	factors[0][GetParam() ? 1 : 5] = GetParam() ? 0.125 : 0.25;
	std::vector<uint16_t> input(Dim.GetSize(), 0), output(Dim.GetSize(), 0), expected(Dim.GetSize(), 0);
	input[Dim.x + 1] = 160;
	// Les coins sont relies lineairement par le splatting : intégration de deux trapezes.
	expected[Dim.x + 1] = GetParam() ? 70 : 60;
	expected[Dim.x + 2] = GetParam() ? 90 : 100;
	Alignment align;
	align.Align(input.data(), output.data(), Dim, factors, 1);
	EXPECT_EQ(output, expected);
	EXPECT_EQ(std::accumulate(output.begin(), output.end(), uint64_t{0}), 160);
}

INSTANTIATE_TEST_SUITE_P(Warps, AlignmentWarpTests, testing::Bool(),
                         [](const testing::TestParamInfo<bool>& info) { return info.param ? "Cubic" : "Quadratic"; });

#pragma region Rendus speciaux
/// <summary>Produit les rendus et vérifie la relecture exacte du TIFF effectivement enregistre.</summary>
TEST_P(AlignmentTests, SavesAndReloadsTransformedStack)
{
	if (!exists(OUTPUT_PATH)) { create_directories(OUTPUT_PATH); }
	const auto& [outputPath, factors, upsampling] = GetParam();

	std::vector<uint16_t> input;
	Vec3I dim{};

	// 1. Charger
	LoadTiff3D(REF_FILE_PATH.string().c_str(), input, dim);

	// 2. Allouer
	std::vector<uint16_t> output(input.size() * upsampling * upsampling, 0);

	// 3. Algorithme
	Alignment align;
	align.Align(input.data(), output.data(), dim, factors, upsampling);

	// 4. Sauvegarder
	const Vec3I Newdim{.x = dim.x * upsampling, .y = dim.y * upsampling, .z = dim.z};

	SaveTiff3D(outputPath.c_str(), output.data(), Newdim);

	if (GetParam().factors == std::array<std::array<double, 10>, 2>{}) {
		EXPECT_TRUE(std::ranges::all_of(output, [](const uint16_t value) { return value == 0; }));
	}
	else { EXPECT_TRUE(std::ranges::any_of(output, [](const uint16_t value) { return value > 0; })); }
	std::vector<uint16_t> reloaded;
	Vec3I reloadedDim{};
	LoadTiff3D(outputPath.c_str(), reloaded, reloadedDim);
	EXPECT_EQ(reloadedDim.x, Newdim.x);
	EXPECT_EQ(reloadedDim.y, Newdim.y);
	EXPECT_EQ(reloadedDim.z, Newdim.z);
	EXPECT_EQ(reloaded, output);
}

#ifndef INSTANTIATE_TEST_SUITE_P
#define INSTANTIATE_TEST_SUITE_P INSTANTIATE_TEST_CASE_P
#endif

INSTANTIATE_TEST_SUITE_P(
	TestSuites, AlignmentTests,
	::testing::Values(
		AlignTestParams{ (OUTPUT_PATH / "Align_1_Blank.tiff").string(), {}, 1 },
		AlignTestParams{ (OUTPUT_PATH / "Align_2_Copy.tiff").string(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[1][8] = 1.0;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_3_translate_0.25.tiff").string(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 0.25;
		factors[1][8] = 1.0;
		factors[1][9] = 0.25;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_3_translate_0.50.tiff").string(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 0.5;
		factors[1][8] = 1.0;
		factors[1][9] = 0.5;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_3_translate_0.75.tiff").string(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 0.75;
		factors[1][8] = 1.0;
		factors[1][9] = 0.75;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_3_translate_1.00.tiff").string(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 1.0;
		factors[1][8] = 1.0;
		factors[1][9] = 1.0;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_4_zoom.tiff").string(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 2.0;
		factors[1][8] = 2.0;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_5_transpose.tiff").string(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][8] = 1.0;
		factors[1][7] = 1.0;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_6_translate_0.50_up_sampling_2.tiff").string(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 0.5;
		factors[1][8] = 1.0;
		factors[1][9] = 0.5;
		return factors;
		}(), 2
		}
	),
	([](const testing::TestParamInfo<AlignTestParams>& info) {
		constexpr std::array Names{"Collapsed", "Identity", "TranslateQuarter", "TranslateHalf", "TranslateThreeQuarters",
		"TranslateOne", "Zoom", "Transpose", "TranslateHalfUpsampled"};
		return Names[info.index];
		}));
#pragma endregion Rendus speciaux
} // namespace palm_tracer::signal_processing
