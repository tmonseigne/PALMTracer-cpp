#include "gtest/gtest.h"

#include <algorithm>
#include <filesystem>

#include "TiffManagement.hpp"
#include "signal_processing//Alignment.hpp"

namespace palm_tracer {
const std::filesystem::path PATH          = std::filesystem::path(__FILE__).parent_path();
const std::filesystem::path REF_FILE_PATH = PATH / "Input/Alignement_Reference.tiff";
const std::filesystem::path OUTPUT_PATH   = PATH / "Output";

//-------------------------------------------------------------------------------------------------
struct AlignTestParams {
	std::string outputPath;
	std::array<std::array<double, 10>, 2> factors;
	int upsampling;
};

//-------------------------------------------------------------------------------------------------
class AlignmentTests : public testing::TestWithParam<AlignTestParams> { };

//-------------------------------------------------------------------------------------------------
TEST(AlignmentTests, BlankCallWithReference)
{
	Alignment align;

	constexpr size_t Height     = 100, Width = 100, Planes = 10;
	constexpr size_t ImSize     = Height * Width;
	constexpr size_t Size       = ImSize * Planes;
	constexpr size_t Upsampling = 1;
	constexpr Vec3I Dim         = {.x = Width, .y = Height, .z = Planes};
	const auto HighSNR          = new uint16_t[ImSize],
	           LowSNR           = new uint16_t[ImSize],
	           Input            = new uint16_t[Size],
	           Output           = new uint16_t[Size];

	// Création et Remplissage des vecteurs
	std::fill_n(HighSNR, ImSize, 0);
	std::fill_n(LowSNR, ImSize, 0);
	std::fill_n(Input, Size, 0);
	std::fill_n(Output, Size, 0);

	EXPECT_THROW(align.Align(HighSNR, LowSNR, Input, Output, Dim, Upsampling), std::runtime_error);


	delete[] HighSNR;
	delete[] LowSNR;
	delete[] Input;
	delete[] Output;
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
	const auto Input            = new uint16_t[Size],
	           Output           = new uint16_t[Size];

	constexpr std::array<std::array<double, 10>, 2> Factors = {};

	// Création et Remplissage des vecteurs
	std::fill_n(Input, Size, 0);
	std::fill_n(Output, Size, 0);

	align.Align(Input, Output, Dim, Factors, Upsampling);

	delete[] Input;
	delete[] Output;
}

TEST_P(AlignmentTests, AlignmentCore)
{
	if (!std::filesystem::exists(OUTPUT_PATH)) { std::filesystem::create_directories(OUTPUT_PATH); }
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

	SUCCEED();
}

INSTANTIATE_TEST_CASE_P(
	TestSuites, AlignmentTests,
	::testing::Values(
		AlignTestParams{ (OUTPUT_PATH / "Align_1_Blank.tiff").string().c_str(), {}, 1 },
		AlignTestParams{ (OUTPUT_PATH / "Align_2_Copy.tiff").string().c_str(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[1][8] = 1.0;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_3_translate_0.25.tiff").string().c_str(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 0.25;
		factors[1][8] = 1.0;
		factors[1][9] = 0.25;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_3_translate_0.50.tiff").string().c_str(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 0.5;
		factors[1][8] = 1.0;
		factors[1][9] = 0.5;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_3_translate_0.75.tiff").string().c_str(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 0.75;
		factors[1][8] = 1.0;
		factors[1][9] = 0.75;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_3_translate_1.00.tiff").string().c_str(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 1.0;
		factors[1][8] = 1.0;
		factors[1][9] = 1.0;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_4_zoom.tiff").string().c_str(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 2.0;
		factors[1][8] = 2.0;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_5_transpose.tiff").string().c_str(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][8] = 1.0;
		factors[1][7] = 1.0;
		return factors;
		}(), 1
		},
		AlignTestParams{ (OUTPUT_PATH / "Align_6_translate_0.50_up_sampling_2.tiff").string().c_str(),
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 0.5;
		factors[1][8] = 1.0;
		factors[1][9] = 0.5;
		return factors;
		}(), 2
		}
	));
} // namespace palm_tracer
