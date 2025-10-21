#include <gtest/gtest-param-test.h>
#include <gtest/gtest.h>

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
	std::vector<uint16_t> highSNR(ImSize, 0), lowSNR(ImSize, 0), input(Size, 0), output(Size, 0);
	EXPECT_THROW(align.Align(highSNR.data(), lowSNR.data(), input.data(), output.data(), Dim, Upsampling), std::runtime_error);
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
	std::vector<uint16_t> input(Size, 0), output(Size, 0);

	constexpr std::array<std::array<double, 10>, 2> Factors = {};

	align.Align(input.data(), output.data(), Dim, Factors, Upsampling);
}

TEST_P(AlignmentTests, AlignmentCore)
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

	SUCCEED();
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
	));
} // namespace palm_tracer
