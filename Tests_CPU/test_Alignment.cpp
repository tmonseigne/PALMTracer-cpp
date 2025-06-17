#include "gtest/gtest.h"

#include <algorithm>

#include "TiffManagement.hpp"
#include "alignment/Alignment.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
struct AlignTestParams {
	std::string outputPath;
	std::array<std::array<double, 10>, 2> factors;
	int upsampling;
};

//-------------------------------------------------------------------------------------------------
class AlignmentTests : public ::testing::TestWithParam<AlignTestParams> { };

//-------------------------------------------------------------------------------------------------
TEST(AlignmentTests, BaseCreation) { Alignment align; }

//-------------------------------------------------------------------------------------------------
TEST(AlignmentTests, BlankCallWithReference)
{
	Alignment align;

	constexpr size_t Height     = 100, Width = 100, Planes = 10;
	constexpr size_t ImSize     = Height * Width;
	constexpr size_t Size       = ImSize * Planes;
	constexpr size_t Upsampling = 1;
	const auto HighSNR          = new unsigned short[ImSize],
	           LowSNR           = new unsigned short[ImSize],
	           Input            = new unsigned short[Size],
	           Output           = new unsigned short[Size];

	// Création et Remplissage des vecteurs
	std::fill_n(HighSNR, ImSize, 0);
	std::fill_n(LowSNR, ImSize, 0);
	std::fill_n(Input, Size, 0);
	std::fill_n(Output, Size, 0);

	EXPECT_THROW(align.Align(HighSNR, LowSNR, Input, Output, Height, Width, Planes, Upsampling),
	             std::runtime_error);


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
	const auto Input            = new unsigned short[Size],
	           Output           = new unsigned short[Size];

	constexpr std::array<std::array<double, 10>, 2> Factors = {};

	// Création et Remplissage des vecteurs
	std::fill_n(Input, Size, 0);
	std::fill_n(Output, Size, 0);

	align.Align(Input, Output, Height, Width, Planes, Factors, Upsampling);

	delete[] Input;
	delete[] Output;
}

TEST_P(AlignmentTests, AlignmentCore)
{
	const auto& [outputPath, factors, upsampling] = GetParam();

	std::vector<unsigned short> input;
	size_t width, height, planes;

	// 1. Charger
	loadTiff3D("C:/Git/PALMTracer-cpp/Tests_CPU/Input/Alignement_Reference.tiff", input, width, height, planes);

	// 2. Allouer
	std::vector<unsigned short> output(input.size() * upsampling * upsampling, 0);

	// 3. Algorithme
	Alignment align;
	align.Align(input.data(), output.data(), height, width, planes, factors, upsampling);

	// 4. Sauvegarder
	saveTiff3D(outputPath.c_str(), output.data(), width * upsampling, height * upsampling, planes);

	SUCCEED();
}

INSTANTIATE_TEST_CASE_P(
	TestSuites, AlignmentTests,
	::testing::Values(
		AlignTestParams{ "C:/Git/PALMTracer-cpp/Tests_CPU/Output/Align_1_Blank.tiff", {}, 1 },
		AlignTestParams{ "C:/Git/PALMTracer-cpp/Tests_CPU/Output/Align_2_Copy.tiff",
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[1][8] = 1.0;
		return factors;
		}(), 1
		},
		AlignTestParams{ "C:/Git/PALMTracer-cpp/Tests_CPU/Output/Align_3_translate_0.25.tiff",
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 0.25;
		factors[1][8] = 1.0;
		factors[1][9] = 0.25;
		return factors;
		}(), 1
		},
		AlignTestParams{ "C:/Git/PALMTracer-cpp/Tests_CPU/Output/Align_3_translate_0.50.tiff",
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 0.5;
		factors[1][8] = 1.0;
		factors[1][9] = 0.5;
		return factors;
		}(), 1
		},
		AlignTestParams{ "C:/Git/PALMTracer-cpp/Tests_CPU/Output/Align_3_translate_0.75.tiff",
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 0.75;
		factors[1][8] = 1.0;
		factors[1][9] = 0.75;
		return factors;
		}(), 1
		},
		AlignTestParams{ "C:/Git/PALMTracer-cpp/Tests_CPU/Output/Align_3_translate_1.00.tiff",
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 1.0;
		factors[0][9] = 1.0;
		factors[1][8] = 1.0;
		factors[1][9] = 1.0;
		return factors;
		}(), 1
		},
		AlignTestParams{ "C:/Git/PALMTracer-cpp/Tests_CPU/Output/Align_4_zoom.tiff",
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][7] = 2.0;
		factors[1][8] = 2.0;
		return factors;
		}(), 1
		},
		AlignTestParams{ "C:/Git/PALMTracer-cpp/Tests_CPU/Output/Align_5_transpose.tiff",
		[] {
		std::array<std::array<double, 10>, 2> factors = {};
		factors[0][8] = 1.0;
		factors[1][7] = 1.0;
		return factors;
		}(), 1
		},
		AlignTestParams{ "C:/Git/PALMTracer-cpp/Tests_CPU/Output/Align_6_translate_0.50_up_sampling_2.tiff",
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
