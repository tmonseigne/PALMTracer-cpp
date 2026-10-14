#include <gtest/gtest.h>

#include <array>
#include <numbers>
#include "TestUtils.hpp"

#include "signal_processing/FFT.hpp"


namespace palm_tracer::signal_processing {
namespace {
/// <summary>Dimensions et variante FFT dont le nom apparait directement dans les rapports.</summary>
struct FFTCase {
	const char* name; ///< Nom explicite du scenario.
	Vec3I dim;        ///< Dimensions du volume.
	bool real3D;      ///< FFT directe 3D ou separable par axes.
};

/// <summary>Affiche le scenario sans representation binaire de la structure.</summary>
void PrintTo(const FFTCase& value, std::ostream* stream) { *stream << value.name; }

/// <summary>Validation independante des dimensions et des deux implementations FFT.</summary>
class FFTVolumeTests : public testing::TestWithParam<FFTCase> { };
}

/// <summary>Compare chaque frequence a une DFT directe, sans reutiliser l'algorithme inverse.</summary>
TEST_P(FFTVolumeTests, MatchesDirectDiscreteTransform)
{
	const auto& [name, dim, real3D] = GetParam();
	std::vector<double> input(dim.GetSize());
	for (size_t i = 0; i < input.size(); ++i) { input[i] = static_cast<double>(i % 7) - 3; }
	std::vector<std::complex<double>> actual;
	FFT3D(input, actual, dim, real3D);
	ASSERT_EQ(actual.size(), input.size());
	for (size_t kz = 0; kz < dim.z; ++kz) {
		for (size_t ky = 0; ky < dim.y; ++ky) {
			for (size_t kx = 0; kx < dim.x; ++kx) {
				std::complex<double> expected{};
				for (size_t z = 0; z < dim.z; ++z) {
					for (size_t y = 0; y < dim.y; ++y) {
						for (size_t x = 0; x < dim.x; ++x) {
							const double Phase = -2 * std::numbers::pi * (static_cast<double>(kx * x) / static_cast<double>(dim.x)
							                                              + static_cast<double>(ky * y) / static_cast<double>(dim.y) + static_cast<double>(
								                                              kz * z) / static_cast<double>(dim.z));
							expected += input[z * dim.GetPlaneSize() + y * dim.x + x] * std::polar(1.0, Phase);
						}
					}
				}
				const auto Value = actual[kz * dim.GetPlaneSize() + ky * dim.x + kx];
				EXPECT_NEAR(Value.real(), expected.real(), 1e-10);
				EXPECT_NEAR(Value.imag(), expected.imag(), 1e-10);
			}
		}
	}
}

/// <summary>Controle l'inversion et la normalisation de l'energie selon Parseval.</summary>
TEST_P(FFTVolumeTests, RoundTripAndParseval)
{
	const auto& [name, dim, real3D] = GetParam();
	std::vector<double> input(dim.GetSize());
	for (size_t i = 0; i < input.size(); ++i) { input[i] = static_cast<double>(i % 7) - 3; }
	std::vector<std::complex<double>> spectrum;
	std::vector<double> reconstructed;
	FFT3D(input, spectrum, dim, real3D);
	IFFT3D(spectrum, reconstructed, dim);
	ASSERT_EQ(reconstructed.size(), input.size());
	double spatialEnergy = 0, frequencyEnergy = 0;
	for (size_t i = 0; i < input.size(); ++i) {
		EXPECT_NEAR(reconstructed[i], input[i], 1e-10);
		spatialEnergy += input[i] * input[i];
		frequencyEnergy += std::norm(spectrum[i]);
	}
	EXPECT_NEAR(spatialEnergy, frequencyEnergy / static_cast<double>(input.size()), 1e-10);
}

INSTANTIATE_TEST_SUITE_P(Volumes, FFTVolumeTests,
                         testing::Values(FFTCase{"SingleVoxelSeparable", {.x = 1, .y = 1, .z = 1}, false},
	                         FFTCase{"SingleVoxel3D", {.x = 1, .y = 1, .z = 1}, true},
	                         FFTCase{"OddRectangleSeparable", {.x = 5, .y = 3, .z = 2}, false},
	                         FFTCase{"OddRectangle3D", {.x = 5, .y = 3, .z = 2}, true},
	                         FFTCase{"CubeSeparable", {.x = 4, .y = 4, .z = 4}, false},
	                         FFTCase{"Cube3D", {.x = 4, .y = 4, .z = 4}, true}),
                         [](const testing::TestParamInfo<FFTCase>& info) { return info.param.name; });

//-------------------------------------------------------------------------------------------------
TEST(FFTTest, CrossPowerSpectrumBasic)
{
	// Données FFT simulées.
	const std::vector<std::complex<double>> Ref    = {{1.0, 0.0}, {2.0, 1.0}, {3.0, -1.0}};
	const std::vector<std::complex<double>> Target = {{0.5, 0.5}, {2.0, -1.0}, {1.0, 0.0}};
	const std::vector Expected                     = {Ref[0] * std::conj(Target[0]), Ref[1] * std::conj(Target[1]), Ref[2] * std::conj(Target[2])};
	std::vector<std::complex<double>> out;

	CrossPowerSpectrum(Ref, Target, out);

	ASSERT_EQ(out.size(), Expected.size());
	for (size_t i = 0; i < out.size(); ++i) {
		EXPECT_NEAR(out[i].real(), Expected[i].real(), 1e-12);
		EXPECT_NEAR(out[i].imag(), Expected[i].imag(), 1e-12);
	}
}
} // namespace palm_tracer::signal_processing
