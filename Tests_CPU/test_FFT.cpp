#include "gtest/gtest.h"

#include "FFT.hpp"


namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
TEST(FFTTest, ForwardInverseIdentity)
{
	constexpr size_t X = 4, Y = 4, Z = 4;
	std::vector<double> input(X * Y * Z);

	// Remplissage avec un pattern simple.
	for (size_t i = 0; i < input.size(); ++i) { input[i] = static_cast<double>(i % 10 + 1); }

	std::vector<std::complex<double>> fft_out;
	FFT3D(input, fft_out, X, Y, Z);

	std::vector<double> ifft_out;
	IFFT3D(fft_out, ifft_out, X, Y, Z);

	// Vérifier que le résultat approxime l'entrée originale.
	for (size_t i = 0; i < input.size(); ++i) { EXPECT_NEAR(input[i], ifft_out[i], 1e-9); }
}

//-------------------------------------------------------------------------------------------------
TEST(FFTTest, CrossPowerSpectrumBasic)
{
	// Données FFT simulées.
	const std::vector<std::complex<double>> Ref      = {{1.0, 0.0}, {2.0, 1.0}, {3.0, -1.0}};
	const std::vector<std::complex<double>> Target   = {{0.5, 0.5}, {2.0, -1.0}, {1.0, 0.0}};
	const std::vector<std::complex<double>> Expected = {Ref[0] * std::conj(Target[0]), Ref[1] * std::conj(Target[1]), Ref[2] * std::conj(Target[2])};
	std::vector<std::complex<double>> out;

	CrossPowerSpectrum(Ref, Target, out);

	ASSERT_EQ(out.size(), Expected.size());
	for (size_t i = 0; i < out.size(); ++i) {
		EXPECT_NEAR(out[i].real(), Expected[i].real(), 1e-12);
		EXPECT_NEAR(out[i].imag(), Expected[i].imag(), 1e-12);
	}
}
} // namespace palm_tracer
