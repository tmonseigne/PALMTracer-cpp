///-------------------------------------------------------------------------------------------------
///
/// \file FFT.hpp
/// \brief Utilitaires pour FFT 3D via FFTW3 : calcul direct et inverse, produit de Fourier croisé.
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 02/07/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
/// \details
/// Cette bibliothèque propose des fonctions inline pour calculer :
/// - la FFT 3D directe d'un volume réel (double),
/// - la FFT 3D inverse d'un volume complexe,
/// - le produit de Fourier croisé entre deux volumes FFTés.
///
/// La FFTW3 est utilisée comme backend.
///-------------------------------------------------------------------------------------------------

#pragma once
#include <complex>
#include <fftw3.h>
#include <vector>


/// <summary>Calcul la transformée de Fourier 3D directe d'un volume réel.</summary>
/// <param name="in">Volume d'entrée réel en 1D, taille z*y*x.</param>
/// <param name="out">Volume de sortie complexe (FFT), taille z*y*x.</param>
/// <param name="width">Largeur (dimension x).</param>
/// <param name="height">Hauteur (dimension y).</param>
/// <param name="planes">Nombre de plans (dimension z).</param>
inline void FFT3D(const std::vector<double>& in, std::vector<std::complex<double>>& out, const size_t width, const size_t height, const size_t planes)
{
	const size_t S = planes * height * width;
	auto* fft_in   = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * S));
	auto* fft_out  = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * S));

	// Remplir fft_in : imaginaire = 0
	for (size_t i = 0; i < S; ++i) {
		fft_in[i][0] = in[i];
		fft_in[i][1] = 0.0;
	}

	fftw_plan plan = fftw_plan_dft_3d(static_cast<int>(planes), static_cast<int>(height), static_cast<int>(width),
	                                  fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	out.resize(S);
	for (size_t i = 0; i < S; ++i) { out[i] = std::complex(fft_out[i][0], fft_out[i][1]); }

	fftw_destroy_plan(plan);
	fftw_free(fft_in);
	fftw_free(fft_out);
}

/// <summary>Calcul la transformée de Fourier 3D inverse d'un volume complexe.</summary>
/// <param name="in">Volume d'entrée complexe (FFT), taille z*y*x.</param>
/// <param name="out">Volume de sortie réel normalisé, taille z*y*x.</param>
/// <param name="width">Largeur (dimension x).</param>
/// <param name="height">Hauteur (dimension y).</param>
/// <param name="planes">Nombre de plans (dimension z).</param>
inline void IFFT3D(const std::vector<std::complex<double>>& in, std::vector<double>& out, const size_t width, const size_t height, const size_t planes)
{
	const size_t S = planes * height * width;
	auto* fft_in   = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * S));
	auto* fft_out  = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * S));

	for (size_t i = 0; i < S; ++i) {
		fft_in[i][0] = in[i].real();
		fft_in[i][1] = in[i].imag();
	}

	fftw_plan plan = fftw_plan_dft_3d(static_cast<int>(planes), static_cast<int>(height), static_cast<int>(width),
	                                  fft_in, fft_out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	out.resize(S);
	// Normaliser iFFT (division par taille)
	const auto Norm = static_cast<double>(S);
	for (size_t i = 0; i < S; ++i) { out[i] = fft_out[i][0] / Norm; }

	fftw_destroy_plan(plan);
	fftw_free(fft_in);
	fftw_free(fft_out);
}

/// <summary>Calcule le produit de Fourier croisé G = FFT(ref) * conj(FFT(target)).</summary>
/// <param name="ref">FFT du volume de référence (complexe).</param>
/// <param name="target">FFT du volume cible (complexe).</param>
/// <param name="g">Produit de Fourier croisé résultant (complexe).</param>
inline void CrossPowerSpectrum(const std::vector<std::complex<double>>& ref, const std::vector<std::complex<double>>& target,
                               std::vector<std::complex<double>>& g)
{
	const size_t N = ref.size();
	g.resize(N);
	for (size_t i = 0; i < N; ++i) { g[i] = ref[i] * std::conj(target[i]); }
}
