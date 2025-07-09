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
#include "Structs.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
/// <summary>Calcul la transformée de Fourier 3D directe d'un volume réel.</summary>
/// <param name="in">Volume d'entrée réel en 1D, taille z*y*x.</param>
/// <param name="out">Volume de sortie complexe (FFT), taille z*y*x.</param>
/// <param name="size">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
/// <param name="real3D"></param>
inline void FFT3D(const std::vector<double>& in, std::vector<std::complex<double>>& out, const Coord3D& size, const bool real3D = false)
{
	const size_t S = size.GetSize();
	out.resize(S);
	if (real3D) {
		auto* fft_in  = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * S));
		auto* fft_out = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * S));

		// Remplir fft_in : imaginaire = 0
		for (size_t i = 0; i < S; ++i) {
			fft_in[i][0] = in[i];
			fft_in[i][1] = 0.0;
		}

		fftw_plan plan = fftw_plan_dft_3d(static_cast<int>(size.z), static_cast<int>(size.y), static_cast<int>(size.x),
		                                  fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);

		for (size_t i = 0; i < S; ++i) { out[i] = std::complex(fft_out[i][0], fft_out[i][1]); }

		fftw_destroy_plan(plan);
		fftw_free(fft_in);
		fftw_free(fft_out);
	}
	else { // Version matlab équivalente à fft(fft2(volume), [], 3) on effectue la FFT plan par plan puis sur Z.
		const size_t Sp = size.GetPlaneSize();
		// 1. FFT2D sur chaque plan (Z)
		auto* temp2D     = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * Sp));
		auto* buffer2D   = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * Sp));
		fftw_plan plan2D = fftw_plan_dft_2d(static_cast<int>(size.y), static_cast<int>(size.x), temp2D, buffer2D,FFTW_FORWARD, FFTW_ESTIMATE);
		for (size_t z = 0; z < size.z; ++z) {
			const size_t Offset = z * Sp;
			// Remplir temp2D
			for (size_t i = 0; i < Sp; ++i) {
				temp2D[i][0] = in[Offset + i];
				temp2D[i][1] = 0.0;
			}

			fftw_execute(plan2D);

			// Sauvegarder résultat FFT2 dans out temporairement
			for (size_t i = 0; i < Sp; ++i) { out[Offset + i] = std::complex(buffer2D[i][0], buffer2D[i][1]); }
		}
		fftw_destroy_plan(plan2D);
		fftw_free(temp2D);
		fftw_free(buffer2D);

		// 2. FFT1D le long de Z pour chaque pixel (x,y)
		auto* tempZ   = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * size.z));
		auto* bufferZ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * size.z));

		fftw_plan plan1D = fftw_plan_dft_1d(static_cast<int>(size.z), tempZ, bufferZ,FFTW_FORWARD, FFTW_ESTIMATE);

		for (size_t y = 0; y < size.y; ++y) {
			for (size_t x = 0; x < size.x; ++x) {
				// Extraire colonne Z pour ce (x, y)
				for (size_t z = 0; z < size.z; ++z) {
					tempZ[z][0] = out[z * Sp + y * size.x + x].real();
					tempZ[z][1] = out[z * Sp + y * size.x + x].imag();
				}

				fftw_execute(plan1D);

				// Écrire le résultat
				for (size_t z = 0; z < size.z; ++z) { out[z * Sp + y * size.x + x] = std::complex(bufferZ[z][0], bufferZ[z][1]); }
			}
		}
		fftw_destroy_plan(plan1D);
		fftw_free(tempZ);
		fftw_free(bufferZ);
	}
}

//-------------------------------------------------------------------------------------------------
/// <summary>Calcul la transformée de Fourier 3D inverse d'un volume complexe.</summary>
/// <param name="in">Volume d'entrée complexe (FFT), taille z*y*x.</param>
/// <param name="out">Volume de sortie réel normalisé, taille z*y*x.</param>
/// <param name="size">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
inline void IFFT3D(const std::vector<std::complex<double>>& in, std::vector<double>& out, const Coord3D& size)
{
	const size_t S = size.GetSize();
	out.resize(S);
	auto* fft_in  = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * S));
	auto* fft_out = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * S));

	for (size_t i = 0; i < S; ++i) {
		fft_in[i][0] = in[i].real();
		fft_in[i][1] = in[i].imag();
	}

	fftw_plan plan = fftw_plan_dft_3d(static_cast<int>(size.z), static_cast<int>(size.y), static_cast<int>(size.x),
	                                  fft_in, fft_out, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan);

	// Normaliser iFFT (division par taille)
	const auto Norm = static_cast<double>(S);
	for (size_t i = 0; i < S; ++i) { out[i] = fft_out[i][0] / Norm; }

	fftw_destroy_plan(plan);
	fftw_free(fft_in);
	fftw_free(fft_out);
}

//-------------------------------------------------------------------------------------------------
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

//-------------------------------------------------------------------------------------------------
/// <summary>Recentre la fréquence nulle d’un volume 3D (fftshift), comme en MATLAB.</summary>
/// <param name="volume">Volume 3D en mémoire linéaire (row-major x, y, z).</param>
/// <param name="size">Dimensions du volume source (width=X, height=Y, depth=Z).</param>
inline void FFTShift3D(std::vector<double>& volume, const Coord3D& size)
{
	std::vector temp(volume.size(), 0.0);
	const Coord3D Mid = size.GetMid();

	Coord3D pos{};
	for (pos.z = 0; pos.z < size.z; ++pos.z) {
		Coord3D shifted{};
		shifted.z = (pos.z + Mid.z) % size.z;
		for (pos.y = 0; pos.y < size.y; ++pos.y) {
			shifted.y = (pos.y + Mid.y) % size.y;
			for (pos.x = 0; pos.x < size.x; ++pos.x) {
				shifted.x                = (pos.x + Mid.x) % size.x;
				temp[pos.GetIndex(size)] = volume[shifted.GetIndex(size)];
			}
		}
	}
	volume.swap(temp); // Remplacer volume par la version recentrée
}
} // namespace palm_tracer
