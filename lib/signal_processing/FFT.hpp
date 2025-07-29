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
/// - le produit de Fourier croisé entre deux volumes FFTés,
/// - et un recentrage (fftshift) des fréquences pour visualisation.
///
/// La bibliothèque FFTW3 est utilisée comme backend.
///-------------------------------------------------------------------------------------------------

#pragma once
#include <complex>
#include <fftw3.h>
#include <vector>
#include "Structs.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
/// <summary>Calcule la transformée de Fourier 3D directe d'un volume réel (FFT).</summary>
/// <param name="in">Volume d'entrée réel en 1D, taille z*y*x.</param>
/// <param name="out">Volume de sortie complexe (FFT) (même taille que `in`).</param>
/// <param name="dim">Dimensions du volume 3D (width=X, height=Y, depth=Z).</param>
/// <param name="real3D">Si vrai, effectue une FFT 3D complète. Sinon, FFT2D par plan puis FFT1D sur Z (comme MATLAB).</param>
/// <remarks>Le volume résultant est non normalisé (i.e. de même échelle que l'entrée).</remarks>
inline void FFT3D(const std::vector<double>& in, std::vector<std::complex<double>>& out, const Vec3I& dim, const bool real3D = false)
{
	const size_t S = dim.GetSize();
	out.resize(S);
	if (real3D) {
		auto* fft_in  = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * S));
		auto* fft_out = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * S));

		// Remplir fft_in : imaginaire = 0
		for (size_t i = 0; i < S; ++i) {
			fft_in[i][0] = in[i];
			fft_in[i][1] = 0.0;
		}

		fftw_plan plan = fftw_plan_dft_3d(static_cast<int>(dim.z), static_cast<int>(dim.y), static_cast<int>(dim.x),
		                                  fft_in, fft_out, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);

		for (size_t i = 0; i < S; ++i) { out[i] = std::complex(fft_out[i][0], fft_out[i][1]); }

		fftw_destroy_plan(plan);
		fftw_free(fft_in);
		fftw_free(fft_out);
	}
	else { // Version matlab équivalente à fft(fft2(volume), [], 3) on effectue la FFT plan par plan puis sur Z.
		const size_t Sp = dim.GetPlaneSize();
		// 1. FFT2D sur chaque plan (Z)
		auto* temp2D     = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * Sp));
		auto* buffer2D   = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * Sp));
		fftw_plan plan2D = fftw_plan_dft_2d(static_cast<int>(dim.y), static_cast<int>(dim.x), temp2D, buffer2D,FFTW_FORWARD, FFTW_ESTIMATE);
		for (size_t z = 0; z < dim.z; ++z) {
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

		// 2. FFT1D le long de Z pour chaque pixel (x, y)
		auto* tempZ   = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * dim.z));
		auto* bufferZ = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * dim.z));

		fftw_plan plan1D = fftw_plan_dft_1d(static_cast<int>(dim.z), tempZ, bufferZ,FFTW_FORWARD, FFTW_ESTIMATE);

		for (size_t y = 0; y < dim.y; ++y) {
			for (size_t x = 0; x < dim.x; ++x) {
				// Extraire colonne Z pour ce (x, y)
				for (size_t z = 0; z < dim.z; ++z) {
					tempZ[z][0] = out[z * Sp + y * dim.x + x].real();
					tempZ[z][1] = out[z * Sp + y * dim.x + x].imag();
				}

				fftw_execute(plan1D);

				// Écrire le résultat
				for (size_t z = 0; z < dim.z; ++z) { out[z * Sp + y * dim.x + x] = std::complex(bufferZ[z][0], bufferZ[z][1]); }
			}
		}
		fftw_destroy_plan(plan1D);
		fftw_free(tempZ);
		fftw_free(bufferZ);
	}
}

//-------------------------------------------------------------------------------------------------

/// <summary>Calcule la transformée de Fourier 3D inverse (iFFT) d'un volume complexe.</summary>
/// <param name="in">Volume complexe en entrée (FFT directe préalablement calculée), taille z*y*x.</param>
/// <param name="out">Volume réel en sortie (normalisé par la taille totale), taille z*y*x.</param>
/// <param name="dim">Dimensions du volume 3D (X, Y, Z).</param>
/// <remarks>Le volume de sortie est toujours réel, les parties imaginaires sont ignorées.</remarks>
inline void IFFT3D(const std::vector<std::complex<double>>& in, std::vector<double>& out, const Vec3I& dim)
{
	const size_t S = dim.GetSize();
	out.resize(S);
	auto* fft_in  = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * S));
	auto* fft_out = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * S));

	for (size_t i = 0; i < S; ++i) {
		fft_in[i][0] = in[i].real();
		fft_in[i][1] = in[i].imag();
	}

	fftw_plan plan = fftw_plan_dft_3d(static_cast<int>(dim.z), static_cast<int>(dim.y), static_cast<int>(dim.x),
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
/// <param name="ref">Transformée de Fourier de la référence (complexe).</param>
/// <param name="target">Transformée de Fourier de la cible (complexe).</param>
/// <param name="g">Résultat du produit croisé (complexe).</param>
inline void CrossPowerSpectrum(const std::vector<std::complex<double>>& ref, const std::vector<std::complex<double>>& target,
                               std::vector<std::complex<double>>& g)
{
	const size_t N = ref.size();
	g.resize(N);
	for (size_t i = 0; i < N; ++i) { g[i] = ref[i] * std::conj(target[i]); }
}

//-------------------------------------------------------------------------------------------------
/// <summary>Recentre le spectre fréquentiel d'un volume 3D FFTé (fftshift).</summary>
/// <param name="volume">Volume réel (typiquement une magnitude spectrale) en mémoire linéaire (row-major x, y, z).</param>
/// <param name="dim">Dimensions du volume 3D (width=X, height=Y, depth=Z).</param>
/// <remarks>Effectue un décalage cyclique de la fréquence nulle vers le centre, comme en MATLAB.</remarks>
inline void FFTShift3D(std::vector<double>& volume, const Vec3I& dim)
{
	std::vector temp(volume.size(), 0.0);
	const Vec3I Mid = dim.GetMid();

	Vec3I pos{};
	for (pos.z = 0; pos.z < dim.z; ++pos.z) {
		Vec3I shifted{};
		shifted.z = (pos.z + Mid.z) % dim.z;
		for (pos.y = 0; pos.y < dim.y; ++pos.y) {
			shifted.y = (pos.y + Mid.y) % dim.y;
			for (pos.x = 0; pos.x < dim.x; ++pos.x) {
				shifted.x               = (pos.x + Mid.x) % dim.x;
				temp[pos.GetIndex(dim)] = volume[shifted.GetIndex(dim)];
			}
		}
	}
	volume.swap(temp); // Remplacer volume par la version recentrée
}
} // namespace palm_tracer
