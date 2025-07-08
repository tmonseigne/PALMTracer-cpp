#include "Matrix.hpp"

#include <algorithm>
#include <cmath>
#include <execution>
#include <vector>

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
void Deriv(double* out, const double i0, const double x, const double y, const double sigma, const bool has_sigma)
{
	const auto X_2 = x * x, Y_2 = y * y, Sigma_2 = sigma * sigma;
	const auto Exp = exp(-(X_2 + Y_2) / (2 * Sigma_2));
	const auto P   = i0 * Exp;

	out[0] = 1;                                                      // offset	: ∂PSF/∂offset
	out[1] = Exp;                                                    // I0		: ∂PSF/∂I0
	out[2] = P * x / Sigma_2;                                        // x0		: ∂PSF/∂x0
	out[3] = P * y / Sigma_2;                                        // y0		: ∂PSF/∂y0
	if (has_sigma) { out[4] = P / (sigma * Sigma_2) * (X_2 + Y_2); } // sigma	: ∂PSF/∂sigma
}

//-------------------------------------------------------------------------------------------------
void Deriv(double* out, const double i0, const double x, const double y, const double sigmaX, const double sigmaY, const double theta, const bool has_theta)
{
	// Calculs préliminaires
	const auto X_2 = x * x, Y_2 = y * y, X_Y = x * y, SigmaX_2 = sigmaX * sigmaX, SigmaY_2 = sigmaY * sigmaY;

	const auto CosT   = cos(theta), SinT    = sin(theta);
	const auto CosT_2 = CosT * CosT, SinT_2 = SinT * SinT, CosSinT = CosT * SinT,
	           XCosT  = x * CosT, XSinT     = x * SinT, YCosT      = y * CosT, YSinT = y * SinT;

	const auto Factor1 = (XCosT - YSinT) / sigmaX,
	           Factor2 = (XSinT + YCosT) / sigmaY,
	           Factor3 = 2 * X_Y * CosSinT,
	           Factor4 = 2 * X_Y / SigmaX_2 - 2 * X_Y / SigmaY_2;
	const auto Exp     = exp(-(Factor1 * Factor1 + Factor2 * Factor2));
	const auto P       = i0 * Exp;

	out[0] = 1;                                                                                                                              // offset
	out[1] = Exp;                                                                                                                            // I0
	out[2] = P * (-((-2 * CosT_2 * x + 2 * y * CosSinT) / SigmaX_2) - (-2 * SinT_2 * x - 2 * y * CosSinT) / SigmaY_2);                       // x0
	out[3] = P * (-((2 * SinT_2 * y - 2 * x * CosSinT) / SigmaX_2) - (-2 * CosT_2 * y - 2 * x * CosSinT) / SigmaY_2);                        // y0
	out[4] = P * (2 * (XCosT * XCosT - Factor3 + YSinT * YSinT) / (sigmaX * SigmaX_2));                                                      // sigma X
	out[5] = P * (2 * (XSinT * XSinT + Factor3 + YCosT * YCosT) / (sigmaY * SigmaY_2));                                                      // sigma Y
	if (has_theta) { out[6] = P * (2 * CosSinT * ((X_2 - Y_2) / SigmaX_2 + (Y_2 - X_2) / SigmaY_2) - SinT_2 * Factor4 + CosT_2 * Factor4); } // theta
}

//-------------------------------------------------------------------------------------------------
void InverseGaussJordan(std::vector<double>& in, std::vector<double>& out, const size_t n)
{
	std::ranges::fill(out, 0.0);
	for (size_t i = 0; i < n; ++i) { out[i * n + i] = 1.0; } // Initialisation de la matrice identité

	// limination de Gauss-Jordan
	for (size_t i = 0; i < n; ++i) {
		const double Pivot_Inv = 1.0 / in[i * n + i];

		// Normalisation de la ligne i
		for (size_t j = 0; j < n; ++j) {
			in[i * n + j] *= Pivot_Inv;
			out[i * n + j] *= Pivot_Inv;
		}

		// Élimination des autres lignes
		for (size_t j = 0; j < n; ++j) {
			if (j != i) {
				const double Factor = in[j * n + i];

				for (size_t k = 0; k < n; ++k) {
					in[j * n + k] -= Factor * in[i * n + k];
					out[j * n + k] -= Factor * out[i * n + k];
				}
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------
void InverseLU(std::vector<double>& in, std::vector<double>& out, const size_t n)
{
	size_t b;
	double tmp1 = 0;
	std::vector<double> y(n * n);

	// Décomposition LU
	for (size_t j = 0; j < n; ++j) {
		//calculate upper matrix
		for (size_t i = 0; i <= j; ++i) {
			//deal with ii-1 in the sum, set sum(kk=0->ii-1) when ii=0 to zero
			if (i > 0) {
				for (size_t k = 0; k <= i - 1; ++k) { tmp1 += in[i + k * n] * in[k + j * n]; }
				in[i + j * n] -= tmp1;
				tmp1 = 0;
			}
		}

		for (size_t i = j + 1; i < n; ++i) {
			if (j > 0) {
				for (size_t k = 0; k <= j - 1; ++k) { tmp1 += in[i + k * n] * in[k + j * n]; }
				in[i + j * n] = (1 / in[j + j * n]) * (in[i + j * n] - tmp1);
				tmp1          = 0;
			}
			else { in[i + j * n] = (1 / in[j + j * n]) * in[i + j * n]; }
		}
	}

	//tmp1 = 0; inutile réinitialisé après chaque utilisation
	// Résolution pour chaque colonne
	for (size_t num = 0; num < n; ++num) {
		// calculate yy
		if (num == 0) { y[0] = 1; }
		else { y[0] = 0; }

		for (size_t i = 1; i < n; ++i) {
			if (i == num) { b = 1; }
			else { b = 0; }
			for (size_t j = 0; j <= i - 1; ++j) { tmp1 += in[i + j * n] * y[j]; }
			y[i] = static_cast<double>(b) - tmp1;
			tmp1 = 0;
		}

		// calculate out
		out[n - 1 + num * n] = y[n - 1] / in[(n - 1) + (n - 1) * n];

		for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
			for (size_t j = i + 1; j < n; ++j) { tmp1 += in[i + j * n] * out[j + num * n]; }
			out[i + num * n] = (1 / in[i + i * n]) * (y[i] - tmp1);
			tmp1             = 0;
		}
	}
}

//-------------------------------------------------------------------------------------------------
void CenterOfMass(const std::vector<Pixel>& data, const size_t n, double& x, double& y)
{
	x          = 0.0;
	y          = 0.0;
	double sum = 0.0;

	size_t idx = 0;
	for (size_t i = 0; i < n; ++i) {            // lignes
		for (size_t j = 0; j < n; ++j, ++idx) { // colonnes
			const double Value = data[idx].v;
			y += Value * static_cast<double>(i);
			x += Value * static_cast<double>(j);
			sum += Value;
		}
	}

	if (sum > 0) {
		x /= sum;
		y /= sum;
	}
	else { // Valeur nulle partout ⇾ centre géométrique
		x = static_cast<double>(n) / 2.0;
		y = static_cast<double>(n) / 2.0;
	}
}

//-------------------------------------------------------------------------------------------------
bool Cholesky(const std::vector<double>& a, const size_t n, std::vector<double>& l, std::vector<double>& u)
{
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < i + 1; ++j) {
			double s = 0;
			// Calcul de la somme des produits scalaires entre lignes/colonnes précédentes
			for (size_t k = 0; k < j; ++k) { s += u[i * n + k] * u[j * n + k]; }

			if (i == j) {                    // Diagonale
				if (a[i * n + i] - s >= 0) { // On vérifie que la racine est réelle (matrice définie positive).
					u[i * n + j] = std::sqrt(a[i * n + i] - s);
					l[j * n + i] = u[i * n + j];
				}
				else { return false; } // La matrice n'est pas définie positive → échec
			}
			else { // Hors diagonale : calcul classique de la valeur
				u[i * n + j] = 1.0 / u[j * n + j] * (a[i * n + j] - s);
				l[j * n + i] = u[i * n + j];
			}
		}
	}
	return true;
}

//-------------------------------------------------------------------------------------------------
void Normalisation(std::vector<double>& a, const bool shift)
{
	if (shift) {
		const double Min = *std::ranges::min_element(a);
		for (auto& v : a) { v -= Min; }
	}
	const double Sum  = std::accumulate(a.begin(), a.end(), 0.0);
	const double Norm = (Sum > 0.0 ? 1.0 / Sum : 1.0);
	for (auto& v : a) { v *= Norm; }
}

void PrintVolume(const std::vector<double>& vol, const Coord3D& size)
{
	size_t idx = 0;
	for (size_t z = 0; z < size.z; ++z) {
		std::cout << "Plane " << z + 1 << " :\n";
		for (size_t y = 0; y < size.y; ++y) {
			std::cout << "|\t";
			for (size_t x = 0; x < size.x; ++x, ++idx) { std::cout << vol[idx] << "\t|\t"; }
			std::cout << "\n";
		}
		std::cout << "\n\n";
	}
}
} // namespace palm_tracer
