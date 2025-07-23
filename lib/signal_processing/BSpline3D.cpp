#include "BSpline3D.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <numeric>
#include <vector>

#include "Interpolation.hpp"
#include "signal_processing/Matrix.hpp"


namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
BSpline3D::BSpline3D(const std::vector<double>& volume, const Coord3D& dim, const size_t degree, const Vec3& lambda)
	: dim_(dim), degree_(degree), lambda_({lambda.x, lambda.y, lambda.z}), coeffs_(volume)
{
	getPolesFilter();                                                                  // Calcul des pôles du filtre B-spline pour une régularisation directe.
	for (size_t d = 0; d < N_DIMENSIONS; ++d) { applyFilter(d, filters_[d], c0_[d]); } // Filtrage par dimension

	// Padding final par réflexion
	const auto Pad = static_cast<size_t>(std::floor(static_cast<double>(degree_) / 2.0));
	dimPadded_     = {.x = dim_.x + 2 * Pad + 1, .y = dim_.y + 2 * Pad + 1, .z = dim_.z + 2 * Pad + 1};

	// Générer les indices réfléchis pour chaque dimension
	const auto IdxX = matrix::ReflectIndex1D(dim_.x, Pad),
	           IdxY = matrix::ReflectIndex1D(dim_.y, Pad),
	           IdxZ = matrix::ReflectIndex1D(dim_.z, Pad);

	// Nouveau buffer paddé
	coeffsPadded_.resize(dimPadded_.GetSize());

	for (size_t z = 0; z < dimPadded_.z; ++z) {
		for (size_t y = 0; y < dimPadded_.y; ++y) {
			for (size_t x = 0; x < dimPadded_.x; ++x) {
				const size_t SrcIndex   = (IdxZ[z] * dim_.y + IdxY[y]) * dim_.x + IdxX[x],
				             DstIndex   = (z * dimPadded_.y + y) * dimPadded_.x + x;
				coeffsPadded_[DstIndex] = coeffs_[SrcIndex];
			}
		}
	}
	interpolateHDR();
	getFinalCoeff();
}

//-------------------------------------------------------------------------------------------------
double BSpline3D::Eval(const double x, const size_t deg) // NOLINT(misc-no-recursion)
{
	switch (deg) {
		case 0: {
			if (x >= -0.5 && x < 0.5) { return 1.0; }
			break;
		}

		case 1: {
			if (x >= -1 && x < 0) { return 1.0 + x; }
			if (x >= 0 && x < 1) { return 1.0 - x; }
			break;
		}
		case 2: {
			constexpr auto Edge   = 1.125, // = 9 / 8 = 3^2 / 2^3
			               Center = 0.75,  // = 3 / 4 = 3   / 2^2
			               A      = 1.5,
			               B      = 0.5;

			const double X2 = x * x;
			if (x >= -A && x < -B) { return Edge + A * x + B * X2; }
			if (x >= -B && x < B) { return Center - x * x; }
			if (x >= B && x < A) { return Edge - A * x + B * X2; }
			break;
		}
		case 3: {
			constexpr auto Edge   = 4.0 / 3.0, // 4 / 3
			               Center = 2.0 / 3.0, // 2 / 3
			               A      = 2.0,
			               //B      = 1.0,
			               C = 1.0 / 6.0,
			               D = 0.5;

			const double X2 = x * x, X3 = X2 * x;
			if (x >= -2 && x < -1) { return Edge + A * x + X2 + C * X3; }
			if (x >= -1 && x < 0) { return Center - X2 - D * X3; }
			if (x >= 0 && x < 1) { return Center - X2 + D * X3; }
			if (x >= 1 && x < 2) { return Edge - A * x + X2 - C * X3; }
			break;
		}
		default: {
			const double HalfDegP1 = (static_cast<double>(deg) + 1.0) / 2.0;
			return ((x + HalfDegP1) * Eval(x + 0.5, deg - 1) + (HalfDegP1 - x) * Eval(x - 0.5, deg - 1)) / static_cast<double>(deg);
		}
	}

	return 0;
}

//-------------------------------------------------------------------------------------------------
void BSpline3D::getFinalCoeff()
{
	constexpr size_t N    = 4;
	constexpr size_t Size = 64;
	std::vector a(Size * Size, 0.0), aInv(Size * Size, 0.0);

	// Construction de la matrice A
	for (size_t i = 0; i < N; ++i) {
		const double Dx = static_cast<double>(i) / 3.0;
		for (size_t j = 0; j < N; ++j) {
			const double Dy = static_cast<double>(j) / 3.0;
			for (size_t k = 0; k < N; ++k) {
				const double Dz  = static_cast<double>(k) / 3.0;
				const size_t Row = i * 16 + j * 4 + k;

				for (size_t l = 0; l < N; ++l) {
					const double Px = std::pow(Dx, l);
					for (size_t m = 0; m < N; ++m) {
						const double Py = std::pow(Dy, m);
						for (size_t n = 0; n < N; ++n) {
							const double Pz     = std::pow(Dz, n);
							const size_t Col    = l * 16 + m * 4 + n;
							a[Row * Size + Col] = Px * Py * Pz;
						}
					}
				}
			}
		}
	}
	matrix::InverseGaussJordan(a, aInv, Size);

	// Upsample l’image (x3)
	const Interpolation3D Interp(interpolation_mode::trilinear, true); // Attention il faut l'interpolation Spline pas trilineaire
	Coord3D upDim{};
	const auto Up = Interp.Upsample(coeffsHDR_, dim_, 3, upDim);

	// Calcul des coefficients finaux
	dimFinal_ = dim_ - 1;
	coeffsFinal_.resize(dimFinal_.GetSize() * Size, 0.0);
	size_t base = 0;
	for (size_t z = 0; z < dimFinal_.z; ++z) {
		for (size_t y = 0; y < dimFinal_.y; ++y) {
			for (size_t x = 0; x < dimFinal_.x; ++x, base += 64) {
				// extrait un cube (4x4x4) échantillonné tous les 1/3 ⇾ soit 64 valeurs
				std::array<double, Size> temp{};

				size_t idx = 0;
				for (size_t k = 0; k < N; ++k) {
					const size_t Zi = z * 3 + k;
					for (size_t j = 0; j < N; ++j) {
						const size_t Yi = y * 3 + j;
						for (size_t i = 0; i < N; ++i, ++idx) {
							const size_t Xi = x * 3 + i;
							temp[idx]       = Up[(Zi * upDim.y + Yi) * upDim.x + Xi];
						}
					}
				}

				// Résolution x = A \ temp
				std::array<double, Size> localCoeff{};
				for (size_t i = 0; i < Size; ++i) { for (size_t j = 0; j < Size; ++j) { localCoeff[i] += aInv[i * Size + j] * temp[j]; } }

				// Stockage dans coeffsFinal_
				std::ranges::copy(localCoeff, coeffsFinal_.begin() + base);
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------
void BSpline3D::getPolesFilter()
{
	for (size_t i = 0; i < N_DIMENSIONS; ++i) {
		filters_[i].clear();
		if (lambda_[i] < 0.0) { throw std::invalid_argument("lambda_[i] must be non-negative"); }

		switch (degree_) {
			case 0: {
				if (lambda_[i] <= std::numeric_limits<double>::epsilon()) { c0_[i] = 1.0; }
				else { throw std::invalid_argument("Even degree 0 with nonzero lambda_[i] not supported"); }
				break;
			}
			case 1: {
				if (lambda_[i] <= std::numeric_limits<double>::epsilon()) { c0_[i] = 1.0; }
				else {
					const double R = 1.0 + 1.0 / (2.0 * lambda_[i]) - std::sqrt(1.0 + 4.0 * lambda_[i]) / (2.0 * lambda_[i]);
					filters_[i].emplace_back(R);
					c0_[i] = -1.0 / lambda_[i];
				}
				break;
			}
			case 2: {
				if (lambda_[i] <= std::numeric_limits<double>::epsilon()) {
					filters_[i].emplace_back(2.0 * std::numbers::sqrt2 - 3.0);
					c0_[i] = 8.0;
				}
				else { throw std::invalid_argument("Even degree 2 with nonzero lambda_[i] not supported"); }
				break;
			}
			case 3: {
				if (lambda_[i] <= std::numeric_limits<double>::epsilon()) {
					filters_[i].emplace_back(std::numbers::sqrt3 - 2.0);
					c0_[i] = 6.0;
				}
				else {
					std::vector<double> p{1, -4, 6, -4, 1};
					p[1] += 1.0 / (6.0 * lambda_[i]);
					p[2] += 2.0 / (3.0 * lambda_[i]);
					p[3] += 1.0 / (6.0 * lambda_[i]);

					auto roots = solvePolynomial(p);
					std::ranges::copy_if(roots, std::back_inserter(filters_[i]), [](auto z) { return std::abs(z) <= 1.0; });

					c0_[i] = 1.0 / lambda_[i];
				}
				break;
			}
			case 4: {
				if (lambda_[i] <= std::numeric_limits<double>::epsilon()) {
					filters_[i] = {{-0.36134122590022}, {-0.0137254292973391}};
					c0_[i]      = 384.0;
				}
				else { throw std::invalid_argument("Even degree 4 with nonzero lambda_[i] not supported"); }
				break;
			}
			case 5: {
				if (lambda_[i] <= std::numeric_limits<double>::epsilon()) {
					filters_[i] = {{-0.430575347099973}, {-0.0430962882032647}};
					c0_[i]      = 120.0;
				}
				else {
					std::vector<double> p{1, -6, 15, -20, 15, -6, 1};
					p[1] -= 1.0 / (120.0 * lambda_[i]);
					p[2] -= 13.0 / (60.0 * lambda_[i]);
					p[3] -= 11.0 / (20.0 * lambda_[i]);
					p[4] -= 13.0 / (60.0 * lambda_[i]);
					p[5] -= 1.0 / (120.0 * lambda_[i]);

					auto roots = solvePolynomial(p);
					std::ranges::copy_if(roots, std::back_inserter(filters_[i]), [](auto z) { return std::abs(z) <= 1.0; });

					c0_[i] = -1.0 / lambda_[i];
				}
				break;
			}
			case 6: {
				if (lambda_[i] <= std::numeric_limits<double>::epsilon()) {
					filters_[i] = {{-0.488294589303046}, {-0.0816792710762375}, {-0.00141415180832582}};
					c0_[i]      = 46080.0;
				}
				else { throw std::invalid_argument("Even degree 6 with nonzero lambda_[i] not supported"); }
				break;
			}
			case 7: {
				if (lambda_[i] <= std::numeric_limits<double>::epsilon()) {
					filters_[i] = {{-0.535280430796439}, {-0.122554615192327}, {-0.00914869480960828}};
					c0_[i]      = 5040.0;
				}
				else {
					std::vector<double> p{1, -8, 28, -56, 70, -56, 28, -8, 1};
					p[1] += 1.0 / (5040.0 * lambda_[i]);
					p[2] += 1.0 / (42.0 * lambda_[i]);
					p[3] += 397.0 / (1680.0 * lambda_[i]);
					p[4] += 151.0 / (315.0 * lambda_[i]);
					p[5] += 397.0 / (1680.0 * lambda_[i]);
					p[6] += 1.0 / (42.0 * lambda_[i]);
					p[7] += 1.0 / (5040.0 * lambda_[i]);

					auto roots = solvePolynomial(p);
					std::ranges::copy_if(roots, std::back_inserter(filters_[i]), [](auto z) { return std::abs(z) <= 1.0; });

					c0_[i] = 1.0 / lambda_[i];
				}
				break;
			}
			default: throw std::invalid_argument("Unsupported B-spline degree");
		}
	}
}

//-------------------------------------------------------------------------------------------------
void BSpline3D::applyFilter(const size_t axis, const std::vector<std::complex<double>>& poles, const double gain)
{
	const double K0Tol = std::log(std::numeric_limits<double>::epsilon()); // log(2.22e-16) ≈ -36
	const size_t Nx    = dim_.x, Ny = dim_.y, Nz = dim_.z;
	const auto N       = axis == 0 ? Nx : axis == 1 ? Ny : Nz;

	std::vector<std::complex<double>> tmpCoeffs, tmpLine(N), filteredLine;
	tmpCoeffs.reserve(coeffs_.size());
	for (auto& c : coeffs_) { tmpCoeffs.emplace_back(c); }

	auto getId = [&](const size_t x, const size_t y, const size_t z) { return (z * Ny + y) * Nx + x; };

	for (const auto& Zi : poles) {
		const auto K0                     = static_cast<size_t>(std::ceil(K0Tol / std::log(std::abs(Zi))));
		const std::vector<size_t> Reflect = matrix::ReflectionIndices(N, K0);
		const std::complex<double> C0     = -Zi / (1.0 - Zi * Zi);
		if (axis == 0) { // X
			for (size_t z = 0; z < Nz; ++z) {
				for (size_t y = 0; y < Ny; ++y) {
					for (size_t x = 0; x < Nx; ++x) { tmpLine[x] = tmpCoeffs[getId(x, y, z)]; }
					filteredLine = symExpFilt(tmpLine, C0, Zi, K0, Reflect);
					for (size_t x = 0; x < Nx; ++x) { tmpCoeffs[getId(x, y, z)] = filteredLine[x]; }
				}
			}
		}
		else if (axis == 1) { // Y
			for (size_t z = 0; z < Nz; ++z) {
				for (size_t x = 0; x < Nx; ++x) {
					for (size_t y = 0; y < Ny; ++y) { tmpLine[y] = tmpCoeffs[getId(x, y, z)]; }
					filteredLine = symExpFilt(tmpLine, C0, Zi, K0, Reflect);
					for (size_t y = 0; y < Ny; ++y) { tmpCoeffs[getId(x, y, z)] = filteredLine[y]; }
				}
			}
		}
		else if (axis == 2) { // Z
			for (size_t y = 0; y < Ny; ++y) {
				for (size_t x = 0; x < Nx; ++x) {
					for (size_t z = 0; z < Nz; ++z) { tmpLine[z] = tmpCoeffs[getId(x, y, z)]; }
					filteredLine = symExpFilt(tmpLine, C0, Zi, K0, Reflect);
					for (size_t z = 0; z < Nz; ++z) { tmpCoeffs[getId(x, y, z)] = filteredLine[z]; }
				}
			}
		}
	}
	// Récupération de la partie Réelle et Gain multiplicatif global
	for (size_t i = 0; i < coeffs_.size(); ++i) { coeffs_[i] = std::real(tmpCoeffs[i]) * gain; }
}

//-------------------------------------------------------------------------------------------------
void BSpline3D::interpolateHDR()
{
	coeffsHDR_.resize(dim_.GetSize(), 0.0);
	for (size_t z = 0; z < dim_.z; ++z) {
		const double Z = static_cast<double>(z) + 1;
		for (size_t y = 0; y < dim_.y; ++y) {
			const double Y = static_cast<double>(y) + 1;
			for (size_t x = 0; x < dim_.x; ++x) {
				const double X = static_cast<double>(x) + 1;

				coeffsHDR_[(z * dim_.y + y) * dim_.x + x] = interpolate_at(X, Y, Z);
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------
double BSpline3D::interpolate_at(const double x, const double y, const double z) const
{
	const auto Pad = static_cast<size_t>(std::floor(static_cast<double>(degree_) / 2.0));

	const int Kx0 = static_cast<int>(std::floor(x)) - static_cast<int>(Pad),
	          Ky0 = static_cast<int>(std::floor(y)) - static_cast<int>(Pad),
	          Kz0 = static_cast<int>(std::floor(z)) - static_cast<int>(Pad);

	double result = 0.0;

	for (size_t q = 0; q <= degree_; ++q) {
		const int Iz = Kz0 + static_cast<int>(q);
		if (Iz < 0 || static_cast<size_t>(Iz) >= dimPadded_.z) { continue; }
		const double Bz = Eval(z - static_cast<double>(Iz), degree_);

		for (size_t i = 0; i <= degree_; ++i) {
			const int Iy = Ky0 + static_cast<int>(i);
			if (Iy < 0 || static_cast<size_t>(Iy) >= dimPadded_.y) { continue; }
			const double Byz = Eval(y - static_cast<double>(Iy), degree_) * Bz;

			for (size_t j = 0; j <= degree_; ++j) {
				const int Ix = Kx0 + static_cast<int>(j);
				if (Ix < 0 || static_cast<size_t>(Ix) >= dimPadded_.x) { continue; }
				const double Bx = Eval(x - static_cast<double>(Ix), degree_);

				const size_t Idx = (static_cast<size_t>(Iz) * dimPadded_.y + static_cast<size_t>(Iy)) * dimPadded_.x + static_cast<size_t>(Ix);
				result += coeffsPadded_[Idx] * Bx * Byz;
			}
		}
	}

	return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::complex<double>> BSpline3D::solvePolynomial(const std::vector<double>& coeffs, const double tol, const int max_iter)
{
	const size_t Degree = coeffs.size() - 1;
	// ---------- Degré 3 (formule de Cardano) ----------
	if (Degree == 3) {
		const double A1 = coeffs[1] / coeffs[0],
		             A2 = coeffs[2] / coeffs[0],
		             A3 = coeffs[3] / coeffs[0];

		const double Q = (3.0 * A2 - A1 * A1) / 9.0,
		             R = (9.0 * A1 * A2 - 27.0 * A3 - 2.0 * A1 * A1 * A1) / 54.0;
		const double D = Q * Q * Q + R * R;

		std::vector<std::complex<double>> roots;

		if (D >= 0) {
			const double S = std::cbrt(R + std::sqrt(D));
			const double T = std::cbrt(R - std::sqrt(D));
			roots.emplace_back(-A1 / 3.0 + (S + T));
			roots.emplace_back(-A1 / 3.0 - (S + T) / 2.0 + std::complex<double>(0, std::numbers::sqrt3 * (S - T) / 2.0));
			roots.emplace_back(-A1 / 3.0 - (S + T) / 2.0 - std::complex<double>(0, std::numbers::sqrt3 * (S - T) / 2.0));
		}
		else {
			const double Theta = std::acos(R / std::sqrt(-Q * Q * Q));
			const double X     = 2.0 * std::sqrt(-Q);
			roots.emplace_back(X * std::cos(Theta / 3.0) - A1 / 3.0);
			roots.emplace_back(X * std::cos((Theta + 2.0 * std::numbers::pi) / 3.0) - A1 / 3.0);
			roots.emplace_back(X * std::cos((Theta + 4.0 * std::numbers::pi) / 3.0) - A1 / 3.0);
		}
		return roots;
	}

	// ---------- Degré 5 et 7 : Durand–Kerner ----------

	// Initialisation : racines uniformément réparties dans le plan complexe
	std::vector<std::complex<double>> roots(Degree);
	for (size_t k = 0; k < Degree; ++k) {
		constexpr double R0 = 0.4;
		const double Theta  = 2.0 * std::numbers::pi * static_cast<double>(k) / static_cast<double>(Degree);
		roots[k]            = R0 * std::complex(std::cos(Theta), std::sin(Theta));
	}

	auto eval = [&](const std::complex<double> x) {
		std::complex val = 0.0;
		for (size_t i = 0; i <= Degree; ++i) { val += coeffs[i] * std::pow(x, static_cast<int>(Degree - i)); }
		return val;
	};

	for (int iter = 0; iter < max_iter; ++iter) {
		bool converged                         = true;
		std::vector<std::complex<double>> next = roots;

		for (size_t i = 0; i < Degree; ++i) {
			std::complex denom = 1.0;
			for (size_t j = 0; j < Degree; ++j) { if (j != i) { denom *= roots[i] - roots[j]; } }
			const auto Delta = eval(roots[i]) / denom;
			next[i]          = roots[i] - Delta;
			if (std::abs(Delta) > tol) { converged = false; }
		}
		roots = next;
		if (converged) { break; }
	}

	return roots;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::complex<double>> BSpline3D::symExpFilt(const std::vector<std::complex<double>>& x, const std::complex<double>& c0,
                                                        const std::complex<double>& zi, const size_t k0, const std::vector<size_t>& kVec)
{
	const size_t N = x.size();
	std::vector<std::complex<double>> y(N, 0.0);

	for (size_t k = 0; k < k0; ++k) { y[0] += std::pow(zi, static_cast<int>(k)) * x[kVec[k]]; } // Initialisation via réflexion
	for (size_t k = 1; k < N; ++k) { y[k] = x[k] + zi * y[k - 1]; }                             // Filtrage causal
	y[N - 1] = (2.0 * y[N - 1] - x[N - 1]) * c0;                                                // Symétrisation (équivalent de l'initialisation anticausale)
	for (int k = static_cast<int>(N) - 2; k >= 0; --k) { y[k] = (y[k + 1] - y[k]) * zi; }       // Filtrage anticausal

	return y;
}
} // namespace palm_tracer
