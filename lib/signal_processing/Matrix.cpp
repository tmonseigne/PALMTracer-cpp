#include "Matrix.hpp"


#include <algorithm>
#include <cassert>
#include <cmath>
#include <execution>
#include <numeric>
#include <ranges>
#include <vector>

namespace palm_tracer::matrix {
//=============================================
//==================== I/O ====================
//=============================================

//-------------------------------------------------------------------------------------------------
void Print(const std::vector<double>& data, const Coord3D& size)
{
	size_t idx = 0;
	for (size_t z = 0; z < size.z; ++z) {
		std::cout << "Plane " << z + 1 << " :\n";
		for (size_t y = 0; y < size.y; ++y) {
			std::cout << "|\t";
			for (size_t x = 0; x < size.x; ++x, ++idx) { std::cout << data[idx] << "\t|\t"; }
			std::cout << "\n";
		}
		std::cout << "\n\n";
	}
}

//========================================================
//==================== Basic Computes ====================
//========================================================
//-------------------------------------------------------------------------------------------------
std::vector<bool> GetValidIndex(const std::vector<double>& data)
{
	std::vector<bool> mask;
	mask.reserve(data.size());

	for (const double V : data) { mask.push_back(std::isfinite(V) && !std::isnan(V)); }

	return mask;
}

//-------------------------------------------------------------------------------------------------
size_t CountValid(const std::vector<double>& data)
{
	auto valid_view = data | std::views::filter([](const double v) { return std::isfinite(v) && !std::isnan(v); });
	return std::ranges::distance(valid_view);
}

//-------------------------------------------------------------------------------------------------
double Max(const std::vector<double>& data)
{
	return *std::ranges::max_element(data, {}, [](const double x) {
		return !std::isfinite(x) || std::isnan(x) ? -std::numeric_limits<double>::infinity() : x;
	});
}

//-------------------------------------------------------------------------------------------------
double SumValid(const std::vector<double>& data)
{
	auto valid_view = data | std::ranges::views::filter([](const double v) { return std::isfinite(v) && !std::isnan(v); });
	return std::accumulate(valid_view.begin(), valid_view.end(), 0.0);
}

//-------------------------------------------------------------------------------------------------
double MeanValid(const std::vector<double>& data)
{
	auto valid       = data | std::views::filter([](const double v) { return std::isfinite(v) && !std::isnan(v); });
	const double Sum = std::accumulate(valid.begin(), valid.end(), 0.0);
	const auto Count = std::ranges::distance(valid);
	return Count > 0 ? Sum / static_cast<double>(Count) : std::numeric_limits<double>::quiet_NaN();
}

//-------------------------------------------------------------------------------------------------
double MeanWithMask(const std::vector<double>& data, const std::vector<bool>& mask)
{
	const size_t N = std::min(data.size(), mask.size());

	// Matérialiser le range filtré
	std::vector<double> valid;
	valid.reserve(N); // Optimisation mémoire

	for (size_t i = 0; i < N; ++i) { if (mask[i] && std::isfinite(data[i]) && !std::isnan(data[i])) { valid.push_back(data[i]); } }

	if (valid.empty()) { return std::numeric_limits<double>::quiet_NaN(); }

	const double Sum = std::accumulate(valid.begin(), valid.end(), 0.0);
	return Sum / static_cast<double>(valid.size());
}

//-------------------------------------------------------------------------------------------------
double Median(const std::vector<double>& data)
{
	std::vector<double> filtered;
	for (const double V : data) { if (std::isfinite(V) && !std::isnan(V)) { filtered.push_back(V); } }
	if (filtered.empty()) { return std::numeric_limits<double>::quiet_NaN(); }

	const size_t Mid = filtered.size() / 2;
	// nth_element assure que les n premiers éléments sont plus petits que les n suivants, mais tout n'est pas trié.
	std::ranges::nth_element(filtered, filtered.begin() + static_cast<int64_t>(Mid));
	if (filtered.size() % 2 == 1) { return filtered[Mid]; } // Si impair, l'index central est le bon, sinon il s'agit du ratio entre les 2 éléments centraux.
	const auto A = *std::max_element(filtered.begin(), filtered.begin() + static_cast<int64_t>(Mid));
	return (A + filtered[Mid]) * 0.5;
}

//-------------------------------------------------------------------------------------------------
double Quantile(const std::vector<double>& data, const double q)
{
	std::vector<double> filtered;
	for (const double V : data) { if (std::isfinite(V) && !std::isnan(V)) { filtered.push_back(V); } }
	if (filtered.empty()) { return std::numeric_limits<double>::quiet_NaN(); }

	std::ranges::sort(filtered); // Tri complet du tableau

	// Index réel (interpolation linéaire)
	const double Pos   = q * static_cast<double>(filtered.size() - 1);
	const auto Low     = static_cast<size_t>(std::floor(Pos)),
	           High    = static_cast<size_t>(std::ceil(Pos));
	const double Alpha = Pos - static_cast<double>(Low);

	// L'index est celui le plus proche du Quantile demandé si l'on veut le quantile précis, il faut interpoler entre les 2 index les plus proches.
	// const size_t Idx = std::clamp(static_cast<size_t>(pos), static_cast<size_t>(0), filtered.size() - 1);
	if (Low == High) { return filtered[Low]; }
	return (1.0 - Alpha) * filtered[Low] + Alpha * filtered[High];
}

//-------------------------------------------------------------------------------------------------
double RmsDifference(const std::vector<double>& a, const std::vector<double>& b)
{
	assert(a.size() == b.size());

	size_t count = 0;
	double sum   = 0.0;

	for (size_t i = 0; i < a.size(); ++i) {
		if (std::isnan(a[i]) || std::isnan(b[i])) { continue; }
		const auto D = a[i] - b[i];
		sum += D * D;
		++count;
	}

	return count > 0 ? std::sqrt(sum / static_cast<double>(count)) : std::numeric_limits<double>::quiet_NaN();
}

//-------------------------------------------------------------------------------------------------
void Normalisation(std::vector<double>& a, const bool shift)
{
	if (shift) {
		const double Min = *std::ranges::min_element(a);
		for (auto& v : a) { v -= Min; }
	}
	const double Sum  = std::accumulate(a.begin(), a.end(), 0.0);
	const double Norm = Sum > 0.0 ? 1.0 / Sum : 1.0;
	for (auto& v : a) { v *= Norm; }
}

//-------------------------------------------------------------------------------------------------
void Crop(std::vector<double>& data, Coord3D& size, const bool cropX, const bool cropY, const bool cropZ)
{
	const Coord3D Offset{
		.x = static_cast<size_t>(cropX && size.x > 2 ? 1 : 0),
		.y = static_cast<size_t>(cropY && size.y > 2 ? 1 : 0),
		.z = static_cast<size_t>(cropZ && size.z > 2 ? 1 : 0)
	};

	const Coord3D NewSize{.x = size.x - 2 * Offset.x, .y = size.y - 2 * Offset.y, .z = size.z - 2 * Offset.z};

	std::vector<double> crop;
	crop.reserve(NewSize.GetSize());
	Coord3D pos{};
	for (pos.z = Offset.z; pos.z < size.z - Offset.z; ++pos.z) {
		for (pos.y = Offset.y; pos.y < size.y - Offset.y; ++pos.y) {
			for (pos.x = Offset.x; pos.x < size.x - Offset.x; ++pos.x) { crop.push_back(data[pos.GetIndex(size)]); }
		}
	}

	data = std::move(crop);
	size = NewSize;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ExtractSubVolume(const std::vector<double>& data, const Coord3D& size, const std::vector<size_t>& xrange,
                                     const std::vector<size_t>& yrange, const std::vector<size_t>& zrange)
{
	std::vector<double> out;
	out.reserve(xrange.size() * yrange.size() * zrange.size());
	for (const size_t Z : zrange) { for (const size_t Y : yrange) { for (const size_t X : xrange) { out.push_back(data[(Z * size.y + Y) * size.x + X]); } } }
	return out;
}

//======================================================
//==================== Manipulation ====================
//======================================================

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
std::vector<size_t> ReflectionIndices(const size_t n, const size_t k0)
{
	std::vector<size_t> base;
	base.reserve(2 * n - 2);
	for (size_t i = 0; i < n; ++i) { base.push_back(i); }
	for (size_t i = n - 2; i > 0; --i) { base.push_back(i); }

	std::vector<size_t> out(k0);
	for (size_t i = 0; i < k0; ++i) { out[i] = base[i % base.size()]; }
	return out;
}

//-------------------------------------------------------------------------------------------------
std::vector<size_t> ReflectIndex1D(const size_t n, const size_t pad)
{
	// Construction du motif réfléchi
	std::vector<size_t> base;
	base.reserve(2 * n - 2);
	for (size_t i = 0; i < n; ++i) { base.push_back(i); }
	for (size_t i = n - 2; i > 0; --i) { base.push_back(i); }

	const size_t Period = 2 * n - 2;
	std::vector<size_t> idx;
	idx.reserve(n + 2 * pad);

	for (int i = -static_cast<int>(pad); i <= static_cast<int>(n + pad); ++i) {
		const size_t Wrapped = (i % Period + Period) % Period;
		idx.push_back(base[Wrapped]);
	}

	return idx;
}
} // namespace palm_tracer::matrix
