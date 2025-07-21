#include "Interpolation.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <numeric>

#include "BSpline3D.hpp"
#include "FFT.hpp"


namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
double Interpolation3D::cubic1D(const std::array<double, 4>& v, const double x)
{
	return v[1] + 0.5 * x * (v[2] - v[0] + x * (2.0 * v[0] - 5.0 * v[1] + 4.0 * v[2] - v[3] + x * (3.0 * (v[1] - v[2]) + v[3] - v[0])));
}

double Interpolation3D::cubic3D(const std::vector<double>& volume, const Coord3D& size, const Vec3& pos)
{
	auto get = [&](const int x, const int y, const int z) -> double {
		if (x < 0 || x >= static_cast<int>(size.x) || y < 0 || y >= static_cast<int>(size.y) || z < 0 || z >= static_cast<int>(size.z)) { return 0.0; }
		return volume[(z * size.y + y) * size.x + x];
	};

	const int Ix = static_cast<int>(std::floor(pos.x));
	const int Iy = static_cast<int>(std::floor(pos.y));
	const int Iz = static_cast<int>(std::floor(pos.z));

	const double Dx = pos.x - Ix;
	const double Dy = pos.y - Iy;
	const double Dz = pos.z - Iz;

	std::array<std::array<double, 4>, 4> temp{};
	std::array<double, 4> line{};
	for (int zz = -1; zz <= 2; ++zz) {
		for (int yy = -1; yy <= 2; ++yy) {
			for (int xx = -1; xx <= 2; ++xx) { line[xx + 1] = get(Ix + xx, Iy + yy, Iz + zz); }
			temp[zz + 1][yy + 1] = cubic1D(line, Dx);
		}
	}

	std::array<double, 4> col{};
	for (int zz = 0; zz < 4; ++zz) { col[zz] = cubic1D(temp[zz], Dy); }

	return cubic1D(col, Dz);
}

//-------------------------------------------------------------------------------------------------
double Interpolation3D::cubicWeight1D(double x, const double a)
{
	x = std::abs(x);
	if (x <= 1.0) { return (a + 2) * x * x * x - (a + 3) * x * x + 1; }
	if (x < 2.0) { return a * x * x * x - 5 * a * x * x + 8 * a * x - 4 * a; }
	return 0.0;
}

//-------------------------------------------------------------------------------------------------
double Interpolation3D::cubicWeight3D(const std::vector<double>& volume, const Coord3D& size, const Vec3& pos)
{
	auto get = [&](const int x, const int y, const int z) -> double {
		if (x < 0 || x >= static_cast<int>(size.x) || y < 0 || y >= static_cast<int>(size.y) || z < 0 || z >= static_cast<int>(size.z)) { return 0.0; }
		return volume[(z * size.y + y) * size.x + x];
	};

	const int Ix = static_cast<int>(std::floor(pos.x));
	const int Iy = static_cast<int>(std::floor(pos.y));
	const int Iz = static_cast<int>(std::floor(pos.z));

	double val = 0.0;

	for (int dz = -1; dz <= 2; ++dz) {
		const double Wz = cubicWeight1D(pos.z - (Iz + dz));
		for (int dy = -1; dy <= 2; ++dy) {
			const double Wy = cubicWeight1D(pos.y - (Iy + dy));
			for (int dx = -1; dx <= 2; ++dx) {
				const double Wx = cubicWeight1D(pos.x - (Ix + dx));
				val += get(Ix + dx, Iy + dy, Iz + dz) * Wx * Wy * Wz;
			}
		}
	}

	return val;
}

//-------------------------------------------------------------------------------------------------
double Interpolation3D::trilinear(const std::vector<double>& volume, const Coord3D& size, const Vec3& pos)
{
	auto get = [&](const int x, const int y, const int z) -> double {
		if (x < 0 || x >= static_cast<int>(size.x) || y < 0 || y >= static_cast<int>(size.y) || z < 0 || z >= static_cast<int>(size.z)) { return 0.0; }
		return volume[(z * size.y + y) * size.x + x];
	};

	const int X0 = static_cast<int>(std::floor(pos.x)),
	          Y0 = static_cast<int>(std::floor(pos.y)),
	          Z0 = static_cast<int>(std::floor(pos.z));
	const int X1 = X0 + 1, Y1 = Y0 + 1, Z1 = Z0 + 1;

	const double Dx = pos.x - X0, Dy = pos.y - Y0, Dz = pos.z - Z0;

	// Trilinear interpolation
	const double C00 = get(X0, Y0, Z0) * (1 - Dx) + get(X1, Y0, Z0) * Dx,
	             C01 = get(X0, Y0, Z1) * (1 - Dx) + get(X1, Y0, Z1) * Dx,
	             C10 = get(X0, Y1, Z0) * (1 - Dx) + get(X1, Y1, Z0) * Dx,
	             C11 = get(X0, Y1, Z1) * (1 - Dx) + get(X1, Y1, Z1) * Dx;

	const double C0 = C00 * (1 - Dy) + C10 * Dy,
	             C1 = C01 * (1 - Dy) + C11 * Dy;

	return C0 * (1 - Dz) + C1 * Dz;
}

Vec3 Interpolation3D::SubpixelMax(const std::vector<double>& volume, const Coord3D& size, const Vec3& pos, const double delta, const double w)
{
	const int N     = static_cast<int>(std::round(2.0 * w / delta)) + 1;
	const int HalfN = N / 2;

	double bestVal = -std::numeric_limits<double>::infinity();
	Vec3 res{};
	for (int dz = -HalfN; dz <= HalfN; ++dz) {
		Vec3 shift{};
		shift.z = dz * delta;
		for (int dy = -HalfN; dy <= HalfN; ++dy) {
			shift.y = dy * delta;
			for (int dx = -HalfN; dx <= HalfN; ++dx) {
				shift.x              = dx * delta;
				const Vec3 SamplePos = pos + shift;
				const double V       = trilinear(volume, size, SamplePos);
				if (!std::isnan(V) && V > bestVal) {
					bestVal = V;
					res     = shift;
				}
			}
		}
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> Interpolation3D::Zoom(const std::vector<double>& volume, const Coord3D& size, const size_t level, Coord3D& newDim) const
{
	if (level == 0) {
		newDim = size;
		return volume;
	}
	return Upsample(volume, size, static_cast<double>(static_cast<size_t>(1) << level), newDim);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> Interpolation3D::Upsample(const std::vector<double>& volume, const Coord3D& size, const double factor, Coord3D& newDim) const
{
	const double Step = 1.0 / factor;

	newDim = {
		.x = static_cast<size_t>(std::round(static_cast<double>(size.x - 1) * factor + 1)),
		.y = static_cast<size_t>(std::round(static_cast<double>(size.y - 1) * factor + 1)),
		.z = static_cast<size_t>(std::round(static_cast<double>(size.z - 1) * factor + 1))
	};

	std::vector out(newDim.GetSize(), 0.0);
	Vec3 pos = {.x = 0.0, .y = 0.0, .z = 0.0};

	// Pointeur vers la fonction d'interpolation finale
	using interp_ptr  = double(*)(const std::vector<double>&, const Coord3D&, const Vec3&);
	interp_ptr interp = nullptr;

	switch (mode_) {
		case interpolation_mode::trilinear: interp = &trilinear;
			break;
		case interpolation_mode::cubic: interp = &cubic3D;
			break;
		case interpolation_mode::cubic_weight: interp = &cubicWeight3D;
			break;
	}

	// Boucle principale
	size_t idx = 0;
	for (size_t z = 0; z < newDim.z; ++z) {
		pos.z = static_cast<double>(z) * Step;
		for (size_t y = 0; y < newDim.y; ++y) {
			pos.y = static_cast<double>(y) * Step;
			for (size_t x = 0; x < newDim.x; ++x, ++idx) {
				pos.x    = static_cast<double>(x) * Step;
				out[idx] = (cropBound_ && isOut(pos, size)) ? 0.0 : interp(volume, size, pos);
			}
		}
	}
	return out;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> Interpolation3D::Shift(const std::vector<double>& volume, const Coord3D& size, const Vec3& shift) const
{
	std::vector out(size.x * size.y * size.z, 0.0);
	Vec3 pos = {.x = 0, .y = 0, .z = 0};

	// Pointeur vers la fonction d'interpolation finale
	using interp_ptr  = double(*)(const std::vector<double>&, const Coord3D&, const Vec3&);
	interp_ptr interp = nullptr;

	switch (mode_) {
		case interpolation_mode::trilinear: interp = &trilinear;
			break;
		case interpolation_mode::cubic: interp = &cubic3D;
			break;
		case interpolation_mode::cubic_weight: interp = &cubicWeight3D;
			break;
	}

	size_t idx = 0;
	for (size_t z = 0; z < size.z; ++z) {
		pos.z = static_cast<double>(z) - shift.z;
		for (size_t y = 0; y < size.y; ++y) {
			pos.y = static_cast<double>(y) - shift.y;
			for (size_t x = 0; x < size.x; ++x, ++idx) {
				pos.x    = static_cast<double>(x) - shift.x;
				out[idx] = cropBound_ && isOut(pos, size) ? 0.0 : interp(volume, size, pos);
			}
		}
	}
	return out;
}

//-------------------------------------------------------------------------------------------------
Vec3 Interpolation3D::GetShift(const std::vector<double>& ref, const std::vector<double>& target, const Coord3D& size)
{
	if (std::accumulate(ref.begin(), ref.end(), 0.0) == 0.0 || std::accumulate(target.begin(), target.end(), 0.0) == 0.0) {
		return {.x = 0.0, .y = 0.0, .z = 0.0};
	}


	// Interpolation de niveau 2
	constexpr size_t ZoomLevel  = 2, // Niveau de Zoom de l'image HD
	                 ZoomFactor = 4; // Diviser le shift final par ce nombre pour revenir à l'échelle initiale (2^ZoomLevel)
	Coord3D size2X{};
	auto refHD = Zoom(ref, size, ZoomLevel, size2X);

	// Changement de mode pour l'élément shifté à cause des 0.
	// Reproduire le comportement matlab est impossible pour le moment. Cette interpolation est la plus proche.
	const auto PreviousMode = mode_;
	mode_                   = interpolation_mode::cubic_weight;
	auto tarHD              = Zoom(target, size, ZoomLevel, size2X);
	mode_                   = PreviousMode;

	// Normalisation (mean)
	const double MuRef = std::accumulate(refHD.begin(), refHD.end(), 0.0) / static_cast<double>(refHD.size());
	const double MuTar = std::accumulate(tarHD.begin(), tarHD.end(), 0.0) / static_cast<double>(tarHD.size());
	for (double& v : refHD) { v /= MuRef; }
	for (double& v : tarHD) { v /= MuTar; }

	// FFT3D & Produit croisé
	std::vector<std::complex<double>> fftRef, fftTar, fftG;
	FFT3D(refHD, fftRef, size2X, true);
	FFT3D(tarHD, fftTar, size2X, true);
	CrossPowerSpectrum(fftRef, fftTar, fftG);

	// iFFT & fftshift & post-traitement
	std::vector<double> spatialG;
	IFFT3D(fftG, spatialG, size2X);

	// fftshift (centre le pic de corrélation)
	FFTShift3D(spatialG, size2X);

	// Supprimer les valeurs < 0
	for (double& v : spatialG) { if (v < 0) { v = NAN; } }

	// Normalisation par le nombre de voxels valides
	const size_t Nvalid = std::ranges::count_if(tarHD, [](const double v) { return !std::isnan(v); });
	for (double& v : spatialG) { v = std::isnan(v) ? NAN : v / static_cast<double>(Nvalid) - 1.0; }


	// Maximum brut + interpolation locale
	const auto It = std::ranges::max_element(spatialG, [](const double a, const double b) { return std::isnan(a) || (!std::isnan(b) && a < b); });
	Vec3 shift{};
	if (It != spatialG.end()) {
		const size_t Idx = std::distance(spatialG.begin(), It);
		const Vec3 P     = CoordToVec(GetCoord3D(Idx, size2X));
		// Subpixel refinement par interpolation
		shift = SubpixelMax(spatialG, size2X, P);
	}

	// Calcul du décalage final à l'échelle originale
	return shift / ZoomFactor;
}
} // namespace palm_tracer
