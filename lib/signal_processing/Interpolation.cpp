#include "Interpolation.hpp"

#include <cmath>


namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
double LinearInterpolation3D(const std::vector<double>& volume, const Coord3D& size, const Vec3& pos)
{
	if (pos.x < 0 || pos.x >= static_cast<int>(size.x) || pos.y < 0 || pos.y >= static_cast<int>(size.y) || pos.z < 0 || pos.z >= static_cast<int>(size.z)) { return 0.0; }

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

//-------------------------------------------------------------------------------------------------
std::vector<double> LinearZoom3D(const std::vector<double>& volume, const Coord3D& size, const size_t level, Coord3D& newSize)
{
	const size_t Factor = static_cast<size_t>(1) << level; // 2 ^ level
	const double Ratio = 1.0 / static_cast<double>(Factor);
	newSize             = {
		.x = (size.x - 1) * Factor + 1,
		.y = (size.y - 1) * Factor + 1,
		.z = (size.z - 1) * Factor + 1
	};
	std::vector out(newSize.x * newSize.y * newSize.z, 0.0);

	Vec3 pos = {.x = 0, .y = 0, .z = 0};

	size_t idx = 0;
	for (size_t z = 0; z < newSize.z; ++z) {
		pos.z = Ratio * static_cast<double>(z);
		for (size_t y = 0; y < newSize.y; ++y) {
			pos.y = Ratio * static_cast<double>(y);
			for (size_t x = 0; x < newSize.x; ++x, ++idx) {
				pos.x = Ratio * static_cast<double>(x);
				out[idx] = LinearInterpolation3D(volume, size, pos);
			}
		}
	}
	return out;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> LinearShift3D(const std::vector<double>& volume, const Coord3D& size, const Vec3& shift)
{
	std::vector out(size.x * size.y * size.z, 0.0);

	size_t idx = 0;
	for (size_t z = 0; z < size.z; ++z) {
		for (size_t y = 0; y < size.y; ++y) {
			for (size_t x = 0; x < size.x; ++x, ++idx) {
				const Vec3 Pos = {
					.x = static_cast<double>(x) - shift.x,
					.y = static_cast<double>(y) - shift.y,
					.z = static_cast<double>(z) - shift.z
				};
				out[idx] = LinearInterpolation3D(volume, size, Pos);
			}
		}
	}
	return out;
}

//-------------------------------------------------------------------------------------------------
double CubicInterpolation1D(const std::array<double, 4>& v, const double x)
{
	return v[1] + 0.5 * x * (v[2] - v[0] + x * (2.0 * v[0] - 5.0 * v[1] + 4.0 * v[2] - v[3] + x * (3.0 * (v[1] - v[2]) + v[3] - v[0])));
}

//-------------------------------------------------------------------------------------------------
double CubicInterpolation3D(const std::vector<double>& volume, const Coord3D& size, const Vec3& pos)
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
			temp[zz + 1][yy + 1] = CubicInterpolation1D(line, Dx);
		}
	}

	std::array<double, 4> col{};
	for (int zz = 0; zz < 4; ++zz) { col[zz] = CubicInterpolation1D(temp[zz], Dy); }

	return CubicInterpolation1D(col, Dz);
}


//-------------------------------------------------------------------------------------------------
std::vector<double> CubicShift3D(const std::vector<double>& volume, const Coord3D& size, const Vec3& shift)
{
	std::vector out(size.x * size.y * size.z, 0.0);

	size_t idx = 0;
	for (size_t z = 0; z < size.z; ++z) {
		for (size_t y = 0; y < size.y; ++y) {
			for (size_t x = 0; x < size.x; ++x, ++idx) {
				const Vec3 Pos = {
					.x = static_cast<double>(x) - shift.x,
					.y = static_cast<double>(y) - shift.y,
					.z = static_cast<double>(z) - shift.z
				};
				out[idx] = LinearInterpolation3D(volume, size, Pos);
			}
		}
	}
	return out;
}
} // namespace palm_tracer
