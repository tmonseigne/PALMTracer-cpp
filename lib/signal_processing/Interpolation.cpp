#include "Interpolation.hpp"


namespace palm_tracer {
//-------------------------------------------------------------------------------------------------

double CubicInterpolation1D(const double v0, const double v1, const double v2, const double v3, const double x)
{
	return v1 + 0.5 * x * (v2 - v0 + x * (2.0 * v0 - 5.0 * v1 + 4.0 * v2 - v3 + x * (3.0 * (v1 - v2) + v3 - v0)));
}

//-------------------------------------------------------------------------------------------------
double CubicInterpolation3D(const std::vector<double>& volume, const size_t Z, const size_t Y, const size_t X, const double z, const double y, const double x)
{
	auto get = [&](const int zz, const int yy, const int xx) -> double {
		if (zz < 0 || zz >= static_cast<int>(Z) || yy < 0 || yy >= static_cast<int>(Y) || xx < 0 || xx >= static_cast<int>(X)) { return 0.0; }
		return volume[(zz * Y + yy) * X + xx];
	};

	const int iz = static_cast<int>(std::floor(z));
	const int iy = static_cast<int>(std::floor(y));
	const int ix = static_cast<int>(std::floor(x));

	const double dz = z - iz;
	const double dy = y - iy;
	const double dx = x - ix;

	double temp[4][4];
	for (int zz = -1; zz <= 2; ++zz) {
		for (int yy = -1; yy <= 2; ++yy) {
			double line[4];
			for (int xx = -1; xx <= 2; ++xx) { line[xx + 1] = get(iz + zz, iy + yy, ix + xx); }
			temp[zz + 1][yy + 1] = CubicInterpolation1D(line[0], line[1], line[2], line[3], dx);
		}
	}

	double col[4];
	for (int zz = 0; zz < 4; ++zz) { col[zz] = CubicInterpolation1D(temp[zz][0], temp[zz][1], temp[zz][2], temp[zz][3], dy); }

	return CubicInterpolation1D(col[0], col[1], col[2], col[3], dz);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> CubicZoom3D(const std::vector<double>& volume, const size_t Z, const size_t Y, const size_t X, const int factor)
{
	const size_t Zout = Z * factor;
	const size_t Yout = Y * factor;
	const size_t Xout = X * factor;
	std::vector out(Zout * Yout * Xout, 0.0);

	for (size_t zz = 0; zz < Zout; ++zz) {
		for (size_t yy = 0; yy < Yout; ++yy) {
			for (size_t xx = 0; xx < Xout; ++xx) {
				const double zf                   = static_cast<double>(zz) / factor;
				const double yf                   = static_cast<double>(yy) / factor;
				const double xf                   = static_cast<double>(xx) / factor;
				out[(zz * Yout + yy) * Xout + xx] = CubicInterpolation3D(volume, Z, Y, X, zf, yf, xf);
			}
		}
	}
	return out;
}
} // namespace palm_tracer
