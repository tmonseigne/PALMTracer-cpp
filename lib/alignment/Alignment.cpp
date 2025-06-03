# include "Alignment.hpp"

#include <algorithm>      // std::min, std::max
#include <cmath> // pour std::lround
#include <initializer_list> // std::initializer_list
#include <iostream>
#include <stdexcept>

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
void Alignment::Align(const unsigned short* highSNR, const unsigned short* lowSNR, const unsigned short* input, unsigned short* output,
                      const size_t height, const size_t width, const size_t planes,
                      const size_t upscale, const size_t upscaleMode, const bool upscaleSource, const bool spread)

{
	highSNR_       = highSNR;
	lowSNR_        = lowSNR;
	input_         = input;
	output_        = output;
	height_        = height;
	width_         = width;
	planes_        = planes;
	size_          = width_ * height_ * planes_;
	upscale_       = upscale;
	upscaleMode_   = upscaleMode;
	upsacleSource_ = upscaleSource;
	spread_        = spread;
	computeFactor();
	align();
}

//-------------------------------------------------------------------------------------------------
void Alignment::Align(const unsigned short* input, unsigned short* output, const size_t height, const size_t width, const size_t planes,
                      const std::array<std::array<double, 10>, 2>& factors, const size_t upscale, const size_t upscaleMode, const bool upscaleSource,
                      const bool spread)

{
	highSNR_       = nullptr;
	lowSNR_        = nullptr;
	input_         = input;
	output_        = output;
	height_        = height;
	width_         = width;
	planes_        = planes;
	size_          = width_ * height_ * planes_;
	upscale_       = upscale;
	upscaleMode_   = upscaleMode;
	upsacleSource_ = upscaleSource;
	spread_        = spread;
	factors_       = factors;
	align();
}

//-------------------------------------------------------------------------------------------------
void Alignment::initOutput() const { for (size_t i = 0; i < size_; ++i) { output_[i] = 0; } }

//-------------------------------------------------------------------------------------------------
void Alignment::computeFactor()
{
	throw std::runtime_error("The method to compute factor isn't implemented use precompute constructor : "
		"Alignment(const unsigned short* input, unsigned short* output, size_t height, size_t width, size_t planes, const std::array<std::array<double, 10>, 2>&factors, size_t upscale, size_t upscaleMode, bool spread)");
}

//-------------------------------------------------------------------------------------------------

void Alignment::computeCoord(const double oldX, const double oldY, double& newX, double& newY) const
{
	const double XY  = oldX * oldY,
	             X2  = oldX * oldX,
	             Y2  = oldY * oldY;
	const double X2Y = X2 * oldY,
	             XY2 = oldX * Y2,
	             X3  = X2 * oldX,
	             Y3  = Y2 * oldY;

	newX = factors_[0][0] * X3 + factors_[0][1] * Y3 + factors_[0][2] * X2Y + factors_[0][3] * XY2 +
	       factors_[0][4] * X2 + factors_[0][5] * Y2 + factors_[0][6] * XY +
	       factors_[0][7] * oldX + factors_[0][8] * oldY + factors_[0][9];

	newY = factors_[1][0] * X3 + factors_[1][1] * Y3 + factors_[1][2] * X2Y + factors_[1][3] * XY2 +
	       factors_[1][4] * X2 + factors_[1][5] * Y2 + factors_[1][6] * XY +
	       factors_[1][7] * oldX + factors_[1][8] * oldY + factors_[1][9];
}

//-------------------------------------------------------------------------------------------------
size_t Alignment::getCoord(const size_t x, const size_t y) const { return y * width_ + x; }

//-------------------------------------------------------------------------------------------------
void Alignment::align() const
{
	initOutput();
	const size_t SizeIm = height_ * width_;

	// Calcul d'une carte de transformation
	struct Coord {
		double x, y;
	};

	std::vector transformMap(height_ + 1, std::vector<std::array<double, 2>>(width_ + 1));

	#pragma omp parallel for //collapse(2) non utile de double parallelisation
	for (int y = 0; y <= static_cast<int>(height_); ++y) {
		for (int x = 0; x <= static_cast<int>(width_); ++x) { computeCoord(x, y, transformMap[y][x][0], transformMap[y][x][1]); }
	}


	for (size_t p = 0; p < planes_; ++p) {
		// Si upscale de l'image de base à faire ici les height et width seront à changer
		const size_t IdxPlane = p * SizeIm;
		std::vector<double> tmp(SizeIm, 0.0);
		for (size_t y = 0; y < height_; ++y) {
			const size_t IdxLine = y * width_;
			for (size_t x = 0; x < width_; ++x) {
				const size_t Idx = IdxLine + x;
				// Récupération des 4 points de ma nouvelle surface
				const std::array P = {transformMap[y][x], transformMap[y][x + 1], transformMap[y + 1][x + 1], transformMap[y + 1][x]};

				// Bounding box du quad
				const auto [x0, x1] = std::minmax({P[0][0], P[1][0], P[2][0], P[3][0]});
				const auto [y0, y1] = std::minmax({P[0][1], P[1][1], P[2][1], P[3][1]});
				const size_t Xmin   = static_cast<size_t>(std::max(0.0, std::floor(x0))),
				             Xmax   = static_cast<size_t>(std::min(static_cast<double>(width_ - 1), std::ceil(x1))),
				             Ymin   = static_cast<size_t>(std::max(0.0, std::floor(y0))),
				             Ymax   = static_cast<size_t>(std::min(static_cast<double>(height_ - 1), std::ceil(y1)));

				// Splatting (distribution de l’intensité sur une grille, pondérée par aire de recouvrement de mon quadrilatere p0,p1,p2,p3)
				for (size_t j = Ymin; j <= Ymax; ++j) {
					for (size_t i = Xmin; i <= Xmax; ++i) {
						if (i < 0 || i >= width_ || j < 0 || j >= height_) { continue; }
						// Récupération de l'intersection entre le polygone et le pixel
						auto clipped = clipPolygon(P, i, i + 1, j, j + 1);
						if (!clipped.empty()) {
							// Récupération de l'aire du pixel (i, j) recouverte par mon polygone p0,p1,p2,p3
							const double Area = getArea(clipped);
							tmp[j * width_ + i] += Area * input_[IdxPlane + Idx];
						}
					}
				}
			}
		}
		std::ranges::transform(tmp, output_ + IdxPlane, [](const double val) { return static_cast<unsigned short>(std::lround(val)); });
	}
}

//-------------------------------------------------------------------------------------------------
std::vector<std::array<double, 2>> Alignment::clipPolygon(const std::array<std::array<double, 2>, 4>& poly,
                                                          const double xmin, const double xmax, const double ymin, const double ymax)
{
	auto inside = [](const double x, const double y, const int edge, const double bound) {
		switch (edge) {
			case 0: return x >= bound; // left
			case 1: return x <= bound; // right
			case 2: return y >= bound; // bottom
			case 3: return y <= bound; // top
			default: return false;
		}
	};

	auto intersect = [](const double x0, const double y0, const double x1, const double y1, const int edge, const double bound) {
		double x = x0, y = y0;
		if (edge < 2) { // vertical clip
			x = bound;
			y = y0 + (y1 - y0) * (bound - x0) / (x1 - x0);
		}
		else { // horizontal clip
			y = bound;
			x = x0 + (x1 - x0) * (bound - y0) / (y1 - y0);
		}
		return std::array{x, y};
	};

	std::vector<std::array<double, 2>> output;
	output.reserve(4);
	for (const auto& P : poly) { output.push_back(P); }
	for (int edge = 0; edge < 4; ++edge) {
		std::vector<std::array<double, 2>> input = output;
		output.clear();

		if (input.empty()) { break; }

		for (size_t i = 0; i < input.size(); ++i) {
			const auto& A  = input[i];
			const auto& B  = input[(i + 1) % input.size()];
			const bool Ain = inside(A[0], A[1], edge, (edge % 2 == 0) ? (edge == 0 ? xmin : ymin) : (edge == 1 ? xmax : ymax));
			const bool Bin = inside(B[0], B[1], edge, (edge % 2 == 0) ? (edge == 0 ? xmin : ymin) : (edge == 1 ? xmax : ymax));

			if (Ain && Bin) { output.push_back(B); }
			else if (Ain && !Bin) {
				output.push_back(intersect(A[0], A[1], B[0], B[1], edge,
				                           (edge % 2 == 0) ? (edge == 0 ? xmin : ymin) : (edge == 1 ? xmax : ymax)));
			}
			else if (!Ain && Bin) {
				output.push_back(intersect(A[0], A[1], B[0], B[1], edge,
				                           (edge % 2 == 0) ? (edge == 0 ? xmin : ymin) : (edge == 1 ? xmax : ymax)));
				output.push_back(B);
			}
		}
	}

	return output;
}

//-------------------------------------------------------------------------------------------------
double Alignment::getArea(const std::vector<std::array<double, 2>>& poly)
{
	double area    = 0.0;
	const size_t N = poly.size();
	for (size_t i = 0; i < N; ++i) {
		const auto& [x0, y0] = poly[i];
		const auto& [x1, y1] = poly[(i + 1) % N];
		area += (x0 * y1 - x1 * y0);
	}
	return 0.5 * std::abs(area);
}
}
