# include "Alignment.hpp"

#include <algorithm>      // std::min, std::max
#include <cmath> // pour std::lround
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
	sizeIm_        = width_ * height_;
	size_          = sizeIm_ * planes_;
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
	sizeIm_        = width_ * height_;
	size_          = sizeIm_ * planes_;
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
void Alignment::align()
{
	initOutput();
	precomputeTransformMap();
	applyTransform();
}

//-------------------------------------------------------------------------------------------------
void Alignment::precomputeTransformMap()
{
	transformMap_.resize(sizeIm_); // chaque pixel cible contiendra une liste de contributions
	std::vector transformGrid(height_ + 1, std::vector<Coord>(width_ + 1));

	// 1. computeCoord pour tous les coins
	#pragma omp parallel for
	for (int y = 0; y <= static_cast<int>(height_); ++y) {
		for (int x = 0; x <= static_cast<int>(width_); ++x) { computeCoord(x, y, transformGrid[y][x].x, transformGrid[y][x].y); }
	}

	// 2. splatting inverse : chaque pixel source distribue à des pixels cibles
	#pragma omp parallel for
	for (int y = 0; y < static_cast<int>(height_); ++y) {
		const size_t RowIdx = y * width_;
		for (int x = 0; x < static_cast<int>(width_); ++x) {
			const size_t SrcIdx = RowIdx + x;

			std::vector P       = {transformGrid[y][x], transformGrid[y][x + 1], transformGrid[y + 1][x + 1], transformGrid[y + 1][x]};
			const auto [x0, x1] = std::minmax({P[0].x, P[1].x, P[2].x, P[3].x});
			const auto [y0, y1] = std::minmax({P[0].y, P[1].y, P[2].y, P[3].y});
			const size_t Xmin   = static_cast<size_t>(std::max(0.0, std::floor(x0)));
			const size_t Xmax   = static_cast<size_t>(std::min(static_cast<double>(width_ - 1), std::ceil(x1)));
			const size_t Ymin   = static_cast<size_t>(std::max(0.0, std::floor(y0)));
			const size_t Ymax   = static_cast<size_t>(std::min(static_cast<double>(height_ - 1), std::ceil(y1)));

			for (size_t j = Ymin; j <= Ymax; ++j) {
				const size_t rowdstIdx = j * width_;
				for (size_t i = Xmin; i <= Xmax; ++i) {
					std::vector<Coord> clipped;
					clipPolygon(P, clipped, static_cast<double>(i), static_cast<double>(j));
					if (!clipped.empty()) {
						#pragma omp critical // Evite la concurrence lors d'écriture simultanée à transformMap_[dstIdx]
						transformMap_[rowdstIdx + i].emplace_back(Contribution{.index = SrcIdx, .weight = getArea(clipped)});
					}
				}
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------
void Alignment::applyTransform() const
{
	std::vector<double> tmp(sizeIm_);

	for (size_t p = 0; p < planes_; ++p) {
		std::ranges::fill(tmp, 0.0);
		const size_t PlaneOffset = p * sizeIm_;

		#pragma omp parallel for
		for (int i = 0; i < static_cast<int>(sizeIm_); ++i) {
			double acc = 0.0;
			for (const auto& [Index, Weight] : transformMap_[i]) { acc += Weight * input_[PlaneOffset + Index]; }
			tmp[i] = acc;
		}
		std::ranges::transform(tmp, output_ + PlaneOffset, [](const double v) { return static_cast<unsigned short>(std::lround(v)); });
	}
}

//-------------------------------------------------------------------------------------------------
void Alignment::clipPolygon(const std::vector<Coord>& poly, std::vector<Coord>& clipped, const double xmin, const double ymin)
{
	const double Xmax = xmin + 1, Ymax = ymin + 1;

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
		if (edge < 2) { return Coord{.x = bound, .y = y0 + (y1 - y0) * (bound - x0) / (x1 - x0)}; } // vertical clip
		return Coord{.x = x0 + (x1 - x0) * (bound - y0) / (y1 - y0), .y = bound};                   // horizontal clip
	};

	clipped.reserve(4);
	for (const auto& P : poly) { clipped.push_back(P); }
	for (int edge = 0; edge < 4; ++edge) {
		std::vector<Coord> input = clipped;
		clipped.clear();

		if (input.empty()) { break; }

		for (size_t i = 0; i < input.size(); ++i) {
			const auto& A  = input[i];
			const auto& B  = input[(i + 1) % input.size()];
			const bool Ain = inside(A.x, A.y, edge, (edge % 2 == 0) ? (edge == 0 ? xmin : ymin) : (edge == 1 ? Xmax : Ymax));
			const bool Bin = inside(B.x, B.y, edge, (edge % 2 == 0) ? (edge == 0 ? xmin : ymin) : (edge == 1 ? Xmax : Ymax));

			if (Ain && Bin) { clipped.push_back(B); }
			else if (Ain && !Bin) {
				clipped.push_back(intersect(A.x, A.y, B.x, B.y, edge,
				                            (edge % 2 == 0) ? (edge == 0 ? xmin : ymin) : (edge == 1 ? Xmax : Ymax)));
			}
			else if (!Ain && Bin) {
				clipped.push_back(intersect(A.x, A.y, B.x, B.y, edge,
				                            (edge % 2 == 0) ? (edge == 0 ? xmin : ymin) : (edge == 1 ? Xmax : Ymax)));
				clipped.push_back(B);
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------
double Alignment::getArea(const std::vector<Coord>& poly)
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
