# include "Alignment.hpp"

#include <algorithm> // std::min, std::max
#include <cmath>	 // pour std::lround
#include <iostream>
#include <omp.h>
#include <stdexcept>


namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
void Alignment::Align(const uint16_t* highSNR, const uint16_t* lowSNR, const uint16_t* input, uint16_t* output, const Vec3I& dim, const size_t upsampling)

{
	highSNR_     = highSNR;
	lowSNR_      = lowSNR;
	input_       = input;
	output_      = output;
	dim_         = dim;
	sizeIm_      = dim_.GetPlaneSize();
	sizeStack_   = dim_.GetSize();
	upsampling_  = upsampling;
	upsampling2_ = upsampling * upsampling;
	computeFactor();
	align();
}

//-------------------------------------------------------------------------------------------------
void Alignment::Align(const uint16_t* input, uint16_t* output, const Vec3I& dim, const std::array<std::array<double, 10>, 2>& factors,
                      const size_t upsampling)

{
	highSNR_     = nullptr;
	lowSNR_      = nullptr;
	input_       = input;
	output_      = output;
	dim_         = dim;
	sizeIm_      = dim_.GetPlaneSize();
	sizeStack_   = dim_.GetSize();
	upsampling_  = upsampling;
	upsampling2_ = upsampling * upsampling;
	factors_     = factors;
	align();
}

//-------------------------------------------------------------------------------------------------
void Alignment::computeFactor()
{
	throw std::runtime_error("The method to compute factor isn't implemented use precompute constructor : "
		"Alignment(const uint16_t* input, uint16_t* output, size_t height, size_t width, size_t planes, const std::array<std::array<double, 10>, 2>&factors, size_t upscale, size_t upscaleMode, bool spread)");
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
	std::vector transformGrid(dim_.y + 1, std::vector<Vec2D>(dim_.x + 1));

	// 1. computeCoord pour tous les coins
	#pragma omp parallel for
	for (int y = 0; y <= static_cast<int>(dim_.y); ++y) {
		for (int x = 0; x <= static_cast<int>(dim_.x); ++x) { computeCoord(x, y, transformGrid[y][x].x, transformGrid[y][x].y); }
	}

	// 2. splatting inverse : chaque pixel source distribue à des pixels cibles
	// Application du suréchantillonnage.
	const size_t UpW      = dim_.x * upsampling_;
	const size_t UpH      = dim_.y * upsampling_;
	const size_t UpSizeIm = UpW * UpH;
	const double Step     = 1.0 / static_cast<double>(upsampling_);
	transformMap_.resize(UpSizeIm); // chaque pixel cible contiendra une liste de contributions
	// Supprime les contributions du precedent appel en conservant les capacites allouees.
	for (auto& contributions : transformMap_) { contributions.clear(); }

	#pragma omp parallel for
	for (int y = 0; y < static_cast<int>(dim_.y); ++y) {
		const size_t RowIdx = y * dim_.x;
		for (int x = 0; x < static_cast<int>(dim_.x); ++x) {
			const size_t SrcIdx = RowIdx + x;

			const std::vector P = {transformGrid[y][x], transformGrid[y][x + 1], transformGrid[y + 1][x + 1], transformGrid[y + 1][x]};
			const auto [x0, x1] = std::minmax({P[0].x, P[1].x, P[2].x, P[3].x});
			const auto [y0, y1] = std::minmax({P[0].y, P[1].y, P[2].y, P[3].y});

			// Récupération de la boite englobante avec le suréchantillonnage
			const size_t Xmin = static_cast<size_t>(std::max(0.0, std::floor(x0))) * upsampling_;
			const size_t Xmax = static_cast<size_t>(std::min(static_cast<double>(dim_.x - 1), std::ceil(x1))) * upsampling_;
			const size_t Ymin = static_cast<size_t>(std::max(0.0, std::floor(y0))) * upsampling_;
			const size_t Ymax = static_cast<size_t>(std::min(static_cast<double>(dim_.y - 1), std::ceil(y1))) * upsampling_;

			for (size_t j = Ymin; j <= Ymax; ++j) {
				const double Ypix      = static_cast<double>(j) * Step;
				const size_t RowdstIdx = j * UpW;
				for (size_t i = Xmin; i <= Xmax; ++i) {
					const double Xpix = static_cast<double>(i) * Step;
					std::vector<Vec2D> clipped;
					clipPolygon(P, clipped, Xpix, Ypix, Step);
					if (!clipped.empty()) {
						#pragma omp critical // Evite la concurrence lors d'écriture simultanée à transformMap_[dstIdx]
						transformMap_[RowdstIdx + i].emplace_back(Contribution{.index = SrcIdx, .weight = getArea(clipped)});
					}
				}
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------
void Alignment::applyTransform() const
{
	const int UpSizeIm = static_cast<int>(sizeIm_ * upsampling2_);
	std::vector tmp(UpSizeIm, 0.0);

	for (size_t p = 0; p < dim_.z; ++p) {
		std::ranges::fill(tmp, 0.0);
		const size_t PlaneOffset = p * sizeIm_; // La source reste inchangé

		#pragma omp parallel for
		for (int i = 0; i < UpSizeIm; ++i) {
			double acc = 0.0;
			for (const auto& [Index, Weight] : transformMap_[i]) { acc += Weight * input_[PlaneOffset + Index]; }
			tmp[i] = acc;
		}
		std::ranges::transform(tmp, output_ + p * UpSizeIm, [](const double v) { return static_cast<uint16_t>(std::lround(v)); });
	}
}

//-------------------------------------------------------------------------------------------------
void Alignment::clipPolygon(const std::vector<Vec2D>& poly, std::vector<Vec2D>& clipped, const double xmin, const double ymin, const double pixelSize)
{
	const double Xmax = xmin + pixelSize, Ymax = ymin + pixelSize;

	auto inside = [](const Vec2D& c, const int edge, const double bound) {
		switch (edge) {
			case 0: return c.x >= bound; // left
			case 1: return c.x <= bound; // right
			case 2: return c.y >= bound; // bottom
			case 3: return c.y <= bound; // top
			default: return false;
		}
	};

	auto intersect = [](const Vec2D& a, const Vec2D& b, const int edge, const double bound) {
		if (edge < 2) { return Vec2D{.x = bound, .y = a.y + (b.y - a.y) * (bound - a.x) / (b.x - a.x)}; } // vertical clip
		return Vec2D{.x = a.x + (b.x - a.x) * (bound - a.y) / (b.y - a.y), .y = bound};                   // horizontal clip
	};

	clipped.reserve(4);
	for (const auto& P : poly) { clipped.push_back(P); }
	for (int edge = 0; edge < 4; ++edge) {
		std::vector<Vec2D> input = clipped;
		clipped.clear();

		if (input.empty()) { break; }

		for (size_t i = 0; i < input.size(); ++i) {
			const auto& A = input[i];
			const auto& B = input[(i + 1) % input.size()];
			double bound;
			switch (edge) {
				case 0: bound = xmin;
					break; // left
				case 1: bound = Xmax;
					break; // right
				case 2: bound = ymin;
					break; // bottom
				case 3:
				default: bound = Ymax;
					break; // top
			}

			const bool Ain = inside(A, edge, bound);
			const bool Bin = inside(B, edge, bound);

			if (Ain && Bin) { clipped.push_back(B); }
			else if (Ain && !Bin) { clipped.push_back(intersect(A, B, edge, bound)); }
			else if (!Ain && Bin) {
				clipped.push_back(intersect(A, B, edge, bound));
				clipped.push_back(B);
			}
		}
	}
}

//-------------------------------------------------------------------------------------------------
double Alignment::getArea(const std::vector<Vec2D>& poly) const
{
	double area    = 0.0;
	const size_t N = poly.size();
	for (size_t i = 0; i < N; ++i) {
		const auto& [x0, y0] = poly[i];
		const auto& [x1, y1] = poly[(i + 1) % N];
		area += x0 * y1 - x1 * y0;
	}
	return 0.5 * static_cast<double>(upsampling2_) * std::abs(area);
}
} // namespace palm_tracer
