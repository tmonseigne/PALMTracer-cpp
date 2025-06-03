# include "Alignment.hpp"

#include <cmath> // pour std::lround
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
size_t Alignment::getCoord(const double x, const double y) const { return static_cast<size_t>(std::lround(y)) * width_ + static_cast<size_t>(std::lround(x)); }

//-------------------------------------------------------------------------------------------------
void Alignment::align() const
{
	initOutput();
	const size_t SizeIm = height_ * width_;

	for (size_t p = 0; p < planes_; ++p) {
		// Si upscale de l'image de base à faire ici les height et width seront à changer
		const size_t IdxPlane = p * SizeIm;
		for (size_t y = 0; y < height_; ++y) {
			const size_t IdxLine = y * width_ + IdxPlane;
			for (size_t x = 0; x < width_; ++x) {
				const double X = static_cast<double>(x * upscale_),
				             Y = static_cast<double>(y * upscale_);
				double newX, newY;
				computeCoord(static_cast<double>(x), static_cast<double>(y), newX, newY);
				// Calculer les éléments à modifier en fonction du spread
				output_[getCoord(newX, newY) + IdxPlane] = input_[IdxLine + x];
			}
		}
	}
}
}
