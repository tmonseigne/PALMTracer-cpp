#include "GaussianModel.hpp"

#include <cmath>
#include <stdexcept>

namespace palm_tracer::fit {
//-------------------------------------------------------------------------------------------------
GaussianModel::GaussianModel(const double offset, const double amplitude, const double sigma, const bool hasSigma)
	: offset_(offset), amplitude_(amplitude), sigmaX_(sigma), sigmaY_(sigma), sigmaX2_(sigma * sigma), sigmaY2_(sigma * sigma),
	  cosTheta_(1), sinTheta_(0), cosTheta2_(1), sinTheta2_(0), cosSinTheta_(0), isotropic_(true), hasSigma_(hasSigma), hasTheta_(false) { }

//-------------------------------------------------------------------------------------------------
GaussianModel::GaussianModel(const double offset, const double amplitude, const double sigmaX, const double sigmaY, const double theta, const bool hasTheta)
	: offset_(offset), amplitude_(amplitude), sigmaX_(sigmaX), sigmaY_(sigmaY), sigmaX2_(sigmaX * sigmaX), sigmaY2_(sigmaY * sigmaY),
	  cosTheta_(std::cos(theta)), sinTheta_(std::sin(theta)), cosTheta2_(cosTheta_ * cosTheta_), sinTheta2_(sinTheta_ * sinTheta_),
	  cosSinTheta_(cosTheta_ * sinTheta_), isotropic_(false), hasSigma_(false), hasTheta_(hasTheta) { }

//-------------------------------------------------------------------------------------------------
double GaussianModel::Evaluate(const double x, const double y, const std::span<double> gradient) const
{
	const size_t Size = isotropic_ ? (hasSigma_ ? 5 : 4) : (hasTheta_ ? 7 : 6);
	if (!gradient.empty() && gradient.size() < Size) { throw std::invalid_argument("GaussianModel : gradient trop petit."); }
	const double X2 = x * x, Y2 = y * y;
	if (isotropic_) {
		const double Exp = std::exp(-(X2 + Y2) / (2 * sigmaX2_));
		const double P   = amplitude_ * Exp;
		if (!gradient.empty()) {
			gradient[0] = 1;
			gradient[1] = Exp;
			gradient[2] = P * x / sigmaX2_;
			gradient[3] = P * y / sigmaX2_;
			if (hasSigma_) { gradient[4] = P / (sigmaX_ * sigmaX2_) * (X2 + Y2); }
		}
		return P + offset_;
	}

	const double XCos = x * cosTheta_, XSin         = x * sinTheta_, YCos = y * cosTheta_, YSin = y * sinTheta_;
	const double F1   = (XCos - YSin) / sigmaX_, F2 = (XSin + YCos) / sigmaY_;
	const double Exp  = std::exp(-((F1 * F1) + (F2 * F2)));
	const double P    = amplitude_ * Exp;
	if (!gradient.empty()) {
		const double XY = x * y, F3 = 2 * XY * cosSinTheta_, F4 = (2 * XY / sigmaX2_) - (2 * XY / sigmaY2_);
		gradient[0]     = 1;
		gradient[1]     = Exp;
		gradient[2]     = P * (-(((-2 * cosTheta2_ * x) + (2 * y * cosSinTheta_)) / sigmaX2_) - (((-2 * sinTheta2_ * x) - (2 * y * cosSinTheta_)) / sigmaY2_));
		gradient[3]     = P * ((((2 * sinTheta2_ * y) - (2 * x * cosSinTheta_)) / sigmaX2_) - (((-2 * cosTheta2_ * y) - (2 * x * cosSinTheta_)) / sigmaY2_));
		gradient[4]     = P * (2 * ((XCos * XCos) - F3 + (YSin * YSin)) / (sigmaX_ * sigmaX2_));
		gradient[5]     = P * (2 * ((XSin * XSin) + F3 + (YCos * YCos)) / (sigmaY_ * sigmaY2_));
		if (hasTheta_) { gradient[6] = P * ((2 * cosSinTheta_ * (((X2 - Y2) / sigmaX2_) + ((Y2 - X2) / sigmaY2_))) - (sinTheta2_ * F4) + (cosTheta2_ * F4)); }
	}
	return P + offset_;
}
} // namespace palm_tracer::fit
