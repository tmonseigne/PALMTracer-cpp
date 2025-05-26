#include "GaussianFit.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <execution>
#include <numbers>
#include <numeric>
#include <ranges>
#include <vector>

namespace palm_tracer {
GaussianFit::GaussianFit(const fit_mode mode)
	: mode_(mode), n_(0), nParams_(4), sigma_(1), theta_(0)
{
	switch (mode_) {
		case fit_mode::gauss_2d: nParams_ = 4; // 4 paramètres : x, y, amplitude, offset
			break;
		case fit_mode::gauss_2d_sigma: nParams_ = 5; // 5 : x, y, amplitude, sigma, offset
			break;
		case fit_mode::gauss_3d: nParams_ = 6; // 6 : x, y, amplitude, sigmaX, sigmaY, offset
			break;
		case fit_mode::gauss_3d_theta: nParams_ = 7; // 7 : x, y, amplitude, sigmaX, sigmaY, theta, offset
			break;
	}
	nParams2_ = nParams_ * nParams_;
	pMin_.resize(nParams_);
	p_.resize(nParams_);
	pMax_.resize(nParams_);
}

//-------------------------------------------------------------------------------------------------
void GaussianFit::Run(const unsigned short* image, size_t width, size_t height, double& x, double& y, double& intensity,
                      double& sigmaX, double& sigmaY, double& theta, double& offset, double& rss, size_t roi)
{
	if (roi % 2 == 0) { roi--; } // Not present in GPU original file why ?
	n_     = roi * roi;
	sigma_ = sigmaX;
	theta_ = theta;

	const size_t Half = roi / 2, Ux = static_cast<size_t>(x), Uy = static_cast<size_t>(y);
	if (Half <= Ux && Ux <= width - Half - 1 && Half <= Uy && Uy <= height - Half - 1) {
		// Création des variables
		double lambda = LAMBDA_0;
		std::vector<double> alpha(nParams2_), alphaR(nParams2_), alphaRi(nParams2_),
		                    beta(nParams_), beta_l_deltaP(nParams_), p_try(nParams_), deltaP(nParams_),
		                    dat(n_), hat(n_), weight(n_), deltaY(n_), deltaY_w(n_);
		std::vector<size_t> t(2 * n_);

		// Initialization of dat & t (A vérifier)
		double minI, maxI; // Not present in GPU original file why ?
		preComputeFit(minI, maxI, image, width, Ux, Uy, dat, weight, t, roi);
		initP(x, y, sigmaX, sigmaY, theta, minI, maxI);

		double x2     = lm_matx(t, p_, dat, weight, alpha, beta, hat);
		double x2_old = x2;

		// Trouver la valeur absolue maximale
		double m = std::ranges::max(std::views::transform(beta, [](const double b) { return std::fabs(b); }));

		bool stop = false;
		if (m < EPSILON_1) { stop = true; }

		size_t iteration = 0;
		while (!stop) { // Main Loop
			iteration++;
			std::ranges::fill(deltaP, 0.0);

			std::copy_n(alpha.begin(), nParams2_, alphaR.begin()); // Copy du alpha
			for (size_t i = 0; i < nParams_; ++i) { alphaR[i * nParams_ + i] += lambda * alpha[i * nParams_ + i]; }
			inverse(alphaR, alphaRi);

			for (size_t i = 0; i < nParams_; ++i) { for (size_t j = 0; j < nParams_; ++j) { deltaP[i] += alphaRi[i * nParams_ + j] * beta[j]; } }

			for (size_t i = 0; i < nParams_; ++i) {
				p_try[i] = p_[i] + deltaP[i];
				p_try[i] = std::clamp(p_try[i], pMin_[i], pMax_[i]);
			}

			for (size_t i = 0; i < n_; ++i) { deltaY[i] = dat[i] - psf(p_try, static_cast<double>(t[i * 2]), static_cast<double>(t[i * 2 + 1])); }

			// deltaY .* weight
			std::ranges::transform(deltaY, weight, deltaY_w.begin(), std::multiplies());
			const double X2_Try = std::transform_reduce(deltaY.begin(), deltaY.end(), deltaY_w.begin(), 0.0, std::plus(), std::multiplies());

			// lambda * deltaP + beta
			std::ranges::transform(deltaP, beta, beta_l_deltaP.begin(), [lambda](const double d, const double b) { return lambda * d + b; });

			m = std::inner_product(deltaP.begin(), deltaP.end(), beta_l_deltaP.begin(), 0.0);

			if ((x2 - X2_Try) / (2 * m) > EPSILON_4) { // it IS significantly better
				x2_old = x2;
				std::copy_n(p_try.begin(), nParams_, p_.begin()); // accept p_try

				x2     = lm_matx(t, p_, dat, weight, alpha, beta, hat);
				lambda = std::max(lambda / LAMBDA_DN, LAMBDA_EPSILON_MIN);
			}
			else {               // it IS NOT better
				x2     = x2_old; // do not accept a_try
				lambda = std::min(lambda * LAMBDA_UP, LAMBDA_EPSILON_MAX);
			}

			m = 0;
			for (size_t i = 0; i < nParams_; ++i) { m = std::max(m, std::fabs(deltaP[i] / p_[i])); }

			if (m < EPSILON_2 && iteration > 2) { stop = true; }
			if (x2 / static_cast<double>(n_) < EPSILON_3 && iteration > 2) { stop = true; }

			m = std::ranges::max(std::views::transform(beta, [](const double b) { return std::fabs(b); }));

			if (m < EPSILON_1 && iteration > 2) { stop = true; }
			if (iteration >= N_ITER) { stop = true; }
		}

		x         = p_[1];
		y         = p_[2];
		intensity = p_[0];
		switch (mode_) {
			case fit_mode::gauss_2d: {
				offset = p_[3];
				break;
			}
			case fit_mode::gauss_2d_sigma: {
				sigmaX = p_[3];
				offset = p_[4];
				break;
			}
			case fit_mode::gauss_3d: {
				sigmaX = p_[3] / std::numbers::sqrt2_v<double>;
				sigmaY = p_[4] / std::numbers::sqrt2_v<double>;
				offset = p_[5];
				break;
			}
			case fit_mode::gauss_3d_theta: {
				sigmaX = p_[3] / std::numbers::sqrt2_v<double>;
				sigmaY = p_[4] / std::numbers::sqrt2_v<double>;
				theta  = p_[5];
				offset = p_[6];
				break;
			}
		}
		rss = computeRSS(dat, hat);
	}
}


// ==================== Run Computes ====================
//-------------------------------------------------------------------------------------------------
void GaussianFit::preComputeFit(double& min, double& max, const unsigned short* image, const size_t width,
                                const size_t x, const size_t y, std::vector<double>& dat, std::vector<double>& weight, std::vector<size_t>& t,
                                const size_t roi)
{
	min = std::numeric_limits<double>::max();    // Plus grand double possible (FLT_MAX) // Not present in GPU original file why ?
	max = std::numeric_limits<double>::lowest(); // Plus petit double possible (-FLT_MAX) // Not present in GPU original file why ?

	const int Half = static_cast<int>(roi / 2);
	for (int i = -Half; i <= Half; ++i) {
		const auto Idx = Half + (i + Half) * roi;
		for (int j = -Half; j <= Half; ++j) {
			dat[j + Idx]    = static_cast<double>(image[x + j + (y + i) * width]);
			const auto TIdx = (j + Idx) * 2;
			t[TIdx]         = x + j;
			t[TIdx + 1]     = y + i;

			max = std::max(max, dat[j + Idx]); // Not present in GPU original file why ?
			min = std::min(min, dat[j + Idx]); // Not present in GPU original file why ?
		}
	}

	std::ranges::fill(weight, 1);
	//for (int i = 0; i < n_; ++i) {
	//	weight[i] = 1; //sqrt((n_-nParams_+1) / (dat[ii]*dat[ii]))
	//	//weight[i] *= weight[i];
	//}
}

//-------------------------------------------------------------------------------------------------
void GaussianFit::initP(const double x, const double y, const double sigmaX, const double sigmaY, const double theta, const double minI, const double maxI)
{
	switch (mode_) {
		case fit_mode::gauss_2d: {
			pMin_ = {0, x - 1, y - 1, 0};
			p_    = {maxI - minI, x, y, minI};
			pMax_ = {65535, x + 1, y + 1, 65535};
			break;
		}
		case fit_mode::gauss_2d_sigma: {
			pMin_ = {0, x - 1, y - 1, 0.1, 0.};
			p_    = {maxI - minI, x, y, sigmaX, minI};
			pMax_ = {65535, x + 1, y + 1, 20, 65535};
			break;
		}
		case fit_mode::gauss_3d: {
			pMin_ = {0, x - 1, y - 1, 0.1, 0.1, 0};
			p_    = {maxI - minI, x, y, std::numbers::sqrt2_v<double> * sigmaX, std::numbers::sqrt2_v<double> * sigmaY, minI};
			pMax_ = {65535, x + 1, y + 1, 20, 20, 65535};
			break;
		}
		case fit_mode::gauss_3d_theta: {
			pMin_ = {0, x - 1, y - 1, 0.1, 0.1, -150, 0};
			p_    = {maxI - minI, x, y, std::numbers::sqrt2_v<double> * sigmaX, std::numbers::sqrt2_v<double> * sigmaY, theta, minI};
			pMax_ = {65535, x + 1, y + 1, 20, 20, 150, 65535};
			break;
		}
	}
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::computeRSS(std::vector<double>& dat, std::vector<double>& hat) const
{
	// computing "dat" mean & "hat" mean
	const double Dat_Mean = std::accumulate(dat.begin(), dat.end(), 0.0) / static_cast<double>(n_),
	             Hat_Mean = std::accumulate(hat.begin(), hat.end(), 0.0) / static_cast<double>(n_);

	// Normalisation en place
	for (size_t i = 0; i < n_; ++i) {
		dat[i] -= Dat_Mean;
		hat[i] -= Hat_Mean;
	}

	// Calcul des t1, t2, t3
	const double T1 = std::inner_product(dat.begin(), dat.end(), hat.begin(), 0.0), // dat * hat
	             T2 = std::inner_product(dat.begin(), dat.end(), dat.begin(), 0.0), // dat * dat
	             T3 = std::inner_product(hat.begin(), hat.end(), hat.begin(), 0.0), // hat * hat
	             R  = T1 / std::sqrt(T2 * T3);                                      // Calcul du R

	return std::isnan(R) ? 0 : R * R;
}

// ==================== PSF ====================
//-------------------------------------------------------------------------------------------------
double GaussianFit::psf(const std::vector<double>& p, const double x, const double y) const
{
	switch (mode_) {
		case fit_mode::gauss_2d: return psf(x - p[1], y - p[2], p[0], sigma_, p[3]);
		case fit_mode::gauss_2d_sigma: return psf(x - p[1], y - p[2], p[0], p[3], p[4]);
		case fit_mode::gauss_3d: return psf(x - p[1], y - p[2], p[0], p[3], p[4], theta_, p[5]);
		case fit_mode::gauss_3d_theta: return psf(x - p[1], y - p[2], p[0], p[3], p[4], p[5], p[6]);
	}
	return 0;
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::psf(const double x, const double y, const double i0, const double sigma, const double offset)
{
	return i0 * exp(-(x * x + y * y) / (2 * (sigma * sigma))) + offset;
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::psf(const double x, const double y, const double i0, const double sigmaX, const double sigmaY, const double theta, const double offset)
{
	const auto CosTheta = cos(theta), SinTheta = sin(theta);
	const auto Factor1  = (x * CosTheta - y * SinTheta) / sigmaX,
	           Factor2  = (x * SinTheta + y * CosTheta) / sigmaY;

	return i0 * exp(-(Factor1 * Factor1 + Factor2 * Factor2)) + offset;
}

// ==================== Matrix ====================
//-------------------------------------------------------------------------------------------------
void GaussianFit::deriv(double* out, const std::vector<double>& p, const double x, const double y) const
{
	switch (mode_) {
		case fit_mode::gauss_2d: return deriv(out, x - p[1], y - p[2], p[0], sigma_, false);
		case fit_mode::gauss_2d_sigma: return deriv(out, x - p[1], y - p[2], p[0], p[3], true);
		case fit_mode::gauss_3d: return deriv(out, x - p[1], y - p[2], p[0], p[3], p[4], theta_, false);
		case fit_mode::gauss_3d_theta: return deriv(out, x - p[1], y - p[2], p[0], p[3], p[4], p[5], true);
	}
}

//-------------------------------------------------------------------------------------------------
void GaussianFit::deriv(double* out, const double x, const double y, const double i0, const double sigma, const bool has_sigma)
{
	const auto X_2 = x * x, Y_2 = y * y, Sigma_2 = sigma * sigma;
	const auto Exp = exp(-(X_2 + Y_2) / (2 * Sigma_2));
	const auto P   = i0 * Exp;

	out[0] = Exp;             // I0
	out[1] = P * x / Sigma_2; // x0
	out[2] = P * y / Sigma_2; // y0

	if (has_sigma) {
		out[3] = P / (sigma * Sigma_2) * (X_2 + Y_2); // sigma
		out[4] = 1;                                   // offset
	}
	else { out[3] = 1; } // offset (quand sigma est fixe)
}

//-------------------------------------------------------------------------------------------------
void GaussianFit::deriv(double* out, const double x, const double y, const double i0, const double sigmaX, const double sigmaY,
                        const double theta, const bool has_theta)
{
	const auto X_2     = x * x, Y_2          = y * y, X_Y = x * y, SigmaX_2 = sigmaX * sigmaX, SigmaY_2 = sigmaY * sigmaY;
	const auto CosT    = cos(theta), SinT    = sin(theta);
	const auto CosT_2  = CosT * CosT, SinT_2 = SinT * SinT, CosSinT = CosT * SinT,
	           XCosT   = x * CosT, XSinT     = x * SinT, YCosT      = y * CosT, YSinT = y * SinT;
	const auto Factor1 = (XCosT - YSinT) / sigmaX,
	           Factor2 = (XSinT + YCosT) / sigmaY,
	           Factor3 = 2 * X_Y * CosSinT,
	           Factor4 = 2 * X_Y / SigmaX_2 - 2 * X_Y / SigmaY_2;
	const auto Exp     = exp(-(Factor1 * Factor1 + Factor2 * Factor2));
	const auto P       = i0 * Exp;

	out[0] = Exp;                                                                                                      // I0
	out[1] = P * (-((-2 * CosT_2 * x + 2 * y * CosSinT) / SigmaX_2) - (-2 * SinT_2 * x - 2 * y * CosSinT) / SigmaY_2); // x0
	out[2] = P * (-((2 * SinT_2 * y - 2 * x * CosSinT) / SigmaX_2) - (-2 * CosT_2 * y - 2 * x * CosSinT) / SigmaY_2);  // y0
	out[3] = P * (2 * (XCosT * XCosT - Factor3 + YSinT * YSinT) / (sigmaX * SigmaX_2));                                // sigma X
	out[4] = P * (2 * (XSinT * XSinT + Factor3 + YCosT * YCosT) / (sigmaY * SigmaY_2));                                // sigma Y

	if (has_theta) {
		out[5] = P * (2 * CosSinT * ((X_2 - Y_2) / SigmaX_2 + (Y_2 - X_2) / SigmaY_2) - SinT_2 * Factor4 + CosT_2 * Factor4); // theta
		out[6] = 1;                                                                                                           // offset
	}
	else { out[5] = 1; } // offset (quand theta est fixe)
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::chi2(std::vector<double>& alpha, std::vector<double>& beta, const std::vector<double>& dydp, const std::vector<double>& weight,
                         const std::vector<double>& deltaY) const
{
	std::vector<double> tmp(n_ * nParams_);
	for (size_t i = 0; i < n_; ++i) { for (size_t j = 0; j < nParams_; ++j) { tmp[i * nParams_ + j] = weight[i] * dydp[i * nParams_ + j]; } }

	// Computing alpha (A v�rifier)
	for (size_t i = 0; i < nParams_; ++i) {
		for (size_t j = 0; j < nParams_; ++j) {
			alpha[nParams_ * i + j] = 0;
			for (size_t k = 0; k < n_; ++k) { alpha[nParams_ * i + j] += dydp[k * nParams_ + i] * tmp[k * nParams_ + j]; }
		}
	}

	// Computing beta
	tmp.resize(n_);
	std::ranges::transform(weight, deltaY, tmp.begin(), std::multiplies{});

	for (size_t i = 0; i < nParams_; ++i) {
		beta[i] = 0;
		for (size_t j = 0; j < n_; ++j) { beta[i] += dydp[j * nParams_ + i] * tmp[j]; }
	}

	// Computing Chi_sq
	#ifdef __APPLE__
// Pas de support de std::execution sous macOS (Clang/libc++)
		return std::transform_reduce(deltaY.begin(), deltaY.end(), tmp.begin(), 0.0);
	#else
	return std::transform_reduce(std::execution::seq, deltaY.begin(), deltaY.end(), tmp.begin(), 0.0);
	#endif
}

//-------------------------------------------------------------------------------------------------
double GaussianFit::lm_matx(const std::vector<size_t>& t, const std::vector<double>& p, const std::vector<double>& dat, const std::vector<double>& weight,
                            std::vector<double>& alpha, std::vector<double>& beta, std::vector<double>& hat) const
{
	std::vector<double> delta(n_);
	for (size_t i = 0; i < n_; ++i) {
		hat[i]   = psf(p, static_cast<double>(t[i * 2]), static_cast<double>(t[i * 2 + 1]));
		delta[i] = dat[i] - hat[i];
	}

	// Computing Jacobian matrix (dydp)
	std::vector<double> dydp(n_ * nParams_);
	for (size_t i = 0; i < n_; ++i) { deriv(&dydp[i * nParams_], p, static_cast<double>(t[i * 2]), static_cast<double>(t[i * 2 + 1])); }

	return chi2(alpha, beta, dydp, weight, delta);
}

//-------------------------------------------------------------------------------------------------
void GaussianFit::inverse(std::vector<double>& in, std::vector<double>& out) const
{
	std::ranges::fill(out, 0.0);
	for (size_t i = 0; i < nParams_; ++i) { out[i * nParams_ + i] = 1.0; } // Initialisation de la matrice identité

	// limination de Gauss-Jordan
	for (size_t i = 0; i < nParams_; ++i) {
		const double Pivot     = in[i * nParams_ + i];
		const double Pivot_Inv = 1.0 / Pivot;

		// Normalisation de la ligne i
		for (size_t j = 0; j < nParams_; ++j) {
			in[i * nParams_ + j] *= Pivot_Inv;
			out[i * nParams_ + j] *= Pivot_Inv;
		}

		// Mise � z�ro des autres lignes
		for (size_t j = 0; j < nParams_; ++j) {
			if (j != i) {
				const double Factor = in[j * nParams_ + i];

				for (size_t k = 0; k < nParams_; ++k) {
					in[j * nParams_ + k] -= Factor * in[i * nParams_ + k];
					out[j * nParams_ + k] -= Factor * out[i * nParams_ + k];
				}
			}
		}
	}
}
}
