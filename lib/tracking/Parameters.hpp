///-------------------------------------------------------------------------------------------------
///
/// \file Parameters.hpp
/// \brief
/// \author Thibaut Monseigne (CNRS).
/// \version 1.0.
/// \date 16/05/2025.
/// \copyright <a href="https://choosealicense.com/licenses/agpl-3.0/">GNU Affero General Public License v3.0</a>.
///
///-------------------------------------------------------------------------------------------------

#pragma once
#include <random>

namespace palm_tracer {
class Parameters {
public:
	Parameters()  = default;
	~Parameters() = default;

	// Static
	static constexpr double DIST_MAX     = 20;
	static constexpr size_t N_STEPS      = 5;
	static constexpr size_t N_ITERATIONS = 300;
	static constexpr size_t MODEL        = 2; //the number of the used model

	// Params
	double maxSpeed_;
	double costBirthDeath_;

	// Computes
	double stdIntensity_, stdMean_, stdDistance_, covMeanIntensity_;
	double maxX_, maxY_, maxZ_;

	size_t RandomModulo(const size_t n) { return rng_() % n; }
	double Randomfloat() { return static_cast<double>(rng_()) / static_cast<double>(std::mt19937::max()); }

	std::mt19937 rng_{42}; // Générateur de nombres pseudo-aléatoires, seedé avec 42
};
}
