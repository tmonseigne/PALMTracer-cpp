#include <algorithm>
#include <array>
#include <random>
#include <gtest/gtest.h>

#include "Structs.hpp"
#include "fit/AstigFit.hpp"

namespace palm_tracer::fit {
constexpr double SAMPLING       = 1.0;
constexpr double PIXEL_SIZE     = 160.0;
constexpr double Z_MAX          = 600.0;
constexpr size_t N_POINTS       = 1000;
constexpr astig_model REF_MODEL = {
	std::array{-50.0, 300.0, 0.06, 0.02, 140.0},
	std::array{70.0, 320.0, -0.04, 0.01, 150.0}
};

//-------------------------------------------------------------------------------------------------
static double SigmaModel(const std::array<double, 5>& p, const double z)
{
	// Accès indirect via Estimate n'est pas pratique ici : on recopie l'idée.
	const double Z0   = p[0], W = p[1], C3 = p[2], C4 = p[3], A = p[4];
	const double U    = (z - Z0) / W;
	const double U2   = U * U;
	const double Poly = 1.0 + U2 + C3 * (U2 * U) + C4 * (U2 * U2);
	const double G    = std::sqrt(std::max(0.0, Poly));
	return (SAMPLING / PIXEL_SIZE) * A * G;
}

//-------------------------------------------------------------------------------------------------
static std::vector<Vec3D> GetDataSet(const astig_model& model = REF_MODEL, const double zMax = Z_MAX, const size_t n = N_POINTS)
{
	std::mt19937 rng(42);
	std::normal_distribution noise(0.0, 0.02); // bruit en px

	std::vector<Vec3D> res(n);

	for (size_t i = 0; i < n; ++i) {
		const double T  = (n > 1) ? static_cast<double>(i) / static_cast<double>(n - 1) : 0.0;
		const double Z  = -zMax + T * (2.0 * zMax);
		const double Sx = SigmaModel(model[0], Z) + noise(rng);
		const double Sy = SigmaModel(model[1], Z) + noise(rng);
		res[i]          = Vec3D{.x = Sx, .y = Sy, .z = Z};
	}

	return res;
}

//-------------------------------------------------------------------------------------------------
TEST(AstigFitTests, Calibrate)
{
	const AstigFit Fit(PIXEL_SIZE, SAMPLING);
	const std::vector<Vec3D> Dataset = GetDataSet();

	const astig_model Model = Fit.Calibrate(Dataset);

	for (size_t ax = 0; ax < 2; ++ax) {
		for (size_t k = 0; k < 5; ++k) {
			std::cout << "Param axis " << ax << " [" << k << "] : fitted = " << Model[ax][k] << ", true = " << REF_MODEL[ax][k] << '\n';
			EXPECT_NEAR(Model[ax][k], REF_MODEL[ax][k], std::abs(0.1 * REF_MODEL[ax][k]) + 5.0); // tol 10% + offset 5 unités Vérifs grossières
		}
	}
}

//-------------------------------------------------------------------------------------------------
TEST(AstigFitTests, Estimate)
{
	const AstigFit Fit(PIXEL_SIZE, SAMPLING);
	// création d'un dataset de 10 éléments entre -200 et 200 nm à partir du modele calibré.
	const std::vector<Vec3D> Datas = GetDataSet(REF_MODEL, 200, 10);
	std::vector<Vec2D> sigmas(Datas.size());
	for (size_t i = 0; i < Datas.size(); ++i) { sigmas[i] = Vec2D{.x = Datas[i].x, .y = Datas[i].y}; }
	const auto Res = Fit.Estimate(REF_MODEL, sigmas);

	for (size_t i = 0; i < Datas.size(); ++i) {
		std::cout << "Point " << Datas[i] << ", Estimated z: " << Res[i] << " nm\n";
		EXPECT_NEAR(Res[i], Datas[i].z, 50.0); // tol 50 nm très très large.
	}
}

//-------------------------------------------------------------------------------------------------
TEST(AstigFitTests, CalibrateBad)
{
	const AstigFit Fit(PIXEL_SIZE, SAMPLING);

	std::mt19937 rng(42);
	std::normal_distribution noise(0.0, 1.0); // bruit en px

	// Dataset trop petit
	std::vector<Vec3D> dataset(2);
	for (auto& p : dataset) { p = Vec3D{.x = noise(rng), .y = noise(rng), .z = noise(rng)}; }

	const astig_model Model1 = Fit.Calibrate(dataset);

	for (size_t ax = 0; ax < 2; ++ax) { for (size_t k = 0; k < 5; ++k) { EXPECT_NEAR(Model1[ax][k], -1, 1e-9); } }


	// Dataset trop aléatoire, l'algo devra sans doute tourner plus longtemps.
	dataset.resize(N_POINTS);
	for (auto& p : dataset) { p = Vec3D{.x = noise(rng), .y = noise(rng), .z = noise(rng)}; }
	[[maybe_unused]] const astig_model Model2 = Fit.Calibrate(dataset);
}

//-------------------------------------------------------------------------------------------------
TEST(AstigFitTests, EstimateBad)
{
	const AstigFit Fit(PIXEL_SIZE, SAMPLING);

	std::mt19937 rng(42);
	std::normal_distribution noise(0.0, 1.0); // bruit en px

	// Dataset trop aléatoire, l'estimation sera plus complexe.
	std::vector<Vec2D> datas(N_POINTS);
	for (auto& p : datas) { p = Vec2D{.x = noise(rng), .y = noise(rng)}; }

	[[maybe_unused]] const auto Res = Fit.Estimate(REF_MODEL, datas);
}
} // namespace palm_tracer::fit
