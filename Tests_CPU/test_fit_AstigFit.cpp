#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "Structs.hpp"
#include "fit/AstigFit.hpp"

namespace palm_tracer::fit {
constexpr double SAMPLING       = 1.0;
constexpr double PIXEL_SIZE     = 160.0;
constexpr double Z_MAX          = 600.0;
constexpr size_t N_POINTS       = 1000;
constexpr astig_model REF_MODEL = {
	std::array{-60.0, 300.0, 0.06, 0.02, 140.0},
	std::array{60.0, 300.0, -0.06, 0.02, 140.0}
};

/// <summary>Modèle dont les deux sigmas varient en sens opposés sur le domaine testé, pour une estimation axiale non ambiguë.</summary>
constexpr astig_model ESTIMATE_MODEL = {
	std::array{-1000.0, 1000.0, 0.0, 0.0, 160.0},
	std::array{1000.0, 1000.0, 0.0, 0.0, 160.0}
};

/// <summary>Scénario d'estimation axiale avec bornes de recherche fixes à ±500 nm.</summary>
struct EstimateTestParams {
	double observedZ; ///< <summary>Position utilisée pour générer les sigmas observés, en nm.</summary>
	double expectedZ; ///< <summary>Position attendue après estimation dans le domaine, en nm.</summary>
	const char* name; ///< <summary>Identifiant anglais du scénario Google Test.</summary>
};

constexpr std::array<EstimateTestParams, 9> ESTIMATE_CASES = {{
	{-600.0, -500.0, "BelowLowerBound"},
	{-500.0, -500.0, "AtLowerBound"},
	{-499.0, -499.0, "NearLowerBound"},
	{-123.4, -123.4, "NegativeBetweenGridPoints"},
	{0.0, 0.0, "AtCenter"},
	{123.4, 123.4, "PositiveBetweenGridPoints"},
	{499.0, 499.0, "NearUpperBound"},
	{500.0, 500.0, "AtUpperBound"},
	{600.0, 500.0, "AboveUpperBound"}
}};

/// <summary>Tests paramétrés des limites et du raffinement de l'estimation axiale.</summary>
class AstigFitEstimateTests : public testing::TestWithParam<EstimateTestParams> { };

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
	std::vector<Vec3D> dataset = GetDataSet();

	const astig_model Model = Fit.Calibrate(dataset, false);

	for (size_t ax = 0; ax < 2; ++ax) {
		for (size_t k = 0; k < 5; ++k) {
			std::cout << "Param axis " << ax << " [" << k << "] : fitted = " << Model[ax][k] << ", true = " << REF_MODEL[ax][k] << '\n';
			EXPECT_NEAR(Model[ax][k], REF_MODEL[ax][k], 5); // offset 5 unités Vérifs grossières
		}
	}

	// Modèle centré : on vérifie que les sigmas sont proches à z=0 (sinon on aurait un décalage systématique dans l'estimation de Z)
	for (auto& d : dataset) { d.z += 50.0; } // on décale le dataset de 50 nm pour éviter que le bruit ne masque un éventuel problème de centrage.
	const astig_model ModelC = Fit.Calibrate(dataset, true);
	for (size_t ax = 0; ax < 2; ++ax) {
		for (size_t k = 0; k < 5; ++k) {
			// tol plus large sur z0, on ne s'attend pas à ce qu'il soit exactement égal, mais il doit être proche.
			if (k == 0) { EXPECT_NEAR(Model[ax][k], ModelC[ax][k], 5); }
			else { EXPECT_NEAR(Model[ax][k], ModelC[ax][k], 0.01); }
		}
	}

	const double Sx0 = SigmaModel(ModelC[0], 0.0);
	const double Sy0 = SigmaModel(ModelC[1], 0.0);

	std::cout << "SigmaX(0) = " << Sx0 << '\n';
	std::cout << "SigmaY(0) = " << Sy0 << '\n';
	std::cout << "|diff| = " << std::abs(Sx0 - Sy0) << '\n';

	EXPECT_NEAR(Sx0, Sy0, 0.01); // ou 0.02 / 0.05 selon le bruit
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
		EXPECT_NEAR(Res[i].x, Datas[i].z, 40.0); // tol 40 nm très très large.
		EXPECT_LT(Res[i].y, 0.002); // Erreur par rapport au modèle très faible.
	}
}

//-------------------------------------------------------------------------------------------------
TEST_P(AstigFitEstimateTests, EstimateWithinLimits)
{
	const AstigFit Fit(PIXEL_SIZE, SAMPLING);
	const auto& [observedZ, expectedZ, name] = GetParam();
	SCOPED_TRACE(name);
	const Vec2D Sigma = {.x = SigmaModel(ESTIMATE_MODEL[0], observedZ), .y = SigmaModel(ESTIMATE_MODEL[1], observedZ)};
	const Vec2D Res = Fit.Estimate(ESTIMATE_MODEL, Sigma, {-500.0, 500.0});

	if (std::abs(expectedZ) == 500.0) { EXPECT_DOUBLE_EQ(Res.x, expectedZ); }
	else {
		EXPECT_NEAR(Res.x, expectedZ, 1.0);
		EXPECT_GT(Res.x, -500.0);
		EXPECT_LT(Res.x, 500.0);
	}

	// L'erreur retournée doit correspondre au Z effectivement retenu, y compris en cas de saturation.
	const double Dx = SigmaModel(ESTIMATE_MODEL[0], Res.x) - Sigma.x;
	const double Dy = SigmaModel(ESTIMATE_MODEL[1], Res.x) - Sigma.y;
	EXPECT_NEAR(Res.y, Dx * Dx + Dy * Dy, 1e-12);
	EXPECT_GE(Res.y, 0.0);
	if (std::abs(observedZ) > 500.0) { EXPECT_GT(Res.y, 0.0); }

	// Pour les positions entre deux points de grille, Brent doit réellement améliorer l'estimation initiale.
	if (std::abs(observedZ) == 123.4) {
		const double GridZ = std::round(observedZ / 5.0) * 5.0; // Grille de 201 points sur [-500, 500].
		const double GridDx = SigmaModel(ESTIMATE_MODEL[0], GridZ) - Sigma.x;
		const double GridDy = SigmaModel(ESTIMATE_MODEL[1], GridZ) - Sigma.y;
		EXPECT_LT(Res.y, GridDx * GridDx + GridDy * GridDy);
	}
}

#ifndef INSTANTIATE_TEST_SUITE_P
#define INSTANTIATE_TEST_SUITE_P INSTANTIATE_TEST_CASE_P
#endif

INSTANTIATE_TEST_SUITE_P(
	AxialLimits, AstigFitEstimateTests, testing::ValuesIn(ESTIMATE_CASES),
	[](const testing::TestParamInfo<EstimateTestParams>& info) { return std::string(info.param.name); });

//-------------------------------------------------------------------------------------------------
TEST(AstigFitTests, EstimateBatchWithinLimits)
{
	const AstigFit Fit(PIXEL_SIZE, SAMPLING);
	std::vector<Vec2D> sigmas;
	sigmas.reserve(ESTIMATE_CASES.size());
	for (const auto& scenario : ESTIMATE_CASES) {
		sigmas.push_back({.x = SigmaModel(ESTIMATE_MODEL[0], scenario.observedZ), .y = SigmaModel(ESTIMATE_MODEL[1], scenario.observedZ)});
	}

	const auto Res = Fit.Estimate(ESTIMATE_MODEL, sigmas, {-500.0, 500.0});
	ASSERT_EQ(Res.size(), ESTIMATE_CASES.size());
	for (size_t i = 0; i < Res.size(); ++i) {
		SCOPED_TRACE(ESTIMATE_CASES[i].name);
		const Vec2D Scalar = Fit.Estimate(ESTIMATE_MODEL, sigmas[i], {-500.0, 500.0});
		EXPECT_DOUBLE_EQ(Res[i].x, Scalar.x);
		EXPECT_DOUBLE_EQ(Res[i].y, Scalar.y);
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

	const astig_model Model1 = Fit.Calibrate(dataset, false);
	for (size_t ax = 0; ax < 2; ++ax) { for (size_t k = 0; k < 5; ++k) { EXPECT_NEAR(Model1[ax][k], -1, 1e-9); } }

	// Dataset trop aléatoire, l'algo devra sans doute tourner plus longtemps.
	dataset.resize(N_POINTS);
	for (auto& p : dataset) { p = Vec3D{.x = noise(rng), .y = noise(rng), .z = noise(rng)}; }
	[[maybe_unused]] const astig_model Model2 = Fit.Calibrate(dataset, false);
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
