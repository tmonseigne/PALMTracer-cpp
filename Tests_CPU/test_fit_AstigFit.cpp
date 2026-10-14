#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
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

constexpr std::array<EstimateTestParams, 9> ESTIMATE_CASES = {
	{
		{-600.0, -500.0, "BelowLowerBound"},
		{-500.0, -500.0, "AtLowerBound"},
		{-499.0, -499.0, "NearLowerBound"},
		{-123.4, -123.4, "NegativeBetweenGridPoints"},
		{0.0, 0.0, "AtCenter"},
		{123.4, 123.4, "PositiveBetweenGridPoints"},
		{499.0, 499.0, "NearUpperBound"},
		{500.0, 500.0, "AtUpperBound"},
		{600.0, 500.0, "AboveUpperBound"}
	}
};

/// <summary>Tests paramétrés des limites et du raffinement de l'estimation axiale.</summary>
class AstigFitEstimateTests : public testing::TestWithParam<EstimateTestParams> { };

/// <summary>Accès aux calculs protégés pour isoler les critères de calibration et de qualité axiale.</summary>
class AstigFitAccess final : public AstigFit {
public:
	using AstigFit::AstigFit;
	using AstigFit::getQuality;
	using AstigFit::fitAxis;
	using AstigFit::lm;
	using AstigFit::computeSSE;
	using AstigFit::solver;
	using AstigFit::LAMBDA_0;
	using AstigFit::TOL_STEP;
	using AstigFit::isValidAxis;
};

/// <summary>Écart relatif identique sur les deux axes et score analytique attendu.</summary>
struct QualityTestParams {
	double relativeError; ///< Écart relatif appliqué aux deux sigmas observés.
	double expectedScore; ///< Score analytique attendu pour cet écart.
	const char* name;     ///< Identifiant anglais du scénario.
};

/// <summary>Vérification de l'échelle conventionnelle du score axial.</summary>
class AstigFitQualityTests : public testing::TestWithParam<QualityTestParams> { };

/// <summary>Validité continue d'un axe, indépendante de l'échantillonnage des données de calibration.</summary>
struct AxisDomainParams {
	std::array<double, 5> coefficients; ///< Paramètres [z0, w, c3, c4, A] de l'axe.
	std::array<double, 2> limits;       ///< Bornes de calibration en nanomètres.
	bool valid;                         ///< Validité attendue sur tout l'intervalle fermé.
	const char* name;                   ///< Identifiant anglais du scénario.
};

/// <summary>Scénarios indépendants de contrôle des bornes et extrema du polynôme.</summary>
class AstigFitDomainTests : public testing::TestWithParam<AxisDomainParams> { };

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
/// <summary>Un premier pas rejeté ne doit pas interrompre la calibration avant les essais plus amortis.</summary>
TEST(AstigFitTests, CalibrateAfterRejectedStep)
{
	const AstigFitAccess Fit(PIXEL_SIZE, SAMPLING);
	constexpr std::array Model{0.0, 100.0, 0.0, 0.0, 160.0};
	std::vector<double> z, sigma;
	for (int i = -10; i <= 10; ++i) {
		const double Z = 40.0 * i;
		z.push_back(Z);
		sigma.push_back(SigmaModel(Model, Z));
	}
	// L'heuristique choisit le minimum à z=0, W=étendue/4=200 et A=160, sans bruit.
	constexpr std::array Initial{0.0, 200.0, 0.0, 0.0, 160.0};
	std::vector<double> residue;
	std::vector<std::array<double, 5>> jacobian;
	const double InitialCost = Fit.lm(sigma, z, Initial, residue, jacobian);
	std::array<std::array<double, 5>, 5> alpha{};
	std::array<double, 5> beta{}, delta{}, trial{};
	for (size_t i = 0; i < z.size(); ++i) {
		for (size_t j = 0; j < 5; ++j) {
			beta[j] -= jacobian[i][j] * residue[i];
			for (size_t k = 0; k < 5; ++k) { alpha[j][k] += jacobian[i][j] * jacobian[i][k]; }
		}
	}
	for (size_t i = 0; i < 5; ++i) { alpha[i][i] *= 1.0 + AstigFitAccess::LAMBDA_0; }
	ASSERT_TRUE(AstigFitAccess::solver(alpha, beta, delta));
	for (size_t i = 0; i < 5; ++i) { trial[i] = Initial[i] + delta[i]; }
	trial[1] = std::max(trial[1], AstigFitAccess::TOL_STEP);
	trial[4] = std::max(trial[4], AstigFitAccess::TOL_STEP);
	// Vérifie la précondition de la régression : le premier candidat doit réellement être rejeté.
	ASSERT_GT(Fit.computeSSE(sigma, z, trial), InitialCost);
	const auto Fitted = Fit.fitAxis(sigma, z);
	for (const double Parameter : Fitted) { EXPECT_TRUE(std::isfinite(Parameter)); }
	EXPECT_LT(Fit.computeSSE(sigma, z, Fitted), 1e-10);
	EXPECT_TRUE(AstigFitAccess::isValidAxis(Fitted, z.front(), z.back()));
}

//-------------------------------------------------------------------------------------------------
/// <summary>Les bornes et extrema détectent aussi les polynômes invalides entre deux mesures positives.</summary>
TEST_P(AstigFitDomainTests, CalibrationInterval)
{
	const auto& [coefficients, limits, valid, name] = GetParam();
	SCOPED_TRACE(name);
	EXPECT_EQ(AstigFitAccess::isValidAxis(coefficients, limits[0], limits[1]), valid);
}

#ifndef INSTANTIATE_TEST_SUITE_P
#define INSTANTIATE_TEST_SUITE_P INSTANTIATE_TEST_CASE_P
#endif

INSTANTIATE_TEST_SUITE_P(
	Domains, AstigFitDomainTests,
	testing::Values(
		AxisDomainParams{{0, 200, 0, 0, 160}, {-400, 400}, true, "Quadratic"},
		AxisDomainParams{{0, 200, 0.1, 0, 160}, {-200, 200}, true, "Cubic"},
		AxisDomainParams{{0, 200, 0, -1, 160}, {-200, 200}, true, "ValidOnlyWithinLimits"},
		AxisDomainParams{{0, 200, 0, -1, 160}, {-400, 400}, false, "NegativeBoundary"},
		// P(0)=1 et P(3)=10, mais P(2)=-3 : les seules bornes ou mesures ne suffisent pas.
		AxisDomainParams{{-100, 100, -3, 1, 160}, {-100, 200}, false, "NegativeBetweenSamples"},
		// P(1)=P'(1)=0 : exclure également une largeur exactement nulle à un minimum intérieur.
		AxisDomainParams{{0, 200, -6, 4, 160}, {0, 400}, false, "ZeroInterior"},
		AxisDomainParams{{0, 200, -6, 4, 160}, {0, 200}, false, "ZeroBoundary"},
		AxisDomainParams{{0, 200, -5.99, 4, 160}, {0, 400}, true, "PositiveInteriorMinimum"}),
	[](const testing::TestParamInfo<AxisDomainParams>& info) { return std::string(info.param.name); });

//-------------------------------------------------------------------------------------------------
TEST(AstigFitTests, Estimate)
{
	const AstigFit Fit(PIXEL_SIZE, SAMPLING);
	// création d'un dataset de 10 éléments entre -200 et 200 nm à partir du modèle calibré.
	const std::vector<Vec3D> Datas = GetDataSet(REF_MODEL, 200, 10);
	std::vector<Vec2D> sigmas(Datas.size());
	for (size_t i = 0; i < Datas.size(); ++i) { sigmas[i] = Vec2D{.x = Datas[i].x, .y = Datas[i].y}; }
	const auto Res = Fit.Estimate(REF_MODEL, sigmas);

	for (size_t i = 0; i < Datas.size(); ++i) {
		std::cout << "Point " << Datas[i] << ", Estimated z: " << Res[i] << " nm\n";
		EXPECT_NEAR(Res[i].x, Datas[i].z, 40.0); // tol 40 nm très très large.
		EXPECT_GT(Res[i].y, 0.95);               // Bonne compatibilité des sigmas avec le modèle.
		EXPECT_LE(Res[i].y, 1.0);
	}
}

//-------------------------------------------------------------------------------------------------
TEST_P(AstigFitEstimateTests, EstimateWithinLimits)
{
	const AstigFit Fit(PIXEL_SIZE, SAMPLING);
	const auto& [observedZ, expectedZ, name] = GetParam();
	SCOPED_TRACE(name);
	const Vec2D Sigma = {.x = SigmaModel(ESTIMATE_MODEL[0], observedZ), .y = SigmaModel(ESTIMATE_MODEL[1], observedZ)};
	const Vec2D Res   = Fit.Estimate(ESTIMATE_MODEL, Sigma, {-500.0, 500.0});

	if (std::abs(expectedZ) == 500.0) { EXPECT_DOUBLE_EQ(Res.x, expectedZ); }
	else {
		EXPECT_NEAR(Res.x, expectedZ, 1.0);
		EXPECT_GT(Res.x, -500.0);
		EXPECT_LT(Res.x, 500.0);
	}

	// Le score doit correspondre au Z retenu, y compris en cas de saturation, et non au meilleur point de grille.
	const double PredictedX = SigmaModel(ESTIMATE_MODEL[0], Res.x), PredictedY = SigmaModel(ESTIMATE_MODEL[1], Res.x);
	const double Dx         = PredictedX - Sigma.x, Dy                         = PredictedY - Sigma.y;
	const double Rx         = Dx / PredictedX, Ry                              = Dy / PredictedY;
	EXPECT_NEAR(Res.y, 1.0 / (1.0 + (Rx * Rx + Ry * Ry) / (2 * 0.2 * 0.2)), 1e-12);
	EXPECT_GE(Res.y, 0.0);
	EXPECT_LE(Res.y, 1.0);
	if (std::abs(observedZ) > 500.0) { EXPECT_LT(Res.y, 1.0); }

	// Pour les positions entre deux points de grille, Brent doit réellement améliorer l'estimation initiale.
	if (std::abs(observedZ) == 123.4) {
		const double GridZ  = std::round(observedZ / 5.0) * 5.0; // Grille de 201 points sur [-500, 500].
		const double GridDx = SigmaModel(ESTIMATE_MODEL[0], GridZ) - Sigma.x;
		const double GridDy = SigmaModel(ESTIMATE_MODEL[1], GridZ) - Sigma.y;
		EXPECT_LT(Dx * Dx + Dy * Dy, GridDx * GridDx + GridDy * GridDy);
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
TEST_P(AstigFitQualityTests, RelativeSigmaScore)
{
	const AstigFitAccess Fit(PIXEL_SIZE, SAMPLING);
	const auto& [relativeError, expectedScore, name] = GetParam();
	SCOPED_TRACE(name);
	constexpr double Z   = 123.4;
	const Vec2D Observed = {
		.x = SigmaModel(ESTIMATE_MODEL[0], Z) * (1 + relativeError),
		.y = SigmaModel(ESTIMATE_MODEL[1], Z) * (1 + relativeError)
	};
	EXPECT_NEAR(Fit.getQuality(ESTIMATE_MODEL, Z, Observed), expectedScore, 1e-14);
}

INSTANTIATE_TEST_SUITE_P(
	RelativeErrors, AstigFitQualityTests,
	testing::Values(QualityTestParams{0.0, 1.0, "Exact"}, QualityTestParams{0.01, 400.0 / 401.0, "OnePercent"},
		QualityTestParams{0.05, 16.0 / 17.0, "FivePercent"}, QualityTestParams{0.10, 0.8, "TenPercent"},
		QualityTestParams{0.20, 0.5, "TwentyPercent"}, QualityTestParams{1.0, 1.0 / 26.0, "HundredPercent"},
		QualityTestParams{-0.20, 0.5, "NegativeTwentyPercent"}),
	[](const testing::TestParamInfo<QualityTestParams>& info) { return std::string(info.param.name); });

//-------------------------------------------------------------------------------------------------
/// <summary>Changer l'unité des sigmas ne change ni Z ni le score ; les deux axes contribuent à leur moyenne quadratique.</summary>
TEST(AstigFitTests, QualityScaleInvariant)
{
	const AstigFit Fit(PIXEL_SIZE, SAMPLING), Rescaled(PIXEL_SIZE / 2, SAMPLING);
	const Vec2D Sigma  = {.x = SigmaModel(ESTIMATE_MODEL[0], 123.4) * 1.05, .y = SigmaModel(ESTIMATE_MODEL[1], 123.4) * 0.98};
	const Vec2D First  = Fit.Estimate(ESTIMATE_MODEL, Sigma, {-500, 500});
	const Vec2D Second = Rescaled.Estimate(ESTIMATE_MODEL, {.x = 2 * Sigma.x, .y = 2 * Sigma.y}, {-500, 500});
	EXPECT_DOUBLE_EQ(First.x, Second.x);
	EXPECT_DOUBLE_EQ(First.y, Second.y);
	const double Rx = Sigma.x / SigmaModel(ESTIMATE_MODEL[0], First.x) - 1,
	             Ry = Sigma.y / SigmaModel(ESTIMATE_MODEL[1], First.x) - 1;
	EXPECT_NEAR(First.y, 1.0 / (1.0 + (Rx * Rx + Ry * Ry) / 0.08), 1e-14);
}

//-------------------------------------------------------------------------------------------------
/// <summary>Une largeur physiquement invalide ne produit jamais un bon score.</summary>
TEST(AstigFitTests, InvalidWidthQuality)
{
	const AstigFitAccess Fit(PIXEL_SIZE, SAMPLING);
	EXPECT_DOUBLE_EQ(Fit.getQuality(ESTIMATE_MODEL, 0, {.x = 0, .y = 1}), 0);
	EXPECT_DOUBLE_EQ(Fit.getQuality(ESTIMATE_MODEL, 0, {.x = 1, .y = -1}), 0);
	EXPECT_DOUBLE_EQ(Fit.getQuality(ESTIMATE_MODEL, 0, {.x = std::numeric_limits<double>::quiet_NaN(), .y = 1}), 0);
	auto invalidModel  = ESTIMATE_MODEL;
	invalidModel[0][4] = 0;
	EXPECT_DOUBLE_EQ(Fit.getQuality(invalidModel, 0, {.x = 1, .y = 1}), 0);
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
