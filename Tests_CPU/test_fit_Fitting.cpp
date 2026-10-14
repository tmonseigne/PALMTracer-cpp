#include "TestUtils.hpp"
#include "fit/Fitting.hpp"
#include "fit/GaussianFit.hpp"

#include <gtest/internal/gtest-port.h>

#include <array>
#include <cmath>
#include <numbers>
#include <numeric>

namespace palm_tracer::fit {
namespace {
/// <summary>Expose uniquement les calculs proteges de la classe de base, sans modifier l'API de production.</summary>
class FittingAccess final : public Fitting {
public:
	using Fitting::Fitting;
	using Fitting::fill;
	using Fitting::preCompute;
};

/// <summary>Paramètres physiques connus d'une tache gaussienne, indépendants de PSF().</summary>
struct GaussianCase {
	const char* name; ///< Nom du scénario dans les rapports.
	fit_mode mode;    ///< Mode d'ajustement.
	double sigmaX;    ///< Écart-type horizontal en pixels.
	double sigmaY;    ///< Écart-type vertical en pixels.
	double theta;     ///< Angle en radians.
};

/// <summary>Affiche un nom lisible avec les anciennes versions de GoogleTest.</summary>
void PrintTo(const GaussianCase& value, std::ostream* stream) { *stream << value.name; }

/// <summary>Évalue une gaussienne par rotation des coordonnées plutot que par les coefficients internes.</summary>
double Model(const double x, const double y, const double x0, const double y0, const double amplitude, const double background,
             const double sigmaX, const double sigmaY, const double theta)
{
	const double Dx = x - x0, Dy = y - y0;
	const double U  = std::cos(theta) * Dx - std::sin(theta) * Dy;
	const double V  = std::sin(theta) * Dx + std::cos(theta) * Dy;
	return background + amplitude * std::exp(-0.5 * (U * U / (sigmaX * sigmaX) + V * V / (sigmaY * sigmaY)));
}

/// <summary>Un seul candidat de centroïde entier permet d'isoler le fitting de la détection.</summary>
localization::Segmentation2D Candidate(const Vec2I& dim)
{
	const std::vector<uint16_t> Image(dim.GetSize(), 0);
	localization::ATrous decomposition(Image.data(), dim);
	auto& potential = decomposition.wDecomp_[PLANE_ID];
	for (size_t y = 6; y <= 8; ++y) { for (size_t x = 6; x <= 8; ++x) { potential[y * dim.x + x] = 100; } }
	return localization::Segmentation2D(decomposition, 0, 0);
}

/// <summary>Image quantifiee uint16 d'une tache sous-pixel, sans bruit aleatoire.</summary>
std::vector<uint16_t> GaussianImage(const GaussianCase& scenario, const Vec2I& dim)
{
	std::vector<uint16_t> image(dim.GetSize());
	for (size_t y = 0; y < dim.y; ++y) {
		for (size_t x = 0; x < dim.x; ++x) {
			image[y * dim.x + x] = static_cast<uint16_t>(std::lround(Model(static_cast<double>(x), static_cast<double>(y),
			                                                               7.25, 6.7, 8000, 100, scenario.sigmaX, scenario.sigmaY, scenario.theta)));
		}
	}
	return image;
}

/// <summary>Mêmes assertions de récupération pour chaque modèle gaussien.</summary>
class GaussianFitTests : public testing::TestWithParam<GaussianCase> { };

/// <summary>Conventions de remplissage pour les modèles sans ajustement et gaussiens.</summary>
class FittingOutputTests : public testing::TestWithParam<fit_mode> { };

/// <summary>Capacités inférieure, égale et supérieure au nombre de candidats.</summary>
class FittingCapacityTests : public testing::TestWithParam<size_t> { };
}

#pragma region Fit
/// <summary>Une troncature produit un seul avertissement ; une sortie suffisante reste silencieuse.</summary>
TEST_P(FittingCapacityTests, WarnsOnlyWhenTruncated)
{
	constexpr Vec2I Dim{.x = 15, .y = 15};
	const std::vector<uint16_t> Image(Dim.GetSize(), 0);
	constexpr std::array Params{7.0};
	const auto Segmentation = Candidate(Dim);
	ASSERT_EQ(Segmentation.Infos().size(), 1);
	const Fitting Fit(fit_mode::none, Image.data(), Dim, GetParam(), Params.data());
	std::vector<std::array<double, N_LOCALIZATIONS>> points;
	testing::internal::CaptureStderr();
	Fit.FitAll(Segmentation, points);
	const auto Warning = testing::internal::GetCapturedStderr();
	if (GetParam() == 0) { EXPECT_EQ(Warning, "[WARNING] Density is too high; some points are missed (output capacity exceeded).\n"); }
	else { EXPECT_TRUE(Warning.empty()); }
}

INSTANTIATE_TEST_SUITE_P(Capacities, FittingCapacityTests, testing::Values(size_t{0}, size_t{1}, size_t{2}),
                         ([](const testing::TestParamInfo<size_t>& info) {
	                         constexpr std::array Names{"Truncated", "Exact", "Spare"};
	                         return Names[info.param];
	                         }));

/// <summary>La capacité limite les lignes effectivement calculées, y compris lors d'une réutilisation.</summary>
TEST(FittingTests, RespectsOutputCapacity)
{
	constexpr Vec2I Dim{.x = 15, .y = 15};
	const std::vector<uint16_t> Image(Dim.GetSize(), 0);
	constexpr std::array Params{7.0};
	const auto Segmentation = Candidate(Dim);
	ASSERT_EQ(Segmentation.Infos().size(), 1);
	std::vector<std::array<double, N_LOCALIZATIONS>> points(3);
	const Fitting Full(fit_mode::none, Image.data(), Dim, 1, Params.data());
	Full.FitAll(Segmentation, points);
	ASSERT_EQ(points.size(), 1);
	EXPECT_DOUBLE_EQ(points[0][4], 7);
	EXPECT_DOUBLE_EQ(points[0][5], 7);
	const Fitting Empty(fit_mode::none, Image.data(), Dim, 0, Params.data());
	Empty.FitAll(Segmentation, points);
	EXPECT_TRUE(points.empty());
}
#pragma endregion Fit

#pragma region preCompute
/// <summary>Vérifie chaque coordonnee et valeur d'une ROI décalée dans une image non carrée.</summary>
TEST(FittingTests, ExtractsRoi)
{
	constexpr Vec2I Dim{.x = 7, .y = 5};
	std::vector<uint16_t> image(Dim.GetSize());
	for (size_t y = 0; y < Dim.y; ++y) { for (size_t x = 0; x < Dim.x; ++x) { image[y * Dim.x + x] = static_cast<uint16_t>(100 * y + x); } }
	constexpr std::array Params{3.0};
	const FittingAccess Fit(fit_mode::none, image.data(), Dim, 100, Params.data());
	std::vector<Pixel> grid(9);
	double min = 0, max = 0;
	Fit.preCompute(4, 2, grid, min, max);
	for (size_t i = 0; i < grid.size(); ++i) {
		const size_t X = 3 + i % 3, Y = 1 + i / 3;
		EXPECT_EQ(grid[i].x, X);
		EXPECT_EQ(grid[i].y, Y);
		EXPECT_DOUBLE_EQ(grid[i].v, static_cast<double>(100 * Y + X));
	}
	EXPECT_DOUBLE_EQ(min, 103);
	EXPECT_DOUBLE_EQ(max, 305);
}
#pragma endregion preCompute

#pragma region fill
/// <summary>Contrôle les champs de sortie, l'integrale analytique et la circularite.</summary>
TEST_P(FittingOutputTests, LocalizationFields)
{
	constexpr Vec2I Dim{.x = 7, .y = 5};
	const std::vector<uint16_t> Image(Dim.GetSize(), 0);
	constexpr std::array Params{3.0};
	const FittingAccess Fit(GetParam(), Image.data(), Dim, 100, Params.data());
	localization::SegmentationInfo info{};
	info.Reset();
	info.intensity_ = 123;
	info.surf_      = 9;
	std::array<double, N_LOCALIZATIONS> actual{};
	Fit.fill(info, actual, 3.25, 2.5, 0, 7, 120, 1.2, 2.4, 0.3, 0.98, -1);
	const std::array<double, N_LOCALIZATIONS> Expected{
		0, 0, 0, -1, 3.25, 2.5, 0,
		GetParam() == fit_mode::none ? 123 : 2 * std::numbers::pi * 120 * 1.2 * 2.4,
		1.2, 2.4, 0.3, 0.98, -1, 120, 7, 123, 9, 0.5
	};
	for (size_t i = 0; i < Expected.size(); ++i) { EXPECT_DOUBLE_EQ(actual[i], Expected[i]) << "Champ " << i; }
}

INSTANTIATE_TEST_SUITE_P(Models, FittingOutputTests,
                         testing::Values(fit_mode::none, fit_mode::gaussian_2d_fixed, fit_mode::gaussian_2d_sigma, fit_mode::gaussian_3d_fixed, fit_mode::
	                         gaussian_3d_theta),
                         ([](const testing::TestParamInfo<fit_mode>& info) {
	                         constexpr std::array Names{"NoFit", "GaussianXY", "GaussianSigma", "GaussianSigmaXY", "GaussianTheta"};
	                         return Names[static_cast<size_t>(info.param)];
	                         }));
#pragma endregion fill

#pragma region GaussianFit
/// <summary>Retrouve une tache sous-pixel connue avec une tolérance tenant compte de la quantification uint16.</summary>
TEST_P(GaussianFitTests, RecoversParameters)
{
	constexpr Vec2I Dim{.x = 15, .y = 15};
	const auto& Scenario = GetParam();
	const auto Image     = GaussianImage(Scenario, Dim);
	const std::array Params{9.0, 1.2, 1.8, Scenario.mode == fit_mode::gaussian_3d_theta ? 0.2 : Scenario.theta};
	const GaussianFit Fit(Scenario.mode, Image.data(), Dim, 100, Params.data());
	std::vector<std::array<double, N_LOCALIZATIONS>> points;
	Fit.FitAll(Candidate(Dim), points);
	ASSERT_EQ(points.size(), 1);
	const auto& Point = points[0];
	EXPECT_NEAR(Point[4], 7.25, 0.01);
	EXPECT_NEAR(Point[5], 6.7, 0.01);
	EXPECT_NEAR(Point[8], Scenario.sigmaX, 0.02);
	EXPECT_NEAR(Point[9], Scenario.sigmaY, 0.02);
	EXPECT_NEAR(Point[13], 8000, 5);
	EXPECT_NEAR(Point[14], 100, 2);
	EXPECT_NEAR(Point[7], 2 * std::numbers::pi * 8000 * Scenario.sigmaX * Scenario.sigmaY, 100);
	EXPECT_GE(Point[11], 0.99999);
	if (Scenario.mode == fit_mode::gaussian_3d_theta) { EXPECT_NEAR(Point[10], Scenario.theta, 0.02); }
}

/// <summary>La colonne MSE XY contient la corrélation au carré, recalculee indépendamment sur la ROI.</summary>
TEST_P(GaussianFitTests, ScoreMatchesSquaredCorrelation)
{
	constexpr Vec2I Dim{.x = 15, .y = 15};
	const auto& Scenario = GetParam();
	const auto Image     = GaussianImage(Scenario, Dim);
	const std::array Params{9.0, 1.2, 1.8, Scenario.mode == fit_mode::gaussian_3d_theta ? 0.2 : Scenario.theta};
	const GaussianFit Fit(Scenario.mode, Image.data(), Dim, 100, Params.data());
	std::vector<std::array<double, N_LOCALIZATIONS>> points;
	Fit.FitAll(Candidate(Dim), points);
	ASSERT_EQ(points.size(), 1);
	const auto& Point = points[0];
	std::vector<double> observed, predicted;
	for (size_t y = 3; y <= 11; ++y) {
		for (size_t x = 3; x <= 11; ++x) {
			observed.push_back(Image[y * Dim.x + x]);
			predicted.push_back(Model(static_cast<double>(x), static_cast<double>(y), Point[4], Point[5], Point[13], Point[14], Point[8], Point[9], Point[10]));
		}
	}
	const double MeanD = std::accumulate(observed.begin(), observed.end(), 0.0) / static_cast<double>(observed.size());
	const double MeanM = std::accumulate(predicted.begin(), predicted.end(), 0.0) / static_cast<double>(predicted.size());
	double covariance  = 0, varianceD = 0, varianceM = 0;
	for (size_t i = 0; i < observed.size(); ++i) {
		const double D = observed[i] - MeanD, M = predicted[i] - MeanM;
		covariance += D * M;
		varianceD += D * D;
		varianceM += M * M;
	}
	ASSERT_GT(varianceD * varianceM, 0);
	EXPECT_NEAR(Point[11], covariance * covariance / (varianceD * varianceM), 1e-10);
}

INSTANTIATE_TEST_SUITE_P(Models, GaussianFitTests,
                         testing::Values(GaussianCase{"FixedSigma", fit_mode::gaussian_2d_fixed, 1.2, 1.2, 0},
	                         GaussianCase{"FreeSigma", fit_mode::gaussian_2d_sigma, 1.35, 1.35, 0},
	                         GaussianCase{"FixedAngle", fit_mode::gaussian_3d_fixed, 1.3, 1.9, 0},
	                         GaussianCase{"OneDegree", fit_mode::gaussian_3d_fixed, 1.3, 1.9, std::numbers::pi / 180},
	                         GaussianCase{"FreeAngle", fit_mode::gaussian_3d_theta, 1.3, 1.9, 0.4}),
                         [](const testing::TestParamInfo<GaussianCase>& info) { return info.param.name; });
#pragma endregion GaussianFit
} // namespace palm_tracer::fit
