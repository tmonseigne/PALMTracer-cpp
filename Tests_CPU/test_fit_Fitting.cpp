#include "TestUtils.hpp"
#include "fit/Fitting.hpp"
#include "fit/GaussianFit.hpp"

#include <gtest/internal/gtest-port.h>

#include <array>
#include <cmath>
#include <limits>
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
	const double U  = (std::cos(theta) * Dx) - (std::sin(theta) * Dy);
	const double V  = (std::sin(theta) * Dx) + (std::cos(theta) * Dy);
	return background + (amplitude * std::exp(-0.5 * ((U * U / (sigmaX * sigmaX)) + (V * V / (sigmaY * sigmaY)))));
}

/// <summary>Un seul candidat de centroïde entier permet d'isoler le fitting de la détection.</summary>
localization::Segmentation2D Candidate(const Vec2I& dim)
{
	const std::vector<uint16_t> Image(dim.GetSize(), 0);
	localization::ATrous decomposition(Image.data(), dim);
	auto& potential = decomposition.wDecomp_[PLANE_ID];
	for (size_t y = 6; y <= 8; ++y) { for (size_t x = 6; x <= 8; ++x) { potential[(y * dim.x) + x] = 100; } }
	return localization::Segmentation2D(decomposition, 0, 0);
}

/// <summary>Image quantifiee uint16 d'une tache sous-pixel, sans bruit aleatoire.</summary>
std::vector<uint16_t> GaussianImage(const GaussianCase& scenario, const Vec2I& dim)
{
	std::vector<uint16_t> image(dim.GetSize());
	for (size_t y = 0; y < dim.y; ++y) {
		for (size_t x = 0; x < dim.x; ++x) {
			image[(y * dim.x) + x] = static_cast<uint16_t>(std::lround(Model(static_cast<double>(x), static_cast<double>(y),
			                                                                 7.25, 6.7, 8000, 100, scenario.sigmaX, scenario.sigmaY, scenario.theta)));
		}
	}
	return image;
}

/// <summary>Mêmes assertions de récupération pour chaque modèle gaussien.</summary>
class GaussianFitTests : public testing::TestWithParam<GaussianCase> { };

/// <summary>Réemploi des buffers et appels concurrents pour chacun des modèles gaussiens.</summary>
class GaussianWorkspaceTests : public testing::TestWithParam<fit_mode> { };

/// <summary>Récupération de la forme après une initialisation circulaire à angle libre.</summary>
class GaussianCircularInitializationTests : public testing::TestWithParam<GaussianCase> { };

/// <summary>Conventions de remplissage pour les modèles sans ajustement et gaussiens.</summary>
class FittingOutputTests : public testing::TestWithParam<fit_mode> { };

/// <summary>Capacités inférieure, égale et supérieure au nombre de candidats.</summary>
class FittingCapacityTests : public testing::TestWithParam<size_t> { };
} // namespace

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
	std::vector<std::array<double, localization_columns::Count>> points;
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
	std::vector<std::array<double, localization_columns::Count>> points(3);
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
	for (size_t y = 0; y < Dim.y; ++y) { for (size_t x = 0; x < Dim.x; ++x) { image[(y * Dim.x) + x] = static_cast<uint16_t>((100 * y) + x); } }
	constexpr std::array Params{3.0};
	const FittingAccess Fit(fit_mode::none, image.data(), Dim, 100, Params.data());
	std::vector<Pixel> grid(9);
	double min = 0, max = 0;
	Fit.preCompute(4, 2, grid, min, max);
	for (size_t i = 0; i < grid.size(); ++i) {
		const size_t X = 3 + (i % 3), Y = 1 + (i / 3);
		EXPECT_EQ(grid[i].x, X);
		EXPECT_EQ(grid[i].y, Y);
		EXPECT_DOUBLE_EQ(grid[i].v, static_cast<double>((100 * Y) + X));
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
	std::array<double, localization_columns::Count> actual{};
	Fit.fill(info, actual, 3.25, 2.5, 0, 7, 120, 1.2, 2.4, 0.3, 0.98, -1);
	const std::array<double, localization_columns::Count> Expected{
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
/// <summary>L'angle initialement non identifiable n'empêche pas la récupération du centre et de la covariance.</summary>
TEST_P(GaussianCircularInitializationTests, RecoversCenterAndCovariance)
{
	constexpr Vec2I Dim{.x = 15, .y = 15};
	const auto& Scenario = GetParam();
	const auto Image     = GaussianImage(Scenario, Dim);
	constexpr std::array Params{9.0, 1.2, 1.2, 0.2};
	const GaussianFit Fit(fit_mode::gaussian_3d_theta, Image.data(), Dim, 1, Params.data());
	localization::SegmentationInfo info{};
	info.Reset();
	info.x_ = info.y_ = 7.0;
	std::array<double, localization_columns::Count> point{};
	Fit.Fit(info, point);
	for (const double Value : point) { EXPECT_TRUE(std::isfinite(Value)); }
	ASSERT_GT(point[8], 0);
	ASSERT_GT(point[9], 0);
	EXPECT_NEAR(point[4], 7.25, 0.01);
	EXPECT_NEAR(point[5], 6.7, 0.01);
	// La covariance est indépendante des échanges d'axes et de l'angle arbitraire à l'isotropie.
	const auto Covariance = [](const double sx, const double sy, const double theta) {
		const double C = std::cos(theta), S = std::sin(theta), X2 = sx * sx, Y2 = sy * sy;
		return std::array{(C * C * X2) + (S * S * Y2), (S * S * X2) + (C * C * Y2), C * S * (Y2 - X2)};
	};
	const auto Actual   = Covariance(point[8], point[9], point[10]);
	const auto Expected = Covariance(Scenario.sigmaX, Scenario.sigmaY, Scenario.theta);
	for (size_t i = 0; i < Actual.size(); ++i) { EXPECT_NEAR(Actual[i], Expected[i], 0.02); }
	EXPECT_NEAR(point[13], 8000, 5);
	EXPECT_NEAR(point[14], 100, 2);
}

INSTANTIATE_TEST_SUITE_P(Shapes, GaussianCircularInitializationTests,
                         testing::Values(GaussianCase{"Circular", fit_mode::gaussian_3d_theta, 1.35, 1.35, 0},
	                         GaussianCase{"AlmostCircular", fit_mode::gaussian_3d_theta, 1.35, 1.37, 0.4},
	                         GaussianCase{"Rotated", fit_mode::gaussian_3d_theta, 1.3, 1.9, 0.4}),
                         [](const testing::TestParamInfo<GaussianCase>& info) { return info.param.name; });

/// <summary>Une initialisation non finie retourne les sentinelles usuelles, sans exception dans OpenMP.</summary>
TEST(GaussianFitFailureTests, RejectsNonFiniteInitialization)
{
	constexpr Vec2I Dim{.x = 15, .y = 15};
	const auto Image = GaussianImage({.name = "Circular", .mode = fit_mode::gaussian_3d_theta, .sigmaX = 1.35, .sigmaY = 1.35, .theta = 0}, Dim);
	constexpr std::array Params{9.0, std::numeric_limits<double>::quiet_NaN(), 1.8, 0.0};
	const GaussianFit Fit(fit_mode::gaussian_3d_theta, Image.data(), Dim, 1, Params.data());
	localization::SegmentationInfo info{};
	info.Reset();
	info.x_ = info.y_ = 7.0;
	std::array<double, localization_columns::Count> point{};
	EXPECT_NO_THROW(Fit.Fit(info, point));
	for (const double Value : point) { EXPECT_TRUE(std::isfinite(Value)); }
	for (const size_t Column : {size_t{7}, size_t{8}, size_t{9}, size_t{10}, size_t{11}, size_t{13}, size_t{14}}) { EXPECT_DOUBLE_EQ(point[Column], -1.0); }
}

/// <summary>Retrouve une tache sous-pixel connue avec une tolérance tenant compte de la quantification uint16.</summary>
TEST_P(GaussianFitTests, RecoversParameters)
{
	constexpr Vec2I Dim{.x = 15, .y = 15};
	const auto& Scenario = GetParam();
	const auto Image     = GaussianImage(Scenario, Dim);
	const std::array Params{9.0, 1.2, 1.8, Scenario.mode == fit_mode::gaussian_3d_theta ? 0.2 : Scenario.theta};
	const GaussianFit Fit(Scenario.mode, Image.data(), Dim, 100, Params.data());
	std::vector<std::array<double, localization_columns::Count>> points;
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
	std::vector<std::array<double, localization_columns::Count>> points;
	Fit.FitAll(Candidate(Dim), points);
	ASSERT_EQ(points.size(), 1);
	const auto& Point = points[0];
	std::vector<double> observed, predicted;
	for (size_t y = 3; y <= 11; ++y) {
		for (size_t x = 3; x <= 11; ++x) {
			observed.push_back(Image[(y * Dim.x) + x]);
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
/// <summary>Une ROI retrouve le même résultat après des échecs, un changement de modèle et de taille de ROI.</summary>
TEST_P(GaussianWorkspaceTests, ReuseAfterFailureAndConfigurationChange)
{
	constexpr Vec2I Dim{.x = 15, .y = 15};
	const auto Mode      = GetParam();
	const bool Isotropic = Mode == fit_mode::gaussian_2d_fixed || Mode == fit_mode::gaussian_2d_sigma;
	const GaussianCase Scenario{
		.name = "Workspace", .mode = Mode, .sigmaX = Mode == fit_mode::gaussian_2d_fixed ? 1.2 : 1.35,
		.sigmaY = Isotropic ? (Mode == fit_mode::gaussian_2d_fixed ? 1.2 : 1.35) : 1.9,
		.theta = Isotropic ? 0.0 : 0.4
	};
	const auto Original = GaussianImage(Scenario, Dim);
	auto image          = Original;
	const std::array Params{9.0, 1.2, 1.8, Mode == fit_mode::gaussian_3d_theta ? 0.2 : Scenario.theta};
	const GaussianFit Fit(Mode, image.data(), Dim, 1, Params.data());
	localization::SegmentationInfo info{};
	info.x_ = info.y_ = 7;
	std::array<double, localization_columns::Count> expected{}, actual{}, discarded{};
	Fit.Fit(info, expected);
	ASSERT_GT(expected[7], 0);

	// Une image plate provoque un échec après calcul du score, qui modifie la grille en place.
	std::fill(image.begin(), image.end(), uint16_t{100});
	Fit.Fit(info, discarded);
	EXPECT_DOUBLE_EQ(discarded[7], -1);
	std::copy(Original.begin(), Original.end(), image.begin());

	// Le même travailleur traverse ensuite un autre nombre de pixels/paramètres et un retour anticipé.
	constexpr std::array OtherParams{13.0, std::numeric_limits<double>::quiet_NaN(), 1.8, 0.2};
	const GaussianFit Other(fit_mode::gaussian_3d_theta, image.data(), Dim, 1, OtherParams.data());
	Other.Fit(info, discarded);
	EXPECT_DOUBLE_EQ(discarded[7], -1);
	auto border = info;
	border.x_   = border.y_ = 0;
	Fit.Fit(border, discarded);
	EXPECT_DOUBLE_EQ(discarded[7], -1);
	Fit.Fit(info, actual);
	for (size_t i = 0; i < actual.size(); ++i) { EXPECT_DOUBLE_EQ(actual[i], expected[i]) << "Colonne " << i; }
}

/// <summary>Les appels concurrents sur le même objet gardent les résultats indépendants des autres ROI.</summary>
TEST_P(GaussianWorkspaceTests, ConcurrentCallsMatchSerialResults)
{
	constexpr Vec2I Dim{.x = 15, .y = 15};
	const auto Mode      = GetParam();
	const bool Isotropic = Mode == fit_mode::gaussian_2d_fixed || Mode == fit_mode::gaussian_2d_sigma;
	const GaussianCase Scenario{
		.name = "Parallel", .mode = Mode, .sigmaX = Mode == fit_mode::gaussian_2d_fixed ? 1.2 : 1.35,
		.sigmaY = Isotropic ? (Mode == fit_mode::gaussian_2d_fixed ? 1.2 : 1.35) : 1.9,
		.theta = Isotropic ? 0.0 : 0.4
	};
	const auto Image = GaussianImage(Scenario, Dim);
	const std::array Params{9.0, 1.2, 1.8, Mode == fit_mode::gaussian_3d_theta ? 0.2 : Scenario.theta};
	const GaussianFit Fit(Mode, Image.data(), Dim, 16, Params.data());
	std::array<localization::SegmentationInfo, 16> infos{};
	std::array<std::array<double, localization_columns::Count>, 16> expected{}, actual{};
	for (size_t i = 0; i < infos.size(); ++i) {
		infos[i].x_         = 7.0 + (0.1 * static_cast<double>(i % 3));
		infos[i].y_         = 7.0 - (0.1 * static_cast<double>(i % 2));
		infos[i].intensity_ = static_cast<double>(100 + i);
		infos[i].surf_      = static_cast<double>(9 + i);
		Fit.Fit(infos[i], expected[i]);
	}
	#pragma omp parallel for num_threads(4)
	for (int i = 0; i < static_cast<int>(infos.size()); ++i) { Fit.Fit(infos[i], actual[i]); }
	for (size_t i = 0; i < infos.size(); ++i) {
		for (size_t j = 0; j < actual[i].size(); ++j) { EXPECT_DOUBLE_EQ(actual[i][j], expected[i][j]) << "ROI " << i << ", colonne " << j; }
	}
}

INSTANTIATE_TEST_SUITE_P(Models, GaussianWorkspaceTests,
                         testing::Values(fit_mode::gaussian_2d_fixed, fit_mode::gaussian_2d_sigma, fit_mode::gaussian_3d_fixed, fit_mode::gaussian_3d_theta),
                         ([](const testing::TestParamInfo<fit_mode>& info) {
	                         constexpr std::array Names{"FixedSigma", "FreeSigma", "FixedAngle", "FreeAngle"};
	                         return Names[static_cast<size_t>(info.param) - 1];
	                         }));

#pragma endregion GaussianFit
} // namespace palm_tracer::fit
