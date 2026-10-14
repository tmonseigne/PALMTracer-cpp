#include "Constants.hpp"
#include "TestUtils.hpp"
#include "tracking/Tracks.hpp"

#include <array>
#include <cmath>
#include <map>
#include <tuple>

namespace palm_tracer::tracking {
namespace {
/// <summary>Colonnes du format de suivi commun aux deux langages.</summary>
const std::vector<size_t> TrackColumns{0, 1, 2, 3, 4, 5, 6, 7};

/// <summary>Charge les trajectoires et annule Z pour les calculs 2D.</summary>
Tracks LoadTracksData(const bool is3D)
{
	auto table = LoadCsv(InputPath / "tracking2.csv");
	if (!is3D) { for (auto& row : table) { row[5] = 0; } }
	const auto Input = Flatten(table, TrackColumns);
	Tracks result;
	result.LoadTracks(Input.data(), table.size(), 1.0);
	return result;
}

using TrackingParams = std::tuple<size_t, bool, bool>;

/// <summary>Variantes du suivi par mode d'ajustement, watershed et ordre d'entree.</summary>
class TrackingTests : public testing::TestWithParam<TrackingParams> { };

/// <summary>Modes de reconnexion independants valides sur fichiers CSV.</summary>
class BlinkingTests : public testing::TestWithParam<size_t> { };

/// <summary>Calculs calibres avec et sans coordonnee axiale.</summary>
class TrackAnalysisTests : public testing::TestWithParam<bool> { };

/// <summary>Ajustements non lineaires, complementaires au modele lineaire teste en 2D et 3D.</summary>
class TrackFitTests : public testing::TestWithParam<size_t> { };
} // namespace

#pragma region Tracking
/// <summary>Verifie les associations lorsque des plans sont absents.</summary>
TEST(TracksTests, TrackingWithGaps)
{
	const std::vector<double> InitialPoints = {
		// Plan 3
		1, 3, 0.0, 0.0,   // track 1 -> 4.
		2, 3, 10.0, 10.0, // seul, pas de correspondance.
		3, 3, 5.0, 5.0,   // conflit vers 6, mais perd contre 4.

		// Plan 4
		4, 4, 0.2, 0.0, // proche de 1 et 3, gagne avec 1.
		5, 4, 0.6, 0.0, // aussi proche de 1, mais 1 déjà utilisé.
		6, 4, 5.2, 5.0, // proche de 3 -> track 3.

		// Plan 6 : discontinuité, plan 5 vide
		7, 6, 0.3, 0.0,   // pas relié à 4/5, car plan 5 vide.
		8, 6, 20.0, 20.0, // seul.

		// Plan 7
		9, 7, 0.4, 0.0,   // relié à 7.
		10, 7, 20.2, 20.0 // relié à 8.
	};

	const std::vector<uint64_t> Ref = {
		1, 1, //		1->4
		2, 2, // 2 seul. |
		3, 3, //		 |	3->6
		4, 1, //		1->4 |
		5, 4, // 5 seul.	 |
		6, 3, //			3->6
		7, 5, //				7->9
		8, 6, //		8->10	 |
		9, 5, //		  |		7->9
		10, 6 //		8->10
	};

	constexpr size_t N = 10;
	std::vector<uint64_t> output(N * 2, 0);
	Tracks::Tracking(InitialPoints.data(), output.data(), N, 1);
	for (size_t i = 0; i < N * 2; ++i) { EXPECT_EQ(output[i], Ref[i]); }
}

/// <summary>Compare les associations Id/Track, y compris lorsque les localisations sont inversees.</summary>
TEST_P(TrackingTests, Links)
{
	const auto [fit, watershed, reverse] = GetParam();
	const std::string Suffix             = "103.6_" + std::string(watershed ? "True_" : "False_") + std::to_string(fit) + "_1.0_0.0_7";
	auto localizations                   = LoadCsv(InputPath / "ref" / ("stack-localizations-" + Suffix + ".csv"));
	// Seules les localisations d'intensite integree positive participent au suivi.
	std::erase_if(localizations, [](const auto& row) { return row[7] <= 0; });
	const auto Expected = LoadCsv(InputPath / "ref" / ("stack-tracking-" + Suffix + "-5.csv"));
	if (reverse) { std::ranges::reverse(localizations); }
	const auto Input = Flatten(localizations, {0, 1, 4, 5});
	std::vector<uint64_t> output(2 * localizations.size() + 1, UINT64_MAX);
	ASSERT_EQ(Tracks::Tracking(Input.data(), output.data(), localizations.size(), 5), 2 * localizations.size());
	EXPECT_EQ(output.back(), UINT64_MAX);
	std::map<uint64_t, uint64_t> expected, actual;
	for (const auto& Row : Expected) { ASSERT_TRUE(expected.emplace(static_cast<uint64_t>(Row[2]), static_cast<uint64_t>(Row[0])).second); }
	for (size_t i = 0; i < 2 * localizations.size(); i += 2) { ASSERT_TRUE(actual.emplace(output[i], output[i + 1]).second); }
	EXPECT_EQ(actual, expected);
}

INSTANTIATE_TEST_SUITE_P(, TrackingTests,
                           testing::Combine(testing::Values(size_t{0}, size_t{1}, size_t{2}, size_t{3}, size_t{4}), testing::Bool(), testing::Bool()),
                           ([](const testing::TestParamInfo<TrackingParams>& info) {
	                           constexpr std::array Names{"NoFit", "GaussianXY", "GaussianSigma", "GaussianSigmaXY", "GaussianTheta"};
	                           return std::string(Names[std::get<0>(info.param)]) + (std::get<1>(info.param) ? "WithWatershed" : "WithoutWatershed")
	                           + (std::get<2>(info.param) ? "Reversed" : "Ordered");
	                           }));

/// <summary>Une entree sans point ne doit rien ecrire dans le tampon de sortie.</summary>
TEST(TracksTests, TrackingEmptyInput)
{
	uint64_t output = 123;
	EXPECT_EQ(Tracks::Tracking(nullptr, &output, 0, 5), 0);
	EXPECT_EQ(output, 123);
}
#pragma endregion Tracking

#pragma region BlinkingReconnection
/// <summary>Compare toutes les trajectoires reconnectees aux valeurs attendues.</summary>
TEST_P(BlinkingTests, Reconnection)
{
	const auto Initial  = LoadCsv(InputPath / "tracking.csv");
	const auto Input    = Flatten(Initial, TrackColumns);
	const auto Expected = LoadCsv(InputPath / "ref" / ("tracking-blinking-" + std::to_string(GetParam()) + ".csv"));
	Tracks tracks;
	tracks.LoadTracks(Input.data(), Initial.size(), 1);
	tracks.BlinkingReconnection(GetParam(), 4, 2);
	std::vector<double> output(Input.size() + 1, -12345);
	ASSERT_EQ(tracks.FillOutput(output.data()), Input.size());
	EXPECT_DOUBLE_EQ(output.back(), -12345);
	Table actual;
	for (size_t i = 0; i < Input.size(); i += N_TRACKING) {
		actual.emplace_back(output.begin() + static_cast<ptrdiff_t>(i), output.begin() + static_cast<ptrdiff_t>(i + N_TRACKING));
	}
	ExpectTableNear(actual, Expected);
}

INSTANTIATE_TEST_SUITE_P(, BlinkingTests, testing::Range(size_t{0}, size_t{3}),
                           ([](const testing::TestParamInfo<size_t>& info) {
	                           constexpr std::array Names{"Stationary", "Diffusive", "Linear"};
	                           return Names[info.param];
	                           }));

/// <summary>La reconnexion d'un objet initialement vide conserve des sorties vides.</summary>
TEST(TracksTests, BlinkingEmptyState)
{
	Tracks tracks;
	tracks.BlinkingReconnection(0, 4, 2);
	EXPECT_EQ(tracks.FillOutput(nullptr), 0);
}
#pragma endregion BlinkingReconnection

#pragma region Analyses
/// <summary>Verifie tous les lags et les valeurs manquantes de trajectoires discontinues.</summary>
TEST(TracksTests, MSDWithGaps)
{
	const std::vector<double> InitialTracks = {
		1, 1, 1, 1, 1, 1, 1, 1,
		1, 5, 2, 2, 1, 1, 2, 2,

		2, 1, 3, 1, 1, 1, 3, 3,
		2, 2, 4, 2, 2, 1, 4, 4,
		2, 3, 5, 3, 0, 3, 7, 7,
		2, 4, 6, 7, 4, 5, 8, 8,

		3, 1, 7, 1, 1, 1, 11, 11,
		3, 2, 8, 3, 1, 1, 12, 12,
		3, 6, 9, 8, 1, 1, 13, 13,
		3, 7, 10, 2, 1, 1, 14, 14,
	};
	const std::vector<std::vector<double>> Ref = {
		{1, -1, -1, -1, 1},
		{2, 15.6666667, 27, 61},
		{3, 20, -1, -1, 25, 25, 1},
	};
	Tracks tracks;
	tracks.LoadTracks(InitialTracks.data(), 10, 1);
	const auto Res = tracks.ComputeAllMSD();
	ExpectTableNear(Res, Ref);
}

/// <summary>Compare tous les lags MSD aux references 2D et 3D.</summary>
TEST_P(TrackAnalysisTests, MSD)
{
	const auto Tracks = LoadTracksData(GetParam());
	ExpectTableNear(Tracks.ComputeAllMSD(), LoadCsv(InputPath / "ref" / (GetParam() ? "tracking2-MSD-True.csv" : "tracking2-MSD-False.csv")));
}

/// <summary>Compare le modele lineaire et les metriques initiales aux references 2D et 3D.</summary>
TEST_P(TrackAnalysisTests, Fit)
{
	const auto Tracks = LoadTracksData(GetParam());
	const ComputesParams Params{.mode = fit_mode::linear, .fitLength = 4, .is3D = GetParam(), .pixelSize = 1, .exposureTime = 1};
	ExpectTableNear(Tracks.ComputeAllFitMetrics(Params), LoadCsv(InputPath / "ref" / (GetParam() ? "tracking2-Fit-True.csv" : "tracking2-Fit-False.csv")), 1e-5,
	                1e-5);
}

INSTANTIATE_TEST_SUITE_P(, TrackAnalysisTests, testing::Bool(),
                           [](const testing::TestParamInfo<bool>& info) { return info.param ? "Spatial3D" : "Planar2D"; });

/// <summary>La reference historique de diffusion instantanee est exprimee en logarithme decimal.</summary>
TEST(TracksTests, InstantDiffusion)
{
	const auto Tracks = LoadTracksData(true);
	const ComputesParams Params{.mode = fit_mode::linear, .fitLength = 4, .is3D = true, .pixelSize = 1, .exposureTime = 1};
	auto expected = LoadCsv(InputPath / "ref/tracking2-InstantD-True.csv");
	for (auto& row : expected) { for (size_t i = 1; i < row.size(); ++i) { row[i] = std::pow(10.0, row[i]); } }
	ExpectTableNear(Tracks.ComputeAllInstantDiffusion(Params), expected, 1e-5);
}

/// <summary>Compare les modeles de puissance et exponentiel sur une trajectoire 2D commune.</summary>
TEST_P(TrackFitTests, Fit)
{
	constexpr std::array Modes{fit_mode::linear, fit_mode::power, fit_mode::exponential};
	const auto Tracks = LoadTracksData(false);
	const ComputesParams Params{.mode = Modes[GetParam()], .fitLength = 4, .is3D = false, .pixelSize = 1, .exposureTime = 1};
	ExpectTableNear(Tracks.ComputeAllFitMetrics(Params), LoadCsv(InputPath / "ref" / ("tracking2-Fit-" + std::to_string(GetParam() + 1) + ".csv")), 1e-5, 1e-3);
}

INSTANTIATE_TEST_SUITE_P(, TrackFitTests, testing::Range(size_t{1}, size_t{3}),
                           ([](const testing::TestParamInfo<size_t>& info) {
	                           constexpr std::array Names{"Linear", "Power", "Exponential"};
	                           return Names[info.param];
	                           }));

/// <summary>Les analyses d'un objet sans trajectoire donnent des tableaux vides.</summary>
TEST(TracksTests, AnalysisEmptyState)
{
	const Tracks Tracks;
	const ComputesParams Params{.fitLength = 4};
	EXPECT_TRUE(Tracks.ComputeAllMSD().empty());
	EXPECT_TRUE(Tracks.ComputeAllInstantDiffusion(Params).empty());
	EXPECT_TRUE(Tracks.ComputeAllFitMetrics(Params).empty());
}
#pragma endregion Analyses
} // namespace palm_tracer::tracking
