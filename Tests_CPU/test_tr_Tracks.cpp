#include "Constants.hpp"
#include "TestUtils.hpp"
#include "tracking/Tracks.hpp"

#include <array>
#include <cmath>
#include <limits>
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

/// <summary>Expose les éléments nécessaires aux vérifications ciblées, sans étendre l'API de production.</summary>
class TracksAccess final : public Tracks {
public:
	using Tracks::t_;
	using Tracks::getCandidates;
	using Tracks::getPlaneRanges;
	using Tracks::loadPoints;
};

/// <summary>Variantes du suivi par mode d'ajustement, watershed et ordre d'entrée.</summary>
class TrackingTests : public testing::TestWithParam<TrackingParams> { };

/// <summary>Filtre numérique des candidats : rayon, distance, admissibilité et nom du cas.</summary>
class CandidateRadiusTests : public testing::TestWithParam<std::tuple<double, double, bool, std::string>> { };

/// <summary>Ex aequo ou second voisin disponible, dans les deux ordres d'entrée.</summary>
class TrackingGreedyTests : public testing::TestWithParam<std::tuple<bool, bool>> { };

/// <summary>Scénario de reconnexion avec chaînes attendues explicitement.</summary>
struct BlinkingScenario {
	std::string name;                              ///< Nom du comportement vérifié.
	std::vector<std::array<double, 3>> rows;        ///< Un segment par ligne : identifiant, plan et X.
	std::vector<std::vector<size_t>> expected;      ///< Identifiants des points, regroupés par chaîne attendue.
	size_t maxDuration = 4;                        ///< Délai maximal admissible.
	size_t emptyIndex  = INVALID_INDEX;            ///< Segment à vider pour vérifier son exclusion, le cas échéant.
};

/// <summary>Bornes temporelles et contraintes gloutonnes absentes des seules références CSV.</summary>
class BlinkingIndexTests : public testing::TestWithParam<BlinkingScenario> { };

/// <summary>Modes de reconnexion indépendants valides sur fichiers CSV.</summary>
class BlinkingTests : public testing::TestWithParam<size_t> { };

/// <summary>Priorité temporelle ou départage spatial, avec permutation de l'ordre des candidats.</summary>
class BlinkingPriorityTests : public testing::TestWithParam<std::tuple<bool, bool>> { };

/// <summary>Calculs calibres avec et sans coordonnee axiale.</summary>
class TrackAnalysisTests : public testing::TestWithParam<bool> { };

/// <summary>Ajustements non linéaires, complémentaires au modèle linéaire teste en 2D et 3D.</summary>
class TrackFitTests : public testing::TestWithParam<size_t> { };
} // namespace

#pragma region LoadTracks
/// <summary>Les identifiants clairsemes ne determinent pas les tailles d'allocation.</summary>
TEST(TracksTests, LoadsSparseGroups)
{
	constexpr std::array<double, 3 * N_TRACKING> Input{
		1000000000000, 1, 1, 2, 3, 4, 5, 6,
		1000000000000, 2, 2, 3, 4, 5, 6, 7,
		2000000000000, 3, 3, 4, 5, 6, 7, 8
	};
	TracksAccess tracks;
	tracks.LoadTracks(Input.data(), 3, 2);
	ASSERT_EQ(tracks.t_.size(), 2);
	EXPECT_LT(tracks.t_.capacity(), 10);
	EXPECT_EQ(tracks.t_[0].points.size(), 2);
	EXPECT_LT(tracks.t_[1].points.capacity(), 4);
	std::array<double, Input.size()> output{};
	EXPECT_EQ(tracks.FillOutput(output.data()), Input.size());
	EXPECT_EQ(output, Input);
	EXPECT_DOUBLE_EQ(tracks.t_[0].points[0].pos_.x, 4);
	EXPECT_DOUBLE_EQ(tracks.t_[0].points[0].pos_.z, 4);
	tracks.LoadTracks(nullptr, 0, 2);
	EXPECT_TRUE(tracks.t_.empty());
	EXPECT_EQ(tracks.FillOutput(nullptr), 0);
}
#pragma endregion LoadTracks

#pragma region Tracking
/// <summary>Le filtre conserve le seuil inclusif et les règles numériques du carré du rayon.</summary>
TEST_P(CandidateRadiusTests, FiltersBySquaredRadius)
{
	const auto& [radius, distance, accepted, name] = GetParam();
	const std::array Input{1.0, 1.0, 0.0, 0.0, 2.0, 2.0, distance, 0.0};
	const auto Points = TracksAccess::loadPoints(Input.data(), 2);
	std::vector<Candidate> candidates;
	TracksAccess::getCandidates(Points, {.begin = 0, .end = 1}, {.begin = 1, .end = 2}, radius, candidates);
	ASSERT_EQ(candidates.size(), accepted ? size_t{1} : size_t{0});
	if (accepted) {
		EXPECT_EQ(candidates[0].from, 0);
		EXPECT_EQ(candidates[0].to, 1);
		EXPECT_EQ(candidates[0].distance2, distance * distance);
	}
}

INSTANTIATE_TEST_SUITE_P(Radii, CandidateRadiusTests,
                         testing::Values(
	                         std::make_tuple(2.0, 2.0, true, "InclusiveBoundary"),
	                         std::make_tuple(2.0, std::nextafter(2.0, 3.0), false, "OutsideBoundary"),
	                         std::make_tuple(0.0, 0.0, true, "ZeroRadius"),
	                         std::make_tuple(-1.0, 0.5, true, "NegativeRadius"),
	                         std::make_tuple(1e-200, 1e-200, true, "SquaredUnderflow"),
	                         std::make_tuple(1e-140, 1e-140, true, "SmallNormalRadius"),
	                         std::make_tuple(1e140, 1e140, true, "LargeNormalRadius"),
	                         std::make_tuple(0.1, 0.1, true, "FractionalRadius"),
	                         std::make_tuple(std::numeric_limits<double>::infinity(), 3.0, true, "InfiniteRadius"),
	                         std::make_tuple(std::numeric_limits<double>::quiet_NaN(), 1.0, false, "NanRadius"),
	                         std::make_tuple(1.0, std::numeric_limits<double>::infinity(), false, "InfiniteCoordinate")),
                         [](const testing::TestParamInfo<CandidateRadiusTests::ParamType>& info) { return std::get<3>(info.param); });

/// <summary>Le même buffer conserve sa capacité et ne garde aucun candidat des transitions précédentes, sans conserver de candidat périmé.</summary>
TEST(TracksTests, ReusesCandidateBufferAcrossTransitions)
{
	const std::array Input{1.0, 1.0, 0.0, 0.0, 2.0, 1.0, 0.2, 0.0,
	                       3.0, 2.0, 0.1, 0.0, 4.0, 2.0, 0.3, 0.0,
	                       5.0, 3.0, 3.0, 0.0, 6.0, 4.0, 3.0, 0.0};
	const auto Points = TracksAccess::loadPoints(Input.data(), Input.size() / 4);
	const auto Planes = TracksAccess::getPlaneRanges(Points);
	std::vector<Candidate> candidates;
	candidates.reserve(64);
	const size_t Capacity = candidates.capacity();
	candidates.emplace_back(INVALID_INDEX, INVALID_INDEX, -1);
	TracksAccess::getCandidates(Points, Planes[0], Planes[1], 1, candidates);
	ASSERT_EQ(candidates.size(), 4);
	const Candidate* Storage = candidates.data();
	EXPECT_EQ(candidates.capacity(), Capacity);

	// Une transition sans voisin efface les quatre anciens candidats sans libérer la capacité.
	TracksAccess::getCandidates(Points, Planes[1], Planes[2], 1, candidates);
	EXPECT_TRUE(candidates.empty());
	EXPECT_EQ(candidates.capacity(), Capacity);

	// Le rayon nul remplit la même allocation avec une seule liaison exacte.
	TracksAccess::getCandidates(Points, Planes[2], Planes[3], 0, candidates);
	ASSERT_EQ(candidates.size(), 1);
	EXPECT_EQ(candidates[0].from, 4);
	EXPECT_EQ(candidates[0].to, 5);
	EXPECT_EQ(candidates.data(), Storage);
	EXPECT_EQ(candidates.capacity(), Capacity);

	const PlaneRange Empty{.begin = 0, .end = 0};
	TracksAccess::getCandidates(Points, Empty, Planes[1], 1, candidates);
	EXPECT_TRUE(candidates.empty());
	EXPECT_EQ(candidates.capacity(), Capacity);
	TracksAccess::getCandidates(Points, Planes[0], Empty, 1, candidates);
	EXPECT_TRUE(candidates.empty());
	EXPECT_EQ(candidates.capacity(), Capacity);
}

/// <summary>Les ex aequo suivent les identifiants ; un conflit laisse au second point son autre voisin admissible.</summary>
TEST_P(TrackingGreedyTests, ResolvesConflicts)
{
	const auto [fallback, reverse] = GetParam();
	std::array<std::array<double, 4>, 4> rows{{
		{1, 1, 0, 0}, {2, 1, fallback ? 2.2 : 0.0, 0},
		{3, 2, 1, 0}, {4, 2, fallback ? 3.7 : -1.0, 0}
	}};
	if (reverse) { std::ranges::reverse(rows); }
	std::vector<double> input;
	for (const auto& Row : rows) { input.insert(input.end(), Row.begin(), Row.end()); }
	std::array<uint64_t, 8> output{};
	const std::array<uint64_t, 8> Expected{1, 1, 2, 2, 3, 1, 4, 2};
	EXPECT_EQ(Tracks::Tracking(input.data(), output.data(), rows.size(), fallback ? 2 : 1), output.size());
	EXPECT_EQ(output, Expected);
}

INSTANTIATE_TEST_SUITE_P(Greedy, TrackingGreedyTests, testing::Combine(testing::Bool(), testing::Bool()),
                         [](const testing::TestParamInfo<TrackingGreedyTests::ParamType>& info) {
	                         return std::string(std::get<0>(info.param) ? "FallbackNeighbor" : "ExactTies")
	                         + (std::get<1>(info.param) ? "Reversed" : "Ordered");
                         });

/// <summary>Deux distances distinctes de moins d'un epsilon ne sont pas des ex aequo.</summary>
TEST(TracksTests, ResolvesCloseDistances)
{
	const double D = std::sqrt(std::numeric_limits<double>::epsilon() / 2);
	const std::array Input{1.0, 1.0, D, 0.0, 2.0, 1.0, 0.0, 0.0, 3.0, 2.0, 0.0, 0.0};
	std::array<uint64_t, 6> output{};
	Tracks::Tracking(Input.data(), output.data(), 3, 1);
	EXPECT_EQ(output[3], output[5]);
	EXPECT_NE(output[1], output[5]);
}

/// <summary>Vérifie les associations lorsque des plans sont absents.</summary>
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

/// <summary>Compare les associations Id/Track, y compris lorsque les localisations sont inversées.</summary>
TEST_P(TrackingTests, Links)
{
	const auto [fit, watershed, reverse] = GetParam();
	const std::string Suffix             = "103.6_" + std::string(watershed ? "True_" : "False_") + std::to_string(fit) + "_1.0_0.0_7";
	auto localizations                   = LoadCsv(InputPath / "ref" / ("stack-localizations-" + Suffix + ".csv"));
	// Seules les localisations d'intensité intégrée positive participent au suivi.
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

/// <summary>Une entrée sans point ne doit rien écrire dans le tampon de sortie.</summary>
TEST(TracksTests, TrackingEmptyInput)
{
	uint64_t output = 123;
	EXPECT_EQ(Tracks::Tracking(nullptr, &output, 0, 5), 0);
	EXPECT_EQ(output, 123);
}
#pragma endregion Tracking

#pragma region BlinkingReconnection

/// <summary>Vérifie les chaînes attendues sans recopier l'algorithme de reconnexion dans les tests.</summary>
TEST_P(BlinkingIndexTests, BuildsExpectedChains)
{
	const auto& Scenario = GetParam();
	std::vector<double> input;
	for (const auto& Row : Scenario.rows) {
		input.insert(input.end(), {Row[0], Row[1], Row[0], Row[2], 0, 0, 10, 1});
	}
	TracksAccess tracks;
	tracks.LoadTracks(input.data(), Scenario.rows.size(), 1);
	if (Scenario.emptyIndex != INVALID_INDEX) { tracks.t_[Scenario.emptyIndex].points.clear(); }
	tracks.BlinkingReconnection(0, Scenario.maxDuration, 2);
	ASSERT_EQ(tracks.t_.size(), Scenario.expected.size());
	for (size_t i = 0; i < Scenario.expected.size(); ++i) {
		EXPECT_EQ(tracks.t_[i].id, i + 1);
		std::vector<size_t> ids;
		for (const auto& Point : tracks.t_[i].points) { ids.push_back(Point.id_); }
		EXPECT_EQ(ids, Scenario.expected[i]);
	}
}

INSTANTIATE_TEST_SUITE_P(TemporalIndex, BlinkingIndexTests,
                         testing::Values(
	                         BlinkingScenario{"ZeroDuration", {{1, 1, 0}, {2, 5, 0}}, {{1}, {2}}, 0},
	                         BlinkingScenario{"OutsideWindow", {{1, 1, 0}, {2, 5, 0}}, {{1}, {2}}, 3},
	                         BlinkingScenario{"InclusiveWindow", {{1, 1, 0}, {2, 5, 0}}, {{1, 2}}},
	                         BlinkingScenario{"SameFrameExcluded", {{1, 1, 0}, {2, 1, 0}}, {{1}, {2}}},
	                         BlinkingScenario{"UnlimitedDuration", {{1, 1, 0}, {2, 1e12, 0}}, {{1, 2}}, std::numeric_limits<size_t>::max()},
	                         BlinkingScenario{"ExactTiesOrdered", {{1, 1, 0}, {2, 2, 1}, {3, 2, -1}}, {{1, 2}, {3}}},
	                         BlinkingScenario{"ExactTiesReversed", {{1, 1, 0}, {3, 2, -1}, {2, 2, 1}}, {{1, 3}, {2}}},
	                         BlinkingScenario{"TakenCandidate", {{1, 1, 0}, {2, 1, 0.1}, {3, 2, 0}}, {{1, 3}, {2}}},
	                         BlinkingScenario{"OriginalHeadOrder", {{1, 3, 0}, {2, 1, 0}, {3, 2, 0}}, {{1}, {2, 3}}},
	                         BlinkingScenario{"EmptySegment", {{1, 1, 0}, {2, 2, 0}, {3, 3, 0}}, {{1, 3}}, 4, 1},
	                         // Une distance NaN ne départage pas les candidats de même délai : le premier reste retenu.
	                         BlinkingScenario{"NanDistance", {{1, 1, 0}, {2, 2, std::numeric_limits<double>::quiet_NaN()}, {3, 2, 1}}, {{1, 2}, {3}}}),
                         [](const testing::TestParamInfo<BlinkingScenario>& info) { return info.param.name; });

/// <summary>La borne temporelle reste valide lorsque son addition au dernier plan déborderait size_t.</summary>
TEST(TracksTests, BlinkingNearPlaneLimit)
{
	const std::array<double, 2 * N_TRACKING> Input{1, 1, 1, 0, 0, 0, 10, 1, 2, 2, 2, 0, 0, 0, 10, 1};
	TracksAccess tracks;
	tracks.LoadTracks(Input.data(), 2, 1);
	// Affectation entière directe : un double ne représente pas exactement ces plans extrêmes.
	tracks.t_[0].points[0].plane_ = std::numeric_limits<size_t>::max() - 2;
	tracks.t_[1].points[0].plane_ = std::numeric_limits<size_t>::max() - 1;
	tracks.BlinkingReconnection(0, 10, 2);
	ASSERT_EQ(tracks.t_.size(), 1);
	ASSERT_EQ(tracks.t_[0].points.size(), 2);
	EXPECT_EQ(tracks.t_[0].points[0].id_, 1);
	EXPECT_EQ(tracks.t_[0].points[1].id_, 2);
}

/// <summary>Le candidat temporellement le plus proche gagne ; a même frame, la distance départage, quel que soit l'ordre de parcours.</summary>
TEST_P(BlinkingPriorityTests, CandidatePriority)
{
	const auto [sameFrame, reverse] = GetParam();
	std::array<std::array<double, N_TRACKING>, 3> input{
		{
			{1, 1, 1, 0, 0, 0, 10, 1},
			{2, sameFrame ? 2.0 : 3.0, 2, sameFrame ? 2.0 : 1.0, 0, 0, 10, 1},
			{3, 2, 3, sameFrame ? 1.0 : 2.0, 0, 0, 10, 1}
		}
	};
	if (reverse) { std::swap(input[1], input[2]); }
	// Les candidats sont deux groupes distincts ; le premier groupe est toujours la tête de chaîne.
	std::vector<double> flat;
	for (const auto& Row : input) { flat.insert(flat.end(), Row.begin(), Row.end()); }
	TracksAccess tracks;
	tracks.LoadTracks(flat.data(), input.size(), 1);
	tracks.BlinkingReconnection(0, 3, 3);
	ASSERT_FALSE(tracks.t_.empty());
	ASSERT_GE(tracks.t_[0].points.size(), 2);
	EXPECT_EQ(tracks.t_[0].points[0].id_, 1);
	EXPECT_EQ(tracks.t_[0].points[1].id_, 3);
}

INSTANTIATE_TEST_SUITE_P(Priorities, BlinkingPriorityTests, testing::Combine(testing::Bool(), testing::Bool()),
                         [](const testing::TestParamInfo<std::tuple<bool, bool>>& info) {
                         return std::string(std::get<0>(info.param) ? "SameFrameNearestDistance" : "EarlierFrameBeforeDistance")
                         + (std::get<1>(info.param) ? "Reversed" : "Ordered");
                         });

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
/// <summary>Vérifie tous les lags et les valeurs manquantes de trajectoires discontinues.</summary>
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

/// <summary>Compare tous les lags MSD aux références 2D et 3D.</summary>
TEST_P(TrackAnalysisTests, MSD)
{
	const auto Tracks = LoadTracksData(GetParam());
	ExpectTableNear(Tracks.ComputeAllMSD(), LoadCsv(InputPath / "ref" / (GetParam() ? "tracking2-MSD-True.csv" : "tracking2-MSD-False.csv")));
}

/// <summary>Compare le modèle linéaire et les métriques initiales aux références 2D et 3D.</summary>
TEST_P(TrackAnalysisTests, Fit)
{
	const auto Tracks = LoadTracksData(GetParam());
	const ComputesParams Params{.mode = fit_mode::linear, .fitLength = 4, .is3D = GetParam(), .pixelSize = 1, .exposureTime = 1};
	ExpectTableNear(Tracks.ComputeAllFitMetrics(Params), LoadCsv(InputPath / "ref" / (GetParam() ? "tracking2-Fit-True.csv" : "tracking2-Fit-False.csv")), 1e-5,
	                1e-5);
}

INSTANTIATE_TEST_SUITE_P(, TrackAnalysisTests, testing::Bool(),
                           [](const testing::TestParamInfo<bool>& info) { return info.param ? "Spatial3D" : "Planar2D"; });

/// <summary>La référence historique de diffusion instantanée est exprimee en logarithme decimal.</summary>
TEST(TracksTests, InstantDiffusion)
{
	const auto Tracks = LoadTracksData(true);
	constexpr ComputesParams Params{.mode = fit_mode::linear, .fitLength = 4, .is3D = true, .pixelSize = 1, .exposureTime = 1};
	auto expected = LoadCsv(InputPath / "ref/tracking2-InstantD-True.csv");
	for (auto& row : expected) { for (size_t i = 1; i < row.size(); ++i) { row[i] = std::pow(10.0, row[i]); } }
	ExpectTableNear(Tracks.ComputeAllInstantDiffusion(Params), expected, 1e-5);
}

/// <summary>Compare les modèles de puissance et exponentiel sur une trajectoire 2D commune.</summary>
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
	constexpr ComputesParams Params{.fitLength = 4};
	EXPECT_TRUE(Tracks.ComputeAllMSD().empty());
	EXPECT_TRUE(Tracks.ComputeAllInstantDiffusion(Params).empty());
	EXPECT_TRUE(Tracks.ComputeAllFitMetrics(Params).empty());
}
#pragma endregion Analyses
} // namespace palm_tracer::tracking
