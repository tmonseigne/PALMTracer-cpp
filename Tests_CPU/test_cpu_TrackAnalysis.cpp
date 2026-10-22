#include "TestUtils.hpp"
#include "../CPU/TrackAnalysis.hpp"

#include <array>
#include <span>
#include <tuple>

namespace palm_tracer::cpu {
namespace {
constexpr double Guard = -12345.0; ///< Valeur protégeant les bornes des sorties.

/// <summary>Décode les lignes irrégulières dans les limites du tampon fourni.</summary>
Table Decode(const std::span<const double> buffer)
{
	Table result;
	size_t offset = 0;
	while (offset < buffer.size() && buffer[offset] >= 0) {
		const auto Length = static_cast<size_t>(buffer[offset++]);
		if (Length > buffer.size() - offset) {
			ADD_FAILURE() << "Ligne irrégulière hors tampon";
			return {};
		}
		result.emplace_back(buffer.begin() + static_cast<ptrdiff_t>(offset), buffer.begin() + static_cast<ptrdiff_t>(offset + Length));
		offset += Length;
	}
	EXPECT_LT(offset, buffer.size());
	if (offset < buffer.size()) { EXPECT_DOUBLE_EQ(buffer[offset], -1.0); }
	return result;
}

/// <summary>Combinaisons indépendantes des analyses activables sans fit global.</summary>
class AnalysisSelectionTests : public testing::TestWithParam<std::tuple<bool, bool>> { };

/// <summary>Modes exposés par la DLL, y compris le repli historique sur le modèle linéaire.</summary>
class AnalysisFitModeTests : public testing::TestWithParam<uint64_t> { };
} // namespace

/// <summary>Vérifie la sélection des sorties, les trous temporels et une diffusion calculable analytiquement.</summary>
TEST_P(AnalysisSelectionTests, SelectsOutputsAndPreservesAnalyticalValues)
{
	const auto [MSD, Instant] = GetParam();
	const std::vector<double> Input{
		7, 1, 1, 0, 0, 0, 1, 1,
		7, 3, 2, 2, 0, 0, 1, 1,
		11, 1, 3, 0, 0, 0, 1, 1,
		11, 2, 4, 1, 0, 0, 1, 1,
		11, 3, 5, 2, 0, 0, 1, 1,
		11, 4, 6, 3, 0, 0, 1, 1,
		11, 5, 7, 4, 0, 0, 1, 1,
	};
	constexpr std::array Params{4.0};
	std::vector<double> msd(13, Guard), instant(6, Guard);
	ASSERT_TRUE(AnalyzeTracks(Input.data(), MSD ? msd.data() + 1 : nullptr, Instant ? instant.data() + 1 : nullptr,
		nullptr, 7, MSD, Instant, false, 1.0, 1.0, 0, Params.data()));
	if (MSD) { ExpectTableNear(Decode(std::span(msd).subspan(1, 11)), {{7, -1, 4}, {11, 1, 4, 9, 16}}, 1e-12); }
	if (Instant) {
		// La régression des MSD 1, 4, 9, 16 sur les lags 1 à 4 a une pente de 5 ; D = 5 / (4 × 1).
		ExpectTableNear(Decode(std::span(instant).subspan(1, 4)), {{11, 1.25}}, 1e-5);
	}
	EXPECT_DOUBLE_EQ(msd.front(), Guard);
	EXPECT_DOUBLE_EQ(msd.back(), Guard);
	EXPECT_DOUBLE_EQ(instant.front(), Guard);
	EXPECT_DOUBLE_EQ(instant.back(), Guard);
}

INSTANTIATE_TEST_SUITE_P(, AnalysisSelectionTests, testing::Combine(testing::Bool(), testing::Bool()),
                           [](const testing::TestParamInfo<std::tuple<bool, bool>>& info) {
                           return std::string(std::get<0>(info.param) ? "MSD" : "NoMSD") + (std::get<1>(info.param) ? "Instant" : "NoInstant");
                           });

/// <summary>Compare le choix du modèle et la sérialisation des métriques aux CSV existants.</summary>
TEST_P(AnalysisFitModeTests, MatchesExistingFitReference)
{
	auto input = LoadCsv(InputPath / "tracking2.csv");
	for (auto& row : input) { row[5] = 0; }
	const auto Values    = Flatten(input, {0, 1, 2, 3, 4, 5, 6, 7});
	const auto Mode      = GetParam();
	const auto Reference = Mode == 99 ? uint64_t{1} : Mode;
	const auto Expected  = LoadCsv(InputPath / "ref" / ("tracking2-Fit-" + std::to_string(Reference) + ".csv"));
	size_t capacity      = 1;
	for (const auto& Row : Expected) { capacity += Row.size() + 1; }
	std::vector<double> output(capacity + 2, Guard);
	constexpr std::array Params{4.0};
	ASSERT_TRUE(AnalyzeTracks(Values.data(), nullptr, nullptr, output.data() + 1, input.size(),
		false, false, false, 1.0, 1.0, Mode, Params.data()));
	ExpectTableNear(Decode(std::span(output).subspan(1, capacity)), Expected, 1e-5, 1e-3);
	EXPECT_DOUBLE_EQ(output.front(), Guard);
	EXPECT_DOUBLE_EQ(output.back(), Guard);
}

INSTANTIATE_TEST_SUITE_P(, AnalysisFitModeTests, testing::Values(uint64_t{1}, uint64_t{2}, uint64_t{3}, uint64_t{99}),
                           ([](const testing::TestParamInfo<uint64_t>& info) {
	                           switch (info.param) {
	                           case 1: return "Linear";
	                           case 2: return "Power";
	                           case 3: return "Exponential";
	                           default: return "UnknownModeUsesLinear";
	                           }
	                           }));

/// <summary>Une entrée vide écrit uniquement les sentinelles des trois sorties activées.</summary>
TEST(TrackAnalysisAdapterTests, WritesEmptySentinels)
{
	std::array msd{Guard, Guard, Guard}, instant{Guard, Guard, Guard}, fit{Guard, Guard, Guard};
	constexpr std::array Params{4.0};
	ASSERT_TRUE(AnalyzeTracks(nullptr, msd.data() + 1, instant.data() + 1, fit.data() + 1, 0,
		true, true, true, 1.0, 1.0, 3, Params.data()));
	for (const auto& Output : {msd, instant, fit}) {
		EXPECT_DOUBLE_EQ(Output.front(), Guard);
		EXPECT_DOUBLE_EQ(Output[1], -1.0);
		EXPECT_DOUBLE_EQ(Output.back(), Guard);
	}
}
} // namespace palm_tracer::cpu
