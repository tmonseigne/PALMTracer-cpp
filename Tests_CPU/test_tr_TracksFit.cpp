#include "TestUtils.hpp"
#include "tracking/TracksFit.hpp"

#include <array>
#include <cmath>

namespace palm_tracer::tracking {
namespace {
/// <summary>Acces de test aux calculs proteges, sans changer l'interface du moteur.</summary>
class TracksFitAccess final : public TracksFit {
public:
	using TracksFit::TracksFit;

	/// <summary>Evalue le modele sur des abscisses connues.</summary>
	std::vector<double> Evaluate(const std::vector<double>& parameters, const std::vector<double>& x)
	{
		x_         = x;
		fitLength_ = x.size();
		f_.resize(fitLength_);
		base(parameters);
		return f_;
	}

	/// <summary>Evalue le Jacobien analytique sur les memes abscisses.</summary>
	std::vector<std::vector<double>> Gradient(const std::vector<double>& parameters, const std::vector<double>& x)
	{
		p_         = parameters;
		x_         = x;
		fitLength_ = x.size();
		df_.assign(fitLength_, std::vector<double>(parameters.size()));
		deriv();
		return df_;
	}

	/// <summary>Calcule la metrique de determination sur des observations explicites.</summary>
	double Score(const std::vector<double>& parameters, const std::vector<double>& x, const std::vector<double>& y)
	{
		x_         = x;
		y_         = y;
		fitLength_ = x.size();
		f_.resize(fitLength_);
		return computeR2(parameters);
	}

	/// <summary>Expose les lags et valeurs conserves apres retrait des trous temporels.</summary>
	std::array<std::vector<double>, 2> Prepare(const std::vector<double>& input)
	{
		fillXY(input);
		return {x_, y_};
	}
};

/// <summary>Valeurs des trois modeles, calculees sans utiliser le moteur de fitting.</summary>
double Model(const fit_mode mode, const std::vector<double>& p, const double x)
{
	switch (mode) {
		case fit_mode::linear: return p[0] * x + p[1];
		case fit_mode::power: return p[1] * std::pow(x, p[0]) + p[2];
		case fit_mode::exponential: return (p[3] - p[1]) * std::exp(-(x - p[0]) / p[2]) + p[1];
	}
	return 0;
}

/// <summary>Parametres synthetiques dont les valeurs restent positives sur la fenetre de test.</summary>
std::vector<double> Parameters(const fit_mode mode)
{
	switch (mode) {
		case fit_mode::linear: return {2, 3};
		case fit_mode::power: return {1.3, 2, 0};
		case fit_mode::exponential: return {0, 20, 3, 0};
	}
	return {};
}

/// <summary>Validation independante de chacun des modeles d'analyse de trajectoire.</summary>
class TracksFitModelTests : public testing::TestWithParam<fit_mode> { };
}

#pragma region Fit
/// <summary>Retrouve une serie synthetique ; l'exponentielle est verifiee en predictions car x0 et b ne sont pas identifiables separement.</summary>
TEST_P(TracksFitModelTests, RecoversSyntheticCurve)
{
	const ComputesParams Params{.mode = GetParam(), .fitLength = 12};
	TracksFit fit(Params);
	const auto Expected = Parameters(GetParam());
	std::vector<double> input(Params.fitLength);
	for (size_t i = 0; i < input.size(); ++i) { input[i] = Model(GetParam(), Expected, static_cast<double>(i + 1)); }
	EXPECT_NEAR(fit.Fit(input), 1, 1e-5);
	const auto Actual = fit.GetP();
	ASSERT_EQ(Actual.size(), Expected.size());
	for (size_t i = 0; i < input.size(); ++i) { EXPECT_NEAR(Model(GetParam(), Actual, static_cast<double>(i + 1)), input[i], 0.01) << "Lag " << i + 1; }
	if (GetParam() != fit_mode::exponential) { for (size_t i = 0; i < Actual.size(); ++i) { EXPECT_NEAR(Actual[i], Expected[i], 0.01); } }
}

/// <summary>Les changements de modele et de longueur sur un objet reutilise donnent le resultat d'un objet neuf.</summary>
TEST(TracksFitTests, ReconfigureAndReuse)
{
	TracksFit reused(ComputesParams{.mode = fit_mode::linear, .fitLength = 6});
	constexpr std::array Modes{fit_mode::linear, fit_mode::power, fit_mode::exponential, fit_mode::linear};
	for (size_t k = 0; k < Modes.size(); ++k) {
		SCOPED_TRACE("Configuration " + std::to_string(k));
		const ComputesParams Params{.mode = Modes[k], .fitLength = 6 + 2 * k};
		reused.SetMode(Params);
		TracksFit fresh(Params);
		std::vector<double> input(Params.fitLength);
		for (size_t i = 0; i < input.size(); ++i) { input[i] = Model(Params.mode, Parameters(Params.mode), static_cast<double>(i + 1)); }
		EXPECT_DOUBLE_EQ(reused.Fit(input), fresh.Fit(input));
		EXPECT_EQ(reused.GetP(), fresh.GetP());
		EXPECT_EQ(reused.HasBlink(), fresh.HasBlink());
	}
}
#pragma endregion Fit

#pragma region fillXY
/// <summary>Les lags d'origine sont conserves apres retrait des valeurs manquantes et le scintillement est reinitialise.</summary>
TEST(TracksFitTests, PreservesMissingLagPositions)
{
	const ComputesParams Params{.mode = fit_mode::linear, .fitLength = 6};
	TracksFitAccess fit(Params);
	const auto Data = fit.Prepare({5, -1, 9, 11, -1, 15});
	EXPECT_EQ(Data[0], (std::vector<double>{1, 3, 4, 6}));
	EXPECT_EQ(Data[1], (std::vector<double>{5, 9, 11, 15}));
	EXPECT_TRUE(fit.HasBlink());
	fit.SetFitLength(Params);
	const auto Full = fit.Prepare({5, 7, 9, 11, 13, 15});
	EXPECT_EQ(Full[0], (std::vector<double>{1, 2, 3, 4, 5, 6}));
	EXPECT_FALSE(fit.HasBlink());
}
#pragma endregion fillXY

#pragma region base
/// <summary>Compare chaque modele aux equations analytiques sur des abscisses non entieres.</summary>
TEST_P(TracksFitModelTests, ModelValues)
{
	TracksFitAccess fit(ComputesParams{.mode = GetParam(), .fitLength = 4});
	const auto Params = Parameters(GetParam());
	const std::vector<double> X{0.5, 1.25, 2.5, 4};
	const auto Values = fit.Evaluate(Params, X);
	ASSERT_EQ(Values.size(), X.size());
	for (size_t i = 0; i < X.size(); ++i) { EXPECT_NEAR(Values[i], Model(GetParam(), Params, X[i]), 1e-12); }
}
#pragma endregion base

#pragma region deriv
/// <summary>Compare toutes les derivees analytiques aux differences finies centrees.</summary>
TEST_P(TracksFitModelTests, GradientMatchesFiniteDifferences)
{
	TracksFitAccess fit(ComputesParams{.mode = GetParam(), .fitLength = 4});
	const auto Params = Parameters(GetParam());
	const std::vector<double> X{0.5, 1.25, 2.5, 4};
	const auto Gradient   = fit.Gradient(Params, X);
	constexpr double Step = 1e-5;
	for (size_t j = 0; j < Params.size(); ++j) {
		auto plus = Params, minus = Params;
		plus[j] += Step;
		minus[j] -= Step;
		for (size_t i = 0; i < X.size(); ++i) {
			const double Expected = (Model(GetParam(), plus, X[i]) - Model(GetParam(), minus, X[i])) / (2 * Step);
			EXPECT_NEAR(Gradient[i][j], Expected, 1e-7) << "Lag " << i << ", parametre " << j;
		}
	}
}

INSTANTIATE_TEST_SUITE_P(Models, TracksFitModelTests,
                         testing::Values(fit_mode::linear, fit_mode::power, fit_mode::exponential),
                         ([](const testing::TestParamInfo<fit_mode>& info) {
	                         constexpr std::array Names{"Linear", "Power", "Exponential"};
	                         return Names[static_cast<size_t>(info.param)];
	                         }));
#pragma endregion deriv

#pragma region computeR2
/// <summary>Le R2 de trajectoire est 1-SSE/SST, distinct de la correlation au carre du fitting gaussien.</summary>
TEST(TracksFitTests, DeterminationScore)
{
	TracksFitAccess fit(ComputesParams{.mode = fit_mode::linear, .fitLength = 4});
	const std::vector<double> X{1, 2, 3, 4}, Y{5, 7, 9, 11};
	EXPECT_DOUBLE_EQ(fit.Score({2, 3}, X, Y), 1);
	EXPECT_DOUBLE_EQ(fit.Score({2, 4}, X, Y), 0.8);
	EXPECT_DOUBLE_EQ(fit.Score({0, 8}, X, Y), 0);
	EXPECT_LT(fit.Score({0, 30}, X, Y), 0);
}
#pragma endregion computeR2
} // namespace palm_tracer::tracking
