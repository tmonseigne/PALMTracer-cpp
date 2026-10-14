#include "TestUtils.hpp"
#include "tracking/TracksFit.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace palm_tracer::tracking {
namespace {
/// <summary>Accès de test aux calculs proteges, sans changer l'interface du moteur.</summary>
class TracksFitAccess final : public TracksFit {
public:
	using TracksFit::TracksFit;

	/// <summary>Observe la longueur configurée et les tailles de stockage après un fit.</summary>
	std::array<size_t, 5> StorageSizes() const { return {fitLength_, x_.size(), y_.size(), f_.size(), df_.size()}; }

	/// <summary>Évalue le modèle sur des abscisses connues.</summary>
	std::vector<double> Evaluate(const std::vector<double>& parameters, const std::vector<double>& x)
	{
		x_         = x;
		fitLength_ = x.size();
		f_.resize(fitLength_);
		base(parameters);
		return f_;
	}

	/// <summary>Évalue le Jacobien analytique sur les mêmes abscisses.</summary>
	std::vector<std::vector<double>> Gradient(const std::vector<double>& parameters, const std::vector<double>& x)
	{
		p_         = parameters;
		x_         = x;
		fitLength_ = x.size();
		df_.assign(fitLength_, std::vector<double>(parameters.size()));
		deriv();
		return df_;
	}

	/// <summary>Calcule la corrélation au carré sur des observations explicites.</summary>
	double Score(const std::vector<double>& parameters, const std::vector<double>& x, const std::vector<double>& y)
	{
		x_         = x;
		y_         = y;
		fitLength_ = x.size();
		f_.resize(fitLength_);
		return computeR2(parameters);
	}

	/// <summary>Calcule le coût d'optimisation sur des observations explicites.</summary>
	double Error(const std::vector<double>& parameters, const std::vector<double>& x, const std::vector<double>& y)
	{
		x_         = x;
		y_         = y;
		fitLength_ = x.size();
		f_.resize(fitLength_);
		return computeError(parameters);
	}

	/// <summary>Expose les lags et valeurs conservés après retrait des trous temporels.</summary>
	std::array<std::vector<double>, 2> Prepare(const std::vector<double>& input)
	{
		fillXY(input);
		return {x_, y_};
	}

	/// <summary>Résout un système augmente connu pour vérifier le pivotage indépendamment de l'optimiseur.</summary>
	bool Solve(const std::vector<std::vector<double>>& augmented, std::vector<double>& solution)
	{
		alpha_             = augmented;
		const bool Success = computeBeta();
		solution           = beta_;
		return Success;
	}
};

/// <summary>Valeurs des trois modèles, calculées sans utiliser le moteur de fitting.</summary>
double Model(const fit_mode mode, const std::vector<double>& p, const double x)
{
	switch (mode) {
		case fit_mode::linear: return p[0] * x + p[1];
		case fit_mode::power: return p[1] * std::pow(x, p[0]) + p[2];
		case fit_mode::exponential: return (p[3] - p[1]) * std::exp(-(x - p[0]) / p[2]) + p[1];
	}
	return 0;
}

/// <summary>Paramètres synthétiques dont les valeurs restent positives sur la fenêtre de test.</summary>
std::vector<double> Parameters(const fit_mode mode)
{
	switch (mode) {
		case fit_mode::linear: return {2, 3};
		case fit_mode::power: return {1.3, 2, 0};
		case fit_mode::exponential: return {0, 20, 3, 0};
	}
	return {};
}

/// <summary>Validation indépendante de chacun des modèles d'analyse de trajectoire.</summary>
class TracksFitModelTests : public testing::TestWithParam<fit_mode> { };

/// <summary>Réutilisation après compactage partiel ou scintillement invalidant le fit.</summary>
class TracksFitReuseTests : public testing::TestWithParam<std::vector<double>> { };
}

#pragma region Fit
/// <summary>Retrouve une série synthétique ; l'exponentielle est vérifiée en prédictions car x0 et b ne sont pas identifiables séparément.</summary>
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

/// <summary>Sur données bruitées, le coût obtenu ne dépasse pas celui des paramètres ayant généré la courbe.</summary>
TEST_P(TracksFitModelTests, MinimizesSquaredResiduals)
{
	const ComputesParams Params{.mode = GetParam(), .fitLength = 12};
	TracksFitAccess fit(Params);
	const auto Expected = Parameters(GetParam());
	constexpr std::array Noise{0.2, -0.3, 0.1, 0.4, -0.2, 0.1, -0.4, 0.2, 0.3, -0.1, 0.2, -0.3};
	std::vector<double> x(Noise.size()), input(Noise.size());
	double noiseCost = 0;
	for (size_t i = 0; i < input.size(); ++i) {
		x[i]     = static_cast<double>(i + 1);
		input[i] = Model(GetParam(), Expected, x[i]) + Noise[i];
		noiseCost += Noise[i] * Noise[i];
	}
	const double Score = fit.Fit(input);
	EXPECT_GE(Score, 0);
	EXPECT_LE(Score, 1);
	EXPECT_LE(fit.Error(fit.GetP(), x, input), noiseCost + 1e-6);
}

/// <summary>Le fit linéaire retrouve la solution des moindres carrés, même si les premiers lags suggerent une autre droite.</summary>
TEST(TracksFitTests, LinearLeastSquares)
{
	TracksFit fit(ComputesParams{.mode = fit_mode::linear, .fitLength = 6});
	const std::vector<double> Input{5, 7, 9, 11, 50, 70};
	// Solution analytique : moyenne X = 3.5, moyenne Y = 152/6, covariance = 228, variance X = 17.5.
	constexpr double Slope = 228.0 / 17.5, Intercept = 152.0 / 6.0 - Slope * 3.5;
	EXPECT_GE(fit.Fit(Input), 0);
	const auto Actual = fit.GetP();
	ASSERT_EQ(Actual.size(), 2);
	EXPECT_NEAR(Actual[0], Slope, 1e-5);
	EXPECT_NEAR(Actual[1], Intercept, 1e-5);
}

/// <summary>Les changements de modèle et de longueur sur un objet réutilisé donnent le résultat d'un objet neuf.</summary>
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

/// <summary>Une fenêtre avec des lags absents ne modifie pas la configuration utilisée par la trajectoire suivante.</summary>
TEST_P(TracksFitReuseTests, RestoresWindowAfterMissingLags)
{
	const ComputesParams Params{.mode = fit_mode::linear, .fitLength = 6};
	TracksFitAccess reused(Params);
	const double FirstScore = reused.Fit(GetParam());
	if (std::ranges::count_if(GetParam(), [](const double value) { return value >= 0; }) < 2) { EXPECT_DOUBLE_EQ(FirstScore, -1); }
	EXPECT_TRUE(reused.HasBlink());
	constexpr std::array<size_t, 5> ExpectedSizes{6, 6, 6, 6, 6};
	EXPECT_EQ(reused.StorageSizes(), ExpectedSizes);
	// Les derniers lags ne suivent pas la droite des premiers : une fenêtre raccourcie donnerait un autre résultat.
	const std::vector<double> Full{5, 7, 9, 11, 50, 70};
	TracksFit fresh(Params);
	EXPECT_DOUBLE_EQ(reused.Fit(Full), fresh.Fit(Full));
	EXPECT_EQ(reused.GetP(), fresh.GetP());
	EXPECT_FALSE(reused.HasBlink());
	EXPECT_EQ(reused.StorageSizes(), ExpectedSizes);
}

INSTANTIATE_TEST_SUITE_P(MissingLags, TracksFitReuseTests,
                         testing::Values(std::vector<double>{5, -1, 9, 11, -1, 15}, std::vector<double>{-1, -1, 9, -1, -1, -1},
	                         std::vector<double>{-1, -1, -1, -1, -1, -1}),
                         ([](const testing::TestParamInfo<std::vector<double>>& info) {
	                         constexpr std::array Names{"PartialLags", "SingleLag", "NoLags"};
	                         return Names[info.index];
	                         }));
#pragma endregion Fit

#pragma region fillXY
/// <summary>Les lags d'origine sont conservés après retrait des valeurs manquantes et le scintillement est réinitialisé.</summary>
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
/// <summary>Compare chaque modèle aux équations analytiques sur des abscisses non entières.</summary>
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
/// <summary>Compare toutes les dérivées analytiques aux différences finies centrées.</summary>
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

#pragma region computeBeta
/// <summary>Un pivot initial nul doit être remplace par une ligne valide, et non provoquer un échec du solveur.</summary>
TEST(TracksFitTests, SolvesZeroLeadingPivot)
{
	TracksFitAccess fit(ComputesParams{.mode = fit_mode::linear, .fitLength = 2});
	std::vector<double> solution;
	ASSERT_TRUE(fit.Solve({{0, 1, 2}, {1, 1, 5}}, solution));
	ASSERT_EQ(solution.size(), 2);
	EXPECT_DOUBLE_EQ(solution[0], 3);
	EXPECT_DOUBLE_EQ(solution[1], 2);
}
#pragma endregion computeBeta

#pragma region computeError
/// <summary>La SSE distingue des modèles dont la corrélation au carré est pourtant identique.</summary>
TEST(TracksFitTests, SquaredResidualCost)
{
	TracksFitAccess fit(ComputesParams{.mode = fit_mode::linear, .fitLength = 4});
	const std::vector<double> X{1, 2, 3, 4}, Y{5, 7, 9, 11};
	EXPECT_DOUBLE_EQ(fit.Error({2, 3}, X, Y), 0);
	EXPECT_DOUBLE_EQ(fit.Error({2, 4}, X, Y), 4);
	EXPECT_DOUBLE_EQ(fit.Error({0, 8}, X, Y), 20);
	EXPECT_DOUBLE_EQ(fit.Error({0, 30}, X, Y), 1956);
}
#pragma endregion computeError

#pragma region computeR2
/// <summary>Le score final partage la convention gaussienne, invariante aux transformations affines.</summary>
TEST(TracksFitTests, CorrelationScore)
{
	TracksFitAccess fit(ComputesParams{.mode = fit_mode::linear, .fitLength = 4});
	const std::vector<double> X{1, 2, 3, 4}, Y{5, 7, 9, 11};
	EXPECT_DOUBLE_EQ(fit.Score({2, 3}, X, Y), 1);
	EXPECT_DOUBLE_EQ(fit.Score({2, 4}, X, Y), 1);
	EXPECT_DOUBLE_EQ(fit.Score({-2, 4}, X, Y), 1);
	EXPECT_DOUBLE_EQ(fit.Score({0, 8}, X, Y), 0);
	EXPECT_DOUBLE_EQ(fit.Score({0, 30}, X, Y), 0);
	EXPECT_DOUBLE_EQ(fit.Score({2, 3}, X, {8, 8, 8, 8}), 0);
}

/// <summary>La corrélation d'une série non affine reste inférieure à 1 et correspond au calcul analytique centre.</summary>
TEST(TracksFitTests, ImperfectCorrelationScore)
{
	TracksFitAccess fit(ComputesParams{.mode = fit_mode::linear, .fitLength = 4});
	EXPECT_NEAR(fit.Score({1, 0}, {1, 2, 3, 4}, {1, 2, 4, 8}), 0.92, 1e-14);
}
#pragma endregion computeR2
} // namespace palm_tracer::tracking
