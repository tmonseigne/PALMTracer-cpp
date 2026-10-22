#include "Fit.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <numeric>
#include <span>
#include <vector>


namespace palm_tracer::tracking {
//-------------------------------------------------------------------------------------------------
Fit::Fit(const ComputesParams& p) { SetMode(p); }

//-------------------------------------------------------------------------------------------------
void Fit::SetMode(const ComputesParams& p)
{
	mode_  = p.mode;
	pSize_ = getPSize();

	// Parsing des paramètres de Fit (pour le moment non implémenté à par fitLength)...
	// Récupération du Fit length et allocation des tableaux dépendant de la taille de la fenêtre
	SetFitLength(p);

	// Allocation des tableaux dépendant du nombre de paramètres, l'utilisation de Assign qui remet à 0 les éléments au lieu du resize qui ne met à 0 que les
	// nouveaux.
	p_.assign(pSize_, 0.0);
	alpha_.assign(pSize_, std::vector(pSize_ + 1, 0.0)); // Attention Alpha possede N lignes et N + 1 colonnes
	alphaR_.assign(pSize_, std::vector(pSize_, 0.0));    // Attention Alpha possede N lignes et N + 1 colonnes
	beta_.assign(pSize_, 0.0);
	pTry_.assign(pSize_, 0.0);
}

//-------------------------------------------------------------------------------------------------
void Fit::SetFitLength(const ComputesParams& p)
{
	fitLength_ = p.fitLength;
	x_.assign(fitLength_, 0.0);
	y_.assign(fitLength_, 0.0);
	f_.assign(fitLength_, 0.0);
	df_.assign(fitLength_, std::vector(pSize_, 0.0));
}

//-------------------------------------------------------------------------------------------------
double Fit::Compute(const std::vector<double>& input)
{
	if (!checkInputValidity(input)) { return -1; }
	const size_t ConfiguredLength = fitLength_;

	// Initialisation
	fillXY(input);
	if (fitLength_ < 2) {
		SetFitLength(ComputesParams{.fitLength = ConfiguredLength});
		return -1; // Sécurité supplémentaire en cas de scintillement massif.
	}
	initP();
	double lambda = 0.001, lambdaOld = lambda, lambdaFactor = LAMBDA_FACTOR_HIGH,
	       errOld = computeError(p_);
	if (!std::isfinite(errOld)) {
		if (fitLength_ != ConfiguredLength) { SetFitLength(ComputesParams{.fitLength = ConfiguredLength}); }
		return -1;
	}

	for (size_t iteration = 0, k = 0; iteration < N_ITER && k < N_K; ++iteration) { // Main Loop
		if (errOld <= 0 || !std::isfinite(lambda)) { break; }
		// Étape LM : construit alpha avec lambda, puis résout pour beta (δp)
		computeAlpha(lambda);
		const bool Ok = computeBeta();
		if (!Ok) { // Échec de la résolution → pas nul, on durcit lambda et on retente
			lambda *= lambdaFactor;
			// Vérification des oscillations (non-convergence), on teste un facteur de lambda plus petit.
			if (iteration % 4 == 0) {
				lambdaFactor = (std::abs(lambdaOld - lambda) <= std::numeric_limits<double>::epsilon()) ? LAMBDA_FACTOR_LOW : LAMBDA_FACTOR_HIGH;
				lambdaOld    = lambda;
			}
			continue;
		}

		// Calcul du candidat (p + δp)
		for (size_t i = 0; i < pSize_; ++i) { pTry_[i] = p_[i] + beta_[i]; }

		// Calcul de l'erreur pour l'essai en cours
		const double ErrTry = computeError(pTry_);

		// Un essai rejete ne modifie ni les paramètres ni le coût de l'état accepté.
		const bool Accepted = std::isfinite(ErrTry) && ErrTry <= errOld;
		if (!Accepted) { lambda *= lambdaFactor; }
		else {
			const double Rho = (errOld - ErrTry) / errOld;
			k                = Rho < EPSILON_RHO ? k + 1 : 0;
			lambda /= lambdaFactor;
			p_     = pTry_;  // accepte le nouveau point
			errOld = ErrTry; // met à jour le coût
		}

		// Vérification des oscillations (non-convergence), on teste un facteur de lambda plus petit.
		if (iteration % 4 == 0) {
			lambdaFactor = (std::abs(lambdaOld - lambda) <= std::numeric_limits<double>::epsilon()) ? LAMBDA_FACTOR_LOW : LAMBDA_FACTOR_HIGH;
			lambdaOld    = lambda;
		}

		// Seul un pas accepté peut déclencher la convergence.
		if (!Accepted) { continue; }
		double m = 0.0;
		for (size_t i = 0; i < pSize_; ++i) { m = std::max(m, std::abs(beta_[i]) / std::max(std::abs(p_[i]), 1e-12)); }
		if (iteration >= 2) {                                                           // Si on est déjà à 2 itérations
			if (m < EPSILON_STEP) { break; }                                            // Si le pas est trop faible
			if (errOld / static_cast<double>(fitLength_) < EPSILON_ERR_STOP) { break; } // Si l'erreur est trop faible
		}
	}
	const double Score = computeR2(p_); // Calculer le R2 sur les seuls lags valides avant de restaurer les buffers.
	if (fitLength_ != ConfiguredLength) { SetFitLength(ComputesParams{.fitLength = ConfiguredLength}); }
	return Score;
}

//-------------------------------------------------------------------------------------------------
size_t Fit::getPSize() const
{
	switch (mode_) {
		case fit_mode::linear: return 2;
		case fit_mode::power: return 3;
		case fit_mode::exponential: return 4;
	}
	return 0;
}

//-------------------------------------------------------------------------------------------------
bool Fit::checkInputValidity(const std::vector<double>& input) const
{
	//Vérification de la taille du Fit.
	if (fitLength_ < 2) {
		//std::cerr << "[Fit] Fit length too short (" << fitLength_ << ") need at least a value of 2" << std::endl;
		return false;
	}

	//Vérification entre la taille de l'entrée et le fitlength.
	if (input.size() < fitLength_) {
		//std::cerr << "[Fit] Not enough sample in input (" << N << ") for a Fit length of " << fitLength_ << std::endl;
		return false;
	}

	// Cas particuliers
	switch (mode_) {
		case fit_mode::linear: { break; }
		case fit_mode::power: {
			// Si les deux premières valeurs sont négatives ou égales à 0.
			if (input[0] <= 0 || input[1] <= 0) {
				//std::cerr << "[Fit] The two first input must be strictly positive (" << input[0] << " and " << input[1] << ")." << std::endl;
				return false;
			}
			break;
		}
		case fit_mode::exponential: { break; }
	}
	return true;
}

//-------------------------------------------------------------------------------------------------
void Fit::fillXY(const std::vector<double>& input)
{
	if (input.size() < fitLength_) { return; } // Déjà vérifier en amont
	hasBlink_ = false;                         // (Re)Initialisation de la détection de scintillement.
	std::ranges::fill(y_, 0.0);                // On initialise à 0, certaines cases peuvent ne pas être remplies.
	size_t idx = 0;
	for (size_t i = 0; i < fitLength_; ++i) {
		if (input[i] < 0) { continue; } // Le msd n'a pas été calculé (blinking fort - aucun élément à un lag donné - sur cette fenêtre).
		x_[idx] = static_cast<double>(i + 1);
		y_[idx] = input[i];
		idx++;
	}
	// Redimensionnement des tableaux x_, y_, f_ et df_ pour ne garder que les éléments valides.
	if (idx != fitLength_) {
		hasBlink_  = true;
		fitLength_ = idx;
		x_.resize(fitLength_);
		y_.resize(fitLength_);
		f_.resize(fitLength_);
		df_.resize(fitLength_);
	}
}

//-------------------------------------------------------------------------------------------------
void Fit::initP()
{
	switch (mode_) {
		case fit_mode::linear: { // p = [a, b]
			p_[0] = y_[1] - y_[0];
			p_[1] = y_[0] - p_[0];
			break;
		}
		case fit_mode::power: { // p = [α, a, b]
			p_[2] = 0.0;
			p_[1] = y_[0];
			p_[0] = std::log(y_[1] / y_[0]) / std::numbers::ln2;
			break;
		}
		case fit_mode::exponential: {        // p = [X0, a, Tau, b]
			p_[0]             = 0.0;         // x0
			p_[1]             = y_.back();   // a (plateau)
			p_[3]             = 0.0;         // b (valeur à x0)
			const double Half = p_[1] / 2.0; // Estimation de la moitié.
			if (y_[0] < p_[1]) {
				size_t i = 0;
				for (; i < fitLength_; ++i) { if (y_[i] >= Half) { break; } }
				p_[2] = static_cast<double>(i) / std::numbers::ln2;
			}
			else { p_[2] = 1.0; }
			break;
		}
	}
}

//-------------------------------------------------------------------------------------------------
void Fit::base(const std::vector<double>& p)
{
	switch (mode_) {
		case fit_mode::linear: { // p = [a, b],  f(x) = a * x + b
			const double A = p[0], B = p[1];
			for (size_t i = 0; i < fitLength_; ++i) { f_[i] = (A * x_[i]) + B; }
			break;
		}
		case fit_mode::power: { // p = [α, a, b],  f(x) = a * x^α + b
			const double P = p[0], A = p[1], B = p[2];
			for (size_t i = 0; i < fitLength_; ++i) { f_[i] = (A * std::pow(x_[i], P)) + B; }
			break;
		}
		case fit_mode::exponential: { // p = [X0, a, Tau, b],  f(x) = (b - a) * exp(-(x - x_0)/(tau)) + a
			const double X0 = p[0], A = p[1], Tau = p[2], Factor = p[3] - A;
			for (size_t i = 0; i < fitLength_; ++i) { f_[i] = (Factor * std::exp(-(x_[i] - X0) / Tau)) + A; }
			break;
		}
	}
}

//-------------------------------------------------------------------------------------------------
void Fit::deriv()
{
	switch (mode_) {
		case fit_mode::linear: { // p = [a, b]
			for (size_t i = 0; i < fitLength_; ++i) {
				df_[i][0] = x_[i]; // ∂f/∂a
				df_[i][1] = 1;     // ∂f/∂b
			}
			break;
		}
		case fit_mode::power: { // p = [α, a, b]
			const double Alpha = p_[0], A = p_[1];
			for (size_t i = 0; i < fitLength_; ++i) {
				const double XAlpha = std::pow(x_[i], Alpha); // x^α
				// ∂/∂p : a * ln(x) * x^α ; cas x≈0 → 0
				df_[i][0] = std::abs(x_[i]) <= std::numeric_limits<double>::epsilon() ? 0 : A * std::log(x_[i]) * XAlpha;
				df_[i][1] = XAlpha; // ∂/∂a
				df_[i][2] = 1;      // ∂/∂b
			}
			break;
		}
		case fit_mode::exponential: { // p = [X0, a, Tau, b]
			const double X0 = p_[0], C = p_[1], Tau = p_[2], Factor = p_[3] - C, D0 = Factor / Tau, D2 = Factor / (Tau * Tau);
			for (size_t i = 0; i < fitLength_; ++i) {
				const double Exp = std::exp(-(x_[i] - X0) / Tau); // E = exp(-(x - X0)/Tau)  (si Tau==0, expression singulière : pas de garde)
				df_[i][0]        = D0 * Exp;                      // ∂/∂X0
				df_[i][1]        = 1 - Exp;                       // ∂/∂a
				df_[i][2]        = D2 * (x_[i] - X0) * Exp;       // ∂/∂Tau
				df_[i][3]        = Exp;                           // ∂/∂b
			}
			break;
		}
	}
}

//-------------------------------------------------------------------------------------------------
void Fit::computeAlpha(const double lambda)
{
	// Initialisation
	for (auto& row : alpha_) { std::ranges::fill(row, 0.0); }
	base(p_);
	deriv();

	// Calcul de l'alpha (Accumulation partie supérieure et RHS)
	for (size_t i = 0; i < fitLength_; ++i) {
		const double Ri = y_[i] - f_[i]; // Résidu réel VS Modèle sur p_
		for (size_t j = 0; j < pSize_; ++j) {
			alpha_[j][pSize_] += Ri * df_[i][j];                                           // RHS (J^T r)
			for (size_t k = j; k < pSize_; ++k) { alpha_[j][k] += df_[i][j] * df_[i][k]; } // partie supérieure (J^T J)
		}
	}
	// Symétrisation (partie inférieure)
	for (size_t i = 1; i < pSize_; ++i) { for (size_t j = 0; j < i; ++j) { alpha_[i][j] = alpha_[j][i]; } }
	// Levenberg–Marquardt scaling diag par 1 + λ (au lieu de λ dans le cas classique)
	for (size_t i = 0; i < pSize_; ++i) { alpha_[i][i] *= 1 + lambda; }
}

//-------------------------------------------------------------------------------------------------
bool Fit::computeBeta()
{
	// Initialisation
	std::ranges::fill(beta_, 0.0);
	constexpr double Eps = std::numeric_limits<double>::epsilon();

	// Triangularisation
	for (size_t i = 0; i < pSize_; ++i) {
		pivot(i); // Pivot partiel échelonné sur la tranche [i... pSize_]
		const double Aii = alpha_[i][i];
		if (std::abs(Aii) <= Eps) { return false; }
		for (size_t j = i + 1; j < pSize_; ++j) {
			const double M = alpha_[j][i] / Aii;
			alphaR_[j][i]  = M;
			alpha_[j][i]   = 0; // On force le zéro sous-diagonal.
			// Mise à jour des colonnes i + 1... pSize_ (inclut RHS à pSize_)
			for (size_t k = i + 1; k <= pSize_; ++k) { alpha_[j][k] -= M * alpha_[i][k]; }
		}
	}

	// Remontée
	for (int i = static_cast<int>(pSize_ - 1); i >= 0; --i) {
		const double Aii = alpha_[i][i];
		if (std::abs(Aii) <= Eps) { return false; }
		double s = alpha_[i][pSize_];
		for (size_t j = i + 1; j < pSize_; ++j) { s -= alpha_[i][j] * beta_[j]; }
		beta_[i] = s / Aii;
	}
	return true;
}

//-------------------------------------------------------------------------------------------------
void Fit::pivot(const size_t index)
{
	// Recherche du maximum échelonné par ligne (m)
	std::vector<double> maxRow(pSize_ - index);
	for (size_t i = index, k = 0; i < pSize_; ++i, ++k) {
		// max sur j = index..pSize_-1 (coefficients, sans RHS)
		const double Vi = std::abs(alpha_[i][index]);
		maxRow[k]       = Vi;
		for (size_t j = index + 1; j < pSize_; ++j) { maxRow[k] = std::max(std::abs(alpha_[i][j]), maxRow[k]); }
		// Normalisation : m = |a(i,index)| / max_row  si max_row ≠ 0
		if (std::abs(maxRow[k]) > std::numeric_limits<double>::epsilon()) { maxRow[k] = Vi / maxRow[k]; }
	}
	double max = maxRow[0];

	// Pivot
	for (size_t i = index + 1, k = 1; i < pSize_; ++i, ++k) {
		if (max >= maxRow[k]) { continue; }
		// Échanger toutes les colonnes de index .. pSize_ (incluant RHS à la colonne pSize_)
		const auto A = std::span(alpha_[index]).subspan(index);
		const auto B = std::span(alpha_[i]).subspan(index);
		std::ranges::swap_ranges(A, B);
		max = maxRow[k];
	}
}

//-------------------------------------------------------------------------------------------------
double Fit::computeError(const std::vector<double>& p)
{
	base(p);
	double error = 0;
	for (size_t i = 0; i < fitLength_; ++i) {
		const double Residual = y_[i] - f_[i];
		error += Residual * Residual;
	}
	return error;
}

//-------------------------------------------------------------------------------------------------
double Fit::computeR2(const std::vector<double>& p)
{
	base(p); // Calcule la prédiction f_ à partir des paramètres
	if (fitLength_ == 0) { return std::numeric_limits<double>::quiet_NaN(); }

	const double MuY = std::accumulate(y_.begin(), y_.end(), 0.0) / static_cast<double>(fitLength_),
	             MuF = std::accumulate(f_.begin(), f_.end(), 0.0) / static_cast<double>(fitLength_);

	// Même score que le fit gaussien : corrélation centrée au carré, pas 1-SSE/SST.
	double covariance = 0, varianceY = 0, varianceF = 0;
	for (size_t i = 0; i < fitLength_; ++i) {
		const double Dy = y_[i] - MuY, Df = f_[i] - MuF;
		covariance += Dy * Df;
		varianceY += Dy * Dy;
		varianceF += Df * Df;
	}
	if (varianceY <= 0 || varianceF <= 0) { return 0; }
	const double R = covariance / std::sqrt(varianceY) / std::sqrt(varianceF);
	return std::isfinite(R) ? std::clamp(R * R, 0.0, 1.0) : 0;
}
} // namespace palm_tracer::tracking
