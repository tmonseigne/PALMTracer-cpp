// TODO(S) A VOIR !!!!
#include "TracksFit.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <numbers>
#include <numeric>
#include <span>
#include <vector>


namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
TracksFit::TracksFit(const fit_mode mode, const std::vector<double>& fitParams) { SetMode(mode, fitParams); }

//-------------------------------------------------------------------------------------------------
void TracksFit::SetMode(const fit_mode mode, const std::vector<double>& fitParams)
{
	mode_  = mode;
	pSize_ = getPSize();

	// Parsing des paramètres de Fit (pour le moment non implémenté à par fitLength)...
	// Récupération du Fit length et allocation des tableaux dépendant de la taille de la fenêtre
	SetFitLength(static_cast<size_t>(fitParams[0]));

	// Allocation des tableaux dépendant du nombre de paramètres
	p_.resize(pSize_);
	alpha_.resize(pSize_, std::vector<double>(pSize_ + 1)); // Attention Alpha possede N lignes et N + 1 colonnes
	alphaR_.resize(pSize_, std::vector<double>(pSize_));    // Attention Alpha possede N lignes et N + 1 colonnes
	beta_.resize(pSize_);
	pTry_.resize(pSize_);
}

//-------------------------------------------------------------------------------------------------
void TracksFit::SetFitLength(const size_t fitLength)
{
	fitLength_ = fitLength;
	x_.resize(fitLength_);
	y_.resize(fitLength_);
	f_.resize(fitLength_);
	df_.resize(fitLength_, std::vector<double>(pSize_));
}

//-------------------------------------------------------------------------------------------------
double TracksFit::Fit(const std::vector<double>& input)
{
	if (!checkInputValidity(input)) { return -1; }

	// Initialisation
	fillXY(input);
	initP(input);
	double mu     = 0,
	       lambda = 0.001, lambdaOld = lambda, lambdaFactor = LAMBDA_FACTOR_HIGH,
	       errOld = computeR2(p_);

	for (size_t iteration = 0, k = 0; iteration < N_ITER && k < N_K; ++iteration) { // Main Loop
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
		const double ErrTry = computeR2(pTry_);

		// Mise à jour de mu et lambda
		mu               = (ErrTry > EPSILON_ERR) ? ErrTry / errOld : 0.0;
		const double Rho = (ErrTry - errOld) / (ErrTry > EPSILON_ERR ? ErrTry : 1);

		if (mu > EPSILON_MU) { lambda *= lambdaFactor; } // Mauvaise iteration, on durcit lambda et on retente
		else {
			if (std::abs(Rho) < EPSILON_RHO) { ++k; } // stagnation détectée
			lambda /= lambdaFactor;
			p_     = pTry_;  // accepte le nouveau point
			errOld = ErrTry; // met à jour le coût
		}

		// Vérification des oscillations (non-convergence), on teste un facteur de lambda plus petit.
		if (iteration % 4 == 0) {
			lambdaFactor = (std::abs(lambdaOld - lambda) <= std::numeric_limits<double>::epsilon()) ? LAMBDA_FACTOR_LOW : LAMBDA_FACTOR_HIGH;
			lambdaOld    = lambda;
		}

		// TODO Ajouts de critères d'arrêts supplémentaires ? Gaussian Fit est bien plus sévère (plus complexe aussi) et là, on a 500 interations par défaut
		// if (...) { break; }
	}
	return computeR2(p_); // R² final
}

//-------------------------------------------------------------------------------------------------
size_t TracksFit::getPSize() const
{
	switch (mode_) {
		case fit_mode::linear: return 2;
		case fit_mode::power: return 3;
		case fit_mode::exponential: return 4;
	}
	return 0;
}

//-------------------------------------------------------------------------------------------------
bool TracksFit::checkInputValidity(const std::vector<double>& input) const
{
	//Vérification de la taille du Fit.
	if (fitLength_ < 2) {
		std::cerr << "[TracksFit] Fit length too short (" << fitLength_ << ") need at least a value of 2" << std::endl;
		return false;
	}

	//Vérification entre la taille de l'entrée et le fitlength.
	const size_t N = input.size();
	if (N < fitLength_) {
		std::cerr << "[TracksFit] Not enough sample in input (" << N << ") for a Fit length of " << fitLength_ << std::endl;
		return false;
	}

	// Cas particuliers
	switch (mode_) {
		case fit_mode::linear: { break; }
		case fit_mode::power: {
			// Si les deux premières valeurs sont négatives ou égales à 0.
			if (!(input[0] > 0 && input[1] > 0)) {
				std::cerr << "[TracksFit] The two first input must be strictly positive (" << input[0] << " and " << input[1] << ")." << std::endl;
				return false;
			}
			break;
		}
		case fit_mode::exponential: { break; }
	}
	return true;
}

//-------------------------------------------------------------------------------------------------
void TracksFit::fillXY(const std::vector<double>& input)
{
	if (input.size() < fitLength_) { return; } // Déjà vérifier en amont
	size_t toRemove = 0;
	std::ranges::fill(y_, 0.0); // On initialise à 0, certaines cases peuvent ne pas être remplies.
	for (size_t i = 0; i < fitLength_; ++i) {
		x_[i] = static_cast<double>(i + 1);
		if (input[i] >= 0) { y_[i - toRemove] = input[i]; }
		else { toRemove++; } // Le msd n'a pas été calculé (blinking fort - aucun élément à un lag donné - sur cette fenêtre).
	}
	nY_ = fitLength_ - toRemove;
	// TODO: on ne redimensionne pas les tableaux, mais on prend en compte les éléments "sauté" pour l'initialisation des paramètres lors du fit.
	// Ce doit être une coquille, mais dans le code VB, on parcourt Y pour calculer le R2 et la matrice Alpha.
	// Mais comme il y a potentiellement des éléments, à la fin, initialisé à 0 et que le tableau conserve une taille de fitlength.
	// Le comportement est-il voulu ? faut-il parcourir plutôt jusqu'à la fin des Y "correct".
	// On utilise dans base et deriv X et sa taille donc `fitLength`, mais c'est juste un compteur de fitLength pour le moment.
	// (la partie histogramme modifie snas doute le comportement à vérifier.)
}

//-------------------------------------------------------------------------------------------------
void TracksFit::initP(const std::vector<double>& input)
{
	switch (mode_) {
		case fit_mode::linear: { // p = [a, b]
			p_[0] = input[1] - input[0];
			p_[1] = input[0] - p_[0];
			break;
		}
		case fit_mode::power: { // p = [p, a, b]
			p_[2] = 0.0;
			p_[1] = input[0];
			p_[0] = std::log(input[1] / input[0]) / std::numbers::ln2;
			break;
		}
		case fit_mode::exponential: { // p = [X0, a, Tau, b]
			// TODO : Attention ici. J'ai un sérieux doute. Je devrais prendre Y selon moi qui possède les inputs moins les -1 et des 0 à la fin car non-redimensionné.
			// Sinon prendre le vrai dernier élément de input
			p_[1]             = input[nY_ - 1]; //p_[1] = y_[nY_ - 1]; ou p_[1] = input.back(); // a (plateau)
			p_[0]             = 0.0;            // x0
			p_[3]             = 0.0;            // b (valeur à x0)
			const double Half = p_[1] / 2.0;    // Estimation de la moitié.
			if (input[0] < p_[1]) {
				//p_[2] = -1.0 / (std::log(1 - input[0] / p_[1])); // inutile, on le réalloue direct après
				size_t i = 0;
				for (; i < fitLength_; ++i) { if (input[i] >= Half) { break; } }
				p_[2] = static_cast<double>(i) / std::numbers::ln2;
			}
			else { p_[2] = 1.0; }
			break;
		}
	}
}

//-------------------------------------------------------------------------------------------------
void TracksFit::base(const std::vector<double>& p)
{
	switch (mode_) {
		case fit_mode::linear: { // p = [a, b],  f(x) = a * x + b
			const double A = p[0], B = p[1];
			for (size_t i = 0; i < fitLength_; ++i) { f_[i] = A * x_[i] + B; }
			break;
		}
		case fit_mode::power: { // p = [p, a, b],  f(x) = a * x^p + b
			const double P = p[0], A = p[1], B = p[2];
			for (size_t i = 0; i < fitLength_; ++i) { f_[i] = A * std::pow(x_[i], P) + B; }
			break;
		}
		case fit_mode::exponential: { // p = [X0, a, Tau, b],  f(x) = (b - a) * exp(-(x - x_0)/(tau)) + a
			const double X0 = p[0], A = p[1], Tau = p[2], Factor = p[3] - A;
			for (size_t i = 0; i < fitLength_; ++i) { f_[i] = Factor * std::exp(-(x_[i] - X0) / Tau) + A; }
			break;
		}
	}
}

//-------------------------------------------------------------------------------------------------
void TracksFit::deriv()
{
	switch (mode_) {
		case fit_mode::linear: { // p = [a, b]
			for (size_t i = 0; i < fitLength_; ++i) {
				df_[i][0] = x_[i]; // ∂f/∂a
				df_[i][1] = 1;     // ∂f/∂b
			}
			break;
		}
		case fit_mode::power: { // p = [p, a, b]
			const double P = p_[0], A = p_[1];
			for (size_t i = 0; i < fitLength_; ++i) {
				const double PowX = std::pow(x_[i], P); // x^p
				// ∂/∂p : a * ln(x) * x^p ; cas x≈0 → 0
				df_[i][0] = std::abs(x_[i]) <= std::numeric_limits<double>::epsilon() ? 0 : A * std::log(x_[i]) * PowX;
				df_[i][1] = PowX; // ∂/∂a
				df_[i][2] = 1;    // ∂/∂b
			}
			break;
		}
		case fit_mode::exponential: { // p = [X0, a, Tau, b]
			const double X0 = p_[0], C = p_[1], Tau = p_[2], Factor = p_[3] - C, D0 = Factor / Tau, D2 = Factor / (Tau * Tau);
			for (size_t i = 0; i < fitLength_; ++i) {
				const double Exp = std::exp(-(x_[i] - X0) / Tau); // E = exp(-(x - X0)/Tau)  (si Tau==0, expr. singulière : pas de garde)
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
void TracksFit::computeAlpha(const double lambda)
{
	// Initialisation
	for (auto& row : alpha_) { std::ranges::fill(row, 0.0); }
	base(p_);
	deriv();

	// Calcul de l'alpha (Accumulation partie supérieure et RHS)
	for (size_t i = 0; i < fitLength_; ++i) { // TODO: longueurs : fitLength_ vs nY_ les 0 de fin de Y_ sont assumé si fitLength_
		const double Ri = y_[i] - f_[i];      // Résidu réel VS Modèle sur p_
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
bool TracksFit::computeBeta()
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
			alpha_[j][i]   = 0; // on force le zéro sous-diagonal
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
void TracksFit::pivot(const size_t index)
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
	if (std::abs(max) <= std::numeric_limits<double>::epsilon()) { return; }

	// Pivot
	for (size_t i = index + 1, k = 1; i < pSize_; ++i, ++k) {
		if (max >= maxRow[k]) { continue; }
		// Échanger toutes les colonnes de index .. pSize_ (incluant RHS à la colonne pSize_)
		auto a = std::span(alpha_[index]).subspan(index);
		auto b = std::span(alpha_[i]).subspan(index);
		std::ranges::swap_ranges(a, b);
		max = maxRow[k];
	}
}

//-------------------------------------------------------------------------------------------------
double TracksFit::computeR2(const std::vector<double>& p)
{
	base(p); // Calcules la prédiction f_ à partir des paramètres
	const size_t N = y_.size();
	if (N == 0) { return std::numeric_limits<double>::quiet_NaN(); }

	// Moyenne de Y.
	const double Mu = std::accumulate(y_.begin(), y_.end(), 0.0) / static_cast<double>(N);

	// Sommes des carrés : résiduels (s1) et totale (s2)
	double s1 = 0.0, s2 = 0.0;
	for (size_t i = 0; i < N; ++i) {
		const double D1 = y_[i] - f_[i],
		             D2 = y_[i] - Mu;
		s1 += D1 * D1;
		s2 += D2 * D2;
	}
	if (s2 == 0.0) { return 0.0; } // Variance nulle : éviter /0 ; comportement choisi : R² = 0.0 (documenté ci-dessus).
	return 1 - (s1 / s2);
}
} // namespace palm_tracer
