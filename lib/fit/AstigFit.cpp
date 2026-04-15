#include "AstigFit.hpp"
#include <limits>
#include <ranges>
#include <tuple> // Pour std::tie

namespace palm_tracer::fit {
//==================== Fonctions Principales ====================
//-------------------------------------------------------------------------------------------------
astig_model AstigFit::Calibrate(const std::vector<Vec3D>& input, const bool center) const
{
	const size_t N = input.size();
	std::vector<double> sx(N), sy(N), z(N);
	// Remplissage des vecteurs
	for (size_t i = 0; i < N; ++i) {
		sx[i] = input[i].x;
		sy[i] = input[i].y;
		z[i]  = input[i].z;
	}

	astig_model res;
	res[0] = fitAxis(sx, z); // Axe X
	res[1] = fitAxis(sy, z); // Axe Y

	if (center) {
		const double Zc = findAstigmaticCenter(res, z);
		res[0][0] -= Zc, res[1][0] -= Zc;
	}

	return res;
}

//-------------------------------------------------------------------------------------------------
Vec2D AstigFit::Estimate(const astig_model& model, const Vec2D& sigma, const std::array<double, 2>& zLimit, const size_t zStep) const
{
	const double DeltaZ = zLimit[1] - zLimit[0];
	double bestZ        = zLimit[0], bestE = std::numeric_limits<double>::infinity();

	// Première passe rapide pour trouver une bonne estimation initiale
	for (size_t i = 0; i < zStep; ++i) {
		const double T = static_cast<double>(i) / static_cast<double>(zStep - 1);
		const double Z = zLimit[0] + T * DeltaZ;
		const double E = getError(model, Z, sigma);
		if (E < bestE) {
			bestE = E;
			bestZ = Z;
		}
	}

	// Utilisation de la méthode de Brent pour affiner l'estimation
	const double Dz = DeltaZ / static_cast<double>(zStep - 1);
	const double A  = std::max(zLimit[0], std::min(bestZ - 2.0 * Dz, zLimit[1]));
	const double B  = std::max(zLimit[0], std::min(bestZ + 2.0 * Dz, zLimit[1]));
	return brent(A, B, model, sigma);
}

//-------------------------------------------------------------------------------------------------
std::vector<Vec2D> AstigFit::Estimate(const astig_model& model, const std::vector<Vec2D>& sigmas, const std::array<double, 2>& zLimit,
                                      const size_t zStep) const
{
	const int N = static_cast<int>(sigmas.size());
	std::vector<Vec2D> res(N);

	#pragma omp parallel for
	for (int i = 0; i < N; ++i) { res[i] = Estimate(model, sigmas[i], zLimit, zStep); }
	return res;
}

//==================== Modèle / coût ====================
//-------------------------------------------------------------------------------------------------
double AstigFit::getSigma(const std::array<double, N_PARAMS>& coeffs, const double z) const
{
	const double Z0 = coeffs[0], W = coeffs[1], C3 = coeffs[2], C4 = coeffs[3], A = coeffs[4]; // Extraction des coefficients

	const double U  = (z - Z0) / W;
	const double U2 = U * U;
	const double U3 = U2 * U;
	const double U4 = U3 * U;

	const double Poly  = 1.0 + U2 + C3 * U3 + C4 * U4;       // Polynôme en U.
	const double Sigma = A * std::sqrt(std::max(0.0, Poly)); // Sigma en nm (ajout d'une sécurité pour éviter les racines négatives)
	return sampling_ / pixelSize_ * Sigma;                   // Sigma en pixel.
}

//-------------------------------------------------------------------------------------------------
double AstigFit::getError(const astig_model& model, const double z, const Vec2D& sigma) const
{
	const Vec2D S     = {.x = getSigma(model[0], z), .y = getSigma(model[1], z)};
	const Vec2D Delta = S - sigma;
	return Delta.x * Delta.x + Delta.y * Delta.y;
}

//==================== Calibration LM ====================
//-------------------------------------------------------------------------------------------------
std::array<double, 5> AstigFit::fitAxis(const std::vector<double>& sigma, const std::vector<double>& z) const
{
	if (sigma.size() < 5 || sigma.size() != z.size()) { // Vérifie si les données sont suffisantes
		std::array<double, N_PARAMS> res{};
		res.fill(-1.0);
		return res; // Renvoie des -1 en cas d'erreur.
	}

	// ---------- Initialisation heuristique ----------
	const auto ItMin  = std::ranges::min_element(sigma);
	const auto IdxMin = static_cast<size_t>(std::distance(sigma.begin(), ItMin));
	const double SMin = std::max(TOL_STEP, sigma[IdxMin]);

	double zMin = z[0], zMax = z[0];
	for (const double Z : z) {
		zMin = std::min(zMin, Z);
		zMax = std::max(zMax, Z);
	}
	const double DeltaZ = std::max(1.0, zMax - zMin);

	// p = [z0, w, c3, c4, A]
	std::array p
	{
		z[IdxMin],                      // z0 ~ plan où sigma est min
		std::max(100.0, DeltaZ / 4.0),  // w  ~ range/4
		0.0, 0.0,                       // C3, C4
		SMin * (pixelSize_ / sampling_) // A (nm) ~ SMin converti en nm : sigma_px = (sampling/pixelSize)*A => A ~ sigma_px * pixelSize / sampling
	};

	// COntraintes pour éviter les divisions par zéro
	p[1] = std::max(p[1], TOL_STEP); // W
	p[4] = std::max(p[4], TOL_STEP); // A

	// ---------- LM ----------
	const size_t N = z.size();
	double lambda  = LAMBDA_0, last_cost = std::numeric_limits<double>::infinity();

	std::vector<double> residue;
	std::vector<std::array<double, N_PARAMS>> jacobian;
	std::array<std::array<double, N_PARAMS>, N_PARAMS> alpha{};
	std::array<double, N_PARAMS> beta{}, delta{}, pTry{};

	for (size_t iteration = 0; iteration < N_ITER; ++iteration) { // Main Loop
		// Calcul du Jacobien et du résidu
		const double Cost = lm(sigma, z, p, residue, jacobian);
		if (last_cost < std::numeric_limits<double>::infinity()) {                    // Si ce n'est pas la première iteration.
			const double Rel = std::abs(Cost - last_cost) / std::max(1.0, last_cost); // Changement relatif du coût
			if (Rel < TOL_COST) { break; }                                            // Critère d'arrêt sur le coût
		}
		last_cost = Cost;

		// Construction des matrices alpha et beta
		for (auto& a : alpha) { a.fill(0.0); }
		beta.fill(0.0);
		delta.fill(0.0);

		for (size_t i = 0; i < N; ++i) {
			for (size_t j = 0; j < N_PARAMS; ++j) {
				beta[j] += jacobian[i][j] * residue[i];
				for (size_t k = 0; k < N_PARAMS; ++k) { alpha[j][k] += jacobian[i][j] * jacobian[i][k]; }
			}
		}

		for (size_t i = 0; i < N_PARAMS; ++i) {
			alpha[i][i] *= 1.0 + lambda;
			beta[i] = -beta[i];
		}

		// Échec de la résolution → pas nul, on durcit lambda et on retente
		if (!solver(alpha, beta, delta)) {
			lambda *= 10.0;
			continue;
		}

		// Calcul du candidat (p + δp)
		for (size_t i = 0; i < N_PARAMS; ++i) { pTry[i] = p[i] + delta[i]; }
		// Contraintes minimales
		pTry[1] = std::max(pTry[1], TOL_STEP); // W
		pTry[4] = std::max(pTry[4], TOL_STEP); // A

		const double CostTry = computeSSE(sigma, z, pTry);

		// Pas relatif
		double stepMax = 0.0;
		for (size_t i = 0; i < N_PARAMS; ++i) {
			const double Denom = std::max(1e-12, std::abs(p[i]));
			stepMax            = std::max(stepMax, std::abs(delta[i]) / Denom);
		}

		if (CostTry < Cost) {
			p      = pTry;
			lambda = std::max(1e-12, lambda / 5.0);
			if (stepMax < TOL_STEP) { break; }
		}
		else { lambda *= 10.0; }
	}
	return p;
}

//-------------------------------------------------------------------------------------------------
double AstigFit::lm(const std::vector<double>& sigma, const std::vector<double>& z, const std::array<double, N_PARAMS>& p,
                    std::vector<double>& r, std::vector<std::array<double, N_PARAMS>>& j) const
{
	const size_t N = z.size();
	r.assign(N, 0.0);
	j.assign(N, {0.0, 0.0, 0.0, 0.0, 0.0});
	double cost        = 0.0;
	const double Z0    = p[0], W = std::max(p[1], TOL_STEP), C3 = p[2], C4 = p[3], A = std::max(p[4], TOL_STEP); // Extraction des coefficients
	const double Scale = sampling_ / pixelSize_;

	for (size_t i = 0; i < N; ++i) {
		const double U  = (z[i] - Z0) / W;
		const double U2 = U * U;
		const double U3 = U2 * U;
		const double U4 = U3 * U;

		const double Poly = 1.0 + U2 + C3 * U3 + C4 * U4; // Polynôme en U.
		const double S    = std::max(1e-18, Poly);        // Ajout d'une sécurité pour éviter les racines négatives ou nulles pour le jacobien
		const double G    = std::sqrt(S);                 // Sigma en nm
		const double Pred = Scale * A * G;                // Sigma en pixel.

		const double Res = Pred - sigma[i];
		r[i]             = Res;
		cost += Res * Res;

		// Calcul du Jacobien numérique
		const double DsDu  = 2.0 * U + 3.0 * C3 * U2 + 4.0 * C4 * U3;
		const double Inv2G = 0.5 / G;

		j[i][0] = Scale * A * (Inv2G * DsDu) * (-1.0 / W); // dσ/dz0
		j[i][1] = Scale * A * (Inv2G * DsDu) * (-U / W);   // dσ/dw
		j[i][2] = Scale * A * Inv2G * U3;                  // dσ/dc3
		j[i][3] = Scale * A * Inv2G * U4;                  // dσ/dc4
		j[i][4] = Scale * G;                               // dσ/dA
	}

	return cost;
}

//-------------------------------------------------------------------------------------------------
double AstigFit::computeSSE(const std::vector<double>& sigma, const std::vector<double>& z, const std::array<double, N_PARAMS>& p) const
{
	const size_t N = z.size();
	double res     = 0;
	for (size_t i = 0; i < N; ++i) {
		const double Pred = getSigma(p, z[i]);
		const double R    = Pred - sigma[i];
		res += R * R;
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
bool AstigFit::solver(std::array<std::array<double, N_PARAMS>, N_PARAMS>& a, std::array<double, N_PARAMS>& b, std::array<double, N_PARAMS>& x)
{
	// Élimination Gauss avec pivot partiel.
	for (size_t i = 0; i < N_PARAMS; ++i) {
		// Pivot
		size_t piv     = i;
		double max_abs = std::abs(a[i][i]);
		for (size_t j = i + 1; j < N_PARAMS; ++j) {
			const double V = std::abs(a[j][i]);
			if (V > max_abs) {
				max_abs = V;
				piv     = j;
			}
		}
		if (max_abs <= std::numeric_limits<double>::epsilon()) { return false; }
		if (piv != i) {
			for (size_t j = i; j < N_PARAMS; ++j) { std::swap(a[i][j], a[piv][j]); }
			std::swap(b[i], b[piv]);
		}

		// Normalisation ligne i
		const double D = a[i][i];
		for (size_t j = i; j < N_PARAMS; ++j) { a[i][j] /= D; }
		b[i] /= D;

		// Élimination
		for (size_t j = 0; j < N_PARAMS; ++j) {
			if (j == i) { continue; }
			const double F = a[j][i];
			if (std::abs(F) <= std::numeric_limits<double>::epsilon()) { continue; }
			for (size_t k = i; k < N_PARAMS; ++k) { a[j][k] -= F * a[i][k]; }
			b[j] -= F * b[i];
		}
	}

	for (size_t i = 0; i < N_PARAMS; ++i) { x[i] = b[i]; }
	return true;
}

//-------------------------------------------------------------------------------------------------
Vec2D AstigFit::brent(double a, double b, const astig_model& model, const Vec2D& sigma) const
{
	if (b < a) { std::swap(a, b); } // Normalisation de l'intervalle : on s'assure que a ≤ b.

	// Constantes de la section dorée (Phi étant le nombre d'or)
	// La section dorée permet de réduire l'intervalle [a, b] avec une convergence garantie, sans hypothèse forte sur la fonction.
	constexpr double InvPhi  = 0.6180339887498948482;
	constexpr double InvPhi2 = 1.0 - InvPhi;

	// Points "mémoire" : x, w, v
	// (x, fx) : meilleur point courant (minimum courant)
	// (w, fw) : second meilleur (historique)
	// (v, fv) : troisième meilleur (historique)
	// Ces 3 points permettent de proposer une interpolation parabolique (rapide) quand la fonction est bien régulière localement.
	double x = a + InvPhi2 * (b - a);
	double w = x, v = x;

	// Évaluations aux points mémorisés
	double fx = getError(model, x, sigma);
	double fw = fx, fv = fx;

	// Pas d / e
	// d : déplacement choisi à l'itération courante (x ⇾ x+d)
	// e : "déplacement précédent" (sert à décider si une parabole est autorisée)
	// Dans Brent, e mémorise la longueur du mouvement qui a précédé, et contrôle le déclenchement de l'étape parabolique.
	double d = 0.0, e = 0.0;


	for (size_t iter = 0; iter < BRENT_ITER; ++iter) {
		// Critère d'arrêt basé sur la taille de l'intervalle et la tolérance
		const double M    = 0.5 * (a + b);
		const double Tol1 = BRENT_TOL * std::abs(x) + 1e-12, Tol2 = 2.0 * Tol1;

		// Condition classique de Brent : si x est suffisamment proche du milieu de l'intervalle au regard de sa demi-largeur, on considère que l'intervalle est "assez petit".
		if (std::abs(x - M) <= Tol2 - 0.5 * (b - a)) { return Vec2D{.x = x, .y = fx}; }

		// ---------- Tentative d'interpolation parabolique. ----------
		// Elle peut accélérer énormément si la fonction est bien lisse près du minimum.
		// Mais elle peut être dangereuse (proposer un point hors intervalle, ou faire un pas trop agressif) : d'où les tests de validité.
		bool parabolic = false;

		// On n'autorise la parabole que si le mouvement précédent "e" est significatif.
		// Au début (e=0), Brent démarre en section dorée (robuste) car l'historique n'est pas encore informatif.
		if (std::abs(e) > Tol1) {
			// Formules "Brent" pour calculer un pas parabolique d = p/q.
			// R, p, q sont des quantités intermédiaires liées à la parabole passant par les 3 points x, w, v. Le minimum de cette parabole tombe en x + p/q.
			const double R = (x - w) * (fx - fv);
			double q       = (x - v) * (fx - fw);
			double p       = (x - v) * q - (x - w) * R;
			q              = 2.0 * (q - R);

			// Convention de signe : on veut q >= 0 pour simplifier les tests.
			if (q > 0.0) { p = -p; }
			q = std::abs(q);

			// ETmp mémorise l'ancien e (longueur de pas "avant-avant"). Ensuite e prend la valeur du dernier pas d.
			const double ETmp = e;
			e                 = d;

			// Conditions de validité du pas parabolique :
			// 1) pas trop grand : |p| < |0.5*q*ETmp| ⇒ le pas proposé doit être raisonnable par rapport au pas précédent.
			// 2) pas dans [a, b] : q * (a-x) < p < q * (b-x) ⇒ le point U = x + p/q doit rester dans l'intervalle.
			// Si tout est OK ⇒ on accepte la parabole.
			if (std::abs(p) < std::abs(0.5 * q * ETmp) && p > q * (a - x) && p < q * (b - x)) {
				parabolic      = true;  // On valide la parabole.
				d              = p / q; // Nouveau pas parabolique.
				const double U = x + d; // Point candidat issu de la parabole.
				// Protection : éviter un U trop proche des bords, ce qui pourrait nuire à la réduction effective de l'intervalle.
				if (U - a < Tol2 || b - U < Tol2) { d = x < M ? Tol1 : -Tol1; }
			}
		}

		// ---------- Étape robuste : section dorée si parabole refusée ----------
		// Si la parabole n'est pas fiable, on retombe sur une étape de section dorée, qui réduit forcément l'intervalle.
		if (!parabolic) {
			e = x < M ? b - x : a - x; // Distance du point x vers le bord "du côté du minimum" (vers le côté où l'intervalle est plus large).
			d = InvPhi2 * e;           // Fraction "dorée" de e.
		}

		// On force un déplacement minimum (Tol1) pour éviter une stagnation due à des arrondis (x + d == x).
		const double U  = x + (std::abs(d) >= Tol1 ? d : d > 0.0 ? Tol1 : -Tol1);
		const double Fu = getError(model, U, sigma);

		if (Fu <= fx) { // U devient le nouveau meilleur point
			// Mise à jour de l'intervalle
			if (U < x) { b = x; }
			else { a = x; }

			// Décalage de l'historique : v <- w <- x <- U
			std::tie(v, fv, w, fw, x, fx) = std::tie(w, fw, x, fx, U, Fu);
		}
		else { // U n'est pas meilleur que x, mais il réduit quand même l'intervalle.
			// Mise à jour de l'intervalle
			if (U < x) { a = U; }
			else { b = U; }

			// Mise à jour de l'historique : on remplace w ou v si U est meilleur qu'eux.
			if (Fu <= fw || std::abs(x - w) <= Z_TOL) { std::tie(v, fv, w, fw) = std::tie(w, fw, U, Fu); }
			else if (Fu <= fv || std::abs(v - x) <= Z_TOL || std::abs(v - w) <= Z_TOL) { std::tie(v, fv) = std::tie(U, Fu); }
		}
	}

	return Vec2D{.x = x, .y = fx};
}

//-------------------------------------------------------------------------------------------------
double AstigFit::findAstigmaticCenter(const astig_model& model, const std::vector<double>& z) const
{
	if (z.empty()) { return 0.0; }

	auto [zMinIt, zMaxIt] = std::ranges::minmax_element(z);
	const double zMin     = *zMinIt;
	const double zMax     = *zMaxIt;

	if (zMax <= zMin) { return zMin; }

	constexpr size_t N_SAMPLES = 2048;
	constexpr size_t N_BISECT  = 64;

	auto delta = [this, &model](const double zz) -> double { return getSigma(model[0], zz) - getSigma(model[1], zz); };

	double prevZ = zMin, prevD            = delta(prevZ);
	double bestD = std::abs(prevD), bestZ = prevZ;

	for (size_t i = 1; i < N_SAMPLES; ++i) {
		const double T  = static_cast<double>(i) / static_cast<double>(N_SAMPLES - 1);
		const double zz = zMin + T * (zMax - zMin);
		const double dd = delta(zz);

		const double absD = std::abs(dd);
		if (absD < bestD) { bestD = absD, bestZ = zz; }

		// Si on a un zéro exact ou un changement de signe, on fait une bissection locale pour affiner la solution.
		if ((prevD == 0.0) || (dd == 0.0) || ((prevD < 0.0) != (dd < 0.0))) {
			double leftZ  = prevZ;
			double rightZ = zz;
			double leftD  = prevD;
			double rightD = dd;

			for (size_t k = 0; k < N_BISECT; ++k) {
				const double midZ = 0.5 * (leftZ + rightZ);
				const double midD = delta(midZ);

				if (std::abs(midD) < bestD) {
					bestD = std::abs(midD);
					bestZ = midZ;
				}

				if (midD == 0.0) { return midZ; }

				if ((leftD < 0.0) != (midD < 0.0)) {
					rightZ = midZ;
					rightD = midD;
				}
				else {
					leftZ = midZ;
					leftD = midD;
				}
			}

			return bestZ;
		}

		prevZ = zz;
		prevD = dd;
	}

	return bestZ;
}
} // namespace palm_tracer::fit
