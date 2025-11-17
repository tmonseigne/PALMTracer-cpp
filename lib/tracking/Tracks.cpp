#include "Tracks.hpp"
#include "Constants.hpp"

#include <cmath>
#include <omp.h>

#include "TracksFit.hpp"

namespace palm_tracer::tracking {
//==================== DLL Input/output ====================
//-------------------------------------------------------------------------------------------------
void Tracks::LoadTracking(const double* input, const size_t nRow, const double pixelSize)
{
	// Remise à zéro.
	t_.clear();

	// Récupération du nombre de trajectoires pour préallouer le tableau
	const auto NTracks = static_cast<size_t>(input[(nRow - 1) * N_TRACKING]);
	t_.reserve(nRow);

	// Parcours des lignes
	track_t track;           // Trajectoire en cours (réutilisée entre groupes)
	track.reserve(64);       // Heuristique: évite quelques reallocs sur petites tracks
	size_t currentTrack = 1; // input[0]
	for (size_t i = 0, r = 0; i < nRow; ++i, r += N_TRACKING) {
		const auto TrackId = static_cast<size_t>(input[r]);
		// Nouveau groupe de trajectoire détecté
		if (TrackId != currentTrack) {
			if (!track.empty()) { t_.push_back(std::move(track)); } // Move, pas de copie, `track` est maintenant vide, mais conserve sa capacité.
			track.clear();                                          // Redondant après move, mais sûr si branché autrement
			currentTrack = TrackId;
		}

		// Ajout du point (emplace_back = construction in-place)
		const Vec3D PosPx = {.x = input[r + 3], .y = input[r + 4], .z = input[r + 5]};
		const Vec3D Pos   = {.x = PosPx.x * pixelSize, .y = PosPx.y * pixelSize, .z = PosPx.z};
		track.emplace_back(Point(static_cast<size_t>(input[r + 2]), // Id
		                         static_cast<size_t>(input[r + 1]), // Plan
		                         PosPx,                             // Position en pixel
		                         Pos,                               // Position
		                         input[r + 6],                      // Intensité
		                         input[r + 7]));                    // Surface
	}
	if (!track.empty()) { t_.push_back(std::move(track)); } // Ne pas oublier d'ajouter la dernière trajectoire
}

//-------------------------------------------------------------------------------------------------
size_t Tracks::FillOutput(double* output) const
{
	size_t idx = 0;
	for (size_t i = 0; i < t_.size(); ++i) {
		for (const auto& P : t_[i]) {
			output[idx++] = static_cast<double>(i + 1);    // Track
			output[idx++] = static_cast<double>(P.plane_); // Plane
			output[idx++] = static_cast<double>(P.id_);    // Id
			output[idx++] = P.posPx_.x;                    // X
			output[idx++] = P.posPx_.y;                    // Y
			output[idx++] = P.posPx_.z;                    // Z
			output[idx++] = P.intensity_;                  // Intensité
			output[idx++] = P.surface_;                    // Surface
		}
	}
	return idx;
}

//==================== Computes ====================
//-------------------------------------------------------------------------------------------------
void Tracks::BlinkingReconnection(const size_t mode, const size_t maxDuration, const double maxSpeed)
{
	if (t_.empty()) { return; }
	const size_t N = t_.size();
	std::vector taken(N, false);             // True si l'index appartient déjà à une chaîne (hors tête)
	std::vector<std::vector<size_t>> chains; // Liste des chaînes d'indices (ex: { [2,4,5], [3], [6,7], ... })
	chains.reserve(N);                       // Taille maximale possible

	// --- Construction des chaînes d'indices de connexion ---
	for (size_t head = 0; head < N; ++head) { // Trajectoire à tester
		if (t_[head].empty()) { continue; }   // Au cas où, mais théoriquement impossible.
		if (taken[head]) { continue; }        // Trajectoire déjà absorbée par une chaîne précédente
		std::vector<size_t> chain;            // Nouvelle Trajectoire
		chain.push_back(head);                // Ajout de la trajectoire initiale

		// On tente d'allonger la chaîne de façon gloutonne
		size_t current = head;
		while (true) {
			if (t_[current].empty()) { break; }                   // Au cas où, mais théoriquement impossible.
			const Point& Lp = t_[current].back();                 // Dernier point de la trajectoire i.
			size_t bJ       = std::numeric_limits<size_t>::max(), // Indice (j) du meilleur point.
			       bT       = std::numeric_limits<size_t>::max(); // Différence temporelle avec le meilleur point.
			double bD       = std::numeric_limits<double>::max(); // Différence spatiale avec le meilleur point.

			// Cherche le meilleur j non encore pris (hors tête) et non vide
			for (size_t j = head + 1; j < N; ++j) {
				if (j == current || taken[j] || t_[j].empty()) { continue; }
				const Point& Fp = t_[j].front();          // Premier point de la trajectoire j.
				if (Fp.plane_ <= Lp.plane_) { continue; } // Vérification de continuité temporelle.
				const size_t Dp = Fp.plane_ - Lp.plane_;  // Calcul du temps (nombre de plans entre les trajectoires).
				// (à vérifier, mais un break est possible, car les trajectoires sont triées par plan de démarrage)
				if (Dp > maxDuration || Dp > bT) { continue; }                               // Vérification de validité temporelle.
				const double D    = Fp.GetDistance(Lp),                                      // Calcul de la distance.
				             MaxD = getMaxDistance(mode, static_cast<double>(Dp), maxSpeed); // Calcul de la distance maximale.
				if (D >= bD || D > MaxD) { continue; }                                       // Vérification de validité Spatiale.

				bJ = j;  // Mise à jour de l'élément sélectionné.
				bT = Dp; // Mise à jour de la durée (nombre de plan).
				bD = D;  // Mise à jour de la distance.
			}

			if (bJ == std::numeric_limits<size_t>::max()) { break; } // Aucune trajectoire de trouvé, on sort du while
			chain.push_back(bJ);                                     // Ajout de la trajectoire à la liste
			taken[bJ] = true;                                        // Ajout d'un token "trajectoire prise"
			current   = bJ;                                          // Mise à jour de la trajectoire courante (qui est la dernière de la liste)
		}
		chains.emplace_back(std::move(chain)); // Ajout de la chaîne de trajectoire au tableau final.
	}

	// --- Recomposition finale dans un nouveau conteneur, par move-concat ---
	std::vector<track_t> new_tracks;   // Création d'une liste de trajectoires temporaire
	new_tracks.reserve(chains.size()); // On a notre nombre de trajectoires finales
	for (const auto& Ch : chains) {
		// Concatène tous les points des indices de la chaîne dans un seul track
		track_t acc;

		// (option : size hint) somme des tailles pour limiter les reallocs
		size_t totalPts = 0;
		for (const size_t Idx : Ch) { totalPts += t_[Idx].size(); }
		acc.reserve(totalPts);

		for (const size_t Idx : Ch) {
			auto& tr = t_[Idx];
			acc.insert(acc.end(), std::make_move_iterator(tr.begin()), std::make_move_iterator(tr.end()));
			tr.clear(); // libère les références, pas la capacité (sans importance ici)
		}
		new_tracks.emplace_back(std::move(acc));
	}
	t_.swap(new_tracks); // t_ devient le résultat reconstruit (complexité : O(1))
}

//-------------------------------------------------------------------------------------------------
std::vector<double> Tracks::ComputeMSD(const track_t& track, const size_t start, size_t end)
{
	if (end == 0 || end > track.size()) { end = track.size(); }   // Borne de End (pour le comportement par défaut)
	if (end < 2 + start) { return std::vector<double>(); }        // Si l'on a moins de 2 éléments à analyser
	const size_t L = track[end - 1].plane_ - track[start].plane_; // Lag maximum possible
	// if (L <= 0) { return std::vector<double>(); } // Impossible : les trajectoires n'ont pas de points sur le même plan et elles sont triées.

	std::vector sum(L, 0.0), res(L, -1.0);
	std::vector<size_t> cnt(L, 0);

	// Accumulation (toutes paires, regroupées par Lag)
	for (size_t i = start; i < end - 1; ++i) {
		const Point& Pi = track[i];
		for (size_t j = i + 1; j < end; ++j) {
			const Point& Pj  = track[j];
			const size_t Lag = (Pj.plane_ - Pi.plane_) - 1;
			if (Lag >= L) { continue; } // Impossible normalement, mais au cas où.
			const double D = Pj.GetDistance(Pi);
			sum[Lag] += D * D;
			cnt[Lag]++;
		}
	}

	// Sortie de taille exacte, Lag vides = -1
	for (size_t i = 0; i < L; ++i) { if (cnt[i] != 0) { res[i] = sum[i] / static_cast<double>(cnt[i]); } }
	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::vector<double>> Tracks::ComputeAllMSD() const
{
	const auto N = t_.size(); // Nombre de trajectoires
	std::vector<std::vector<double>> res(N);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(N); ++i) { res[static_cast<size_t>(i)] = (t_[i].size() < 2) ? std::vector<double>() : ComputeMSD(t_[i]); }
	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> Tracks::ComputeInstantDiffusion(const track_t& track, const ComputesParams& p)
{
	if (p.fitLength < 1 || p.fitLength > track.size()) { return std::vector<double>{}; }
	const size_t NPoint     = p.fitLength + 1; // Il faut prendre un point de plus pour avoir un vecteur de fitLength à la fin du calcul du MSD
	const double SlopeRatio = (p.is3D ? 6.0 : 4.0) * p.exposureTime; // 2*d*Δt
	const auto N            = track.size() - (p.fitLength); // Nombre de points moins la taille de la fenêtre (-1) si Nt = 3, Fit = 2, Np = 3 ⇒ N = 1
	std::vector res(N, -1.0);
	auto fitP = p;
	fitP.mode = fit_mode::linear;

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(N); ++i) {
		TracksFit fit(fitP);
		auto msd = ComputeMSD(track, i, i + NPoint);               // Il faut prendre un point de plus pour avoir un vecteur de fitLength à la fin
		if (msd.size() > p.fitLength) { msd.resize(p.fitLength); } // En cas de blinking, on supprime les lags plus longs que fitLength
		if (msd.size() < p.fitLength) { res[i] = -1; }             // En cas de fenêtre plus courte, mais impossible géré en amont
		else {
			const double Err = fit.Fit(msd);

			if ((Err < 0) || fit.HasBlink()) { res[i] = -1; }
			// Une erreur négative indique un problème durant le fitting, un blinking "fort" invalide la fenêtre
			else {
				const double S = fit.GetP()[0] / SlopeRatio; // Pente
				res[i]         = S < INSTANT_DIFFUSION_MIN ? INSTANT_DIFFUSION_MIN : S;
			}
		}
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::vector<double>> Tracks::ComputeAllInstantDiffusion(const ComputesParams& p) const
{
	const auto N = t_.size(); // Nombre de trajectoires
	std::vector<std::vector<double>> res(N);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(N); ++i) { res[static_cast<size_t>(i)] = ComputeInstantDiffusion(t_[i], p); }
	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> Tracks::ComputeFitMetrics(const track_t& track, const ComputesParams& p)
{
	std::vector res(getFitResSize(p.mode), -1.0);
	const double L = static_cast<double>(track.back().plane_ - track.front().plane_);
	res[0]         = L + 1; // Delta de longueur de la track (en nombre de plans)
	res[1]         = 0.0;
	for (const auto& P : track) { res[1] += P.intensity_; } // Calcul de l'intensité totale

	const auto MSD          = ComputeMSD(track);                     // Calcul du MSD
	const double SlopeRatio = (p.is3D ? 6.0 : 4.0) * p.exposureTime; // 2*d*Δt

	// Fit linéaire des N premiers (avec fitLength p[0])
	auto fitP = p;
	fitP.mode = fit_mode::linear;
	TracksFit fit(fitP);
	double err               = fit.Fit(MSD);
	std::vector<double> pFit = fit.GetP();
	// Une erreur négative indique un problème durant le fitting, un blinking invalide la fenêtre de base
	if ((err < 0) || fit.HasBlink()) { return std::vector<double>{}; }

	double s = pFit[0] / SlopeRatio; // Pente
	res[2]   = s < INSTANT_DIFFUSION_MIN ? INSTANT_DIFFUSION_MIN : s;
	res[3]   = pFit[1];
	res[4]   = err;

	// Fit sur l'ensemble des MSD
	fitP.mode      = p.mode;
	fitP.fitLength = MSD.size();

	fit.SetMode(fitP);
	err  = fit.Fit(MSD);
	pFit = fit.GetP();
	if (err < 0) { return std::vector<double>{}; } // Si échec, on arrête.
	// (ici le blinking est autorisé, on fait toute la trajectoire, il y en aura forcément en cas de reconnexion)

	// Remplissage suivant le mode
	switch (p.mode) {
		case fit_mode::linear: {
			s      = pFit[0] / SlopeRatio; // Pente
			res[5] = s < INSTANT_DIFFUSION_MIN ? INSTANT_DIFFUSION_MIN : s;
			res[6] = pFit[1];
			res[7] = err;

			break;
		}
		case fit_mode::power: {
			res[5] = pFit[0]; // Coefficient de puissance
			res[6] = pFit[1]; // Coefficient multiplicateur
			res[7] = err;

			const auto &Tb = track.back(), Tf = track.front();
			const auto D   = Tb.GetDistance(Tf) / static_cast<double>(Tb.plane_ - Tf.plane_); // Distance moyenne par frame
			res[8]         = D / p.exposureTime;                                              // Vitesse moyenne
			break;
		}
		case fit_mode::exponential: {
			if (pFit[1] + pFit[3] <= 0) { return std::vector<double>{}; } // Contrainte de positivité.
			res[5] = pFit[1];                                             // Plateau
			res[6] = pFit[2] / p.exposureTime;                            // Tau // Divisé par temps d'exposition
			res[7] = pFit[3];                                             // b
			res[8] = err;
			res[9] = std::sqrt(0.75 * (pFit[1] + pFit[3]));
			break;
		}
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::vector<double>> Tracks::ComputeAllFitMetrics(const ComputesParams& p) const
{
	const auto N = t_.size(); // Nombre de trajectoires
	std::vector<std::vector<double>> res(N);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(N); ++i) { res[static_cast<size_t>(i)] = ComputeFitMetrics(t_[i], p); }
	return res;
}

//-------------------------------------------------------------------------------------------------
double Tracks::getMaxDistance(const size_t mode, const double duration, const double maxSpeed)
{
	switch (mode) {
		case 0: return maxSpeed;                                // Immobile
		default: case 1: return maxSpeed * std::sqrt(duration); // Diffusive
		case 2: return maxSpeed * duration;                     // Linear
	}
}

//-------------------------------------------------------------------------------------------------
size_t Tracks::getFitResSize(const fit_mode mode)
{
	switch (mode) {
		case fit_mode::linear: return 8;       // [L, I, Da, Db, Derr, a, b, err]
		case fit_mode::power: return 9;        // [L, I, Da, Db, Derr, α, a, err, μ(dist)]
		case fit_mode::exponential: return 10; // [L, I, Da, Db, Derr, a, τ, b, err, sqrt(3/4*(b+a))]
	}
	return 0;
}

//-------------------------------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const Tracks& tracks)
{
	os << "[\n";
	for (size_t i = 0; i < tracks.t_.size(); ++i) {
		os << "  Track " << i + 1 << ": | ";
		for (const auto& P : tracks.t_[i]) { os << P << " | "; }
		os << "\n";
	}
	os << "]";
	return os;
}
} // namespace palm_tracer::tracking
