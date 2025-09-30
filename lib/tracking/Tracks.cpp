#include "Tracks.hpp"

#include "Constants.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
void Tracks::LoadTracking(const double* input, const size_t nRow)
{
	// Remise à zéro.
	t_.clear();

	//Récupération du nombre de trajectoires pour préallouer le tableau
	const auto NTracks = static_cast<size_t>(input[(nRow - 1) * N_TRACKING]);
	t_.reserve(NTracks);

	// Parcours des lignes
	track_t track;           // Trajectoire en cours (réutilisée entre groupes)
	track.reserve(16);       // Heuristique: évite quelques reallocs sur petites tracks
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
		track.emplace_back(Point(static_cast<size_t>(input[r + 2]),                         // Id
		                         static_cast<size_t>(input[r + 1]),                         // Plan
		                         {.x = input[r + 3], .y = input[r + 4], .z = input[r + 5]}, // Position
		                         input[r + 6],                                              // Intensité
		                         input[r + 7]));                                            // Surface
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
			output[idx++] = P.pos_.x;                      // X
			output[idx++] = P.pos_.y;                      // Y
			output[idx++] = P.pos_.z;                      // Z
			output[idx++] = P.intensity_;                  // Intensity
			output[idx++] = P.surface_;                    // Surface
		}
	}
	return idx;
}

//-------------------------------------------------------------------------------------------------
void Tracks::BlinkingReconnection(const size_t mode, const size_t maxDuration, const double maxSpeed)
{
	if (t_.empty()) { return; }
	const size_t N = t_.size();
	std::vector taken(N, false);             // True si l’index appartient déjà à une chaîne (hors tête)
	std::vector<std::vector<size_t>> chains; // Liste des chaînes d’indices (ex: { [2,4,5], [3], [6,7], ... })
	chains.reserve(N);                       // Taille maximale possible

	// --- Construction des chaînes d’indices de connexion ---
	for (size_t head = 0; head < N; ++head) { // Trajectoire à tester
		if (t_[head].empty()) { continue; }   // Au cas où, mais théoriquement impossible.
		if (taken[head]) { continue; }        // Trajectoire déjà absorbée par une chaîne précédente
		std::vector<size_t> chain;            // Nouvelle Trajectoire
		chain.push_back(head);                // Ajout de la trajectoire initiale

		// On tente d’allonger la chaîne de façon gloutonne
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
				const double D    = (Fp.pos_ - Lp.pos_).Norm(),                              // Calcul de la distance.
				             MaxD = getMaxDistance(mode, static_cast<double>(Dp), maxSpeed); // Calcul de la distance maximale.
				if (D >= bD || D > MaxD) { continue; }                                       // Vérification de validité Spatiale.
				bJ = j;                                                                      // Mise à jour de l'élément sélectionné.
				bT = Dp;                                                                     // Mise à jour de la durée (nombre de plan).
				bD = D;                                                                      // Mise à jour de la distance.
			}
			if (bJ == std::numeric_limits<size_t>::max()) { break; } // Aucune trajectoire de trouvé, on sort du while
			chain.push_back(bJ);                                     // Ajout de la trajectoire à la list
			taken[bJ] = true;                                        // Ajout d'un token "trajectoire prise"
			current   = bJ;                                          // Mise à jour de la trajectoire courrante (qui est la dernière de la liste)
		}
		chains.emplace_back(std::move(chain)); // Ajout de la chaine de trajectoire au tableau final.
	}

	// --- Recomposition finale dans un nouveau conteneur, par move-concat ---
	std::vector<track_t> new_tracks;   // Création d'une liste de trajectoire temporaire
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
double Tracks::getMaxDistance(const size_t mode, const double duration, const double maxSpeed)
{
	switch (mode) {
		case 0: return maxSpeed;                                // Immobile
		default: case 1: return maxSpeed * std::sqrt(duration); // Diffusive
		case 2: return maxSpeed * duration;                     // Linear
	}
}

//-------------------------------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const Tracks& tracks)
{
	os << "[\n";
	for (size_t i = 0; i < tracks.t_.size(); ++i) {
		os << "  Track " << i + 1 << ": [";
		const auto& T = tracks.t_[i];
		for (size_t j = 0; j < T.size(); ++j) {
			os << "(" << T[j].id_ << ", " << T[j].plane_ << ", " << T[j].pos_ << ")";
			if (j < T.size() - 1) { os << ", "; }
		}
		os << "]\n";
	}
	os << "]";
	return os;
}
} // namespace palm_tracer
