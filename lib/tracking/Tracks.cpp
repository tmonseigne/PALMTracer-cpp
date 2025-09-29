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
	// Parcours de toutes les trajectoires
	size_t i = 0;
	while (i < t_.size()) {
		if (t_[i].empty()) {
			t_.erase(t_.begin() + static_cast<std::ptrdiff_t>(i)); // ne pas incrémenter i : après l’effacement, l’élément suivant glisse à la position i
			continue;
		}                                                     // Au cas où, mais théoriquement impossible.
		const Point& Lp = t_[i].back();                       // Dernier point de la trajectoire i.
		size_t cI       = std::numeric_limits<size_t>::max(), // Indice courrant.
		       cT       = std::numeric_limits<size_t>::max(); // Différence temporelle courrante.
		double cD       = std::numeric_limits<double>::max(); // Différence spatiale courrante.
		// Parcours de toutes les trajectoires suivantes (pas la peine de tester les trajectoires antérieures, elles sont triées par plan du premier point)
		for (size_t j = i + 1; j < t_.size(); ++j) {
			const Point& Fp = t_[j].front();          // Premier point de la trajectoire j.
			if (Fp.plane_ <= Lp.plane_) { continue; } // Vérification de continuité temporelle.
			const size_t Dp = Fp.plane_ - Lp.plane_;  // Calcul du temps (nombre de plans entre les trajectoires).
			// (à vérifier, mais un break est possible, car les trajectoires sont triées par plan de démarrage)
			if (Dp > maxDuration || Dp > cT) { continue; }                                  // Vérification de validité temporelle.
			const double Dist    = (Fp.pos_ - Lp.pos_).Norm(),                              // Calcul de la distance.
			             MaxDist = getMaxDistance(mode, static_cast<double>(Dp), maxSpeed); // Calcul de la distance maximale.
			if (Dist >= cD || Dist > MaxDist) { continue; }                                 // Vérification de validité Spatiale.
			cI = j;                                                                         // Mise à jour de l'élément sélectionné.
			cT = Dp;                                                                        // Mise à jour de la durée (nombre de plan).
			cD = Dist;                                                                      // Mise à jour de la distance.
		}
		if (cI < t_.size()) { connect(i, cI); } // Connexion des trajectoires
		else { ++i; }                           // On passe à la trajectoire suivbante
	}
	// Autre possibilité, noter toutes les trajectoires à connecter et recopier à la toute fin dans un nouveau tableau
}

//-------------------------------------------------------------------------------------------------
double Tracks::getMaxDistance(const size_t mode, const double duration, const double maxSpeed)
{
	switch (mode) {
		case 0: return maxSpeed;
		default: case 1: return maxSpeed * std::sqrt(duration);
		case 3: return maxSpeed * duration;
	}
}

//-------------------------------------------------------------------------------------------------
void Tracks::connect(const size_t track1, const size_t track2)
{
	std::cout << "Connect track " << track1 << " and " << track2 << "\n";
	if (track1 >= t_.size() || track2 >= t_.size() || track1 == track2) { return; }
	auto& a = t_[track1];
	auto& b = t_[track2];
	// Concat par move (aucune copie de Point)
	a.reserve(a.size() + b.size()); // reservation de l'espace
	a.insert(a.end(), std::make_move_iterator(b.begin()), std::make_move_iterator(b.end()));
	b.clear();                                                  // b est désormais vide
	t_.erase(t_.begin() + static_cast<std::ptrdiff_t>(track2)); // Retire l'entrée track2 du conteneur (rétrécit t_ et décale les indices > track2)
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
