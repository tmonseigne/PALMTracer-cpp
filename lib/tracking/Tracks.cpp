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
	std::cout << "Blinking Reconnection Begin\n";
	std::cout << *this << "\n";

	std::cout << "Blinking Reconnection Finish\n";
	std::cout << *this << "\n";
}

//-------------------------------------------------------------------------------------------------
void Tracks::connect(const size_t track1, const size_t track2) { std::cout << "Connect track " << track1 << " and " << track2 << "\n"; }

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
