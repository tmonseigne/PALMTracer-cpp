#include "Tracks.hpp"
#include "Constants.hpp"

#include <algorithm>
#include <ostream>
#include <ranges>
#include <utility>

namespace palm_tracer::tracking {
namespace col = tracking_columns;
namespace input_col = tracking_input_columns;
namespace output_col = tracking_output_columns;

//==================== DLL Input/output ====================
//-------------------------------------------------------------------------------------------------
void Tracks::LoadTracks(const double* input, const size_t nRow, const double pixelSize)
{
	// Remise à zéro.
	t_.clear();

	if (nRow == 0) { return; }

	// Les identifiants peuvent être clairsemes : reserver selon les groupes réels.
	size_t groups = 1;
	for (size_t i = 1; i < nRow; ++i) {
		if (static_cast<size_t>(input[i * col::Count + col::Track]) != static_cast<size_t>(input[(i - 1) * col::Count + col::Track])) { ++groups; }
	}
	t_.reserve(groups);
	for (size_t begin = 0; begin < nRow;) {
		const auto Id = static_cast<size_t>(input[begin * col::Count + col::Track]);
		size_t end    = begin + 1;
		while (end < nRow && static_cast<size_t>(input[end * col::Count + col::Track]) == Id) { ++end; }
		Track currentTrack;
		currentTrack.id = Id;
		currentTrack.points.reserve(end - begin);
		for (size_t i = begin; i < end; ++i) {
			const size_t R    = i * col::Count;
			const Vec3D PosPx = {.x = input[R + col::X], .y = input[R + col::Y], .z = input[R + col::Z]};
			const Vec3D Pos   = {.x = PosPx.x * pixelSize, .y = PosPx.y * pixelSize, .z = PosPx.z};
			// Ajout du point (emplace_back = construction in-place)
			currentTrack.points.emplace_back(static_cast<size_t>(input[R + col::Id]), static_cast<size_t>(input[R + col::Plane]),
			                                 PosPx, Pos, input[R + col::Intensity], input[R + col::Surface]);
		}
		t_.emplace_back(std::move(currentTrack));
		begin = end;
	}
}

//-------------------------------------------------------------------------------------------------
size_t Tracks::FillOutput(double* output) const
{
	size_t idx = 0;
	for (const auto& T : t_) {
		for (const auto& P : T.points) {
			output[idx + col::Track]     = static_cast<double>(T.id);     // Track
			output[idx + col::Plane]     = static_cast<double>(P.plane_); // Plane
			output[idx + col::Id]        = static_cast<double>(P.id_);    // Id
			output[idx + col::X]         = P.posPx_.x;                    // X
			output[idx + col::Y]         = P.posPx_.y;                    // Y
			output[idx + col::Z]         = P.posPx_.z;                    // Z
			output[idx + col::Intensity] = P.intensity_;                  // Intensité
			output[idx + col::Surface]   = P.surface_;                    // Surface
			idx += col::Count;
		}
	}
	return idx;
}

//-------------------------------------------------------------------------------------------------
uint64_t Tracks::Tracking(const double* input, uint64_t* output, const size_t nPoints, const double maxDistance)
{
	auto trackingPoints = loadPoints(input, nPoints); // Liste des input de mon tableau input, à construire
	const auto Planes   = getPlaneRanges(trackingPoints);
	track(trackingPoints, Planes, maxDistance);
	return fillTrackingIdOutput(trackingPoints, output);
}

//==================== Tracking ====================
//-------------------------------------------------------------------------------------------------
std::vector<TrackingPoint> Tracks::loadPoints(const double* points, const size_t nPoints)
{
	std::vector<TrackingPoint> res(nPoints);

	for (size_t i = 0, idx = 0; i < nPoints; ++i, idx += input_col::Count) {
		auto& p          = res[i];
		p.point.id_      = static_cast<size_t>(points[idx + input_col::Id]);
		p.point.plane_   = static_cast<size_t>(points[idx + input_col::Plane]);
		p.point.posPx_.x = points[idx + input_col::X];
		p.point.posPx_.y = points[idx + input_col::Y];
	}

	std::ranges::sort(res, [](const TrackingPoint& a, const TrackingPoint& b) {
		if (a.point.plane_ != b.point.plane_) { return a.point.plane_ < b.point.plane_; }
		return a.point.id_ < b.point.id_;
	});

	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<PlaneRange> Tracks::getPlaneRanges(const std::vector<TrackingPoint>& points)
{
	std::vector<PlaneRange> res;
	if (points.empty()) { return res; }

	const size_t MinPlane = points.front().point.plane_;
	const size_t MaxPlane = points.back().point.plane_;

	res.resize(MaxPlane - MinPlane + 1);

	size_t begin = 0;
	size_t plane = MinPlane;

	for (size_t i = 0; i < points.size(); ++i) {
		const size_t CurrentPlane = points[i].point.plane_;

		while (plane < CurrentPlane) {
			res[plane - MinPlane] = {.begin = begin, .end = begin};
			++plane;
		}

		begin = i;
		while (i + 1 < points.size() && points[i + 1].point.plane_ == CurrentPlane) { ++i; }

		res[CurrentPlane - MinPlane] = {.begin = begin, .end = i + 1};
		begin                        = i + 1;
		plane                        = CurrentPlane + 1;
	}

	return res;
}

//-------------------------------------------------------------------------------------------------
void Tracks::getCandidates(const std::vector<TrackingPoint>& points, const PlaneRange current, const PlaneRange next,
                           const double maxDistance, std::vector<Candidate>& candidates)
{
	candidates.clear();
	candidates.reserve(current.end - current.begin);
	if (current.begin == current.end || next.begin == next.end) { return; }
	const double Radius2 = maxDistance * maxDistance;
	for (size_t i = current.begin; i < current.end; ++i) {
		for (size_t j = next.begin; j < next.end; ++j) {
			const double D2 = points[i].GetSqDistance(points[j]);
			if (D2 <= Radius2) { candidates.emplace_back(i, j, D2); }
		}
	}
}

//-------------------------------------------------------------------------------------------------
void Tracks::track(std::vector<TrackingPoint>& points, const std::vector<PlaneRange>& planes, const double maxDistance)
{
	const size_t Transitions = planes.empty() ? 0 : planes.size() - 1;
	// Crée une équipe fixe : chaque thread réutilise son buffer de candidats entre les transitions.
	#pragma omp parallel
	{
		std::vector<Candidate> candidates;
		// Répartit les transitions entre les threads déjà créés ; la séparation des deux pragmas évite de recréer le buffer à chaque itération.
		#pragma omp for schedule(dynamic)
		for (int plane = 0; plane < static_cast<int>(Transitions); ++plane) {
			const PlaneRange& Current = planes[plane];
			const PlaneRange& Next    = planes[plane + 1];

			if (Current.begin == Current.end || Next.begin == Next.end) { continue; }

			// Examine toutes les paires et conserve tous les candidats admissibles.
			getCandidates(points, Current, Next, maxDistance, candidates);

			// Tri des candidats par distance croissante, puis par index croissant pour garantir la stabilité
			std::ranges::sort(candidates, [](const Candidate& a, const Candidate& b) {
				// Comparaison exacte nécessaire pour garantir un ordre faible strict ; une tolérance n'est pas transitive.
				// ReSharper disable once CppClangTidyClangDiagnosticFloatEqual
				if (a.distance2 != b.distance2) { return a.distance2 < b.distance2; }
				if (a.from != b.from) { return a.from < b.from; }
				return a.to < b.to;
			});

			// Pour chaque candidat, on relie les points si aucun des deux n'est déjà relié
			for (const Candidate& C : candidates) {
				// les points ne sont accedés pour le plan P uniquement pour next et pour le plan P+1 pour previous, donc pas de data race possible
				TrackingPoint& from = points[C.from];
				TrackingPoint& to   = points[C.to];

				// Si l'un des deux est déjà relié, on ignore ce candidat. Utile, car on a trié par distance croissante, donc on garde le plus proche.
				if (from.next != INVALID_INDEX || to.previous != INVALID_INDEX) { continue; }

				from.next   = C.to;
				to.previous = C.from;
			}
		}
	}

	// Attribution des identifiants de trajectoire
	size_t trackId = 1;
	for (TrackingPoint& point : points) {
		if (point.previous != INVALID_INDEX) { continue; }           // On ne commence que par les points qui n'ont pas de précédent (têtes de trajectoire).
		for (TrackingPoint* current = &point; current != nullptr;) { // Parcours de la chaîne de points reliés
			current->trackId = trackId;                              // Attribution de l'identifiant de trajectoire
			if (current->next == INVALID_INDEX) { break; }           // Fin de la chaîne
			current = &points[current->next];                        // Passage au point suivant
		}

		++trackId;
	}
}

//-------------------------------------------------------------------------------------------------
uint64_t Tracks::fillTrackingIdOutput(const std::vector<TrackingPoint>& points, uint64_t* output)
{
	size_t idx = 0;
	for (const auto& P : points) {
		output[idx + output_col::Id]    = P.point.id_; // Id
		output[idx + output_col::Track] = P.trackId;   // Track
		idx += output_col::Count;
	}
	return idx;
}


//-------------------------------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const Tracks& tracks)
{
	os << "[\n";
	for (const auto& T : tracks.t_) {
		os << "  Track " << T.id << ": | ";
		for (const auto& P : T.points) { os << P << " | "; }
		os << "\n";
	}
	os << "]";
	return os;
}
} // namespace palm_tracer::tracking
