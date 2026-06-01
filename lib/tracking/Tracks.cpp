#include "Tracks.hpp"
#include "Constants.hpp"

#include <cmath>
#include <omp.h>

#include "TracksFit.hpp"

namespace palm_tracer::tracking {
//==================== DLL Input/output ====================
//-------------------------------------------------------------------------------------------------
void Tracks::LoadTracks(const double* input, const size_t nRow, const double pixelSize)
{
	// Remise à zéro.
	t_.clear();

	// Récupération du nombre de trajectoires pour préallouer le tableau
	const auto NTracks = static_cast<size_t>(input[(nRow - 1) * N_TRACKING]); // Uniquement le numéro de trajectoire de la dernière ligne
	t_.reserve(NTracks);

	Track currentTrack;                              // Trajectoire en cours (réutilisée entre groupes)
	currentTrack.id = static_cast<size_t>(input[0]); // Identifiant original non forcément égal à 0.
	currentTrack.points.reserve(64);                 // Heuristique: évite quelques reallocs sur petites tracks

	for (size_t i = 0, r = 0; i < nRow; ++i, r += N_TRACKING) {
		const auto Id = static_cast<size_t>(input[r]); // Identifiant de la trajectoire du point actuel

		// Nouveau groupe de trajectoire détecté
		if (Id != currentTrack.id) {
			if (!currentTrack.points.empty()) { t_.emplace_back(std::move(currentTrack)); }

			currentTrack    = Track{};
			currentTrack.id = Id;
			currentTrack.points.reserve(64);
		}

		const Vec3D PosPx = {.x = input[r + 3], .y = input[r + 4], .z = input[r + 5]};
		const Vec3D Pos   = {.x = PosPx.x * pixelSize, .y = PosPx.y * pixelSize, .z = PosPx.z};

		// Ajout du point (emplace_back = construction in-place)
		currentTrack.points.emplace_back(static_cast<size_t>(input[r + 2]), // Id
		                                 static_cast<size_t>(input[r + 1]), // Plan
		                                 PosPx,                             // Position en pixel
		                                 Pos,                               // Position
		                                 input[r + 6],                      // Intensité
		                                 input[r + 7]);                     // Surface
	}

	if (!currentTrack.points.empty()) { t_.emplace_back(std::move(currentTrack)); } // Ne pas oublier d'ajouter la dernière trajectoire
}

//-------------------------------------------------------------------------------------------------
size_t Tracks::FillOutput(double* output) const
{
	size_t idx = 0;
	for (const auto& T : t_) {
		for (const auto& P : T.points) {
			output[idx++] = static_cast<double>(T.id);     // Track
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

//-------------------------------------------------------------------------------------------------
uint64_t Tracks::Tracking(const double* input, uint64_t* output, const size_t nPoints, const double maxDistance)
{
	auto trackingPoints = loadPoints(input, nPoints); // Liste des input de mon tableau input, à construire
	const auto Planes   = getPlaneRanges(trackingPoints);
	track(trackingPoints, Planes, maxDistance);
	return fillTrackingIdOutput(trackingPoints, output);
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
	for (size_t head = 0; head < N; ++head) {      // Trajectoire à tester
		if (t_[head].points.empty()) { continue; } // Au cas où, mais théoriquement impossible.
		if (taken[head]) { continue; }             // Trajectoire déjà absorbée par une chaîne précédente
		std::vector<size_t> chain;                 // Nouvelle Trajectoire
		chain.push_back(head);                     // Ajout de la trajectoire initiale

		// On tente d'allonger la chaîne de façon gloutonne
		size_t current = head;
		while (true) {
			if (t_[current].points.empty()) { break; }            // Au cas où, mais théoriquement impossible.
			const Point& Lp = t_[current].points.back();          // Dernier point de la trajectoire i.
			size_t bJ       = std::numeric_limits<size_t>::max(), // Indice (j) du meilleur point.
			       bT       = std::numeric_limits<size_t>::max(); // Différence temporelle avec le meilleur point.
			double bD       = std::numeric_limits<double>::max(); // Différence spatiale avec le meilleur point.

			// Cherche le meilleur j non encore pris (hors tête) et non vide
			for (size_t j = head + 1; j < N; ++j) {
				if (j == current || taken[j] || t_[j].points.empty()) { continue; }
				const Point& Fp = t_[j].points.front();   // Premier point de la trajectoire j.
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
	std::vector<Track> new_tracks;     // Création d'une liste de trajectoires temporaire
	new_tracks.reserve(chains.size()); // On a notre nombre de trajectoires finales
	size_t track_id = 1;
	for (const auto& Ch : chains) {
		// Concatène tous les points des indices de la chaîne dans un seul track
		Track acc;
		acc.id = track_id++; // acc.id = t_[Ch.front()].id; si l'on veut garder les numéro précédents, mais il y aura des trou dnas les numéros de trajectoires.

		// (option : size hint) somme des tailles pour limiter les reallocs
		size_t totalPts = 0;
		for (const size_t Idx : Ch) { totalPts += t_[Idx].points.size(); }
		acc.points.reserve(totalPts);

		for (const size_t Idx : Ch) {
			auto& tr = t_[Idx];
			acc.points.insert(acc.points.end(), std::make_move_iterator(tr.points.begin()), std::make_move_iterator(tr.points.end()));
			tr.points.clear(); // libère les références, pas la capacité (sans importance ici)
		}
		new_tracks.emplace_back(std::move(acc));
	}
	t_.swap(new_tracks); // t_ devient le résultat reconstruit (complexité : O(1))
}

//-------------------------------------------------------------------------------------------------
std::vector<double> Tracks::ComputeMSD(const Track& track, const size_t start, size_t end)
{
	const auto& P = track.points;
	if (end == 0 || end > P.size()) { end = P.size(); }    // Borne de End (pour le comportement par défaut)
	if (end < 2 + start) { return std::vector<double>(); } // Si l'on a moins de 2 éléments à analyser
	const size_t L = P[end - 1].plane_ - P[start].plane_;  // Lag maximum possible
	// if (L ≤ 0) { return std::vector<double>(); } // Impossible : les trajectoires n'ont pas de points sur le même plan et elles sont triées.

	std::vector sum(L, 0.0), res(L, -1.0);
	std::vector<size_t> cnt(L, 0);

	// Accumulation (toutes paires, regroupées par Lag)
	for (size_t i = start; i < end - 1; ++i) {
		const Point& Pi = P[i];
		for (size_t j = i + 1; j < end; ++j) {
			const Point& Pj  = P[j];
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
	for (int i = 0; i < static_cast<int>(N); ++i) {
		const auto Idx     = static_cast<size_t>(i);
		const Track& Track = t_[Idx];

		if (Track.points.size() < 2) { res[Idx] = std::vector<double>(); }
		else {
			std::vector<double> msd = ComputeMSD(Track);
			msd.insert(msd.begin(), static_cast<double>(Track.id)); // Ajout de l'identifiant de trajectoire en première colonne
			res[Idx] = std::move(msd);
		}
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> Tracks::ComputeInstantDiffusion(const Track& track, const ComputesParams& p)
{
	const auto& P = track.points;
	if (p.fitLength < 1 || p.fitLength > P.size()) { return std::vector<double>{}; }
	const size_t NPoint     = p.fitLength + 1; // Il faut prendre un point de plus pour avoir un vecteur de fitLength à la fin du calcul du MSD
	const double SlopeRatio = (p.is3D ? 6.0 : 4.0) * p.exposureTime; // 2*d*Δt
	const auto N            = P.size() - (p.fitLength); // Nombre de points moins la taille de la fenêtre (-1) si Nt = 3, Fit = 2, Np = 3 ⇒ N = 1
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
	for (int i = 0; i < static_cast<int>(N); ++i) {
		const size_t Idx   = static_cast<size_t>(i);
		const Track& Track = t_[Idx];

		auto diff = ComputeInstantDiffusion(Track, p);
		diff.insert(diff.begin(), static_cast<double>(Track.id)); // Ajout de l'identifiant de trajectoire en première colonne
		res[Idx] = std::move(diff);
	}
	return res;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> Tracks::ComputeFitMetrics(const Track& track, const ComputesParams& p)
{
	const auto& P = track.points;
	std::vector res(getFitResSize(p.mode), -1.0);
	const auto L = static_cast<double>(P.back().plane_ - P.front().plane_);
	size_t idx   = 0;
	res[idx++]   = L + 1; // Delta de longueur de la track (en nombre de plans)
	res[idx]     = 0.0;
	for (const auto& Point : P) { res[idx] += Point.intensity_; } // Calcul de l'intensité totale
	idx++;

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

	double s   = pFit[0] / SlopeRatio; // Pente
	res[idx++] = s < INSTANT_DIFFUSION_MIN ? INSTANT_DIFFUSION_MIN : s;
	res[idx++] = pFit[1];
	res[idx++] = err;

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
			s          = pFit[0] / SlopeRatio; // Pente
			res[idx++] = s < INSTANT_DIFFUSION_MIN ? INSTANT_DIFFUSION_MIN : s;
			res[idx++] = pFit[1];
			res[idx++] = err;

			break;
		}
		case fit_mode::power: {
			res[idx++] = pFit[0]; // Coefficient de puissance
			res[idx++] = pFit[1]; // Coefficient multiplicateur
			res[idx++] = err;

			const auto &Tb = P.back(), Tf = P.front();
			const auto D   = Tb.GetDistance(Tf) / static_cast<double>(Tb.plane_ - Tf.plane_); // Distance moyenne par frame
			res[idx++]     = D / p.exposureTime;                                              // Vitesse moyenne
			break;
		}
		case fit_mode::exponential: {
			if (pFit[1] + pFit[3] <= 0) { return std::vector<double>{}; } // Contrainte de positivité.
			res[idx++] = pFit[1];                                         // Plateau
			res[idx++] = pFit[2] / p.exposureTime;                        // Tau // Divisé par temps d'exposition
			res[idx++] = pFit[3];                                         // b
			res[idx++] = err;
			res[idx++] = std::sqrt(0.75 * (pFit[1] + pFit[3]));
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
	for (int i = 0; i < static_cast<int>(N); ++i) {
		const size_t Idx   = static_cast<size_t>(i);
		const Track& Track = t_[Idx];

		auto metrics = ComputeFitMetrics(Track, p);
		metrics.insert(metrics.begin(), static_cast<double>(Track.id)); // Ajout de l'identifiant de trajectoire en première colonne
		res[Idx] = std::move(metrics);
	}
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

//==================== Tracking ====================
//-------------------------------------------------------------------------------------------------
std::vector<TrackingPoint> Tracks::loadPoints(const double* points, const size_t nPoints)
{
	std::vector<TrackingPoint> res(nPoints);

	for (size_t i = 0, idx = 0; i < nPoints; ++i) {
		auto& p          = res[i];
		p.point.id_      = static_cast<size_t>(points[idx++]);
		p.point.plane_   = static_cast<size_t>(points[idx++]);
		p.point.posPx_.x = points[idx++];
		p.point.posPx_.y = points[idx++];
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
void Tracks::track(std::vector<TrackingPoint>& points, const std::vector<PlaneRange>& planes, const double maxDistance)
{
	const double MaxDistance2 = maxDistance * maxDistance; // Pré-calcul du carré de la distance maximale pour éviter les racines carrées répétées

	#pragma omp parallel for schedule(dynamic)
	for (int plane = 0; plane < static_cast<int>(planes.size() - 1); ++plane) {
		const PlaneRange& Current = planes[plane];
		const PlaneRange& Next    = planes[plane + 1];

		if (Current.begin == Current.end || Next.begin == Next.end) { continue; }

		std::vector<Candidate> candidates;
		candidates.reserve((Current.end - Current.begin) * (Next.end - Next.begin));

		for (size_t i = Current.begin; i < Current.end; ++i) {
			for (size_t j = Next.begin; j < Next.end; ++j) {
				const double D2 = points[i].GetSqDistance(points[j]);
				if (D2 <= MaxDistance2) { candidates.emplace_back(i, j, D2); }
			}
		}

		std::ranges::sort(candidates, [](const Candidate& a, const Candidate& b) {
			if (std::abs(a.distance2 - b.distance2) > std::numeric_limits<double>::epsilon()) { return a.distance2 < b.distance2; }
			if (a.from != b.from) { return a.from < b.from; }
			return a.to < b.to;
		});

		for (const Candidate& C : candidates) {
			TrackingPoint& from = points[C.from];
			TrackingPoint& to   = points[C.to];

			if (from.next != INVALID_INDEX || to.previous != INVALID_INDEX) { continue; }

			from.next   = C.to;
			to.previous = C.from;
		}
	}

	size_t trackId = 1;
	for (TrackingPoint& point : points) {
		if (point.previous != INVALID_INDEX) { continue; }
		for (TrackingPoint* current = &point; current != nullptr;) {
			current->trackId = trackId;
			if (current->next == INVALID_INDEX) { break; }
			current = &points[current->next];
		}

		++trackId;
	}
}

//-------------------------------------------------------------------------------------------------
uint64_t Tracks::fillTrackingIdOutput(const std::vector<TrackingPoint>& points, uint64_t* output)
{
	size_t idx = 0;
	for (const auto& P : points) {
		output[idx++] = P.point.id_; // Id
		output[idx++] = P.trackId;   // Track
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
