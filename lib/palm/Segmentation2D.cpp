#include "Segmentation2D.hpp"
#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <queue>
#include <utility> // std::pair

#include "Constants.hpp"
#include "Structs.hpp"

namespace palm_tracer {
//-------------------------------------------------------------------------------------------------
namespace {
struct Filter {
	bool enabled;
	double value;
};

// Création d'un tableau constexpr avec tous à false + 0.0 sauf
constexpr std::array<Filter, N_LOCALIZATIONS> MakeDefaultFilter()
{
	std::array<Filter, N_LOCALIZATIONS> params{};
	for (auto& p : params) { p = {.enabled = false, .value = 0.0}; }
	// Paramètres à filtrer
	params[0] = {.enabled = true, .value = VOLUME_MIN};
	return params;
}

constexpr std::array<Filter, N_LOCALIZATIONS> SEGMENTATION_FILTER = MakeDefaultFilter();
///< Tableau de paramètres de filtre (min/max de surface, intensité, morphologie...)

/// <summary>Vérifie si le filtre est actif et valide.</summary>
/// <param name="f">Filtre.</param>
/// <param name="value">Valeur.</param>
/// <param name="lessThan">Comparatif (True si plus petit que, False sinon).</param>
/// <returns></returns>
bool CheckFilter(const Filter& f, const double value, const bool lessThan) { return f.enabled && (lessThan ? value <= f.value : value >= f.value); }
} // namespace

//-------------------------------------------------------------------------------------------------
Segmentation2D::Segmentation2D(const ATrous& decomposition, const double threshold, const double ratio)
	: dim_(decomposition.dim_), size_(decomposition.size_), ratio_(ratio)
{
	// Initialization
	map_.resize(size_);
	output_.resize(size_);

	#pragma omp parallel for
	for (int i = 0; i < static_cast<int>(size_); ++i) {
		if (decomposition.wDecomp_[PLANE_ID][i] > threshold) {
			output_[i] = decomposition.wDecomp_[PLANE_ID][i];
			map_[i]    = 1;
		}
		else {
			map_[i]    = 0;
			output_[i] = 0;
		}
	}
	WatershedIntelligent();
	ComputeStatistics();
	FilterData();
}

//-------------------------------------------------------------------------------------------------
void Segmentation2D::WatershedIntelligent()
{
	getLocalMax();                                 // first, find the local maximum to set the different germinates
	watershed();                                   // first watershed
	if (static_cast<int>(ratio_) == 0) { return; } // Ratio test
	filterStructureWS();                           // filter the germinates to minimize the oversegmentation
	watershed();                                   // second watershed
}

//-------------------------------------------------------------------------------------------------
void Segmentation2D::FilterData()
{
	const size_t OldSize = infos_.size();

	// Marquage des régions sur les bords
	std::vector bounds(infos_.size(), false);
	if (SEGMENTATION_FILTER[12].enabled) {
		for (size_t i = 0; i < dim_.y; ++i) {
			const auto Idx = i * dim_.x;
			if (map_[Idx] != 0) { bounds[map_[Idx]] = true; }
			if (map_[Idx + dim_.x - 1] != 0) { bounds[map_[Idx + dim_.x - 1]] = true; }
		}
		for (size_t i = 1; i < dim_.x - 1; ++i) {
			if (map_[i] != 0) { bounds[map_[i]] = true; }
			if (map_[(dim_.y - 1) * dim_.x + i] != 0) { bounds[map_[(dim_.y - 1) * dim_.x + i]] = true; }
		}
	}

	// Première passe : filtrage par critère
	auto isFiltered = [&](const SegmentationInfo& info, const size_t index) {
		const double Mean_Surf = info.mean_ * info.surf_;

		return
			CheckFilter(SEGMENTATION_FILTER[0], info.surf_, true) ||  // min Surf
			CheckFilter(SEGMENTATION_FILTER[1], info.surf_, false) || // max Surf
			CheckFilter(SEGMENTATION_FILTER[2], Mean_Surf, true) ||   // min Inten
			CheckFilter(SEGMENTATION_FILTER[3], Mean_Surf, false) ||  // max Inten
			(SEGMENTATION_FILTER[12].enabled && bounds[index]);       // no struct on boundaries
	};

	std::vector<SegmentationInfo> filtered;
	filtered.reserve(OldSize);
	std::vector<size_t> oldToNew(OldSize, 0);
	size_t newIdx = 0;

	for (size_t i = 0; i < OldSize; ++i) {
		if (!isFiltered(infos_[i], i)) {
			oldToNew[i]      = newIdx + 1;
			filtered.emplace_back(infos_[i]);
			++newIdx;
		}
	}

	// Réaffectation
	infos_ = std::move(filtered);

	// Mise à jour de l'image segmentée
	for (size_t i = 0; i < size_; ++i) {
		if (map_[i] > 0) {
			const size_t OldIndex = map_[i] - 1;
			map_[i]               = oldToNew[OldIndex] > 0 ? static_cast<uint16_t>(oldToNew[OldIndex]) : 0;
		}
	}
}

//-------------------------------------------------------------------------------------------------
void Segmentation2D::ComputeStatistics()
{
	for (auto& info : infos_) { info.Reset(); }

	// Accumulation des mesures ne pas paralléliser
	for (size_t i = 0, idx = 0; i < dim_.y; ++i) {
		for (size_t j = 0; j < dim_.x; ++j, ++idx) {
			const uint16_t S = map_[idx];
			if (S == 0) { continue; }
			auto& info       = infos_[S - 1];
			const double Val = output_[idx];

			// Mesures
			const auto Di = static_cast<double>(i), Dj = static_cast<double>(j);
			info.surf_++;
			info.x_ += Val * Dj;
			info.y_ += Val * Di;
			info.mean_ += Val;
			info.stdev_ += Val * Val;
		}
	}

	// Post-traitement
	for (auto& info : infos_) {
		if (static_cast<int>(info.surf_) == 0) { continue; }
		info.x_ /= info.mean_;
		info.y_ /= info.mean_;
		info.intensity_ = info.mean_;
		info.stdev_     = std::sqrt((info.stdev_ - info.mean_) / info.surf_);
		info.mean_ /= info.surf_;
	}
}


//-------------------------------------------------------------------------------------------------
void Segmentation2D::alloc(const size_t nbStruct)
{
	infos_.resize(nbStruct);
}

//-------------------------------------------------------------------------------------------------
void Segmentation2D::filterStructureWS()
{
	const size_t NbStruct = infos_.size();
	if (NbStruct == 0) { return; }

	std::vector<WSFusion> wsFusions(NbStruct);

	bool change = true;
	int label   = 0;

	while (change) {
		change = false;
		// Initialisation des structures
		for (size_t i = 0; i < NbStruct; ++i) {
			auto& s                    = wsFusions[i];
			s.maxPotentialWithNeighbor = -1e10;
			s.minPotential             = +1e10;
			s.maxPotential             = -1e10;
			s.structInContact          = 0;
			s.locked                   = false;
			s.structToBecome           = static_cast<int>(i);
			s.newIndex                 = 0;
		}

		// Analyse des structures
		for (size_t i = 0, idx = 0; i < dim_.y; ++i) {
			for (size_t j = 0; j < dim_.x; ++j, ++idx) {
				const uint16_t Sid = map_[idx];
				if (Sid == 0) { continue; }
				auto& s = wsFusions[Sid - 1];

				s.minPotential = std::min(s.minPotential, output_[idx]);

				if (output_[idx] > s.maxPotential) {
					s.maxPotential = output_[idx];
					s.x            = static_cast<int>(i);
					s.y            = static_cast<int>(j);
				}

				// Recherche des voisins
				for (size_t ii = std::max(i, static_cast<size_t>(1)) - 1; ii <= std::min(i + 1, dim_.y - 1); ++ii) {
					for (size_t jj = std::max(j, static_cast<size_t>(1)) - 1; jj <= std::min(j + 1, dim_.x - 1); ++jj) {
						if (ii == i && jj == j) { continue; }

						const uint16_t NeighborSid = map_[ii * dim_.x + jj];
						if (NeighborSid != 0 && NeighborSid != Sid) {
							if (output_[idx] > s.maxPotentialWithNeighbor) {
								s.maxPotentialWithNeighbor = output_[idx];
								s.structInContact          = NeighborSid - 1;
							}
						}
					}
				}
			}
		}

		// Fusion des structures si conditions remplies
		for (size_t i = 0; i < NbStruct; ++i) {
			auto& s = wsFusions[i];
			if (s.structInContact <= 0) { continue; }

			const auto& Neighbor = wsFusions[s.structInContact];

			if (!s.locked && !Neighbor.locked) {
				const double Delta = s.maxPotential - s.minPotential;
				if (std::fabs(Delta) <= DBL_EPSILON) { continue; }

				const double Ratio = (s.maxPotential - s.maxPotentialWithNeighbor) / Delta;
				if (Ratio < ratio_) {
					if (s.maxPotential < Neighbor.maxPotential) {
						s.structToBecome                    = s.structInContact;
						wsFusions[s.structInContact].locked = true;
						s.locked                            = true;
						change                              = true;
					}
				}
			}
		}

		// Réattribution des nouvelles étiquettes
		label = 0;
		for (int i = 0; i <= static_cast<int>(wsFusions.size()); ++i) { if (wsFusions[i].structToBecome == i) { wsFusions[i].newIndex = ++label; } }

		for (uint16_t& s : map_) {
			if (s != 0) {
				const auto& Sid = wsFusions[s - 1];
				s               = static_cast<uint16_t>(wsFusions[Sid.structToBecome].newIndex);
			}
		}
	}

	// Nettoyage infos_ une fois les fusions terminées
	std::vector<SegmentationInfo> updatedInfos(static_cast<size_t>(label));

	for (size_t oldIdx = 0; oldIdx < wsFusions.size(); ++oldIdx) {
		const auto& Fusion = wsFusions[oldIdx];
		if (Fusion.newIndex > 0) {
			updatedInfos[Fusion.newIndex - 1]        = infos_[oldIdx];
		}
	}
	infos_ = std::move(updatedInfos);

	// Marquage d'un seul point par structure pour Watershed final
	std::ranges::fill(map_, 0);
	for (int i = 0; i < static_cast<int>(wsFusions.size()); ++i) {
		if (wsFusions[i].structToBecome == i) { map_[wsFusions[i].x * dim_.x + wsFusions[i].y] = static_cast<uint16_t>(wsFusions[i].newIndex); }
	}
}

//-------------------------------------------------------------------------------------------------
void Segmentation2D::getLocalMax()
{
	std::ranges::fill(map_, 0);
	uint16_t label = 1;

	for (size_t i = 0, idx = 0; i < dim_.y; ++i) {
		for (size_t j = 0; j < dim_.x; ++j, ++idx) {
			if (output_[idx] > THRESHOLD) {
				bool err = false;
				//test if it is a new maximum from the vicinity
				for (size_t ii = i == 0 ? 0 : i - 1; ii <= std::min(i + 1, dim_.y - 1) && !err; ++ii) {
					for (size_t jj = j == 0 ? 0 : j - 1; jj <= std::min(j + 1, dim_.x - 1) && !err; ++jj) {
						if (ii != i || jj != j) { if (output_[idx] < output_[ii * dim_.x + jj]) { err = true; } }
					}
				}
				if (!err) { //new maximum found
					uint16_t maxi = 0;
					err           = false;
					for (size_t ii = i == 0 ? 0 : i - 1; ii <= std::min(i + 1, dim_.y - 1) && !err; ++ii) {
						const size_t Idx1 = ii * dim_.x;
						for (size_t jj = j == 0 ? 0 : j - 1; jj <= std::min(j + 1, dim_.x - 1) && !err; ++jj) {
							if (ii != i || jj != j) {
								if (maxi < map_[Idx1 + jj]) {
									maxi = map_[Idx1 + jj];
									err  = true;
								}
							}
						}
					}
					if (maxi > 0) { map_[idx] = maxi; }
					else {
						map_[idx] = label;
						label++;
					}
				}
			}
		}
	}
	alloc(label - 1);
}

//-------------------------------------------------------------------------------------------------
void Segmentation2D::watershed()
{
	std::priority_queue<Pixel> heap;

	for (size_t x = 0, idx = 0; x < dim_.y; ++x) {
		for (size_t y = 0; y < dim_.x; ++y, ++idx) {
			if (map_[idx] == 0) { continue; }

			bool vois = false;
			for (int dx = -1; dx <= 1 && !vois; ++dx) {
				const int Nx = static_cast<int>(x) + dx;
				if (Nx < 0 || Nx >= static_cast<int>(dim_.y)) { continue; }
				const size_t RowOffset = static_cast<size_t>(Nx) * dim_.x;

				for (int dy = -1; dy <= 1 && !vois; ++dy) {
					const int Ny = static_cast<int>(y) + dy;
					if ((dx == 0 && dy == 0) || Ny < 0 || Ny >= static_cast<int>(dim_.x)) { continue; }

					if (map_[RowOffset + Ny] == 0) { vois = true; }
				}
			}

			if (vois) { heap.emplace(Pixel{.x = x, .y = y, .v = -output_[idx]}); }
		}
	}

	int it = 0;
	while (!heap.empty()) {
		const Pixel P = heap.top();
		heap.pop();
		++it;

		for (int dx = -1; dx <= 1; ++dx) {
			const int Nx = static_cast<int>(P.x) + dx;
			if (Nx < 0 || Nx >= static_cast<int>(dim_.y)) { continue; }

			for (int dy = -1; dy <= 1; ++dy) {
				const int Ny = static_cast<int>(P.y) + dy;
				if ((dx == 0 && dy == 0) || Ny < 0 || Ny >= static_cast<int>(dim_.x)) { continue; }

				const size_t Nid = Nx * dim_.x + Ny;
				if (output_[Nid] > THRESHOLD && map_[Nid] == 0) {
					const double Cle = -output_[Nid] + static_cast<double>(it) * EPSILON;
					heap.emplace(Pixel{.x = static_cast<size_t>(Nx), .y = static_cast<size_t>(Ny), .v = Cle});
					map_[Nid] = map_[P.x * dim_.x + P.y];
				}
			}
		}
	}
}
} // namespace palm_tracer
