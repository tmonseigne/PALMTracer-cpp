#include "Segmentation2D.hpp"
#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <queue>
#include <ranges>
#include <unordered_map>
#include <utility> // std::pair

#include "Constants.hpp"

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

	bool CheckFilter(const Filter& p, const double value, const bool lessThan) { return p.enabled && (lessThan ? (value <= p.value) : (value >= p.value)); }
}

//-------------------------------------------------------------------------------------------------
Segmentation2D::Segmentation2D(const ATrous& decomposition, const double threshold, const double ratio)
	: height_(decomposition.height_), width_(decomposition.width_), size_(decomposition.size_), ratio_(ratio)
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
	if (SEGMENTATION_FILTER[12].enabled == 1) {
		for (size_t i = 0; i < height_; ++i) {
			const auto Idx = i * width_;
			if (map_[Idx] != 0) { bounds[map_[Idx]] = true; }
			if (map_[Idx + width_ - 1] != 0) { bounds[map_[Idx + width_ - 1]] = true; }
		}
		for (size_t i = 1; i < width_ - 1; ++i) {
			if (map_[i] != 0) { bounds[map_[i]] = true; }
			if (map_[(height_ - 1) * width_ + i] != 0) { bounds[map_[(height_ - 1) * width_ + i]] = true; }
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
			infos_[i].index_ = static_cast<int>(newIdx);
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
			map_[i]               = oldToNew[OldIndex] > 0 ? static_cast<unsigned short>(oldToNew[OldIndex]) : 0;
		}
	}
}

//-------------------------------------------------------------------------------------------------
void Segmentation2D::ComputeStatistics()
{
	for (auto& info : infos_) { info.Reset(); }

	// Accumulation des mesures ne pas paralléliser
	for (size_t i = 0, idx = 0; i < height_; ++i) {
		for (size_t j = 0; j < width_; ++j, ++idx) {
			const unsigned short s = map_[idx];
			if (s == 0) { continue; }
			auto& info       = infos_[s - 1];
			const double Val = output_[idx];

			// Mesures
			const double Di = static_cast<double>(i), Dj = static_cast<double>(j);
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
	for (size_t i = 0; i < nbStruct; ++i) { infos_[i].index_ = i; }
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
		for (size_t i = 0, idx = 0; i < height_; ++i) {
			for (size_t j = 0; j < width_; ++j, ++idx) {
				const unsigned short Sid = map_[idx];
				if (Sid == 0) { continue; }
				auto& s = wsFusions[Sid - 1];

				s.minPotential = std::min(s.minPotential, output_[idx]);

				if (output_[idx] > s.maxPotential) {
					s.maxPotential = output_[idx];
					s.x            = static_cast<int>(i);
					s.y            = static_cast<int>(j);
				}

				// Recherche des voisins
				for (size_t ii = std::max(i, static_cast<size_t>(1)) - 1; ii <= std::min(i + 1, height_ - 1); ++ii) {
					for (size_t jj = std::max(j, static_cast<size_t>(1)) - 1; jj <= std::min(j + 1, width_ - 1); ++jj) {
						if (ii == i && jj == j) { continue; }

						const unsigned short NeighborSid = map_[ii * width_ + jj];
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

		for (unsigned short& s : map_) {
			if (s != 0) {
				const auto& Sid = wsFusions[s - 1];
				s               = static_cast<unsigned short>(wsFusions[Sid.structToBecome].newIndex);
			}
		}
	}

	// Nettoyage infos_ une fois les fusions terminées
	std::vector<SegmentationInfo> updatedInfos(static_cast<size_t>(label));

	for (size_t oldIdx = 0; oldIdx < wsFusions.size(); ++oldIdx) {
		const auto& Fusion = wsFusions[oldIdx];
		if (Fusion.newIndex > 0) {
			updatedInfos[Fusion.newIndex - 1]        = infos_[oldIdx];
			updatedInfos[Fusion.newIndex - 1].index_ = Fusion.newIndex - 1;
		}
	}
	infos_ = std::move(updatedInfos);

	// Marquage d'un seul point par structure pour Watershed final
	std::ranges::fill(map_, 0);
	for (int i = 0; i <= static_cast<int>(wsFusions.size()); ++i) {
		if (wsFusions[i].structToBecome == i) { map_[wsFusions[i].x * width_ + wsFusions[i].y] = static_cast<unsigned short>(wsFusions[i].newIndex); }
	}
}

//-------------------------------------------------------------------------------------------------
void Segmentation2D::getLocalMax()
{
	std::ranges::fill(map_, 0);
	unsigned short label = 1;

	for (size_t i = 0, idx = 0; i < height_; ++i) {
		for (size_t j = 0; j < width_; ++j, ++idx) {
			if (output_[idx] > THRESHOLD) {
				bool err = false;
				//test if it is a new maximum from the vicinity
				for (size_t ii = i == 0 ? 0 : i - 1; ii <= std::min(i + 1, height_ - 1) && !err; ++ii) {
					for (size_t jj = j == 0 ? 0 : j - 1; jj <= std::min(j + 1, width_ - 1) && !err; ++jj) {
						if (ii != i || jj != j) { if (output_[idx] < output_[ii * width_ + jj]) { err = true; } }
					}
				}
				if (!err) { //new maximum found
					unsigned short maxi = 0;
					err                 = false;
					for (size_t ii = i == 0 ? 0 : i - 1; ii <= std::min(i + 1, height_ - 1) && !err; ++ii) {
						const size_t Idx1 = ii * width_;
						for (size_t jj = j == 0 ? 0 : j - 1; jj <= std::min(j + 1, width_ - 1) && !err; ++jj) {
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
using coord_2d = std::pair<int, int>;

struct Compare {
	bool operator()(const std::pair<double, coord_2d>& a, const std::pair<double, coord_2d>& b) const { return a.first > b.first; }
};


//-------------------------------------------------------------------------------------------------
void Segmentation2D::watershed()
{
	std::priority_queue<std::pair<double, coord_2d>, std::vector<std::pair<double, coord_2d>>, Compare> heap;

	for (size_t x = 0, idx = 0; x < height_; ++x) {
		for (size_t y = 0; y < width_; ++y, ++idx) {
			if (map_[idx] == 0) { continue; }

			bool vois = false;
			for (int dx = -1; dx <= 1 && !vois; ++dx) {
				const int Nx = static_cast<int>(x) + dx;
				if (Nx < 0 || Nx >= static_cast<int>(height_)) { continue; }
				const size_t RowOffset = static_cast<size_t>(Nx) * width_;

				for (int dy = -1; dy <= 1 && !vois; ++dy) {
					const int Ny = static_cast<int>(y) + dy;
					if ((dx == 0 && dy == 0) || Ny < 0 || Ny >= static_cast<int>(width_)) { continue; }

					if (map_[RowOffset + Ny] == 0) { vois = true; }
				}
			}

			if (vois) {
				double cle = -output_[idx];
				heap.emplace(cle, coord_2d{x, y});
			}
		}
	}

	int it = 0;
	while (!heap.empty()) {
		auto [cle, p] = heap.top();
		heap.pop();
		++it;

		for (int dx = -1; dx <= 1; ++dx) {
			int nx = p.first + dx;
			if (nx < 0 || nx >= static_cast<int>(height_)) { continue; }

			for (int dy = -1; dy <= 1; ++dy) {
				int ny = p.second + dy;
				if ((dx == 0 && dy == 0) || ny < 0 || ny >= static_cast<int>(width_)) { continue; }

				const size_t Nid = nx * width_ + ny;
				if (output_[Nid] > THRESHOLD && map_[Nid] == 0) {
					double newCle = -output_[Nid] + static_cast<double>(it) * EPSILON;
					heap.emplace(newCle, coord_2d{nx, ny});
					map_[Nid] = map_[p.first * width_ + p.second];
				}
			}
		}
	}
}
}
