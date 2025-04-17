#include "Segmentation2D.h"
#include <algorithm>
#include <queue>
#include <unordered_map>
#include <utility> // std::pair

#define EPSILON 0.0000001f

void Segmentation2D::alloc(const size_t nbStruct)
{
	desalloc();
	infos_.resize(nbStruct);
	for (size_t i = 0; i < nbStruct; ++i) {
		infos_[i].index_ = i;
	}
}

void Segmentation2D::desalloc() { infos_.clear(); }


size_t Segmentation2D::FilterData(std::vector<unsigned short>& seg, const double* param, const unsigned short* useParam)
{
	const size_t OldSize = infos_.size();

	// Marquage des régions sur les bords
	std::vector bounds(infos_.size(), false);
	if (useParam[12] == 1) {
		for (size_t i = 0; i < height_; ++i) {
			if (seg[i * width_] != 0) { bounds[seg[i * width_]] = true; }
			if (seg[i * width_ + width_ - 1] != 0) { bounds[seg[i * width_ + width_ - 1]] = true; }
		}
		for (size_t i = 1; i < width_ - 1; ++i) {
			if (seg[i] != 0) { bounds[seg[i]] = true; }
			if (seg[(height_ - 1) * width_ + i] != 0) { bounds[seg[(height_ - 1) * width_ + i]] = true; }
		}
	}

	// Première passe : filtrage par critère
	auto isFiltered = [&](const Info2D& info, const size_t index)
	{
		const double Mean_Surf = info.mean_ * info.surf_;
		const double Morph_Ratio = std::min(info.sigma_.x, info.sigma_.y) / std::max(info.sigma_.x, info.sigma_.y);

		return
			(useParam[0] && info.surf_ <= param[0]) ||    // minSurf
			(useParam[1] && info.surf_ >= param[1]) ||    // maxSurf
			(useParam[2] && Mean_Surf <= param[2]) ||    // minInten
			(useParam[3] && Mean_Surf >= param[3]) ||    // maxInten
			(useParam[4] && info.max_.val <= param[4]) || // minMaxi
			(useParam[5] && info.max_.val >= param[5]) || // maxMaxi
			(useParam[7] && Morph_Ratio >= param[7]) ||  // minMorphParam
			(useParam[8] && Morph_Ratio <= param[8]) ||  // maxMorphParam
			(useParam[12] && bounds[index]);             // no struct on boundaries
	};

	std::vector<Info2D> filtered;
	filtered.reserve(OldSize);
	std::vector<size_t> oldToNew(OldSize, 0);
	size_t newIdx = 0;

	for (size_t i = 0; i < OldSize; ++i) {
		if (!isFiltered(infos_[i], i)) {
			oldToNew[i] = newIdx + 1;
			infos_[i].index_ = static_cast<int>(newIdx);
			filtered.emplace_back(infos_[i]);
			++newIdx;
		}
	}

	// Deuxième passe : filtrage des structures proches de plus brillantes (param[6])
	if (useParam[6]) {
		const size_t NewSize = filtered.size();
		std::vector<Info2D> finalFiltered;
		filtered.reserve(NewSize);
		for (size_t i = 0; i < NewSize; ++i) {
			bool unselected = false;

			for (size_t j = 0; j < NewSize && !unselected; ++j) {
				if (i == j) { continue; }
				const double Dx = static_cast<double>(filtered[i].max_.x - filtered[j].max_.x),
				             Dy = static_cast<double>(filtered[i].max_.y - filtered[j].max_.y);
				const double Dist = std::sqrt(Dx * Dx + Dy * Dy);

				if (Dist > 1.0 && Dist < param[6] && filtered[j].max_.val > filtered[i].max_.val) { unselected = true; }
			}

			if (!unselected) { finalFiltered.emplace_back(filtered[i]); }
		}

		filtered = std::move(finalFiltered);
	}

	// Réaffectation
	infos_ = std::move(filtered);

	// Mise à jour de l’image segmentée
	for (size_t i = 0; i < size_; ++i) {
		if (seg[i] > 0) {
			const size_t OldIndex = seg[i] - 1;
			seg[i] = oldToNew[OldIndex] > 0 ? static_cast<unsigned short>(oldToNew[OldIndex]) : 0;
		}
	}
	return infos_.size();
}

int Segmentation2D::ComputeStatistics(const std::vector<double>& data, const std::vector<unsigned short>& seg, const std::vector<double>& waveData,
                                      const size_t gaussFit)
{
	double val;
	size_t cptPoints = 0;

	for (auto& info : infos_) { info.Reset(height_, width_); }

	// Accumulation des mesures
	for (size_t i = 0, idx = 0; i < height_; ++i) {
		for (size_t j = 0; j < width_; ++j, ++idx) {
			const unsigned short s = seg[idx];
			if (s == 0) { continue; }
			cptPoints++;
			auto& info = infos_[s - 1];
			val = data[idx];

			// Mesures
			const double Di = static_cast<double>(i), Dj = static_cast<double>(j);
			info.surf_++;
			info.moment_.x += val * Di;
			info.moment_.y += val * Dj;
			info.momentSurf_.x += Di;
			info.momentSurf_.y += Dj;
			info.mean_ += val;
			info.var_.xx += val * Di * Di;
			info.var_.yy += val * Dj * Dj;
			info.var_.xy += val * Di * Dj;
			info.stdev_ += val * val;

			info.boundingBox_.xMax = std::max(info.boundingBox_.xMax, i);
			info.boundingBox_.yMax = std::max(info.boundingBox_.yMax, j);
			info.boundingBox_.xMin = std::min(info.boundingBox_.xMin, i);
			info.boundingBox_.yMin = std::min(info.boundingBox_.yMin, j);

			if (val > info.max_.val || static_cast<int>(info.surf_) == 1) {
				info.max_.val = val;
				info.max_.x = i;
				info.max_.y = j;
			}
			if (val < info.min_.val || static_cast<int>(info.surf_) == 1) {
				info.min_.val = val;
				info.min_.x = i;
				info.min_.y = j;
			}

			// Bordure
			if (i == 0 || i == height_ - 1 || j == 0 || j == width_ - 1) { info.peri_++; }
			else {
				if (seg[(i - 1) * width_ + j] != s) { info.peri_++; }
				if (seg[(i + 1) * width_ + j] != s) { info.peri_++; }
				if (seg[idx - 1] != s) { info.peri_++; }
				if (seg[idx + 1] != s) { info.peri_++; }
			}
		}
	}

	// Post-traitement
	for (auto& info : infos_) {
		if (static_cast<int>(info.surf_) == 0) { continue; }

		if (gaussFit == 5) {
			info.moment_ = {.x = 0., .y = 0.};
			double mean = 0.;

			for (int ii = -2; ii <= 2; ++ii) {
				for (int jj = -2; jj <= 2; ++jj) {
					const int X = info.max_.x + ii, Y = info.max_.y + jj;
					if (X >= 0 && X < static_cast<int>(height_) && Y >= 0 && Y < static_cast<int>(width_)) {
						val = waveData[X * width_ + Y];
						info.moment_.x += val * X;
						info.moment_.y += val * Y;
						mean += val;
					}
				}
			}

			info.moment_.x /= mean;
			info.moment_.y /= mean;
			info.mean_ /= info.surf_;
			info.intensity_ = info.mean_;
		}
		else {
			info.moment_.x /= info.mean_;
			info.moment_.y /= info.mean_;
			info.momentSurf_.x /= info.surf_;
			info.momentSurf_.y /= info.surf_;
			info.var_.xx = info.var_.xx / info.mean_ - info.moment_.x * info.moment_.x;
			info.var_.yy = info.var_.yy / info.mean_ - info.moment_.y * info.moment_.y;
			info.var_.xy = info.var_.xy / info.mean_ - info.moment_.x * info.moment_.y;
			info.intensity_ = info.mean_;
			info.stdev_ = std::sqrt((info.stdev_ - info.mean_) / info.surf_);
			info.mean_ /= info.surf_;
		}
	}

	return static_cast<int>(cptPoints);
}

/********************************************************************************************************************
int Segmentation2D::filterStructureWS(...) filter the structure of the watershed.
std::vector<double> &potential: input image
std::vector<unsigned short> &localMax: local maximum
double thresholdRatio: value to distinguish if 2 strucutres in contact have to be considered as a single structure (
							oversegmentation)
********************************************************************************************************************/
int Segmentation2D::filterStructureWS(const std::vector<double>& potential, std::vector<unsigned short>& seg, const double thresholdRatio)
{
	const size_t NbStruct = infos_.size();
	if (NbStruct == 0) { return 0; }

	std::vector<ws_fusion_2d> wsFusions(NbStruct);

	bool change = true;
	int label = 0;

	while (change) {
		change = false;
		// Initialisation des structures
		for (size_t i = 0; i < NbStruct; ++i) {
			auto& s = wsFusions[i];
			s.maxPotentialWithNeighbor = -1e10;
			s.minPotential = +1e10;
			s.maxPotential = -1e10;
			s.structInContact = 0;
			s.kept = true;
			s.locked = false;
			s.structToBecome = static_cast<int>(i);
			s.newIndex = 0;
		}

		// Analyse des structures
		for (size_t i = 0, idx = 0; i < height_; ++i) {
			for (size_t j = 0; j < width_; ++j, ++idx) {
				const unsigned short Sid = seg[idx];
				if (Sid == 0) { continue; }
				auto& s = wsFusions[Sid - 1];

				s.minPotential = std::min(s.minPotential, potential[idx]);

				if (potential[idx] > s.maxPotential) {
					s.maxPotential = potential[idx];
					s.x = static_cast<int>(i);
					s.y = static_cast<int>(j);
				}

				// Recherche voisin
				for (size_t ii = std::max(i, static_cast<size_t>(1)) - 1; ii <= std::min(i + 1, height_ - 1); ++ii) {
					for (size_t jj = std::max(j, static_cast<size_t>(1)) - 1; jj <= std::min(j + 1, width_ - 1); ++jj) {
						if (ii == i && jj == j) { continue; }

						const unsigned short NeighborSid = seg[ii * width_ + jj];
						if (NeighborSid != 0 && NeighborSid != Sid) {
							if (potential[idx] > s.maxPotentialWithNeighbor) {
								s.maxPotentialWithNeighbor = potential[idx];
								s.structInContact = NeighborSid - 1;
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
				if (std::abs(Delta) <= DBL_EPSILON) { continue; }

				const double Ratio = (s.maxPotential - s.maxPotentialWithNeighbor) / Delta;
				if (Ratio < thresholdRatio) {
					if (s.maxPotential < Neighbor.maxPotential) {
						s.structToBecome = s.structInContact;
						wsFusions[s.structInContact].locked = true;
						s.locked = true;
						change = true;
					}
				}
			}
		}

		// Réattribution des nouvelles étiquettes
		label = 0;
		for (int i = 0; i <= static_cast<int>(wsFusions.size()); ++i) { if (wsFusions[i].structToBecome == i) { wsFusions[i].newIndex = ++label; } }

		for (unsigned short& s : seg) {
			if (s != 0) {
				const auto& Sid = wsFusions[s - 1];
				s = static_cast<unsigned short>(wsFusions[Sid.structToBecome].newIndex);
			}
		}
	}

	// Nettoyage infos_ une fois les fusions terminées
	std::vector<Info2D> updatedInfos(static_cast<size_t>(label));

	for (size_t oldIdx = 0; oldIdx < wsFusions.size(); ++oldIdx) {
		const auto& Fusion = wsFusions[oldIdx];
		if (Fusion.newIndex > 0) {
			updatedInfos[Fusion.newIndex - 1] = infos_[oldIdx];
			updatedInfos[Fusion.newIndex - 1].index_ = Fusion.newIndex - 1;
		}
	}
	infos_ = std::move(updatedInfos);

	// Marquage d’un seul point par structure pour Watershed final
	std::ranges::fill(seg, 0);
	for (int i = 0; i <= static_cast<int>(wsFusions.size()); ++i) {
		if (wsFusions[i].structToBecome == i) { seg[wsFusions[i].x * width_ + wsFusions[i].y] = static_cast<unsigned short>(wsFusions[i].newIndex); }
	}

	return static_cast<int>(infos_.size());
}


/********************************************************************************************************************
void Segmentation2D::getLocalMax(...) compute local maximum of the image an mark them in the image seg.
std::vector<double> &image: input image
std::vector<unsigned short> &localMax: local maximum
double threshold: a point to be considered as a local maximum have to be > to threshold.
********************************************************************************************************************/
void Segmentation2D::getLocalMax(const std::vector<double>& potential, std::vector<unsigned short>& seg, const double threshold)
{
	std::ranges::fill(seg, 0);
	unsigned short label = 1;

	for (size_t i = 0, idx = 0; i < height_; ++i) {
		for (size_t j = 0; j < width_; ++j, ++idx) {
			if (potential[idx] > threshold) {
				bool err = false;
				//test if it is a new maximum from the vicinity
				for (size_t ii = i == 0 ? 0 : i - 1; ii <= std::min(i + 1, height_ - 1) && !err; ++ii) {
					for (size_t jj = j == 0 ? 0 : j - 1; jj <= std::min(j + 1, width_ - 1) && !err; ++jj) {
						if (ii != i || jj != j) { if (potential[idx] < potential[ii * width_ + jj]) { err = true; } }
					}
				}
				if (!err) { //new maximum found
					unsigned short maxi = 0;
					err = false;
					for (size_t ii = i == 0 ? 0 : i - 1; ii <= std::min(i + 1, height_ - 1) && !err; ++ii) {
						const size_t Idx1 = ii * width_;
						for (size_t jj = j == 0 ? 0 : j - 1; jj <= std::min(j + 1, width_ - 1) && !err; ++jj) {
							if (ii != i || jj != j) {
								if (maxi < seg[Idx1 + jj]) {
									maxi = seg[Idx1 + jj];
									err = true;
								}
							}
						}
					}
					if (maxi > 0) { seg[idx] = maxi; }
					else {
						seg[idx] = label;
						label++;
					}
				}
			}
		}
	}
	alloc(label - 1);
}

using coord = std::pair<int, int>;

struct Compare
{
	// min-heap (priorit� la plus faible d'abord)
	bool operator()(const std::pair<double, coord>& a, const std::pair<double, coord>& b) const { return a.first > b.first; }
};

void Segmentation2D::watershed(const std::vector<double>& potential, std::vector<unsigned short>& seg, const double threshold) const
{
	std::priority_queue<std::pair<double, coord>, std::vector<std::pair<double, coord>>, Compare> heap;

	for (size_t x = 0, idx = 0; x < height_; ++x) {
		for (size_t y = 0; y < width_; ++y, ++idx) {
			if (seg[idx] == 0) { continue; }

			bool vois = false;
			for (int dx = -1; dx <= 1 && !vois; ++dx) {
				const int Nx = x + dx;
				if (Nx < 0 || Nx >= height_) { continue; }
				const size_t RowOffset = static_cast<size_t>(Nx) * width_;

				for (int dy = -1; dy <= 1 && !vois; ++dy) {
					const int Ny = y + dy;
					if ((dx == 0 && dy == 0) || Ny < 0 || Ny >= width_) { continue; }

					if (seg[RowOffset + Ny] == 0) { vois = true; }
				}
			}

			if (vois) {
				double cle = -potential[idx];
				heap.emplace(cle, coord{x, y});
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
			if (nx < 0 || nx >= height_) { continue; }

			for (int dy = -1; dy <= 1; ++dy) {
				int ny = p.second + dy;
				if ((dx == 0 && dy == 0) || ny < 0 || ny >= width_) { continue; }

				const size_t Nid = nx * width_ + ny;
				if (potential[Nid] > threshold && seg[Nid] == 0) {
					double newCle = -potential[Nid] + static_cast<double>(it) * EPSILON;
					heap.emplace(newCle, coord{nx, ny});
					seg[Nid] = seg[p.first * width_ + p.second];
				}
			}
		}
	}
}

/********************************************************************************************************************
int Segmentation2D::watershedIntelligent(...) compute watershed segmentation insied the reagion upper the threshold value.
It tries to minimize the oversegmentation.
std::vector<double> &potential: defines the potential function the makes the watershed segmented area growing
std::vector<unsigned short> &seg: it contains the result of the watershed segmentation.
double threshold: define the threshold value to distinguish the segmented structures
double thresholdRatio: value 0<= X <=1. Makes a fusion to minimize the oversegmentation
						if X=0 simple watershed
						if X=1 no watershed
return value: number of detected structures
********************************************************************************************************************/
int Segmentation2D::WatershedIntelligent(const std::vector<double>& potential, std::vector<unsigned short>& seg, const double threshold,
                                         const double thresholdRatio)
{
	int ret;

	getLocalMax(potential, seg, threshold); //first, find the local maximum to set the different germinates
	watershed(potential, seg, threshold);   //first watershed

	if (static_cast<int>(thresholdRatio) == 0) { return static_cast<int>(infos_.size()); }

	//filter the germinates to minimize the oversegmentation
	if ((ret = filterStructureWS(potential, seg, thresholdRatio)) < 0) { return ret; }

	watershed(potential, seg, threshold); //second watershed
	return static_cast<int>(infos_.size());
}

#undef EPSILON
