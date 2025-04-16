#include "Segmentation2D.h"
#include <algorithm>
#include <limits>
#include <queue>
#include <tuple>   // std::tie
#include <utility> // std::pair
#include "tools/NoError.h"

#define EPSILON 0.0000001f

void Segmentation2D::desalloc()
{
	InfoStruct2D* ptr = headInfo_;
	while (ptr) {
		headInfo_ = ptr;
		ptr = ptr->next_;
		delete headInfo_;
	}
	tabInfo_.clear();
	nbInfoStruct_ = 0;
	headInfo_ = nullptr;
}

int Segmentation2D::alloc(const size_t nbStruct)
{
	desalloc();

	tabInfo_.resize(nbStruct + 1);

	if (nbStruct > 0) {
		headInfo_ = new InfoStruct2D(1);

		tabInfo_[0] = nullptr;
		tabInfo_[1] = headInfo_;
		InfoStruct2D* ptr = headInfo_;
		for (size_t i = 2; i <= nbStruct; ++i) {
			ptr->next_ = new InfoStruct2D(i);
			tabInfo_[i] = ptr->next_;
			if (ptr->next_ == nullptr) { return S_ALLOC_ERR; }
			ptr = ptr->next_;
		}
	}

	nbInfoStruct_ = nbStruct;
	return RET_OK;
}


size_t Segmentation2D::FilterData(std::vector<unsigned short>& seg, const double* param, const unsigned short* useParam)
{
	InfoStruct2D *ptr = headInfo_, *ptrMax, *ptrOld = nullptr;
	size_t ind = 0;
	bool unselected;
	double valMax;
	const size_t NbInfoStructOld = nbInfoStruct_;

	char* tabBound = nullptr;
	if (useParam[12] == 1) {
		tabBound = new char[nbInfoStruct_ + 1];
		for (size_t i = 0; i < nbInfoStruct_ + 1; ++i) { tabBound[i] = 0; }
		for (size_t i = 0; i < height_; ++i) {
			if (seg[i * width_] != 0) { tabBound[seg[i * width_]] = 1; }
			if (seg[i * width_ + width_ - 1] != 0) { tabBound[seg[i * width_ + width_ - 1]] = 1; }
		}
		for (size_t i = 1; i < width_ - 1; ++i) {
			if (seg[i] != 0) { tabBound[seg[i]] = 1; }
			if (seg[(height_ - 1) * width_ + i] != 0) { tabBound[seg[(height_ - 1) * width_ + i]] = 1; }
		}
	}

	while (ptr) {
		ind++;
		unselected = false;
		const auto& Data = ptr->data_;
		const auto MorphRatio = std::min(Data.sigma.x, Data.sigma.y) / std::max(Data.sigma.x, Data.sigma.y);
		if ((useParam[0] && Data.surf <= param[0]) ||                    //minSurf
		    (useParam[1] && Data.surf >= param[1]) ||                    //maxSurf
		    (useParam[2] && Data.mean * Data.surf <= param[2]) ||        //minInten
		    (useParam[3] && Data.mean * Data.surf >= param[3]) ||        //maxInten
		    (useParam[4] && Data.max.val <= param[4]) ||                 //minMaxi
		    (useParam[5] && Data.max.val >= param[5]) ||                 //maxMaxi
		    (useParam[7] && MorphRatio >= param[7]) ||                   //minMorphParam
		    (useParam[8] && MorphRatio <= param[8]) ||                   //maxMorphParam
		    (useParam[12] && tabBound[ind] == 1)) { unselected = true; } //no struct on boundaries


		if (unselected) {
			if (ptrOld == nullptr) {
				headInfo_ = ptr->next_;
				delete ptr;
				ptr = headInfo_;
			}
			else {
				ptrOld->next_ = ptr->next_;
				delete ptr;
				ptr = ptrOld->next_;
			}
			tabInfo_[ind] = nullptr;
			nbInfoStruct_--;
		}
		else {
			ptrOld = ptr;
			ptr = ptr->next_;
		}
	}
	if (useParam[12] == 1) { delete [] tabBound; }


	if (useParam[6] == 1) {
		size_t oldIndex;
		ptr = headInfo_;
		ind = 0;
		ptrOld = nullptr;
		while (ptr) {
			ind++;
			unselected = false;
			for (size_t i = 1; i <= NbInfoStructOld && !unselected; i++) {
				if (tabInfo_[i]) {
					const auto& Data = ptr->data_;
					const auto& TabData = tabInfo_[i]->data_;
					const double D = sqrt(static_cast<double>((Data.max.x - TabData.max.x) * (Data.max.x - TabData.max.x)
					                                          + (Data.max.y - TabData.max.y) * (Data.max.y - TabData.max.y)));
					if (D > 1 && D < param[6] && TabData.max.val > Data.max.val) { unselected = true; }
				}
			}

			if (unselected) {
				if (ptrOld == nullptr) {
					headInfo_ = ptr->next_;
					oldIndex = ptr->data_.index;
					delete ptr;
					ptr = headInfo_;
				}
				else {
					ptrOld->next_ = ptr->next_;
					oldIndex = ptr->data_.index;
					delete ptr;
					ptr = ptrOld->next_;
				}
				tabInfo_[oldIndex] = nullptr;
				nbInfoStruct_--;
			}
			else {
				ptrOld = ptr;
				ptr = ptr->next_;
			}
		}
	}


	for (size_t indexFilter = 9; indexFilter <= 11; ++indexFilter) {
		if (useParam[indexFilter] == 1 && nbInfoStruct_ > static_cast<size_t>(param[indexFilter]) && static_cast<int>(param[indexFilter]) > 0) {
			//N first structure with the biggest surface (indexFilter=9), the biggest surface (indexFilter=10) or the biggest total intensity (indexFilter=11)
			ptrOld = nullptr;
			for (size_t i = 0; i < static_cast<size_t>(param[indexFilter]); ++i) {
				if (ptrOld == nullptr) { ptr = headInfo_; }
				else { ptr = ptrOld->next_; }

				switch (indexFilter) {
					case 9: valMax = ptr->data_.surf;
						break;
					case 10: valMax = ptr->data_.surf * ptr->data_.mean;
						break;
					case 11:
					default: valMax = ptr->data_.max.val;
						break;
				}

				ptrMax = nullptr;

				while (ptr->next_) {
					switch (indexFilter) {
						case 9: if (valMax <= ptr->next_->data_.surf) {
								valMax = ptr->next_->data_.surf;
								ptrMax = ptr;
							}
							break;
						case 10: if (valMax <= ptr->next_->data_.surf * ptr->next_->data_.mean) {
								valMax = ptr->next_->data_.surf * ptr->next_->data_.mean;
								ptrMax = ptr;
							}
							break;
						case 11: if (valMax <= ptr->next_->data_.max.val) {
								valMax = ptr->next_->data_.max.val;
								ptrMax = ptr;
							}
							break;
						default: ;
					}
					ptr = ptr->next_;
				}
				if (ptrMax) {
					if (ptrOld) {
						ptr = ptrOld->next_;
						ptrOld->next_ = ptrMax->next_;
						ptrMax->next_ = ptrMax->next_->next_;
						ptrOld->next_->next_ = ptr;
						ptrOld = ptrOld->next_;
					}
					else {
						ptr = headInfo_;
						headInfo_ = ptrMax->next_;
						ptrMax->next_ = headInfo_->next_;
						headInfo_->next_ = ptr;
						ptrOld = headInfo_;
					}
				}
				else if (ptrOld == nullptr) { ptrOld = headInfo_; }
				else { ptrOld = ptrOld->next_; }
			}
			ptr = ptrOld->next_;
			ptrOld->next_ = nullptr;
			while (ptr) {
				ptrMax = ptr->next_;
				tabInfo_[ptr->data_.index] = nullptr;
				delete ptr;
				nbInfoStruct_--;
				ptr = ptrMax;
			}
		}
	}


	ptr = headInfo_;
	for (size_t i = 1; i <= nbInfoStruct_; ++i) {
		ptr->data_.index = i;
		ptr = ptr->next_;
	}


	for (size_t i = 0; i < size_; ++i) {
		if (seg[i] != 0) {
			if (tabInfo_[seg[i]] == nullptr) { seg[i] = 0; }
			else { seg[i] = static_cast<unsigned short>(tabInfo_[seg[i]]->data_.index); }
		}
	}


	tabInfo_.clear();
	tabInfo_.resize(nbInfoStruct_ + 1);
	tabInfo_[0] = nullptr;
	tabInfo_[1] = headInfo_;
	ptr = headInfo_;
	for (ind = 2; ind <= nbInfoStruct_; ind++) {
		tabInfo_[ind] = ptr->next_;
		ptr = ptr->next_;
	}

	return nbInfoStruct_;
}

int Segmentation2D::ComputeStatistics(const std::vector<double>& data, const std::vector<unsigned short>& seg, const std::vector<double>& waveData,
                                      const size_t gaussFit, const bool doRealoc)
{
	double val;
	size_t cptPoints = 0;
	int ret;

	if (doRealoc) { if ((ret = alloc(nbInfoStruct_)) < 0) { return ret; } }

	//initalisation of array

	for (size_t ind = 1; ind <= nbInfoStruct_; ++ind) {
		tabInfo_[ind]->data_.peri = 0;
		tabInfo_[ind]->data_.surf = 0;
		tabInfo_[ind]->data_.morpho = 0.;
		tabInfo_[ind]->data_.max.x = 0;
		tabInfo_[ind]->data_.max.y = 0;
		tabInfo_[ind]->data_.max.val = 0.;
		tabInfo_[ind]->data_.min.x = 0;
		tabInfo_[ind]->data_.min.y = 0;
		tabInfo_[ind]->data_.min.val = 0;
		tabInfo_[ind]->data_.stdev = 0.;
		tabInfo_[ind]->data_.moment.x = 0.;
		tabInfo_[ind]->data_.moment.y = 0.;
		tabInfo_[ind]->data_.momentSurf.x = 0.;
		tabInfo_[ind]->data_.momentSurf.y = 0.;
		tabInfo_[ind]->data_.mean = 0.;
		tabInfo_[ind]->data_.intensity = 0.;
		tabInfo_[ind]->data_.var.xx = 0.;
		tabInfo_[ind]->data_.var.yy = 0.;
		tabInfo_[ind]->data_.var.xy = 0.;
		tabInfo_[ind]->data_.sigma.x = 0.;
		tabInfo_[ind]->data_.sigma.y = 0.;
		tabInfo_[ind]->data_.angle = 0.;
		tabInfo_[ind]->data_.boundingBox.xMax = 0;
		tabInfo_[ind]->data_.boundingBox.yMax = 0;
		tabInfo_[ind]->data_.boundingBox.xMin = height_;
		tabInfo_[ind]->data_.boundingBox.yMin = width_;
		tabInfo_[ind]->data_.nbSubstructures = 0;
	}

	for (size_t i = 0, idx = 0; i < height_; ++i) {
		for (size_t j = 0; j < width_; ++j, ++idx) {
			if (seg[idx] > 0) {
				cptPoints++;
				val = data[idx];
				auto& currData = tabInfo_[seg[idx]]->data_;

				const double Di = static_cast<double>(i), Dj = static_cast<double>(j);

				currData.surf++;
				currData.moment.x += val * Di;
				currData.moment.y += val * Dj;
				currData.momentSurf.x += Di;
				currData.momentSurf.y += Dj;
				currData.mean += val;
				currData.var.xx += val * Di * Di;
				currData.var.yy += val * Dj * Dj;
				currData.var.xy += val * Di * Dj;
				currData.stdev += val * val;

				currData.boundingBox.xMax = std::max<size_t>(currData.boundingBox.xMax, i);
				currData.boundingBox.yMax = std::max<size_t>(currData.boundingBox.yMax, j);
				currData.boundingBox.xMin = std::min<size_t>(currData.boundingBox.xMin, i);
				currData.boundingBox.yMin = std::min<size_t>(currData.boundingBox.yMin, j);


				if (val > currData.max.val || static_cast<int>(currData.surf) == 1) {
					currData.max.val = val;
					currData.max.x = i;
					currData.max.y = j;
				}

				if (val < currData.min.val || static_cast<int>(currData.surf) == 1) {
					currData.min.val = val;
					currData.min.x = i;
					currData.min.y = j;
				}

				if (i == 0) { currData.peri++; }
				if (i == height_ - 1) { currData.peri++; }
				if (j == 0) { currData.peri++; }
				if (j == width_ - 1) { currData.peri++; }
				if (i != 0 && i != height_ - 1 && j != 0 && j != width_ - 1) {
					if (seg[(i - 1) * width_ + j] != seg[idx]) { currData.peri++; }
					if (seg[(i + 1) * width_ + j] != seg[idx]) { currData.peri++; }
					if (seg[idx - 1] != seg[idx]) { currData.peri++; }
					if (seg[idx + 1] != seg[idx]) { currData.peri++; }
				}
			}
		}
	}


	for (size_t ind = 1; ind <= nbInfoStruct_; ind++) {
		auto& currData = tabInfo_[ind]->data_;
		if (currData.surf > 0) {
			if (gaussFit == 5) {
				currData.moment.x = 0.;
				currData.moment.y = 0.;
				double mean = 0.;

				for (auto ii = -2; ii <= 2; ii++) {     // height_
					for (auto jj = -2; jj <= 2; jj++) { // width_
						const int Refx = ii + currData.max.x, Refy = jj + currData.max.y;

						if (Refx >= 0 && Refx < height_ && Refy < width_ && Refy >= 0) {
							val = static_cast<double>(waveData[(Refx) * width_ + (Refy)]);
							currData.moment.x += val * static_cast<double>(Refx);
							currData.moment.y += val * static_cast<double>(Refy);
							mean += val;
						}
					}
				}

				currData.moment.x /= mean;
				currData.moment.y /= mean;
				currData.mean /= currData.surf;
				currData.intensity = currData.mean;

				/***********************************/
			}
			else {
				currData.moment.x /= currData.mean;
				currData.moment.y /= currData.mean;

				currData.momentSurf.x /= currData.surf;
				currData.momentSurf.y /= currData.surf;
				currData.var.xx /= currData.mean;
				currData.var.yy /= currData.mean;
				currData.var.xy /= currData.mean;
				currData.var.xx -= currData.moment.x * currData.moment.x;
				currData.var.yy -= currData.moment.y * currData.moment.y;
				currData.var.xy -= currData.moment.x * currData.moment.y;
				currData.intensity = currData.mean;
				currData.stdev -= currData.mean;
				currData.mean /= currData.surf;
				currData.stdev /= currData.surf;
				currData.stdev = sqrt(currData.stdev);
			}
		}
	}
	return cptPoints;
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
	std::vector< wsFusion2D> ptrStr(nbInfoStruct_ + 1);
	int label;
	bool change = true;
	while (change) {
		const size_t OldNbInfoStruct = nbInfoStruct_;

		for (size_t i = 1; i <= nbInfoStruct_; i++) {
			ptrStr[i].maxPotentialWithNeighbor = -1e10;
			ptrStr[i].minPotentialOnTheStructure = +1e10;
			ptrStr[i].maxPotentialOnTheStructure = -1e10;
			ptrStr[i].structureInContact = 0;
			ptrStr[i].kept = true;
			ptrStr[i].locked = false;
			ptrStr[i].structureToBecaume = i;
			ptrStr[i].newIndex = 0;
		}

		size_t idx = 0;
		for (size_t i = 0; i < height_; ++i) {
			for (size_t j = 0; j < width_; ++j, ++idx) {
				if (seg[idx] > 0) {
					ptrStr[seg[idx]].minPotentialOnTheStructure = std::min(ptrStr[seg[idx]].minPotentialOnTheStructure, potential[idx]);
					if (ptrStr[seg[idx]].maxPotentialOnTheStructure < potential[idx]) {
						ptrStr[seg[idx]].x = i;
						ptrStr[seg[idx]].y = j;
						ptrStr[seg[idx]].maxPotentialOnTheStructure = potential[idx];
					}
					for (size_t ii = i == 0 ? 0 : i - 1; ii <= std::min(i + 1, height_ - 1); ii++) {
						for (size_t jj = j == 0 ? 0 : j - 1; jj <= std::min(j + 1, width_ - 1); jj++) {
							if (ii != i || jj != j) {
								if (seg[ii * width_ + jj] != 0) {
									if (seg[idx] != seg[ii * width_ + jj]) {
										if (potential[idx] > ptrStr[seg[idx]].maxPotentialWithNeighbor) {
											ptrStr[seg[idx]].maxPotentialWithNeighbor = potential[idx];
											ptrStr[seg[idx]].structureInContact = seg[ii * width_ + jj];
										}
									}
								}
							}
						}
					}
				}
			}
		}


		change = false;
		for (size_t i = 1; i <= OldNbInfoStruct; ++i) {
			if (ptrStr[i].structureInContact > 0) {
				if (ptrStr[i].locked == false && ptrStr[ptrStr[i].structureInContact].locked == false) {
					const double Ratio = (ptrStr[i].maxPotentialOnTheStructure - ptrStr[i].maxPotentialWithNeighbor) / (
						                     ptrStr[i].maxPotentialOnTheStructure - ptrStr[i].minPotentialOnTheStructure);

					if (Ratio < thresholdRatio) {
						if (ptrStr[i].maxPotentialOnTheStructure < ptrStr[ptrStr[i].structureInContact].maxPotentialOnTheStructure) {
							ptrStr[i].structureToBecaume = ptrStr[i].structureInContact;
							nbInfoStruct_--;
							ptrStr[ptrStr[i].structureInContact].locked = true;
							ptrStr[i].locked = true;
							change = true;
						}
					}
				}
			}
		}

		label = 0;

		for (size_t i = 1; i <= OldNbInfoStruct; ++i) {
			if (ptrStr[i].structureToBecaume == i) {
				label++;
				ptrStr[i].newIndex = label;
			}
		}

		for (size_t i = 0; i < size_; ++i) { if (seg[i] != 0) { seg[i] = static_cast<unsigned short>(ptrStr[ptrStr[seg[i]].structureToBecaume].newIndex); } }
	}

	std::ranges::fill(seg, 0);

	label = 0;
	if (nbInfoStruct_ == 0) { return 0; }
	for (size_t i = 1; i <= nbInfoStruct_; ++i) {
		if (ptrStr[i].structureToBecaume == i) {
			label++;
			seg[ptrStr[i].x * width_ + ptrStr[i].y] = static_cast<unsigned short>(label);
		}
	}
	return label - 1;
}

/********************************************************************************************************************
int Segmentation2D::getLocalMax(...) compute local maximum of the image an mark them in the image seg.
std::vector<double> &image: input image
std::vector<unsigned short> &localMax: local maximum
double threshold: a point to be considered as a local maximum have to be > to threshold.
********************************************************************************************************************/
unsigned short Segmentation2D::getLocalMax(const std::vector<double>& potential, std::vector<unsigned short>& seg, const double threshold) const
{
	std::ranges::fill(seg, 0);
	unsigned short label = 1;

	for (size_t i = 0, idx = 0; i < height_; ++i) {
		for (size_t j = 0; j < width_; ++j, ++idx) {
			if (potential[idx] > threshold) {
				bool err = false;
				//test if it is a new maximum from the vicinity
				for (size_t ii = i == 0 ? 0 : i - 1; ii <= std::min(i + 1, height_ - 1) && !err; ii++) {
					for (size_t jj = j == 0 ? 0 : j - 1; jj <= std::min(j + 1, width_ - 1) && !err; jj++) {
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
	return label - 1;
}

using coord = std::pair<int, int>;

struct Compare
{
	// min-heap (priorit� la plus faible d'abord)
	bool operator()(const std::pair<double, coord>& a, const std::pair<double, coord>& b) const { return a.first > b.first; }
};

size_t Segmentation2D::watershed(const std::vector<double>& potential, std::vector<unsigned short>& seg, const double threshold) const
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

	return nbInfoStruct_;
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

	//first, find the local maximum to set the different germinates
	if ((nbInfoStruct_ = getLocalMax(potential, seg, threshold)) < 0) { return nbInfoStruct_; }

	//first watershed
	if ((ret = watershed(potential, seg, threshold)) < 0) { return ret; }

	if (static_cast<int>(thresholdRatio) == 0) { return nbInfoStruct_; }

	//filter the germinates to minimize the oversegmentation
	if ((ret = filterStructureWS(potential, seg, thresholdRatio)) < 0) { return ret; }

	//second watershed
	if ((ret = watershed(potential, seg, threshold)) < 0) { return ret; }
	return nbInfoStruct_;
}

#undef EPSILON
