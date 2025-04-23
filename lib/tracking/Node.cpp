/*******************************************************************************************************************
									Node Class
								  Node representing one structure
*******************************************************************************************************************/

/********************************************************************************************************************
Class Node contains the caracteristic off the node, and the chain of the other nodes (time after and time before) 
links to the current objet. 
********************************************************************************************************************/

#include "Node.h"

#include <algorithm>
#include <numbers>

constexpr double WEIGHT_DISTANCE = 0.6;
constexpr double WEIGHT_INTENSITY = 0.2;
constexpr double WEIGHT_MEAN_INTENSITY = 0.2;
constexpr double MIN_ENERGY_FUSION = 0;

//=================================================================
//==================== Constructor, Destructor ====================
//=================================================================
/********************************************************************************************************************
Constructuor(size_t time) of Node class.
The current objets is at the time "time".
********************************************************************************************************************/
Node::Node(const size_t time, Parameters* params)
	: params_(params), time_(time),
	  future_(nullptr), blinkFuture_(nullptr), past_(nullptr), blinkPast_(nullptr),
	  fission1Future_(nullptr), fission2Future_(nullptr), fissionPast_(nullptr),
	  fusion1Past_(nullptr), fusion2Past_(nullptr), fusionFuture_(nullptr),
	  trackedIndex_(0), colorIndex_{},
	  neighborsPast_(nullptr), neighborsFuture_(nullptr),
	  neighborsPastCount_(0), neighborsFutureCount_(0),
	  neighborsBlinkPast_(nullptr), neighborsBlinkFuture_(nullptr),
	  neighborsBlinkPastCount_(0), neighborsBlinkFutureCount_(0) { }

/********************************************************************************************************************
Destructuor() of Node class.
Delete all the linkage objets.
********************************************************************************************************************/
Node::~Node()
{
	ChainNeighbor* l;
	while (neighborsFuture_ != nullptr) {
		l = neighborsFuture_->next_;
		delete neighborsFuture_;
		neighborsFuture_ = l;
	}

	while (neighborsPast_ != nullptr) {
		l = neighborsPast_->next_;
		delete neighborsPast_;
		neighborsPast_ = l;
	}

	trackedIndex_ = 0;
}

/*********************************************************************************************
double Node::makeInitializationConnection()
  Make the initial connections. For each strucutres, chose the one with the highest likeness
  if it is avaliable. If it is not possible no connection.

  It returns the likeness.

  New version from 2005/01/12
*********************************************************************************************/
void Node::MakeInitializationConnection()
{
	Node* nMax = nullptr;
	double max = -1e10;

	for (const ChainNeighbor* c = neighborsFuture_; c != nullptr; c = c->next_) {
		if (c->neighbor_->past_ == nullptr) {
			const double Val = Likeness(c->neighbor_);

			if (Val > max) {
				max = Val;
				nMax = c->neighbor_;
			}
		}
	}
	if (nMax != nullptr) {
		nMax->past_ = this;
		future_ = nMax;
	}
}

/********************************************************************************************************************
sortByIntersityNeighbors(ChainNeighbor* neighborsHead) sort the chained list begining at neighborsHead by the
intensity.
********************************************************************************************************************/
void Node::SortByIntersityNeighbors(ChainNeighbor** neighborsHead)
{
	if (*neighborsHead != nullptr) {
		if ((*neighborsHead)->next_ != nullptr) {
			bool change;
			do {
				change = false;
				ChainNeighbor* l = *neighborsHead;
				ChainNeighbor* lOld = nullptr;
				while (l->next_ != nullptr) {
					if (l->next_->neighbor_->info_.intensity_ > l->neighbor_->info_.intensity_) {
						ChainNeighbor* l2 = l->next_;
						if (lOld == nullptr) { *neighborsHead = l2; }
						else { lOld->next_ = l2; }
						l->next_ = l2->next_;
						l2->next_ = l;
						change = true;
						lOld = l2;
					}
					else {
						l = l->next_;
						if (lOld == nullptr) { lOld = *neighborsHead; }
					}
				}
			}
			while (change);
		}
	}
}


//==================================================
//==================== Computes ====================
//==================================================
/********************************************************************************************************************
bool metropolis(double deltaE, double temperature)
Metropolis algorithm. metropolis returns a boolean variable that issues a verdict on whether
to accept a reconfiguration that leads to a change deltaE in the objective function. If deltaE<0,
metropolis = true, while if deltaE>0, metropolis is only true with probability exp(-deltaE/temperature), where
temperature is a temperature determined by the annealing schedule.
********************************************************************************************************************/
bool Node::metropolis(const double deltaE, const double temperature) const
{
	if (deltaE >= 0.) { return true; }
	if (params_->Randomfloat() > exp(deltaE / temperature)) { return false; }
	return true;
}

/*********************************************************************************************
double distanceToBorder()
  get the distance from this to the closer border
*********************************************************************************************/
double Node::distanceToBorder() const
{
	const double X = info_.x_, Y = info_.y_;
	const double MaxX = params_->maxX_, MaxY = params_->maxY_;
	return std::min({X, MaxX - X, Y, MaxY - Y});
}

/*********************************************************************************************
double SphericEffect()
  Compute the mean of radius of the node this
*********************************************************************************************/
double Node::sphericEffect() const
{
	return sqrt(info_.surf_ / std::numbers::pi);
	//volume of a sphere is 4/3 * PI *r^3, surface is PI * r^2
	//if constexpr (Parameters::N_DIM == 2) { return sqrt(info_.surf_ / std::numbers::pi); }
	//else { return pow(info_.GetVolume() * 3.0 / 4.0 / std::numbers::pi, 1.0 / 3.0); }
}

/*********************************************************************************************
double Node::likenessDistance(Node* n)
  Compute the likeness from the distance between this and n.
  New version from 2005/01/12
*********************************************************************************************/
double Node::likenessDistance(const Node* n) const
{
	// attention je rajoute le *3 pour les donn�es d'aqi(autoquant)
	if (n == nullptr) { return getLikeness(std::min(distanceToBorder(), params_->stdDistance_), params_->stdDistance_); }
	return getLikeness(Distance(n), params_->stdDistance_);
}

/*********************************************************************************************
double Node::likenessDistanceFusion(Node* n1, Node* n2)
  Compute the likeness from the distance of a fusion between this and n1 and n2.
  New version from 2005/01/12
*********************************************************************************************/
double Node::likenessDistanceFusion(const Node* n1, const Node* n2) const
{
	if (n1 == nullptr) { return 0.; }
	if (n2 == nullptr) { return 0.; }

	return getLikeness(fabs(distance(n1, n2) + n1->Distance(n2) - sphericEffect(n1, n2)), params_->stdDistance_);
}

/*********************************************************************************************
double Node::likenessMeanIntensity(Node* n)
  Compute the likeness from the mean of intensity between this and n.
  New version from 2005/01/12

  We use the hypothesis that
	sigma(mean intensity) proportional to surface^-1/6 simga_noise / Mean in 2D.
*********************************************************************************************/
double Node::likenessMeanIntensity(const Node* n) const
{
	if (n == nullptr) { return getLikeness(info_.mean_, params_->stdMean_ * std::max(1., 3. * params_->stdDistance_ / distanceToBorder())); }
	return getLikeness(info_.mean_ - n->info_.mean_, params_->stdMean_);
}

/*********************************************************************************************
double Node::likenessMeanIntensityFusion(Node* n1, Node* n2)
  Compute the likeness from the mean of intensity in fusion between this and n1 and n2 .
  New version from 2005/01/12

*********************************************************************************************/
double Node::likenessMeanIntensityFusion(const Node* n1, const Node* n2) const
{
	if (n1 == nullptr || n2 == nullptr) { return 0.; } //Error ???
	return getLikeness(info_.mean_ - (n1->info_.intensity_ + n2->info_.intensity_)
	                   / (n1->info_.surf_ + n2->info_.surf_), params_->stdMean_);
	//if (Parameters::N_DIM == 2) { return getLikeness(info_.mean_ - (n1->info_.intensity_ + n2->info_.intensity_) / (n1->info_.surf_ + n2->info_.surf_), params_->stdMean_);}
	//else { return getLikeness(info_.mean_ - (n1->info_.intensity_ + n2->info_.intensity_) / (n1->info_.GetVolume() + n2->info_.GetVolume()), params_->stdMean_);}
}

/*********************************************************************************************
double Node::likenessIntensity(Node* n)
  Compute the likeness from the total intensity between this and n.
  New version from 2005/01/12

  We use the hypothesis that sigma(intensity) proportional to simga_noise/Mean.
*********************************************************************************************/
double Node::likenessIntensity(const Node* n) const
{
	if (n == nullptr) { return getLikeness(info_.intensity_, params_->stdIntensity_ * std::max(1., 3. * params_->stdDistance_ / distanceToBorder())); }
	return getLikeness(info_.intensity_ - n->info_.intensity_, params_->stdIntensity_);
}

/*********************************************************************************************
double Node::likenessIntensityFusion(Node* n1, Node* n2
  Compute the likeness for the intensity in a fusion between this and n1 and n2.
  New version from 2005/01/12

  We use the hypothesis that sigma(intensity) proportional to simga_noise/Mean.
*********************************************************************************************/
double Node::likenessIntensityFusion(const Node* n1, const Node* n2) const
{
	if (n1 == nullptr) { return 0.; } //inutil
	if (n2 == nullptr) { return 0.; } //inutil
	return getLikeness(info_.intensity_ - n1->info_.intensity_ - n2->info_.intensity_, params_->stdIntensity_);
}

/*********************************************************************************************
double Node::likeness(Node* n)
  Compute the likeness between this and n.
  New version from 2005/01/12
*********************************************************************************************/
double Node::Likeness(const Node* n) const
{
	// the /2 is because the sum of the likeness of 2 unconnected nodes has to be compared to the likeness of the connection between those two nodes.
	if (n == nullptr) {
		return (WEIGHT_DISTANCE * likenessDistance(n) + WEIGHT_INTENSITY * likenessIntensity(n)
		        + WEIGHT_MEAN_INTENSITY * likenessMeanIntensity(n)) / 2. - params_->costBirthDeath_;
	}
	return WEIGHT_DISTANCE * likenessDistance(n) + WEIGHT_INTENSITY * likenessIntensity(n) + WEIGHT_MEAN_INTENSITY * likenessMeanIntensity(n);
}

/*********************************************************************************************
double Node::likeness(Node* n)
  Compute the likeness between this and n.
  New version from 2005/01/12
*********************************************************************************************/
double Node::likenessBlink(const Node* n) const
{
	constexpr double Alpha = 0.7;
	return Likeness(n) * Alpha + (Likeness(nullptr) + n->Likeness(nullptr)) * (1. - Alpha / 2.);
}

/*********************************************************************************************
double Node::likenessFusion(Node* n1, Node* n2)
  Compute the likeness in fusion between this and n1 and n2.
  New version from 2005/01/12
*********************************************************************************************/
double Node::likenessFusion(const Node* n1, const Node* n2) const
{
	if (n1 == nullptr) { return Likeness(n2); }
	if (n2 == nullptr) { return Likeness(n1); }
	return (WEIGHT_DISTANCE * likenessDistanceFusion(n1, n2) + WEIGHT_INTENSITY * likenessIntensityFusion(n1, n2)
	        + WEIGHT_MEAN_INTENSITY * likenessMeanIntensityFusion(n1, n2) + MIN_ENERGY_FUSION) * 3. / 2.;
	// the *3/2 is because the sum of the likeness of 1 unconnected node and 2 linked nodes has to be compared to the likeness of the fusion connection between those 3 nodes.
}

/*********************************************************************************************
void Node::computeConstantSigmas(double &intensity,
	double &mean, double &distance, long &countStd)
  Find in the vicinity of each node the one with the biggest likeness to compute the constant
  kInt and kMean and kDist
*********************************************************************************************/
void Node::ComputeConstantSigmas(double& intensity, size_t& nInt, double& mean, size_t& nMean, double& distance, size_t& nDist,
                                 const bool clipping, long* hist) const
{
	const Node* nMax = nullptr;
	double likenessMax = -1e10;

	for (const ChainNeighbor* c = neighborsFuture_; c != nullptr; c = c->next_) {
		const double LikenessVal = Likeness(c->neighbor_);
		if (LikenessVal > likenessMax) {
			likenessMax = LikenessVal;
			nMax = c->neighbor_;
		}
	}


	if (nMax != nullptr) {
		if (hist != nullptr) { hist[static_cast<size_t>(Distance(nMax))]++; }
		else {
			if (clipping) {
				if (fabs(info_.intensity_ - nMax->info_.intensity_) < 3. * params_->stdIntensity_) {
					intensity += (info_.intensity_ - nMax->info_.intensity_) * (info_.intensity_ - nMax->info_.intensity_);
					nInt++;
				}
				if (fabs(info_.mean_ - nMax->info_.mean_) < 3. * params_->stdMean_) {
					mean += (info_.mean_ - nMax->info_.mean_) * (info_.mean_ - nMax->info_.mean_);
					nMean++;
				}

				if (Distance(nMax) < 3. * params_->stdDistance_) {
					distance += Distance(nMax) * Distance(nMax);
					nDist++;
				}
			}
			else {
				intensity += (info_.intensity_ - nMax->info_.intensity_) * (info_.intensity_ - nMax->info_.intensity_);
				mean += (info_.mean_ - nMax->info_.mean_) * (info_.mean_ - nMax->info_.mean_);
				distance += Distance(nMax) * Distance(nMax);
				nInt++;
				nMean++;
				nDist++;
			}
		}
	}
}

//=====================================================
//==================== Make Change ====================
//=====================================================
double Node::MakeChange(const double limit)
{
	//choose an operation between:
	//simple linkage
	//strucutre fusion
	//strucutre dissociation

	if (params_->allowFusionAndFission_) {
		if (future_ != nullptr) {
			switch (params_->RandomModulo(4)) {
				case 0: return makeChangeFusion(limit);
				case 1: return makeChangeFission(limit);
				default: return makeChangeSimpleLinkage(limit);
			}
		}

		if constexpr (Parameters::N_TEMPORAL_WINDOWS > 1 && blinkFuture_ != nullptr) { return makeChangeBlinkLinkage(limit); }
		if (fission1Future_ != nullptr) { return makeChangeFission(limit); }
		if (fusionFuture_ != nullptr) { return makeChangeFusion(limit); }
		if constexpr (Parameters::N_TEMPORAL_WINDOWS > 1 && params_->RandomModulo(2) == 0) { return makeChangeBlinkLinkage(limit); }

		return makeChangeSimpleLinkage(limit);
	}

	if (future_ != nullptr) { return makeChangeSimpleLinkage(limit); }
	if (blinkFuture_ != nullptr) { return makeChangeBlinkLinkage(limit); }
	if constexpr (Parameters::N_TEMPORAL_WINDOWS > 2 && params_->RandomModulo(2) == 0) { return makeChangeBlinkLinkage(limit); }

	return makeChangeSimpleLinkage(limit);
}

double Node::makeChangeFission0(const double limit)
{
	//no fission in old configuration
	//try to find 2 structures to make fission:
	//From
	//1(this)-------2
	//4-------------3(cRand)
	//To
	//1(this)-------2
	//        \-----3(cRand)
	//4-------------NULL

	if (blinkFuture_ != nullptr) { return -1; }   //ERROR
	if (neighborsFutureCount_ < 2) { return 0.; } //no neighbors like 3

	ChainNeighbor* cRand;
	size_t i;
	const size_t IdRandNeighbor = params_->RandomModulo(neighborsFutureCount_);
	//set the node cRand
	for (i = 0, cRand = neighborsFuture_; i < IdRandNeighbor && cRand != nullptr; cRand = cRand->next_, ++i) { }
	//test if the node is good
	if (cRand == nullptr) { return -1.; }                         //ERROR 
	if (cRand->neighbor_ == future_) { return 0.; }               // nothing to do, the new random link is the same as future_.
	if (cRand->neighbor_->fusion1Past_ != nullptr) { return 0.; } // crand is already implicated in a fusion.
	if (cRand->neighbor_->fissionPast_ != nullptr) { return 0.; } // crand is already implicated in a fission.
	if (cRand->neighbor_->blinkPast_ != nullptr) { return 0.; }   // crand is already implicated in a past blink.

	if (Distance(cRand->neighbor_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = likenessFusion(cRand->neighbor_, future_)
	                      + (cRand->neighbor_->past_ == nullptr ? 0. : cRand->neighbor_->past_->Likeness(nullptr))
	                      - Likeness(future_) - cRand->neighbor_->Likeness(cRand->neighbor_->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness < limit // if the new solution is too bad
	// change the configuration
	fission1Future_ = future_;
	fission2Future_ = cRand->neighbor_;
	future_->fissionPast_ = this;
	cRand->neighbor_->fissionPast_ = this;
	future_->past_ = nullptr;
	future_ = nullptr;

	if (cRand->neighbor_->past_ != nullptr) { cRand->neighbor_->past_->future_ = nullptr; }
	cRand->neighbor_->past_ = nullptr;

	return DeltaE;
}

double Node::makeChangeFission1(const double limit, Node* node)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5
	//possiblilty 1
	//To
	//1(this)---------2(fission1Future_)
	//4(cRand)--------3(fission2Future_)
	//NULL------------5

	if (node->Distance(fission2Future_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected
	if (node->fission1Future_ != nullptr) { return 0; }                       // cRand is already implicated in a fission.
	if (node->fusionFuture_ != nullptr) { return 0; }                         // cRand is already implicated in a fusion.
	if (node->blinkFuture_ != nullptr) { return 0.; }                         // crand is already implicated in a future blink.

	const double DeltaE = Likeness(fission1Future_) + node->Likeness(fission2Future_)
	                      + (node->future_ == nullptr ? 0 : node->future_->Likeness(nullptr))
	                      - likenessFusion(fission1Future_, fission2Future_) - node->Likeness(node->future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	future_ = fission1Future_;
	future_->past_ = this;
	if (node->future_ != nullptr) { node->future_->past_ = nullptr; }
	node->future_ = fission2Future_;
	fission2Future_->past_ = node;
	fission1Future_->fissionPast_ = nullptr;
	fission2Future_->fissionPast_ = nullptr;
	fission1Future_ = nullptr;
	fission2Future_ = nullptr;
	return DeltaE;
}

double Node::makeChangeFission2(const double limit, Node* node)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5
	//possiblilty 2
	//To
	//1---------------3(fission2Future_)
	//4(cRand)--------2(fission1Future_)
	//NULL------------5

	if (node->Distance(fission1Future_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	if (node->fission2Future_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (node->fusionFuture_ != nullptr) { return 0; }   // cRand is already implicated in a fusion.
	if (node->blinkFuture_ != nullptr) { return 0.; }   // crand is already implicated in a future blink.

	const double DeltaE = Likeness(fission2Future_) + node->Likeness(fission1Future_)
	                      + (node->future_ == nullptr ? 0 : node->future_->Likeness(nullptr))
	                      - likenessFusion(fission2Future_, fission1Future_) - node->Likeness(node->future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	future_ = fission2Future_;
	future_->past_ = this;
	if (node->future_ != nullptr) { node->future_->past_ = nullptr; }
	node->future_ = fission1Future_;
	fission1Future_->past_ = node;
	fission2Future_->fissionPast_ = nullptr;
	fission1Future_->fissionPast_ = nullptr;
	fission2Future_ = nullptr;
	fission1Future_ = nullptr;

	return DeltaE;
}

double Node::makeChangeFission3(const double limit)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)

	//possiblilty 3
	//To
	//1(this)---------2(fission1Future_)
	//NULL------------3(fission2Future_)

	const double DeltaE = Likeness(fission1Future_) + fission2Future_->Likeness(nullptr) - likenessFusion(fission1Future_, fission2Future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration

	future_ = fission1Future_;
	future_->past_ = this;
	future_->fissionPast_ = nullptr;
	fission2Future_->fissionPast_ = nullptr;
	fission1Future_ = nullptr;
	fission2Future_ = nullptr;

	return DeltaE;
}

double Node::makeChangeFission4(const double limit)
{
	//possiblilty 4
	//To
	//1(this)---------3(fission2Future_)
	//NULL------------2(fission1Future_)
	const double DeltaE = Likeness(fission2Future_) + fission1Future_->Likeness(nullptr) - likenessFusion(fission1Future_, fission2Future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	future_ = fission2Future_;
	future_->past_ = this;
	future_->fissionPast_ = nullptr;
	fission1Future_->fissionPast_ = nullptr;
	fission1Future_ = nullptr;
	fission2Future_ = nullptr;
	return DeltaE;
}

double Node::makeChangeFission5(const double limit, Node* node)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5

	//possibility 5
	//1(this)---------2(fission1Future_)
	//           /----3(fission2Future_)
	//4(cRand)--/-----5

	if (node->Distance(fission2Future_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected
	if (node->blinkFuture_ != nullptr) { return 0.; }                         // crand is already implicated in a future blink.

	const double DeltaE = Likeness(fission1Future_) + node->likenessFusion(fission2Future_, node->future_) -
	                      likenessFusion(fission1Future_, fission2Future_) - node->Likeness(node->future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration

	future_ = fission1Future_;
	future_->past_ = this;
	node->fission1Future_ = fission2Future_;
	node->fission1Future_->fissionPast_ = node;
	node->fission2Future_ = node->future_;
	node->fission2Future_->fissionPast_ = node;

	fission1Future_->fissionPast_ = nullptr;
	fission1Future_ = nullptr;
	fission2Future_ = nullptr;
	node->future_->past_ = nullptr;
	node->future_ = nullptr;

	return DeltaE;
}

double Node::makeChangeFission8(const double limit, Node* node)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5
	//possibility 8
	//1(this)---------3(fission2Future_)
	//           /----2(fission1Future_)
	//4(cRand)--/-----5

	if (node->Distance(fission1Future_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected
	if (node->blinkFuture_ != nullptr) { return 0.; }                         // crand is already implicated in a future blink.

	const double DeltaE = Likeness(fission2Future_) + node->likenessFusion(fission1Future_, node->future_)
	                      - likenessFusion(fission1Future_, fission2Future_) - node->Likeness(node->future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration

	future_ = fission2Future_;
	future_->past_ = this;
	node->fission2Future_ = fission1Future_;
	node->fission2Future_->fissionPast_ = node;
	node->fission1Future_ = node->future_;
	node->fission1Future_->fissionPast_ = node;

	fission2Future_->fissionPast_ = nullptr;
	fission2Future_ = nullptr;
	fission1Future_ = nullptr;
	node->future_->past_ = nullptr;
	node->future_ = nullptr;

	return DeltaE;
}

double Node::makeChangeFission9(const double limit, Node* node)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5
	//possibility 8
	//1(this)---------5
	//           /----2(fission1Future_)
	//4(cRand)--/-----3(fission2Future_)
	if (node->blinkFuture_ != nullptr) { return 0.; } // crand is already implicated in a future blink.
	if (node->future_ != nullptr) { if (Distance(node->future_) >= params_->maxSpeed_) { return 0.; } } // Distance is not respected
	if (node->Distance(fission1Future_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected
	if (node->Distance(fission2Future_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected
	if (node->fission1Future_ != nullptr) { return 0.; }
	if (node->fusionFuture_ != nullptr) { return 0.; }

	const double DeltaE = Likeness(node->future_) + node->likenessFusion(fission1Future_, fission2Future_)
	                      - likenessFusion(fission1Future_, fission2Future_) - node->Likeness(node->future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration


	node->fission1Future_ = fission1Future_;
	fission1Future_->fissionPast_ = node;
	node->fission2Future_ = fission2Future_;
	fission2Future_->fissionPast_ = node;
	future_ = node->future_;
	if (future_ != nullptr) { future_->past_ = this; }
	fission1Future_ = nullptr;
	fission2Future_ = nullptr;
	node->future_ = nullptr;
	return DeltaE;
}

double Node::makeChangeFission6(const double limit, Node* node)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5(node->future_)
	//To
	//possibility 6
	//1(this)--\------2(fission1Future_)
	//          \-----5
	//4(cRand)--------3(fission2Future_)

	if (node->blinkFuture_ != nullptr) { return 0.; } // crand is already implicated in a future blink.
	if (node->future_ != nullptr) { if (Distance(node->future_) >= params_->maxSpeed_) { return 0.; } } // Distance is not respected
	if (node->Distance(fission2Future_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = likenessFusion(fission1Future_, node->future_) + node->Likeness(fission2Future_)
	                      - likenessFusion(fission1Future_, fission2Future_) - node->Likeness(node->future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration

	fission2Future_->past_ = node;
	node->future_->fissionPast_ = this;
	Node* localN = node->future_;
	node->future_ = fission2Future_;
	fission2Future_->fissionPast_ = nullptr;
	fission2Future_ = localN;
	fission2Future_->past_ = nullptr;

	return DeltaE;
}

double Node::makeChangeFission7(const double limit, Node* node)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5(node->future_)
	//To
	//possibility 7
	//1(this)-\-------3(fission2Future_)
	//          \-----5
	//4(cRand)--------2(fission1Future_)

	if (node->blinkFuture_ != nullptr) { return 0.; } // crand is already implicated in a future blink.
	if (node->future_ != nullptr) { if (Distance(node->future_) >= params_->maxSpeed_) { return 0.; } } // Distance is not respected
	if (node->Distance(fission1Future_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = likenessFusion(fission2Future_, node->future_) + node->Likeness(fission1Future_)
	                      - likenessFusion(fission2Future_, fission1Future_) - node->Likeness(node->future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration

	fission1Future_->past_ = node;
	node->future_->fissionPast_ = this;
	Node* localN = node->future_;
	node->future_ = fission1Future_;
	fission1Future_->fissionPast_ = nullptr;
	fission1Future_ = localN;
	fission1Future_->past_ = nullptr;

	return DeltaE;
}

double Node::makeChangeFission(const double limit)
{
	//is the node already makes a Fission with another structure
	if (future_ == nullptr && fission1Future_ == nullptr) { return 0.; } //nothing connected in the future!
	if (blinkFuture_ != nullptr) { return 0.; }                          // crand is already implicated in a future blink.
	if (fusionFuture_ != nullptr) { return 0.; }                         //already implicated in a Fusion

	if (fission1Future_ == nullptr) {
		//no fission in old configuration
		//try to find 2 structures to make fission:
		//From
		//1(this)-------2
		//4-------------3(cRand)

		//To
		//1(this)-------2
		//        \-----3(cRand)
		//4-------------NULL

		return makeChangeFission0(limit);
	}
	//there is a fission in old configuration
	//Try to break or change the fission:
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5

	//To
	//possiblilty 1
	//1---------------2(fission1Future_)
	//4(cRand)--------3(fission2Future_)
	//NULL------------5
	//or
	//possiblilty 2
	//1---------------3(fission2Future_)
	//4(cRand)--------2(fission1Future_)
	//NULL------------5
	//or 
	//possibility 3
	//1(this)---------2(fission1Future_)
	//NULL------------3(fission2Future_)
	//or
	//possibility 4
	//1(this)---------3(fission2Future_)
	//NULL------------2(fission1Future_)

	//possibility 5
	//1(this)---------2(fission1Future_)
	//           /----3(fission2Future_)
	//4(cRand)--/-----5

	//possibility 6
	//1(this)--\------2(fission1Future_)
	//          \-----5
	//4(cRand)--------3(fission2Future_)

	//possibility 7
	//1(this)-\-------3(fission2Future_)
	//          \-----5
	//4(cRand)--------2(fission1Future_)

	//possibility 8
	//1(this)---------3(fission2Future_)
	//           /----2(fission1Future_)
	//4(cRand)--/-----5

	//possibility 9
	//1(this)---------5
	//           /----2(fission1Future_)
	//4(cRand)--/-----3(fission2Future_)


	switch (params_->RandomModulo(8)) {
		case 0: return makeChangeFission3(limit);
		case 1: return makeChangeFission4(limit);
		default: ;
	}

	if (fission2Future_->neighborsPastCount_ < 2) { //no neighbors like "4"
		switch (params_->RandomModulo(2)) {
			case 0: return makeChangeFission3(limit);
			case 1: return makeChangeFission4(limit);
			default: ;
		}
	}

	//try to find a neighbor_ "4" in the neighborhood of "3"
	const size_t IdRandNeighbor = params_->RandomModulo(fission2Future_->neighborsPastCount_);

	ChainNeighbor* cRand;
	size_t i;
	//set the node cRand
	for (i = 0, cRand = fission2Future_->neighborsPast_; i < IdRandNeighbor && cRand != nullptr; cRand = cRand->next_, ++i) { }
	//test if the node is good

	if (cRand->neighbor_ == this) { return 0; }                   // nothing to do, the new random link is the same as this.
	if (cRand->neighbor_->blinkFuture_ != nullptr) { return 0.; } // crand is already implicated in a future blink.

	if (cRand->neighbor_->future_ != nullptr) {
		switch (params_->RandomModulo(7)) {
			case 0: return makeChangeFission5(limit, cRand->neighbor_);
			case 1: return makeChangeFission6(limit, cRand->neighbor_);
			case 2: return makeChangeFission7(limit, cRand->neighbor_);
			case 3: return makeChangeFission8(limit, cRand->neighbor_);
			default: ;
		}
	}

	switch (params_->RandomModulo(3)) {
		case 0: return makeChangeFission1(limit, cRand->neighbor_);
		case 1: return makeChangeFission2(limit, cRand->neighbor_);
		default: return makeChangeFission9(limit, cRand->neighbor_);
	}
}

double Node::makeChangeFusion0(const double limit)
{
	//no fusion in old configuration
	//try to find 2 structures to make fusion:
	//From
	//1(this)-------2(future_)
	//3(cRand)------4(can be NULL)
	//To
	//1-------------2
	//3(cRand)--/
	//NULL----------4

	if (future_->neighborsPastCount_ < 2) { return 0.; } //no neighbors like 3
	const size_t IdRandNeighbor = params_->RandomModulo(future_->neighborsPastCount_);
	size_t i;
	ChainNeighbor* cRand;
	//set the node cRand
	for (i = 0, cRand = future_->neighborsPast_; i < IdRandNeighbor && cRand != nullptr; cRand = cRand->next_, ++i) { }
	//test if the node is good
	if (cRand == nullptr) { return -1.; }                           //ERROR 
	if (cRand->neighbor_ == this) { return 0; }                     // nothing to do, the new random link is the same as this.
	if (cRand->neighbor_->fission1Future_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (cRand->neighbor_->fusionFuture_ != nullptr) { return 0; }   // cRand is already implicated in a fusion.
	if (cRand->neighbor_->blinkFuture_ != nullptr) { return 0; }    // cRand is already implicated in a blink connection

	if (future_->Distance(cRand->neighbor_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = future_->likenessFusion(cRand->neighbor_, this)
	                      + (cRand->neighbor_->future_ != nullptr ? cRand->neighbor_->future_->Likeness(nullptr) : 0.)
	                      - Likeness(future_) - cRand->neighbor_->Likeness(cRand->neighbor_->future_);

	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	//unlink future_

	fusionFuture_ = future_;
	future_->fusion1Past_ = this;
	cRand->neighbor_->fusionFuture_ = future_;
	future_->fusion2Past_ = cRand->neighbor_;

	if (cRand->neighbor_->future_ != nullptr) { cRand->neighbor_->future_->past_ = nullptr; }

	cRand->neighbor_->future_ = nullptr;
	future_ = nullptr;
	fusionFuture_->past_ = nullptr;

	return DeltaE;
}

double Node::makeChangeFusion1(const double limit, Node* node)
{
	Node* tmp;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { tmp = fusionFuture_->fusion2Past_; }
	else { tmp = fusionFuture_->fusion1Past_; }
	//there is a fusion in old configuration
	//Try to break the fusion:
	//From
	//1(this)---------2(fusionFuture_)
	//3(node)----/
	//5---------------4(cRand)
	//To
	//possiblilty 1
	//1---------------2
	//3---------------4(cRand)
	//5---------------NULL


	if (node->fissionPast_ != nullptr) { return 0; }               // cRand is already implicated in a fission.
	if (node->fusion1Past_ != nullptr) { return 0; }               // cRand is already implicated in a fusion.
	if (node->Distance(tmp) >= params_->maxSpeed_) { return 0.; } // Distance is not respected


	const double DeltaE = Likeness(fusionFuture_) + tmp->Likeness(node)
	                      + (node->past_ == nullptr ? 0 : node->past_->Likeness(nullptr))
	                      - fusionFuture_->likenessFusion(tmp, this) - node->Likeness(node->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	fusionFuture_->past_ = this;
	fusionFuture_->fusion1Past_ = nullptr;
	fusionFuture_->fusion2Past_ = nullptr;
	tmp->future_ = node;
	tmp->fusionFuture_ = nullptr;
	future_ = fusionFuture_;
	fusionFuture_ = nullptr;

	if (node->past_ != nullptr) { node->past_->future_ = nullptr; }
	node->past_ = tmp;

	return DeltaE;
}

double Node::makeChangeFusion2(const double limit, Node* node)
{
	Node* tmp;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { tmp = fusionFuture_->fusion2Past_; }
	else { tmp = fusionFuture_->fusion1Past_; }

	//From
	//1(this)---------2(fusionFuture_)
	//3(node)----/
	//5---------------4(cRand)
	//possiblilty 2
	//To
	//1---------------4(cRand)
	//3(node)---------2(fusionFuture_)
	//5---------------NULL
	//try to find a neighbor_ "4" in the neighborhood of "1"

	if (node->fissionPast_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (node->fusion1Past_ != nullptr) { return 0; } // cRand is already implicated in a fusion.

	if (Distance(node) >= params_->maxSpeed_) { return 0.; } // Distance is not respected


	const double DeltaE = Likeness(node) + tmp->Likeness(fusionFuture_)
	                      + (node->past_ == nullptr ? 0 : node->past_->Likeness(nullptr))
	                      - fusionFuture_->likenessFusion(tmp, this) - node->Likeness(node->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	//unlink future_

	fusionFuture_->past_ = tmp;
	fusionFuture_->fusion1Past_ = nullptr;
	fusionFuture_->fusion2Past_ = nullptr;
	tmp->fusionFuture_ = nullptr;
	tmp->future_ = fusionFuture_;
	fusionFuture_ = nullptr;
	future_ = node;
	if (node->past_ != nullptr) { node->past_->future_ = nullptr; }
	node->past_ = this;

	return DeltaE;
}

double Node::makeChangeFusion3(const double limit)
{
	Node* tmp;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { tmp = fusionFuture_->fusion2Past_; }
	else { tmp = fusionFuture_->fusion1Past_; }
	//From
	//1(this)---------2(fusionFuture_)
	//3(node)----/
	//possiblilty 3
	//To
	//1(this)---------2(fusionFuture_)
	//3(node)---------NULL

	const double DeltaE = Likeness(fusionFuture_) + tmp->Likeness(nullptr) - fusionFuture_->likenessFusion(tmp, this);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration

	fusionFuture_->past_ = this;
	fusionFuture_->fusion1Past_ = nullptr;
	fusionFuture_->fusion2Past_ = nullptr;
	future_ = fusionFuture_;
	fusionFuture_ = nullptr;
	tmp->fusionFuture_ = nullptr;
	if (tmp->future_ != nullptr) { return -1.; } //impossible inutile

	return DeltaE;
}

double Node::makeChangeFusion4(const double limit)
{
	Node* tmp;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { tmp = fusionFuture_->fusion2Past_; }
	else { tmp = fusionFuture_->fusion1Past_; }
	//From
	//1(this)---------2(fusionFuture_)
	//3(node)----/

	//possiblilty 4
	//1(this)---------NULL
	//3(node)---------2(fusionFuture_)
	const double DeltaE = tmp->Likeness(fusionFuture_) + Likeness(nullptr) - fusionFuture_->likenessFusion(tmp, this);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	fusionFuture_->past_ = tmp;
	fusionFuture_->fusion1Past_ = nullptr;
	fusionFuture_->fusion2Past_ = nullptr;
	tmp->future_ = fusionFuture_;
	fusionFuture_ = nullptr;
	tmp->fusionFuture_ = nullptr;

	if (future_ != nullptr) { return -1.; } //impossible inutile
	return DeltaE;
}

double Node::makeChangeFusion5(const double limit, Node* node)
{
	Node* tmp;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { tmp = fusionFuture_->fusion2Past_; }
	else { tmp = fusionFuture_->fusion1Past_; }

	//From
	//1(this)-----/---2(fusionFuture_)
	//3(node)----/
	//5---------------4(cRand)


	//possibility 5
	//1(this)---------2(fusionFuture_)
	//3(node)----\
	//5-----------\---4(cRand)

	if (tmp->Distance(node) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	if (node->fissionPast_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (node->fusion1Past_ != nullptr) { return 0; } // cRand is already implicated in a fusion.

	if (node->past_ == nullptr) { return 0; } // not possible to reconnection the fusion.

	const double DeltaE = Likeness(fusionFuture_) + node->likenessFusion(tmp, node->past_)
	                      - fusionFuture_->likenessFusion(tmp, this) - node->Likeness(node->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration

	future_ = fusionFuture_;
	future_->past_ = this;

	tmp->fusionFuture_ = node;
	node->fusion1Past_ = tmp;

	node->past_->fusionFuture_ = node;
	node->fusion2Past_ = node->past_;

	future_->fusion1Past_ = nullptr;
	fusionFuture_ = nullptr;

	future_->fusion2Past_ = nullptr;

	node->past_->future_ = nullptr;
	node->past_ = nullptr;

	return DeltaE;
}

double Node::makeChangeFusion7(const double limit, Node* node)
{
	Node* tmp;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { tmp = fusionFuture_->fusion2Past_; }
	else { tmp = fusionFuture_->fusion1Past_; }

	//From
	//1(this)-----/---2(fusionFuture_)
	//3(node)----/
	//5---------------4(cRand)


	//possibility 7
	//3(node)---------2(fusionFuture_)
	//1(this)----\
	//5-----------\---4(cRand)

	if (node->fissionPast_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (node->fusion1Past_ != nullptr) { return 0; } // cRand is already implicated in a fusion.

	if (node->past_ == nullptr) { return 0; } // not possible to reconnection the fusion.

	if (Distance(node) >= params_->maxSpeed_) { return 0.; } // Distance is not respected
	const double DeltaE = tmp->Likeness(fusionFuture_) + node->likenessFusion(this, node->past_)
	                      - fusionFuture_->likenessFusion(tmp, this) - node->Likeness(node->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration

	tmp->future_ = tmp->fusionFuture_;
	tmp->future_->past_ = tmp;

	fusionFuture_ = node;
	node->fusion1Past_ = this;

	node->past_->fusionFuture_ = node;
	node->fusion2Past_ = node->past_;

	tmp->future_->fusion1Past_ = nullptr;
	tmp->fusionFuture_ = nullptr;

	tmp->future_->fusion2Past_ = nullptr;

	node->past_->future_ = nullptr;
	node->past_ = nullptr;

	return DeltaE;
}

double Node::makeChangeFusion8(const double limit, Node* node)
{
	Node* tmp;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { tmp = fusionFuture_->fusion2Past_; }
	else { tmp = fusionFuture_->fusion1Past_; }
	//From
	//1(this)---------2(fusionFuture_)
	//3(node)----/
	//5---------------4(cRand)

	//possibility 8
	//5---------------2(fusionFuture_)
	//1(this)----\
	//3(node)-----\---4(cRand)

	if (node->fissionPast_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (node->fusion1Past_ != nullptr) { return 0; } // cRand is already implicated in a fusion.

	if (node->past_ == nullptr) { return 0; } // not possible to reconnection the fusion.

	if (Distance(node) >= params_->maxSpeed_) { return 0.; }                       // Distance is not respected
	if (tmp->Distance(node) >= params_->maxSpeed_) { return 0.; }                 // Distance is not respected
	if (node->past_->Distance(fusionFuture_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = node->past_->Likeness(fusionFuture_) + node->likenessFusion(this, tmp)
	                      - fusionFuture_->likenessFusion(tmp, this) - node->Likeness(node->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	fusionFuture_->past_ = node->past_;
	fusionFuture_->past_->future_ = fusionFuture_;

	node->fusion1Past_ = fusionFuture_->fusion1Past_;
	node->fusion2Past_ = fusionFuture_->fusion2Past_;
	node->past_ = nullptr;
	//fusionFuture_->past_->future_=NULL;
	tmp->fusionFuture_ = node;
	fusionFuture_->fusion1Past_ = nullptr;
	fusionFuture_->fusion2Past_ = nullptr;
	fusionFuture_ = node;

	return DeltaE;
}

double Node::makeChangeFusion6(const double limit, Node* node) const
{
	Node* tmp;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { tmp = fusionFuture_->fusion2Past_; }
	else { tmp = fusionFuture_->fusion1Past_; }

	//From
	//1(this)-----/---2(fusionFuture_)
	//3(node)----/
	//5---------------4(cRand)


	//possibility 6
	//1(this)-----/---2(fusionFuture_)
	//5----------/
	//3(node)---------4(cRand)

	if (tmp->Distance(node) >= params_->maxSpeed_) { return 0.; } // Distance is not respected
	if (node->past_ != nullptr) {
		if (node->past_->Distance(fusionFuture_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected
	}


	const double DeltaE = fusionFuture_->likenessFusion(this, node->past_) + node->Likeness(tmp)
	                      - fusionFuture_->likenessFusion(tmp, this) - node->Likeness(node->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration

	if (fusionFuture_->fusion1Past_ == this) {
		node->past_->fusionFuture_ = fusionFuture_;
		node->past_->future_ = nullptr;
		fusionFuture_->fusion2Past_ = node->past_;
		tmp->fusionFuture_ = nullptr;
		tmp->future_ = node;
		node->past_ = tmp;
	}
	else {
		node->past_->fusionFuture_ = fusionFuture_;
		node->past_->future_ = nullptr;
		fusionFuture_->fusion1Past_ = node->past_;
		tmp->fusionFuture_ = nullptr;
		tmp->future_ = node;
		node->past_ = tmp;
	}

	return DeltaE;
}

double Node::makeChangeFusion9(const double limit, Node* node)
{
	Node* tmp;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { tmp = fusionFuture_->fusion2Past_; }
	else { tmp = fusionFuture_->fusion1Past_; }
	Node* tmp2;
	//find "node2"
	if (node->fissionPast_->fission1Future_ == node) { tmp2 = node->fissionPast_->fission2Future_; }
	else { tmp2 = node->fissionPast_->fission1Future_; }

	//possibility 9
	//From 
	//1(this)-----/---2(fusionFuture_)
	//3(node)----/
	//5---------\-----4(cRand)
	//           \----6(node2)
	//To
	//1(this)---------2(fusionFuture_)
	//3(node)---------4(cRand)
	//5---------------6(node2)

	if (tmp->Distance(node) >= params_->maxSpeed_) { return 0.; } // Distance is not respected


	const double DeltaE = Likeness(fusionFuture_) + tmp->Likeness(node) + node->fissionPast_->Likeness(tmp2)
	                      - fusionFuture_->likenessFusion(tmp, this) - node->fissionPast_->likenessFusion(node, tmp2);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit if the new solution is too bad change the configuration

	fusionFuture_->past_ = this;
	future_ = fusionFuture_;
	fusionFuture_->fusion1Past_ = nullptr;
	fusionFuture_->fusion2Past_ = nullptr;
	tmp->fusionFuture_ = nullptr;
	fusionFuture_ = nullptr;
	tmp->future_ = node;
	node->past_ = tmp;
	node->fissionPast_ = nullptr;
	tmp2->fissionPast_->fission1Future_ = nullptr;
	tmp2->fissionPast_->fission2Future_ = nullptr;
	tmp2->fissionPast_->future_ = tmp2;
	tmp2->past_ = tmp2->fissionPast_;
	tmp2->fissionPast_ = nullptr;

	return DeltaE;
}

double Node::makeChangeFusion(const double limit)
{
	//is the node already making a Fusion with another structure
	if (future_ == nullptr && fusionFuture_ == nullptr) { return 0.; } //nothing connected in the future!

	if (fission1Future_ != nullptr) { return 0.; } //already implicated in a Fission

	if (fusionFuture_ == nullptr) {
		//no fusion in old configuration
		//try to find 2 structures to make fusion:
		//From
		//1(this)-------2
		//3-------------4(cRand)

		//To
		//1-------------2
		//3-------/
		//NULL----------4(cRand)

		return makeChangeFusion0(limit);
	}
	//there is a fusion in old configuration
	//Try to break the fusion:
	//From
	//1(this)---------2(fusionFuture_)
	//3(node)----/
	//5---------------4(cRand)

	//To
	//possiblilty 1
	//1---------------2
	//3---------------4(cRand)
	//5---------------NULL
	//or
	//possiblilty 2
	//1---------------4(cRand)
	//3(node)---------2(fusionFuture_)
	//5---------------NULL
	//or 
	//possibility 3
	//1(this)---------2(fusionFuture_)
	//3(node)---------NULL
	//or
	//possibility 4
	//1(this)---------NULL
	//3(node)---------2(fusionFuture_)

	//possibility 5
	//1(this)---------2(fusionFuture_)
	//3(node)----\
		//5-----------\---4(cRand)

	//possibility 6
	//1(this)----/----2(fusionFuture_)
	//5---------/
	//3(node)---------4(cRand)

	//possibility 7
	//3(node)---------2(fusionFuture_)
	//1(this)----\
		//5-----------\---4(cRand)

	//possibility 8
	//5---------------2(fusionFuture_)
	//1(this)----\
		//3(node)-----\---4(cRand)

	//possibility 9
	//From 
	//1(this)---------2(fusionFuture_)
	//3(node)----/
	//5---------------4(cRand)
	//           \----6
	//To
	//1(this)---------2(fusionFuture_)
	//3(node)---------4(cRand)
	//5---------------6


	//just unlink the fusion

	switch (params_->RandomModulo(8)) {
		case 0: return makeChangeFusion3(limit);
		case 1: return makeChangeFusion4(limit);
		default: ;
	}


	ChainNeighbor* cRand;
	size_t i;

	//get the cRand
	if (neighborsFutureCount_ < 2) { //no neighbors like "4", no cRand
		switch (params_->RandomModulo(2)) {
			case 0: return makeChangeFusion3(limit);
			case 1: return makeChangeFusion4(limit);
			default: ;
		}
	}

	const size_t IdRandNeighbor = params_->RandomModulo(neighborsFutureCount_);
	//set the node cRand
	for (i = 0, cRand = neighborsFuture_; i < IdRandNeighbor && cRand != nullptr; cRand = cRand->next_, ++i) { }
	//test if the node is good
	if (cRand == nullptr) { return -1.; } //ERROR

	if (cRand->neighbor_->fissionPast_ != nullptr) { return makeChangeFusion9(limit, cRand->neighbor_); }
	if (cRand->neighbor_->blinkPast_ != nullptr) { return 0; } // cRand is already implicated in a blink connection

	if (cRand->neighbor_->past_ != nullptr) {
		switch (params_->RandomModulo(7)) {
			case 0: return makeChangeFusion5(limit, cRand->neighbor_);
			case 1:
			case 2: return makeChangeFusion6(limit, cRand->neighbor_);
			case 3: return makeChangeFusion7(limit, cRand->neighbor_);
			case 4: return makeChangeFusion8(limit, cRand->neighbor_);
			default: ;
		}
	}

	switch (params_->RandomModulo(2)) {
		case 0: return makeChangeFusion1(limit, cRand->neighbor_);
		default: return makeChangeFusion2(limit, cRand->neighbor_);
	}
}

double Node::makeChangeSimpleLinkage0(const double limit) const
{ //simple deconnection

	//From 
	//Present	Future
	//1(this)-----2(future_)

	//To
	//1(this)-----NULL
	//NULL--------2(future_)

	const double DeltaE = Likeness(nullptr) + future_->Likeness(nullptr) - Likeness(future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	return 0.;
}

double Node::makeChangeSimpleLinkage1(const double limit, Node* node)
{
	//	From
	//	Present		Future
	//	1(this)-----2
	//	NULL--------3(cRand)

	//	To
	//	Present		Future
	//	1(this)-----3(cRand)
	//	NULL--------2

	if (node->blinkPast_ != nullptr) { return 0; } // cRand is already implicated in a blink connection

	if (Distance(node) >= params_->maxSpeed_) { return 0.; } // Distance is not respected


	const double DeltaE = Likeness(node) + (future_ != nullptr ? future_->Likeness(nullptr) : 0)
	                      - node->Likeness(nullptr) - Likeness(future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	//unlink future_
	if (future_ != nullptr) { future_->past_ = nullptr; }
	future_ = node;
	future_->past_ = this;

	return DeltaE;
	//return(0.);
}

double Node::makeChangeSimpleLinkage2(const double limit, Node* node)
{
	// try to break the link
	//	From
	//	Present		Future
	//	1(this)-----2
	//	4-----------3(cRand)
	//	To
	//	Present		Future
	//	1-----------3(cRand)
	//	NULL--------2
	//	4-----------NULL

	if (Distance(node) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = Likeness(node) + (future_ == nullptr ? 0 : future_->Likeness(nullptr))
	                      + (node->past_ == nullptr ? 0 : node->past_->Likeness(nullptr))
	                      - Likeness(future_) - node->Likeness(node->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	//unlink future_


	if (future_ != nullptr) { future_->past_ = nullptr; }
	future_ = node;
	future_->past_->future_ = nullptr;
	future_->past_ = this;

	return DeltaE;
	//return(0.);
}

double Node::makeChangeSimpleLinkage3(const double limit, Node* node, const size_t  count)
{
	//	From
	//	Present		Future
	//	1-----------2
	//	4-----------3(cRand)
	//	NULL--------5(c)

	//find a random one
	const size_t IdRandNeighbor = params_->RandomModulo(count);
	size_t  iCount = 0;

	const ChainNeighbor* c = node->past_->neighborsFuture_;
	do {
		if (c->neighbor_->fusion1Past_ == nullptr && c->neighbor_->fissionPast_ == nullptr) {
			if (c->neighbor_->past_ == nullptr) { iCount++; }
			else { if (c->neighbor_ == future_) { iCount++; } }
		}
		if (iCount <= IdRandNeighbor) { c = c->next_; }
	}
	while (iCount <= IdRandNeighbor);


	//	To
	//	Present		Future
	//	1-----------3(cRand)
	//	4-----------5(c)
	//	NULL--------2
	if (c == nullptr) { return -1; } //ERROR 

	if (c->neighbor_->fusion1Past_ != nullptr) { return -1; }                           //ERROR 
	if (c->neighbor_->fissionPast_ != nullptr) { return -1; }                           //ERROR 
	if (c->neighbor_ != future_) { if (c->neighbor_->past_ != nullptr) { return -1; } } //ERROR 	

	if (Distance(node) >= params_->maxSpeed_) { return 0.; } // Distance is not respected
	if (node->blinkPast_ != nullptr) { return 0; }           // cRand is already implicated in a blink connection
	if (c->neighbor_->blinkPast_ != nullptr) { return 0; }               // c is already implicated in a blink connection
	if (node->past_ != nullptr) {
		if (node->past_->Distance(c->neighbor_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected
	}

	const double DeltaE = Likeness(node) + node->past_->Likeness(c->neighbor_) + (future_ == nullptr ? 0 : future_->Likeness(nullptr))
	                      - Likeness(future_) - node->Likeness(node->past_) - c->neighbor_->Likeness(nullptr);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	//unlink future_

	if (future_ != nullptr) {
		future_->past_ = nullptr; //????
	}
	future_ = node;
	future_->past_->future_ = c->neighbor_;
	future_->past_->future_->past_ = future_->past_;
	future_->past_ = this;

	return DeltaE;
	//return(0.);
}

double Node::makeChangeSimpleLinkage4(const double limit, Node* node)
{
	//	From
	//	Present		Future
	//	1(this)-----2
	//	4-----------3(cRand)
	//	To
	//	Present		Future
	//	1-----------3(cRand)
	//	4-----------2


	if (node->past_ != nullptr) {
		if (node->past_->Distance(future_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected
	}
	if (Distance(node) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = Likeness(node) + future_->Likeness(node->past_)
	                      - node->Likeness(node->past_) - Likeness(future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration

	node->past_->future_ = future_;
	future_->past_ = node->past_;
	future_ = node;
	node->past_ = this;

	return DeltaE;
	//return(0.);
}

double Node::makeChangeBlinkLinkage0(const double limit)
{
	//transform : 
	//	From
	//	t			t+1		t+2
	//	1(this)-------------NULL
	//	NULL----------------2


	//	To
	//	t			t+1		t+2
	//	1(this)-------------2(blinkFuture_)

	//const size_t Dt = params_->RandomModulo((Parameters::N_TEMPORAL_WINDOWS - 1) + 2); Temporal = 1 donc c'est modulo 0
	constexpr size_t Dt = 2;
	size_t nbNodeAtGoodTime = 0;
	ChainNeighbor* cRand;

	for (cRand = neighborsBlinkFuture_; cRand != nullptr; cRand = cRand->next_) {
		if (cRand->neighbor_->time_ - time_ == Dt && cRand->neighbor_->past_ == nullptr && cRand->neighbor_->blinkPast_ == nullptr && cRand->neighbor_->
		    fissionPast_
		    == nullptr && cRand->neighbor_->fusion1Past_ == nullptr) { nbNodeAtGoodTime++; }
	}
	if (nbNodeAtGoodTime == 0) { return 0; }

	const size_t IdRandNeighbor = params_->RandomModulo(nbNodeAtGoodTime);
	size_t i = 0;
	bool exitFor = false;
	for (cRand = neighborsBlinkFuture_; !exitFor && cRand != nullptr;) {
		if (cRand->neighbor_->time_ - time_ == Dt && cRand->neighbor_->past_ == nullptr && cRand->neighbor_->blinkPast_ == nullptr && cRand->neighbor_->
		    fissionPast_
		    == nullptr && cRand->neighbor_->fusion1Past_ == nullptr) {
			if (i == IdRandNeighbor) { exitFor = true; }
			i++;
		}
		if (!exitFor) { cRand = cRand->next_; }
	}


	if (cRand == nullptr) { return -1.; } //ERROR 

	const double DeltaE = likenessBlink(cRand->neighbor_) - cRand->neighbor_->Likeness(nullptr) - Likeness(nullptr);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	cRand->neighbor_->blinkPast_ = this;
	blinkFuture_ = cRand->neighbor_;

	return DeltaE;
}

double Node::makeChangeBlinkLinkage1(const double limit)
{
	//transform : 
	//	From
	//	t			t+1		t+2
	//	1(this)-------------2(blinkFuture_)

	//	To
	//	t			t+1		t+2
	//	1(this)-------------NULL
	//	NULL----------------2
	const double DeltaE = blinkFuture_->Likeness(nullptr) + Likeness(nullptr) - likenessBlink(blinkFuture_);
	if (!metropolis(DeltaE, limit)) { return 0.; } // deltaLikeness<limit // if the new solution is too bad
	// change the configuration


	blinkFuture_->blinkPast_ = nullptr;
	blinkFuture_ = nullptr;

	return DeltaE;
}

double Node::makeChangeBlinkLinkage2(const double limit, Node* node)
{
	//transform : 
	//	From
	//	t			t+1					t+n
	//	1(this)-------------------------2(blinkFuture_)
	//	4---------3(cRand)----...
	//  1=this
	//  2=blinkFuture_
	//	3=node
	//	4=node->past_


	//	To
	//	t			t+1					t+n
	//	4-------------------------------2(blinkFuture_)
	//	1(this)---3(cRand)----...


	const double DeltaE = Likeness(node) + node->past_->likenessBlink(blinkFuture_) - likenessBlink(blinkFuture_)
	                      - node->past_->Likeness(node);
	if (!metropolis(DeltaE, limit)) { return 0.; } // deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	blinkFuture_->blinkPast_ = node->past_;
	node->past_->blinkFuture_ = blinkFuture_;
	node->past_->future_ = nullptr;
	node->past_ = this;
	future_ = node;
	blinkFuture_ = nullptr;

	return DeltaE;
}

double Node::makeChangeBlinkLinkage3(const double limit, Node* node)
{
	//transform : 
	//	From
	//	t					 t+n-1		t+n
	//	1(this)-------------------------2(blinkFuture_)
	//				...------3(cRand)---4
	//  1=this
	//  2=blinkFuture_
	//	3=node
	//	4=node->future_


	//	To
	//	t					 t+n-1		t+n
	//	1(this)-------------------------4
	//				...------3(cRand)---2(blinkFuture_)
	//  makeChangeBlinkLinkage3

	const double DeltaE = likenessBlink(node->future_) + node->Likeness(blinkFuture_)
	                      - likenessBlink(blinkFuture_) - node->future_->Likeness(node);
	if (!metropolis(DeltaE, limit)) { return 0.; } // deltaLikeness<limit // if the new solution is too bad
	// change the configuration

	node->future_->blinkPast_ = this;
	node->future_->past_ = nullptr;
	blinkFuture_->blinkPast_ = nullptr;
	blinkFuture_->past_ = node;
	Node* tmp = node->future_;
	node->future_ = blinkFuture_;
	blinkFuture_ = tmp;

	return DeltaE;
}

double Node::makeChangeBlinkLinkage4(const double limit, Node* node)
{
	//transform : 
	//	From
	//	t			t+1					t+n
	//	1(this)-------------------------2(blinkFuture_)
	//	NULL--------3(cRand2)----...
	//  1=this
	//  2=blinkFuture_
	//	3=cRand2->neighbor_


	//	To
	//	t			t+1					t+n
	//	NULL----------------------------2(blinkFuture_)
	//	1(this)---3(cRand2)----...
	//  makeChangeBlinkLinkage4

	const double DeltaE = Likeness(node) + blinkFuture_->Likeness(nullptr) - likenessBlink(blinkFuture_) - node->Likeness(nullptr);
	if (!metropolis(DeltaE, limit)) { return 0.; } // deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	blinkFuture_->blinkPast_ = nullptr;
	blinkFuture_ = nullptr;
	node->past_ = this;
	future_ = node;

	return DeltaE;
}

double Node::makeChangeBlinkLinkage5(const double limit, Node* node)
{
	//transform : 
	//	From
	//	t					 t+n-1		t+n
	//	1(this)-------------------------2(blinkFuture_)
	//				...------3(cRand2)---NULL
	//  1=this
	//  2=blinkFuture_
	//	3=cRand2->neighbor_


	//	To
	//	t					 t+n-1		t+n
	//	1(this)-------------------------NULL
	//				...------3(cRand2)---2(blinkFuture_)
	//  makeChangeBlinkLinkage5


	const double DeltaE = Likeness(nullptr) + node->Likeness(blinkFuture_) - likenessBlink(blinkFuture_) - node->Likeness(nullptr);
	if (!metropolis(DeltaE, limit)) { return 0.; } // deltaLikeness<limit // if the new solution is too bad
	// change the configuration


	blinkFuture_->blinkPast_ = nullptr;
	blinkFuture_->past_ = node;
	node->future_ = blinkFuture_;
	blinkFuture_ = nullptr;
	return DeltaE;
}

double Node::makeChangeBlinkLinkage(const double limit)
{
	size_t i;
	if (blinkFuture_ == nullptr) {
		//nothing connected in the future
		//transform : 
		//	From
		//	t			t+1		t+n
		//	1(this)-------------NULL
		//	NULL----------------2
		//	To
		//	t			t+1		t+n
		//	1(this)-------------2(blinkFuture_)
		return makeChangeBlinkLinkage0(limit);
	}
	//blinkFuture_!=NULL
	//transform : 
	//	From
	//	t			t+1		t+n
	//	1(this)-------------2(blinkFuture_)

	//	To
	//	t			t+1		t+n
	//	1(this)-------------NULL
	//	NULL----------------2
	//  makeChangeBlinkLinkage1

	//transform : 
	//	From
	//	t			t+1					t+n
	//	1(this)-------------------------2(blinkFuture_)
	//	4---------3(cRand)----...
	//  1=this
	//  2=blinkFuture_
	//	3=cRand->neighbor_
	//	4=cRand->neighbor_->past_


	//	To
	//	t			t+1					t+n
	//	4-------------------------------2(blinkFuture_)
	//	1(this)---3(cRand)----...
	//  makeChangeBlinkLinkage2

	//OR

	//transform : 
	//	From
	//	t					 t+n-1		t+n
	//	1(this)-------------------------2(blinkFuture_)
	//				...------3(cRand)---4
	//  1=this
	//  2=blinkFuture_
	//	3=cRand->neighbor_
	//	4=cRand->neighbor_->future_


	//	To
	//	t					 t+n-1		t+n
	//	1(this)-------------------------4
	//				...------3(cRand)---2(blinkFuture_)
	//  makeChangeBlinkLinkage3


	//OR

	//transform : 
	//	From
	//	t			t+1					t+n
	//	1(this)-------------------------2(blinkFuture_)
	//	NULL--------3(cRand2)----...
	//  1=this
	//  2=blinkFuture_
	//	3=cRand2->neighbor_


	//	To
	//	t			t+1					t+n
	//	NULL----------------------------2(blinkFuture_)
	//	1(this)---3(cRand2)----...
	//  makeChangeBlinkLinkage4

	//OR

	//transform : 
	//	From
	//	t					 t+n-1		t+n
	//	1(this)-------------------------2(blinkFuture_)
	//				...------3(cRand2)---NULL
	//  1=this
	//  2=blinkFuture_
	//	3=cRand2->neighbor_


	//	To
	//	t					 t+n-1		t+n
	//	1(this)-------------------------NULL
	//				...------3(cRand2)---2(blinkFuture_)
	//  makeChangeBlinkLinkage5

	if (params_->RandomModulo(2)) { return makeChangeBlinkLinkage1(limit); }

	//find a potential cRand
	//cRand has to be either just after this or just before blinkFuture_
	if (params_->RandomModulo(2)) { //cRand just after this
		size_t nbNodeAtGoodTime = 0, nbNodeAtGoodTime2 = 0, idRandNeighbor;
		ChainNeighbor* cRand;

		for (cRand = neighborsFuture_; cRand != nullptr; cRand = cRand->next_) {
			if (cRand->neighbor_->time_ - time_ == 1) {
				if (cRand->neighbor_->past_ != nullptr) { nbNodeAtGoodTime++; }
				else {
					if (cRand->neighbor_->fusion1Past_ == nullptr && cRand->neighbor_->fissionPast_ == nullptr && cRand->neighbor_->blinkPast_ == nullptr) {
						nbNodeAtGoodTime2++;
					}
				}
			}
		}
		if (nbNodeAtGoodTime == 0 && nbNodeAtGoodTime2 == 0) { return makeChangeBlinkLinkage1(limit); }

		if ((params_->RandomModulo(2) && nbNodeAtGoodTime != 0) || nbNodeAtGoodTime2 == 0) {
			idRandNeighbor = params_->RandomModulo(nbNodeAtGoodTime);
			i = 0;
			bool exitFor = false;
			for (cRand = neighborsFuture_; !exitFor && cRand != nullptr;) {
				if (cRand->neighbor_->time_ - time_ == 1 && cRand->neighbor_->past_ != nullptr) {
					if (i == idRandNeighbor) { exitFor = true; }
					i++;
				}
				if (!exitFor) { cRand = cRand->next_; }
			}
			if (cRand == nullptr) { return -1.; } //ERROR 

			if (Distance(cRand->neighbor_) >= params_->maxSpeed_) { return 0; } //cRand not linked to this!
			if (blinkFuture_->Distance(cRand->neighbor_->past_) >= params_->maxSpeed_) { return 0; }
			//blinkFuture_ not linked to cRand->neighbor_->past_
			return makeChangeBlinkLinkage2(limit, cRand->neighbor_);
		}
		ChainNeighbor* cRand2;
		idRandNeighbor = params_->RandomModulo(nbNodeAtGoodTime2);
		i = 0;
		bool exitFor = false;
		for (cRand2 = neighborsFuture_; !exitFor && cRand2 != nullptr;) {
			if (cRand2->neighbor_->time_ - time_ == 1 && cRand2->neighbor_->past_ == nullptr && cRand2->neighbor_->fusion1Past_ == nullptr && cRand2->
			    neighbor_->fissionPast_ == nullptr && cRand2->neighbor_->blinkPast_ == nullptr) {
				if (i == idRandNeighbor) { exitFor = true; }
				i++;
			}
			if (!exitFor) { cRand2 = cRand2->next_; }
		}
		if (cRand2 == nullptr) { return -1.; } //ERROR 

		if (Distance(cRand2->neighbor_) >= params_->maxSpeed_) { return 0; } //cRand2 not linked to this!

		return makeChangeBlinkLinkage4(limit, cRand2->neighbor_);
	}
	// cRand is just before blinkFuture_
	size_t nbNodeAtGoodTime = 0, nbNodeAtGoodTime2 = 0, idRandNeighbor;
	ChainNeighbor* cRand;

	for (cRand = blinkFuture_->neighborsPast_; cRand != nullptr; cRand = cRand->next_) {
		if (blinkFuture_->time_ - cRand->neighbor_->time_ == 1) {
			if (cRand->neighbor_->future_ != nullptr) { nbNodeAtGoodTime++; }
			else {
				if (cRand->neighbor_->fission1Future_ == nullptr && cRand->neighbor_->fusionFuture_ == nullptr && cRand->neighbor_->blinkFuture_ ==
				    nullptr) { nbNodeAtGoodTime2++; }
			}
		}
	}
	if (nbNodeAtGoodTime == 0 && nbNodeAtGoodTime2 == 0) { return makeChangeBlinkLinkage1(limit); }

	if ((params_->RandomModulo(2) && nbNodeAtGoodTime != 0) || nbNodeAtGoodTime2 == 0) {
		idRandNeighbor = params_->RandomModulo(nbNodeAtGoodTime);
		i = 0;
		bool exitFor = false;
		for (cRand = blinkFuture_->neighborsPast_; !exitFor && cRand != nullptr;) {
			if (blinkFuture_->time_ - cRand->neighbor_->time_ == 1 && cRand->neighbor_->future_ != nullptr) {
				if (i == idRandNeighbor) { exitFor = true; }
				i++;
			}
			if (!exitFor) { cRand = cRand->next_; }
		}
		if (cRand == nullptr) { return -1.; }                                             //ERROR 
		if (blinkFuture_->Distance(cRand->neighbor_) >= params_->maxSpeed_) { return 0; } //cRand not linked to blinkFuture_!
		if (Distance(cRand->neighbor_->future_) >= params_->maxSpeed_) { return 0; }      //this is not linked to cRand->neighbor_->future_!


		return makeChangeBlinkLinkage3(limit, cRand->neighbor_);
	}
	ChainNeighbor* cRand2;
	idRandNeighbor = params_->RandomModulo(nbNodeAtGoodTime2);
	i = 0;
	bool exitFor = false;
	for (cRand2 = blinkFuture_->neighborsPast_; !exitFor && cRand2 != nullptr;) {
		if (blinkFuture_->time_ - cRand2->neighbor_->time_ == 1 && cRand2->neighbor_->future_ == nullptr && cRand2->neighbor_->fission1Future_ ==
		    nullptr
		    && cRand2->neighbor_->fusionFuture_ == nullptr && cRand2->neighbor_->blinkFuture_ == nullptr) {
			if (i == idRandNeighbor) { exitFor = true; }
			i++;
		}
		if (!exitFor) { cRand2 = cRand2->next_; }
	}
	if (cRand2 == nullptr) { return -1.; }                                             //ERROR 
	if (blinkFuture_->Distance(cRand2->neighbor_) >= params_->maxSpeed_) { return 0; } //cRand2 not linked to blinkFuture_!
	return makeChangeBlinkLinkage5(limit, cRand2->neighbor_);
}

double Node::makeChangeSimpleLinkage(const double limit)
{
	if (future_ == nullptr || params_->RandomModulo(2) == 0) {
		//cRand is a random node in the vicinity of "this" in the future.

		//compute index of cRand
		if (neighborsFutureCount_ == 0) { return 0.; } //no neighbors
		const size_t IdRandNeighbor = params_->RandomModulo(neighborsFutureCount_);

		//set the node cRand
		size_t i;
		ChainNeighbor* cRand;
		for (i = 0, cRand = neighborsFuture_; i < IdRandNeighbor && cRand != nullptr; cRand = cRand->next_, ++i) { }
		//test if the node is good
		if (cRand == nullptr) { return -1.; } //ERROR 
		// the new random link is the same as future_.
		if (cRand->neighbor_ == future_) { return makeChangeSimpleLinkage0(limit); } //try simple deconnection


		if (cRand->neighbor_->fissionPast_ != nullptr) { return 0; } // cRand already implicated in a fission
		if (cRand->neighbor_->fusion1Past_ != nullptr) { return 0; } // cRand already implicated in a fusion
		if (cRand->neighbor_->blinkPast_ != nullptr) { return 0; }   // cRand is already implicated in a blink connection


		//is there a conflict between this and the pastConnextion of cRand->neighbor_?
		if (cRand->neighbor_->past_ == nullptr) { //no conflict cRand->neighbor_ can be linked with this.
			//	From
			//	Present		Future
			//	1(this)-----2
			//	NULL--------3(cRand)

			//	To
			//	Present		Future
			//	1(this)-----3(cRand)
			//	NULL--------2

			return makeChangeSimpleLinkage1(limit, cRand->neighbor_);
		}
		//there is a conflict, the node cRand->neighbor_ is already connected to a node at the same time_ than this
		if (params_->RandomModulo(2) && future_ != nullptr) { //just flip the 2 nodes
			//	From
			//	Present		Future
			//	1(this)-----2
			//	4-----------3(cRand)
			//	To
			//	Present		Future
			//	1-----------3(cRand)
			//	4-----------2
			return makeChangeSimpleLinkage4(limit, cRand->neighbor_);
		}
		//find a new node to be connected
		// try to find a free node to connect the past_ of the conflict node.
		//count the number of free node
		//	From
		//	Present		Future
		//	1(this)-----2
		//	4-----------3(cRand)
		int count = 0;
		for (const ChainNeighbor* c = cRand->neighbor_->past_->neighborsFuture_; c != nullptr; c = c->next_) {
			if (c->neighbor_->fusion1Past_ == nullptr && c->neighbor_->fissionPast_ == nullptr) {
				if (c->neighbor_->past_ == nullptr) { count++; }
				else { if (c->neighbor_ == future_) { count++; } }
			}
		}
		if (count == 0) { //no node available to relink the conflicted node
			// try to break the link
			//	From
			//	Present		Future
			//	1(this)-----2
			//	4-----------3(cRand)
			//	To
			//	Present		Future
			//	1-----------3(cRand)
			//	NULL--------2
			//	4-----------NULL
			return makeChangeSimpleLinkage2(limit, cRand->neighbor_);
		}
		//there is several available nodes to relink the conflicted node
		//	From
		//	Present		Future
		//	1-----------2
		//	4-----------3(cRand)
		//	NULL--------5(c)

		//	To
		//	Present		Future
		//	1-----------3(cRand)
		//	4-----------5(c)
		//	NULL--------2
		return makeChangeSimpleLinkage3(limit, cRand->neighbor_, count);
	}
	return makeChangeSimpleLinkage0(limit);
}
