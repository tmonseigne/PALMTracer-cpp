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
#include <iostream>
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
	: x_(0), y_(0), z_(0), intensity_(0), mean_(0), surf_(0), color_(0), molId_(0),
	  params_(params), time_(time),
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
	while (neighborsFuture_) {
		l = neighborsFuture_->next_;
		delete neighborsFuture_;
		neighborsFuture_ = l;
	}

	while (neighborsPast_) {
		l = neighborsPast_->next_;
		delete neighborsPast_;
		neighborsPast_ = l;
	}
	trackedIndex_ = 0;
}

double Node::Distance(const Node& node) const
{
	const double Dx = x_ - node.x_, Dy = y_ - node.y_;
	return sqrt(Dx * Dx + Dy * Dy);
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

	for (const ChainNeighbor* c = neighborsFuture_; c; c = c->next_) {
		if (!c->neighbor_->past_) {
			const double Val = Likeness(c->neighbor_);

			if (Val > max) {
				max = Val;
				nMax = c->neighbor_;
			}
		}
	}
	if (nMax) {
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
	if (*neighborsHead) {
		if ((*neighborsHead)->next_) {
			bool change;
			do {
				change = false;
				ChainNeighbor* l = *neighborsHead;
				ChainNeighbor* lOld = nullptr;
				while (l->next_) {
					if (l->next_->neighbor_->intensity_ > l->neighbor_->intensity_) {
						ChainNeighbor* l2 = l->next_;
						if (!lOld) { *neighborsHead = l2; }
						else { lOld->next_ = l2; }
						l->next_ = l2->next_;
						l2->next_ = l;
						change = true;
						lOld = l2;
					}
					else {
						l = l->next_;
						if (!lOld) { lOld = *neighborsHead; }
					}
				}
			}
			while (change);
		}
	}
}

//==================== Computes ====================
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

double Node::distance(const Node& n1, const Node& n2) const
{
	const double N1N2X = (n1.x_ * n1.intensity_ + n2.x_ * n2.intensity_) / (n1.intensity_ + n2.intensity_);
	const double N1N2Y = (n1.y_ * n1.intensity_ + n2.y_ * n2.intensity_) / (n1.intensity_ + n2.intensity_);
	const double Dx = x_ - N1N2X, Dy = y_ - N1N2Y;
	return sqrt(Dx * Dx + Dy * Dy);
}

/*********************************************************************************************
double distanceToBorder()
  get the distance from this to the closer border
*********************************************************************************************/
double Node::distanceToBorder() const
{
	const double X = x_, Y = y_;
	const double MaxX = params_->maxX_, MaxY = params_->maxY_;
	return std::min({X, MaxX - X, Y, MaxY - Y});
}

ChainNeighbor* Node::getRandomNeighbor(ChainNeighbor* head, const size_t size) const
{
	const size_t Index = params_->RandomModulo(size);
	size_t i = 0;
	while (head && i < Index) {
		head = head->next_;
		++i;
	}
	return head;
}

/*********************************************************************************************
double SphericEffect()
  Compute the mean of radius of the node this
*********************************************************************************************/
double Node::sphericEffect() const
{
	return sqrt(surf_ / std::numbers::pi);
	//volume of a sphere is 4/3 * PI *r^3, surface is PI * r^2
	//if constexpr (Parameters::N_DIM == 2) { return sqrt(surf_ / std::numbers::pi); }
	//else { return pow(GetVolume() * 3.0 / 4.0 / std::numbers::pi, 1.0 / 3.0); }
}

//==================== likeness ====================
/*********************************************************************************************
double Node::likenessDistance(Node* n)
  Compute the likeness from the distance between this and n.
  New version from 2005/01/12
*********************************************************************************************/
double Node::likenessDistance(const Node* n) const
{
	return !n ? getLikeness(std::min(distanceToBorder_, params_->stdDistance_), params_->stdDistance_) : getLikeness(Distance(*n), params_->stdDistance_);
}

/*********************************************************************************************
double Node::likenessDistanceFusion(Node* n1, Node* n2)
  Compute the likeness from the distance of a fusion between this and n1 and n2.
  New version from 2005/01/12
*********************************************************************************************/
double Node::likenessDistanceFusion(const Node* n1, const Node* n2) const
{
	if (!n1 || !n2) { return 0.; }
	return getLikeness(fabs(distance(*n1, *n2) + n1->Distance(*n2) - sphericEffect(*n1, *n2)), params_->stdDistance_);
}

/*********************************************************************************************
double Node::likenessMean(Node* n)
  Compute the likeness from the mean of intensity between this and n.
  New version from 2005/01/12

  We use the hypothesis that
	sigma(mean intensity) proportional to surface^-1/6 simga_noise / Mean in 2D.
*********************************************************************************************/
double Node::likenessMean(const Node* n) const { return !n ? likenessMean_ : getLikeness(mean_ - n->mean_, params_->stdMean_); }

/*********************************************************************************************
double Node::likenessMeanFusion(Node* n1, Node* n2)
  Compute the likeness from the mean of intensity in fusion between this and n1 and n2 .
  New version from 2005/01/12

*********************************************************************************************/
double Node::likenessMeanFusion(const Node* n1, const Node* n2) const
{
	if (!n1 || !n2) { return 0.; } //Error ???
	return getLikeness(mean_ - (n1->intensity_ + n2->intensity_) / (n1->surf_ + n2->surf_), params_->stdMean_);
	//if (Parameters::N_DIM == 2) { return getLikeness(mean_ - (n1->intensity_ + n2->intensity_) / (n1->surf_ + n2->surf_), params_->stdMean_);}
	//else { return getLikeness(mean_ - (n1->intensity_ + n2->intensity_) / (n1->GetVolume() + n2->GetVolume()), params_->stdMean_);}
}

/*********************************************************************************************
double Node::likenessIntensity(Node* n)
  Compute the likeness from the total intensity between this and n.
  New version from 2005/01/12

  We use the hypothesis that sigma(intensity) proportional to simga_noise/Mean.
*********************************************************************************************/
double Node::likenessIntensity(const Node* n) const { return !n ? likenessIntensity_ : getLikeness(intensity_ - n->intensity_, params_->stdIntensity_); }

/*********************************************************************************************
double Node::likenessIntensityFusion(Node* n1, Node* n2
  Compute the likeness for the intensity in a fusion between this and n1 and n2.
  New version from 2005/01/12

  We use the hypothesis that sigma(intensity) proportional to simga_noise/Mean.
*********************************************************************************************/
double Node::likenessIntensityFusion(const Node* n1, const Node* n2) const
{
	if (!n1 || !n2) { return 0.; } //inutil
	return getLikeness(intensity_ - n1->intensity_ - n2->intensity_, params_->stdIntensity_);
}

/*********************************************************************************************
double Node::likeness(Node* n)
  Compute the likeness between this and n.
  New version from 2005/01/12
*********************************************************************************************/
double Node::Likeness(const Node* n) const
{
	return !n
		       ? (WEIGHT_DISTANCE * likenessDistance() + WEIGHT_INTENSITY * likenessIntensity()
		          + WEIGHT_MEAN_INTENSITY * likenessMean()) / 2. - params_->costBirthDeath_
		       : WEIGHT_DISTANCE * likenessDistance(n) + WEIGHT_INTENSITY * likenessIntensity(n) + WEIGHT_MEAN_INTENSITY * likenessMean(n);
}

void Node::UpdateConstants()
{
	// compute n->likeness, distance to border....
	likenessDistance_ = getLikeness(std::min(distanceToBorder_, params_->stdDistance_), params_->stdDistance_);
	likenessMean_ = getLikeness(mean_, params_->stdMean_ * std::max(1., 3. * params_->stdDistance_ / distanceToBorder_));
	likenessIntensity_ = getLikeness(intensity_, params_->stdIntensity_ * std::max(1., 3. * params_->stdDistance_ / distanceToBorder_));
	likeness_ = (WEIGHT_DISTANCE * likenessDistance() + WEIGHT_INTENSITY * likenessIntensity()
	             + WEIGHT_MEAN_INTENSITY * likenessMean()) / 2. - params_->costBirthDeath_;
}

/*********************************************************************************************
double Node::likeness(Node* n)
  Compute the likeness between this and n.
  New version from 2005/01/12
*********************************************************************************************/
double Node::likenessBlink(const Node* n) const
{
	constexpr double Alpha = 0.7;
	return Likeness(n) * Alpha + (Likeness() + n->Likeness()) * (1. - Alpha / 2.);
}

/*********************************************************************************************
double Node::likenessFusion(Node* n1, Node* n2)
  Compute the likeness in fusion between this and n1 and n2.
  New version from 2005/01/12
*********************************************************************************************/
double Node::likenessFusion(const Node* n1, const Node* n2) const
{
	if (!n1) { return Likeness(n2); }
	if (!n2) { return Likeness(n1); }
	return (WEIGHT_DISTANCE * likenessDistanceFusion(n1, n2) + WEIGHT_INTENSITY * likenessIntensityFusion(n1, n2)
	        + WEIGHT_MEAN_INTENSITY * likenessMeanFusion(n1, n2) + MIN_ENERGY_FUSION) * 3. / 2.;
	// the *3/2 is because the sum of the likeness of 1 unconnected node and 2 linked nodes has to be compared to the likeness of the fusion connection between those 3 nodes.
}

/*********************************************************************************************
void Node::computeConstantSigmas(double &intensity,
	double &mean, double &distance, long &countStd)
  Find in the vicinity of each node the one with the biggest likeness to compute the constant
  kInt and kMean and kDist
*********************************************************************************************/
void Node::ComputeConstantSigmas(double& intensity, size_t& nInt, double& mean, size_t& nMean, double& distance, size_t& nDist) const
{
	const Node* nMax = nullptr;
	double likenessMax = -1e10;

	for (const ChainNeighbor* c = neighborsFuture_; c; c = c->next_) {
		const double LikenessVal = Likeness(c->neighbor_);
		if (LikenessVal > likenessMax) {
			likenessMax = LikenessVal;
			nMax = c->neighbor_;
		}
	}

	if (nMax) {
		intensity += (intensity_ - nMax->intensity_) * (intensity_ - nMax->intensity_);
		mean += (mean_ - nMax->mean_) * (mean_ - nMax->mean_);
		distance += Distance(*nMax) * Distance(*nMax);
		nInt++;
		nMean++;
		nDist++;
	}
}

void Node::ComputeBorderDistance()
{
	const double X = x_, Y = y_;
	const double MaxX = params_->maxX_, MaxY = params_->maxY_;
	distanceToBorder_ = std::min({X, MaxX - X, Y, MaxY - Y});
}

//==================== Make Change ====================
double Node::MakeChange(const double limit)
{
	//choose an operation between:
	//simple linkage
	//strucutre fusion
	//strucutre dissociation

	if (params_->allowFusionAndFission_) {
		if (future_) {
			switch (params_->RandomModulo(4)) {
				case 0: return makeChangeFusion(limit);
				case 1: return makeChangeFission(limit);
				default: return makeChangeSimpleLinkage(limit);
			}
		}

		if constexpr (Parameters::N_TEMPORAL_WINDOWS > 1 && blinkFuture_) { return makeChangeBlinkLinkage(limit); }
		if (fission1Future_) { return makeChangeFission(limit); }
		if (fusionFuture_) { return makeChangeFusion(limit); }
		if constexpr (Parameters::N_TEMPORAL_WINDOWS > 1 && params_->RandomModulo(2) == 0) { return makeChangeBlinkLinkage(limit); }

		return makeChangeSimpleLinkage(limit);
	}

	if (future_) { return makeChangeSimpleLinkage(limit); }
	if (blinkFuture_) { return makeChangeBlinkLinkage(limit); }
	if constexpr (Parameters::N_TEMPORAL_WINDOWS > 2 && params_->RandomModulo(2) == 0) { return makeChangeBlinkLinkage(limit); }

	return makeChangeSimpleLinkage(limit);
}

//==================== Fission ====================
double Node::makeChangeFission(const double limit)
{
	//is the node already makes a Fission with another structure
	if (!future_ && !fission1Future_) { return 0.; } //nothing connected in the future!
	if (blinkFuture_) { return 0.; }                 // crand is already implicated in a future blink.
	if (fusionFuture_) { return 0.; }                //already implicated in a Fusion

	//no fission in old configuration
	//try to find 2 structures to make fission:
	//From
	//1(this)-------2
	//4-------------3(cRand)
	//To
	//1(this)-------2
	//        \-----3(cRand)
	//4-------------NULL
	if (!fission1Future_) { return fissionCreate(limit); }

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
		case 0: return fissionSymmetric(limit, fission1Future_);
		case 1: return fissionSymmetric(limit, fission2Future_);
		default: ;
	}

	if (fission2Future_->neighborsPastCount_ < 2) { //no neighbors like "4"
		switch (params_->RandomModulo(2)) {
			case 0: return fissionSymmetric(limit, fission1Future_);
			case 1: return fissionSymmetric(limit, fission2Future_);
			default: ;
		}
	}

	//try to find a neighbor_ "4" in the neighborhood of "3"
	//set the node cRand
	const ChainNeighbor* cRand = getRandomNeighbor(fission2Future_->neighborsPast_, fission2Future_->neighborsPastCount_);
	//test if the node is good

	if (cRand->neighbor_ == this) { return 0; }        // nothing to do, the new random link is the same as this.
	if (cRand->neighbor_->blinkFuture_) { return 0.; } // crand is already implicated in a future blink.

	if (cRand->neighbor_->future_) {
		switch (params_->RandomModulo(7)) {
			case 0: return fissionReassignBoth(limit, cRand->neighbor_, fission1Future_);
			case 1: return fissionSwapFuture(limit, cRand->neighbor_, fission2Future_);
			case 2: return fissionSwapFuture(limit, cRand->neighbor_, fission1Future_);
			case 3: return fissionReassignBoth(limit, cRand->neighbor_, fission2Future_);
			default: ;
		}
	}

	switch (params_->RandomModulo(3)) {
		case 0: return fissionTransferToNode(limit, cRand->neighbor_, transfer_type::keep_first);
		case 1: return fissionTransferToNode(limit, cRand->neighbor_, transfer_type::keep_second);
		default: return fissionTransferToNode(limit, cRand->neighbor_, transfer_type::give_all);
	}
}

double Node::fissionCreate(const double limit)
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

	if (blinkFuture_) { return -1; }              //ERROR
	if (neighborsFutureCount_ < 2) { return 0.; } //no neighbors like 3

	//set the node cRand
	const ChainNeighbor* cRand = getRandomNeighbor(neighborsFuture_, neighborsFutureCount_);
	//test if the node is good
	if (!cRand) { return -1.; } //ERROR
	Node* n = cRand->neighbor_;
	// n is the same or already implicated in a fusion / fission / past blink.
	if (n == future_ || !isBlinkLinkageAvailable(n)) { return 0.; }

	if (Distance(*n) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = likenessFusion(n, future_)
	                      + (!n->past_ ? 0. : n->past_->Likeness())
	                      - Likeness(future_) - n->Likeness(n->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness < limit // if the new solution is too bad
	// change the configuration
	fission1Future_ = future_;
	fission2Future_ = n;
	future_->fissionPast_ = this;
	n->fissionPast_ = this;
	future_->past_ = nullptr;
	future_ = nullptr;

	if (n->past_) { n->past_->future_ = nullptr; }
	n->past_ = nullptr;

	return DeltaE;
}

double Node::fissionSymmetric(const double limit, Node* keep)
{
	// Old makeChangeFission3 & makeChangeFission4
	Node* tmp = (keep == fission1Future_) ? fission2Future_ : fission1Future_;

	const double DeltaE = Likeness(keep) + tmp->Likeness() - likenessFusion(fission1Future_, fission2Future_);
	if (!metropolis(DeltaE, limit)) { return 0.; }

	future_ = keep;
	tmp->fissionPast_ = nullptr;
	future_->past_ = this;
	future_->fissionPast_ = nullptr;
	fission1Future_ = nullptr;
	fission2Future_ = nullptr;

	return DeltaE;
}

double Node::fissionReassignBoth(const double limit, Node* node, Node* keep)
{
	// Old makeChangeFission5 & makeChangeFission8
	Node* give = (keep == fission1Future_) ? fission2Future_ : fission1Future_;
	if (node->Distance(*give) >= params_->maxSpeed_ || node->blinkFuture_) { return 0.; }

	const double DeltaE = Likeness(keep) + node->likenessFusion(give, node->future_)
	                      - likenessFusion(fission1Future_, fission2Future_) - node->Likeness(node->future_);

	if (!metropolis(DeltaE, limit)) { return 0.; }

	// Changes
	future_ = keep;
	future_->past_ = this;
	node->fission1Future_ = give;
	give->fissionPast_ = node;
	node->fission2Future_ = node->future_;
	node->fission2Future_->fissionPast_ = node;

	// Clean
	keep->fissionPast_ = nullptr;
	fission1Future_ = nullptr;
	fission2Future_ = nullptr;
	node->future_->past_ = nullptr;
	node->future_ = nullptr;

	return DeltaE;
}

double Node::fissionSwapFuture(const double limit, Node* node, Node* swap)
{
	// Old makeChangeFission6 & makeChangeFission7
	if (node->blinkFuture_ ||
	    (node->future_ && Distance(*node->future_) >= params_->maxSpeed_) ||
	    node->Distance(*swap) >= params_->maxSpeed_) { return 0.; }

	const double DeltaE = likenessFusion(swap, node->future_)
	                      + node->Likeness((swap == fission1Future_ ? fission2Future_ : fission1Future_))
	                      - likenessFusion(fission1Future_, fission2Future_) - node->Likeness(node->future_);
	if (!metropolis(DeltaE, limit)) { return 0.; }

	swap->past_ = node;
	node->future_->fissionPast_ = this;

	Node* tmp = node->future_;
	node->future_ = swap;
	swap->fissionPast_ = nullptr;

	if (swap == fission1Future_) {
		fission1Future_ = tmp;
		fission1Future_->past_ = nullptr;
	}
	else {
		fission2Future_ = tmp;
		fission2Future_->past_ = nullptr;
	}

	return DeltaE;
}

double Node::fissionTransferToNode(const double limit, Node* node, const transfer_type type)
{
	// Old makeChangeFission1 & makeChangeFission2 & makeChangeFission9 
	Node* f1 = fission1Future_;
	Node* f2 = fission2Future_;
	if (!f1 || !f2) { return 0.; }
	if (node->blinkFuture_ || node->fusionFuture_) { return 0.; }
	if ((type == transfer_type::keep_first && node->fission1Future_) ||
	    (type == transfer_type::keep_second && node->fission2Future_)) { return 0.; }

	if (type == transfer_type::give_all) {
		if ((node->future_ && Distance(*node->future_) >= params_->maxSpeed_) ||
		    node->Distance(*f1) >= params_->maxSpeed_ ||
		    node->Distance(*f2) >= params_->maxSpeed_) { return 0.; }
	}
	else {
		Node* keep = (type == transfer_type::keep_first ? f1 : f2);
		const Node* give = (type == transfer_type::keep_first ? f2 : f1);
		if (node->Distance(*give) >= params_->maxSpeed_) { return 0.; }
	}

	double DeltaE;
	if (type == transfer_type::give_all) {
		DeltaE = Likeness(node->future_) + node->likenessFusion(f1, f2)
		         - likenessFusion(f1, f2) - node->Likeness(node->future_);
	}
	else {
		const Node* keep = (type == transfer_type::keep_first ? f1 : f2);
		const Node* give = (type == transfer_type::keep_first ? f2 : f1);
		DeltaE = Likeness(keep) + node->Likeness(give)
		         + (node->future_ ? node->future_->Likeness() : 0.)
		         - likenessFusion(f1, f2) - node->Likeness(node->future_);
	}

	if (!metropolis(DeltaE, limit)) { return 0.; }

	// Apply change
	if (type == transfer_type::give_all) {
		node->fission1Future_ = f1;
		node->fission2Future_ = f2;
		f1->fissionPast_ = node;
		f2->fissionPast_ = node;

		future_ = node->future_;
		if (future_) { future_->past_ = this; }
		node->future_ = nullptr;
	}
	else {
		Node* keep = (type == transfer_type::keep_first ? f1 : f2);
		Node* give = (type == transfer_type::keep_first ? f2 : f1);

		future_ = keep;
		future_->past_ = this;
		if (node->future_) { node->future_->past_ = nullptr; }
		node->future_ = give;
		give->past_ = node;
		give->fissionPast_ = node;
		keep->fissionPast_ = nullptr;
	}

	fission1Future_ = nullptr;
	fission2Future_ = nullptr;

	return DeltaE;
}

//==================== Fusion ====================
double Node::makeChangeFusion(const double limit)
{
	//is the node already making a Fusion with another structure
	if (!future_ && !fusionFuture_) { return 0.; } //nothing connected in the future!
	if (fission1Future_) { return 0.; }            //already implicated in a Fission
	//no fusion in old configuration
	//try to find 2 structures to make fusion:
	//From
	//1(this)-------2
	//3-------------4(cRand)
	//To
	//1-------------2
	//3-------/
	//NULL----------4(cRand)
	if (!fusionFuture_) { return fusionCreate(limit); }
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
		case 0: return fusionSymmetric(limit, fusionFuture_->fusion1Past_);
		case 1: return fusionSymmetric(limit, fusionFuture_->fusion2Past_);
		default: ;
	}


	//get the cRand
	if (neighborsFutureCount_ < 2) { //no neighbors like "4", no cRand
		switch (params_->RandomModulo(2)) {
			case 0: return fusionSymmetric(limit, fusionFuture_->fusion1Past_);
			case 1: return fusionSymmetric(limit, fusionFuture_->fusion2Past_);
			default: ;
		}
	}

	const ChainNeighbor* cRand = getRandomNeighbor(neighborsFuture_, neighborsFutureCount_);
	//test if the node is good
	if (!cRand) { return -1.; } //ERROR
	Node* n = cRand->neighbor_;

	if (n->fissionPast_) { return fusionTransferFull(limit, n); }
	if (n->blinkPast_) { return 0; } // cRand is already implicated in a blink connection

	if (n->past_) {
		switch (params_->RandomModulo(7)) {
			case 0: return fusionReassignBoth(limit, n, fusionFuture_->fusion1Past_);
			case 1:
			case 2: return fusionSwapFuture(limit, n, fusionFuture_->fusion2Past_);
			case 3: return fusionSwapFuture(limit, n, fusionFuture_->fusion1Past_);
			case 4: return fusionReassignBoth(limit, n, fusionFuture_->fusion2Past_);
			default: ;
		}
	}

	switch (params_->RandomModulo(2)) {
		case 0: return fusionTransferToNode(limit, n, transfer_type::keep_first);
		default: return fusionTransferToNode(limit, n, transfer_type::keep_second);
	}
}

double Node::fusionCreate(const double limit)
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
	//set the node cRand
	const ChainNeighbor* cRand = getRandomNeighbor(future_->neighborsPast_, future_->neighborsPastCount_);
	//test if the node is good
	if (!cRand) { return -1.; } //ERROR
	Node* n = cRand->neighbor_;
	// n is the same or already implicated in a fission/fusion/blink.
	if (n == this || n->fission1Future_ || n->fusionFuture_ || n->blinkFuture_) { return 0; }

	if (future_->Distance(*n) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = future_->likenessFusion(n, this)
	                      + (n->future_ ? n->future_->Likeness() : 0.)
	                      - Likeness(future_) - n->Likeness(n->future_);

	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	//unlink future_

	fusionFuture_ = future_;
	future_->fusion1Past_ = this;
	n->fusionFuture_ = future_;
	future_->fusion2Past_ = n;

	if (n->future_) { n->future_->past_ = nullptr; }

	n->future_ = nullptr;
	future_ = nullptr;
	fusionFuture_->past_ = nullptr;

	return DeltaE;
}

double Node::fusionSymmetric(const double limit, const Node* keep)
{
	// Old makeChangeFusion3 & makeChangeFusion4
	Node* tmp = (keep == fusionFuture_->fusion1Past_) ? fusionFuture_->fusion2Past_ : fusionFuture_->fusion1Past_;
	const double DeltaE = Likeness(keep) + tmp->Likeness() - fusionFuture_->likenessFusion(tmp, this);
	if (!metropolis(DeltaE, limit)) { return 0.; }

	fusionFuture_->past_ = (keep == this) ? this : tmp;
	fusionFuture_->fusion1Past_ = nullptr;
	fusionFuture_->fusion2Past_ = nullptr;
	(keep == this ? future_ : tmp->future_) = fusionFuture_;
	fusionFuture_ = nullptr;
	(tmp == this ? future_ : tmp->fusionFuture_) = nullptr;

	return DeltaE;
}

double Node::fusionReassignBoth(const double limit, Node* node, const Node* keep)
{
	// Old makeChangeFusion5 & makeChangeFusion8
	Node* tmp = (keep == fusionFuture_->fusion1Past_) ? fusionFuture_->fusion2Past_ : fusionFuture_->fusion1Past_;
	if (tmp->Distance(*node) >= params_->maxSpeed_ || node->fissionPast_ || node->fusion1Past_ || !node->past_) { return 0.; }

	const double DeltaE = Likeness(fusionFuture_) + node->likenessFusion(tmp, node->past_)
	                      - fusionFuture_->likenessFusion(tmp, this) - node->Likeness(node->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; }

	future_ = fusionFuture_;
	future_->past_ = this;
	tmp->fusionFuture_ = node;
	node->fusion1Past_ = tmp;
	node->past_->fusionFuture_ = node;
	node->fusion2Past_ = node->past_;
	node->past_->future_ = nullptr;
	node->past_ = nullptr;
	future_->fusion1Past_ = nullptr;
	future_->fusion2Past_ = nullptr;
	fusionFuture_ = nullptr;

	return DeltaE;
}

double Node::fusionSwapFuture(const double limit, Node* node, const Node* swap) const
{
	// Old makeChangeFusion6 & makeChangeFusion7
	Node* tmp = (swap == fusionFuture_->fusion1Past_) ? fusionFuture_->fusion2Past_ : fusionFuture_->fusion1Past_;
	if (tmp->Distance(*node) >= params_->maxSpeed_) { return 0.; }
	if (node->past_ && node->past_->Distance(*fusionFuture_) >= params_->maxSpeed_) { return 0.; }

	const double DeltaE = fusionFuture_->likenessFusion(this, node->past_) + node->Likeness(tmp)
	                      - fusionFuture_->likenessFusion(tmp, this) - node->Likeness(node->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; }

	node->past_->fusionFuture_ = fusionFuture_;
	node->past_->future_ = nullptr;
	if (swap == fusionFuture_->fusion1Past_) { fusionFuture_->fusion2Past_ = node->past_; }
	else { fusionFuture_->fusion1Past_ = node->past_; }
	tmp->fusionFuture_ = nullptr;
	tmp->future_ = node;
	node->past_ = tmp;

	return DeltaE;
}

double Node::fusionTransferToNode(const double limit, Node* node, const transfer_type type)
{
	// Old makeChangeFusion1 & makeChangeFusion2
	Node* f1 = fusionFuture_->fusion1Past_;
	Node* f2 = fusionFuture_->fusion2Past_;
	Node* keep = (type == transfer_type::keep_first) ? f1 : f2;
	Node* give = (keep == f1) ? f2 : f1;

	if (node->fissionPast_ || node->fusion1Past_) { return 0.; }
	if (node->Distance(*give) >= params_->maxSpeed_) { return 0.; }

	const double DeltaE = Likeness(keep) + node->Likeness(give) + (node->past_ ? node->past_->Likeness() : 0.)
	                      - fusionFuture_->likenessFusion(f1, f2) - node->Likeness(node->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; }

	fusionFuture_->past_ = this;
	fusionFuture_->fusion1Past_ = nullptr;
	fusionFuture_->fusion2Past_ = nullptr;
	give->fissionPast_ = node;
	node->future_ = give;
	give->past_ = node;
	future_ = fusionFuture_;
	fusionFuture_ = nullptr;
	if (node->past_) { node->past_->future_ = nullptr; }
	node->past_ = keep;

	return DeltaE;
}

double Node::fusionTransferFull(const double limit, Node* node)
{
	// Old makeChangeFusion9
	Node* tmp = (fusionFuture_->fusion1Past_ == this) ? fusionFuture_->fusion2Past_ : fusionFuture_->fusion1Past_;
	Node* tmp2 = (node->fissionPast_->fission1Future_ == node) ? node->fissionPast_->fission2Future_ : node->fissionPast_->fission1Future_;

	if (tmp->Distance(*node) >= params_->maxSpeed_) { return 0.; }

	const double DeltaE = Likeness(fusionFuture_) + tmp->Likeness(node) + node->fissionPast_->Likeness(tmp2)
	                      - fusionFuture_->likenessFusion(tmp, this)
	                      - node->fissionPast_->likenessFusion(node, tmp2);
	if (!metropolis(DeltaE, limit)) { return 0.; }

	fusionFuture_->past_ = this;
	future_ = fusionFuture_;
	tmp->fusionFuture_ = nullptr;
	tmp->future_ = node;
	node->past_ = tmp;
	node->fissionPast_ = nullptr;

	tmp2->fissionPast_->fission1Future_ = nullptr;
	tmp2->fissionPast_->fission2Future_ = nullptr;
	tmp2->fissionPast_->future_ = tmp2;
	tmp2->past_ = tmp2->fissionPast_;
	tmp2->fissionPast_ = nullptr;
	fusionFuture_->fusion1Past_ = nullptr;
	fusionFuture_->fusion2Past_ = nullptr;
	fusionFuture_ = nullptr;

	return DeltaE;
}

//==================== Linkage ====================
double Node::makeChangeSimpleLinkage(const double limit)
{
	if (!future_ || params_->RandomModulo(2) == 0) {
		if (neighborsFutureCount_ == 0) { return 0.; } //no neighbors
		//cRand is a random node in the vicinity of "this" in the future.
		const ChainNeighbor* cRand = getRandomNeighbor(neighborsFuture_, neighborsFutureCount_);
		//test if the node is good
		if (!cRand) { return -1.; } //ERROR
		Node* target = cRand->neighbor_;
		// the new random link is the same as future_.
		if (target == future_) { return simpleLinkageDisconnection(limit); } //try simple deconnection
		// target already implicated in a fission / fusion / blink connection
		if (!isBlinkLinkageAvailable(target)) { return 0; }

		//is there a conflict between this and the pastConnextion of target?
		//no conflict target can be linked with this.
		//	From
		//	Present		Future
		//	1(this)-----2
		//	NULL--------3(cRand)
		//	To
		//	Present		Future
		//	1(this)-----3(cRand)
		//	NULL--------2
		if (!target->past_) { return simpleLinkageReassignment(limit, target); }

		//there is a conflict, the node target is already connected to a node at the same time_ than this
		//just flip the 2 nodes
		//	From
		//	Present		Future
		//	1(this)-----2
		//	4-----------3(cRand)
		//	To
		//	Present		Future
		//	1-----------3(cRand)
		//	4-----------2
		if (params_->RandomModulo(2) && future_) { return simpleLinkageSwap(limit, target); }

		//find a new node to be connected
		// try to find a free node to connect the past_ of the conflict node.
		//count the number of free node
		//	From
		//	Present		Future
		//	1(this)-----2
		//	4-----------3(cRand)
		int count = 0;
		for (const ChainNeighbor* c = target->past_->neighborsFuture_; c; c = c->next_) { if (isSimpleLinkageAvailable(c->neighbor_)) { count++; } }

		//no node available to relink the conflicted node
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
		if (count == 0) { return simpleLinkageBreakAndReconnect(limit, target); }

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
		return simpleLinkageDoubleSwap(limit, target, count);
	}
	return simpleLinkageDisconnection(limit);
}

bool Node::isSimpleLinkageAvailable(const Node* node) const { return !node->fusion1Past_ && !node->fissionPast_ && (!node->past_ || node == future_); }

double Node::simpleLinkageDisconnection(const double limit) const
{
	//From 
	//Present	Future
	//1(this)-----2(future_)
	//To
	//1(this)-----NULL
	//NULL--------2(future_)

	const double DeltaE = Likeness() + future_->Likeness() - Likeness(future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	return 0.;
}

double Node::simpleLinkageReassignment(const double limit, Node* node)
{
	//	From
	//	Present		Future
	//	1(this)-----2
	//	NULL--------3(cRand)
	//	To
	//	Present		Future
	//	1(this)-----3(cRand)
	//	NULL--------2

	if (node->blinkPast_) { return 0; }                       // cRand is already implicated in a blink connection
	if (Distance(*node) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = Likeness(node) + (future_ ? future_->Likeness() : 0) - node->Likeness() - Likeness(future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	//unlink future_
	if (future_) { future_->past_ = nullptr; }
	future_ = node;
	future_->past_ = this;

	return DeltaE;
	//return(0.);
}

double Node::simpleLinkageBreakAndReconnect(const double limit, Node* node)
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

	if (Distance(*node) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = Likeness(node) + (!future_ ? 0 : future_->Likeness())
	                      + (!node->past_ ? 0 : node->past_->Likeness()) - Likeness(future_) - node->Likeness(node->past_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	//unlink future_


	if (future_) { future_->past_ = nullptr; }
	future_ = node;
	future_->past_->future_ = nullptr;
	future_->past_ = this;

	return DeltaE;
	//return(0.);
}

double Node::simpleLinkageDoubleSwap(const double limit, Node* node, const size_t count)
{
	//	From
	//	Present		Future
	//	1-----------2
	//	4-----------3(cRand)
	//	NULL--------5(c)

	//find a random one
	const size_t IdRandNeighbor = params_->RandomModulo(count);
	size_t iCount = 0;

	const ChainNeighbor* c = node->past_->neighborsFuture_;
	do {
		if (isSimpleLinkageAvailable(c->neighbor_)) { iCount++; }
		if (iCount <= IdRandNeighbor) { c = c->next_; }
	}
	while (iCount <= IdRandNeighbor);


	//	To
	//	Present		Future
	//	1-----------3(cRand)
	//	4-----------5(c)
	//	NULL--------2
	if (!c || !isSimpleLinkageAvailable(c->neighbor_)) { return -1; }                             //ERROR 	
	if (Distance(*node) >= params_->maxSpeed_) { return 0.; }                                     // Distance is not respected
	if (node->blinkPast_ || c->neighbor_->blinkPast_) { return 0; }                               // node and c is already implicated in a blink connection
	if (node->past_ && node->past_->Distance(*c->neighbor_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = Likeness(node) + node->past_->Likeness(c->neighbor_) + (!future_ ? 0 : future_->Likeness())
	                      - Likeness(future_) - node->Likeness(node->past_) - c->neighbor_->Likeness();
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	//unlink future_

	if (future_) { future_->past_ = nullptr; } //????
	future_ = node;
	future_->past_->future_ = c->neighbor_;
	future_->past_->future_->past_ = future_->past_;
	future_->past_ = this;

	return DeltaE;
	//return(0.);
}

double Node::simpleLinkageSwap(const double limit, Node* node)
{
	//	From
	//	Present		Future
	//	1(this)-----2
	//	4-----------3(cRand)
	//	To
	//	Present		Future
	//	1-----------3(cRand)
	//	4-----------2

	if (node->past_ && node->past_->Distance(*future_) >= params_->maxSpeed_) { return 0.; } // Distance is not respected
	if (Distance(*node) >= params_->maxSpeed_) { return 0.; }                                // Distance is not respected

	const double DeltaE = Likeness(node) + future_->Likeness(node->past_) - node->Likeness(node->past_) - Likeness(future_);
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration

	node->past_->future_ = future_;
	future_->past_ = node->past_;
	future_ = node;
	node->past_ = this;

	return DeltaE;
	//return(0.);
}

//==================== Blink Linkage ====================
double Node::makeChangeBlinkLinkage(const double limit)
{
	//nothing connected in the future
	//transform : 
	//	From
	//	t			t+1		t+n
	//	1(this)-------------NULL
	//	NULL----------------2
	//	To
	//	t			t+1		t+n
	//	1(this)-------------2(blinkFuture_)
	if (!blinkFuture_) { return blinkLinkageConnection(limit); }

	//blinkFuture_!=NULL
	//transform : 
	//	From
	//	t			t+1		t+n
	//	1(this)-------------2(blinkFuture_)
	//	To
	//	t			t+1		t+n
	//	1(this)-------------NULL
	//	NULL----------------2
	//  blinkLinkageDisconnection

	//transform : 
	//	From
	//	t			t+1					t+n
	//	1(this)-------------------------2(blinkFuture_)
	//	4---------3(cRand)----...
	//  1=this, 2=blinkFuture_, 3=cRand->neighbor_, 4=cRand->neighbor_->past_
	//	To
	//	t			t+1					t+n
	//	4-------------------------------2(blinkFuture_)
	//	1(this)---3(cRand)----...
	//  blinkLinkageSwapBegin

	//OR

	//transform : 
	//	From
	//	t					 t+n-1		t+n
	//	1(this)-------------------------2(blinkFuture_)
	//				...------3(cRand)---4
	//  1=this, 2=blinkFuture_, 3=cRand->neighbor_, 4=cRand->neighbor_->future_
	//	To
	//	t					 t+n-1		t+n
	//	1(this)-------------------------4
	//				...------3(cRand)---2(blinkFuture_)
	//  blinkLinkageSwapEnd

	//OR

	//transform : 
	//	From
	//	t			t+1					t+n
	//	1(this)-------------------------2(blinkFuture_)
	//	NULL--------3(cRand2)----...
	//  1=this, 2=blinkFuture_, 3=cRand2->neighbor_
	//	To
	//	t			t+1					t+n
	//	NULL----------------------------2(blinkFuture_)
	//	1(this)---3(cRand2)----...
	//  blinkLinkageSwapStart

	//OR

	//transform : 
	//	From
	//	t					 t+n-1		t+n
	//	1(this)-------------------------2(blinkFuture_)
	//				...------3(cRand2)---NULL
	//  1=this, 2=blinkFuture_, 3=cRand2->neighbor_
	//	To
	//	t					 t+n-1		t+n
	//	1(this)-------------------------NULL
	//				...------3(cRand2)---2(blinkFuture_)
	//  blinkLinkageSwapFinish

	if (params_->RandomModulo(2)) { return blinkLinkageDisconnection(limit); }

	//find a potential cRand
	//cRand has to be either just after this or just before blinkFuture_
	if (params_->RandomModulo(2)) { //cRand just after this
		size_t nbNodeAtGoodTime = 0, nbNodeAtGoodTime2 = 0, idRandNeighbor;
		ChainNeighbor* cRand;

		for (cRand = neighborsFuture_; cRand; cRand = cRand->next_) {
			if (cRand->neighbor_->time_ - time_ == 1) {
				if (cRand->neighbor_->past_) { nbNodeAtGoodTime++; }
				else if (isBlinkLinkageAvailable(cRand->neighbor_)) { nbNodeAtGoodTime2++; }
			}
		}
		if (nbNodeAtGoodTime == 0 && nbNodeAtGoodTime2 == 0) { return blinkLinkageDisconnection(limit); }

		if ((params_->RandomModulo(2) && nbNodeAtGoodTime != 0) || nbNodeAtGoodTime2 == 0) {
			idRandNeighbor = params_->RandomModulo(nbNodeAtGoodTime);
			size_t i = 0;
			bool exitFor = false;
			for (cRand = neighborsFuture_; !exitFor && cRand;) {
				if (cRand->neighbor_->time_ - time_ == 1 && cRand->neighbor_->past_) {
					if (i == idRandNeighbor) { exitFor = true; }
					i++;
				}
				if (!exitFor) { cRand = cRand->next_; }
			}
			if (!cRand) { return -1.; } //ERROR 

			if (Distance(*cRand->neighbor_) >= params_->maxSpeed_) { return 0; } //cRand not linked to this!
			if (blinkFuture_->Distance(*cRand->neighbor_->past_) >= params_->maxSpeed_) { return 0; }
			//blinkFuture_ not linked to cRand->neighbor_->past_
			return blinkLinkageSwapBegin(limit, cRand->neighbor_);
		}
		ChainNeighbor* cRand2;
		idRandNeighbor = params_->RandomModulo(nbNodeAtGoodTime2);
		size_t i = 0;
		bool exitFor = false;
		for (cRand2 = neighborsFuture_; !exitFor && cRand2;) {
			if (cRand2->neighbor_->time_ - time_ == 1 && !cRand2->neighbor_->past_ && isBlinkLinkageAvailable(cRand2->neighbor_)) {
				if (i == idRandNeighbor) { exitFor = true; }
				i++;
			}
			if (!exitFor) { cRand2 = cRand2->next_; }
		}
		if (!cRand2) { return -1.; } //ERROR 

		if (Distance(*cRand2->neighbor_) >= params_->maxSpeed_) { return 0; } //cRand2 not linked to this!

		return blinkLinkageSwapStart(limit, cRand2->neighbor_);
	}
	// cRand is just before blinkFuture_
	size_t nbNodeAtGoodTime = 0, nbNodeAtGoodTime2 = 0, idRandNeighbor;
	ChainNeighbor* cRand;

	for (cRand = blinkFuture_->neighborsPast_; cRand; cRand = cRand->next_) {
		if (blinkFuture_->time_ - cRand->neighbor_->time_ == 1) {
			if (cRand->neighbor_->future_) { nbNodeAtGoodTime++; }
			else if (!cRand->neighbor_->fission1Future_ && !cRand->neighbor_->fusionFuture_ && !cRand->neighbor_->blinkFuture_) { nbNodeAtGoodTime2++; }
		}
	}
	if (nbNodeAtGoodTime == 0 && nbNodeAtGoodTime2 == 0) { return blinkLinkageDisconnection(limit); }

	if ((params_->RandomModulo(2) && nbNodeAtGoodTime != 0) || nbNodeAtGoodTime2 == 0) {
		idRandNeighbor = params_->RandomModulo(nbNodeAtGoodTime);
		size_t i = 0;
		bool exitFor = false;
		for (cRand = blinkFuture_->neighborsPast_; !exitFor && cRand;) {
			if (blinkFuture_->time_ - cRand->neighbor_->time_ == 1 && cRand->neighbor_->future_) {
				if (i == idRandNeighbor) { exitFor = true; }
				i++;
			}
			if (!exitFor) { cRand = cRand->next_; }
		}
		if (!cRand) { return -1.; }                                                        //ERROR 
		if (blinkFuture_->Distance(*cRand->neighbor_) >= params_->maxSpeed_) { return 0; } //cRand not linked to blinkFuture_!
		if (Distance(*cRand->neighbor_->future_) >= params_->maxSpeed_) { return 0; }      //this is not linked to cRand->neighbor_->future_!


		return blinkLinkageSwapEnd(limit, cRand->neighbor_);
	}
	ChainNeighbor* cRand2;
	idRandNeighbor = params_->RandomModulo(nbNodeAtGoodTime2);
	size_t i = 0;
	bool exitFor = false;
	for (cRand2 = blinkFuture_->neighborsPast_; !exitFor && cRand2;) {
		if (blinkFuture_->time_ - cRand2->neighbor_->time_ == 1 && !cRand2->neighbor_->future_
		    && !cRand2->neighbor_->fission1Future_ && !cRand2->neighbor_->fusionFuture_ && !cRand2->neighbor_->blinkFuture_) {
			if (i == idRandNeighbor) { exitFor = true; }
			i++;
		}
		if (!exitFor) { cRand2 = cRand2->next_; }
	}
	if (!cRand2) { return -1.; }                                                        //ERROR 
	if (blinkFuture_->Distance(*cRand2->neighbor_) >= params_->maxSpeed_) { return 0; } //cRand2 not linked to blinkFuture_!
	return blinkLinkageSwapFinish(limit, cRand2->neighbor_);
}

bool Node::isBlinkLinkageAvailable(const Node* node) { return !node->fusion1Past_ && !node->fissionPast_ && !node->blinkPast_; }

double Node::blinkLinkageConnection(const double limit)
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

	for (cRand = neighborsBlinkFuture_; cRand; cRand = cRand->next_) {
		if (cRand->neighbor_->time_ - time_ == Dt && !cRand->neighbor_->past_ && isBlinkLinkageAvailable(cRand->neighbor_)) { nbNodeAtGoodTime++; }
	}
	if (nbNodeAtGoodTime == 0) { return 0; }

	const size_t IdRandNeighbor = params_->RandomModulo(nbNodeAtGoodTime);
	size_t i = 0;
	bool exitFor = false;
	for (cRand = neighborsBlinkFuture_; !exitFor && cRand;) {
		if (cRand->neighbor_->time_ - time_ == Dt && !cRand->neighbor_->past_ && isBlinkLinkageAvailable(cRand->neighbor_)) {
			if (i == IdRandNeighbor) { exitFor = true; }
			i++;
		}
		if (!exitFor) { cRand = cRand->next_; }
	}


	if (!cRand) { return -1.; } //ERROR 

	const double DeltaE = likenessBlink(cRand->neighbor_) - cRand->neighbor_->Likeness() - Likeness();
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad

	// change the configuration
	cRand->neighbor_->blinkPast_ = this;
	blinkFuture_ = cRand->neighbor_;

	return DeltaE;
}

double Node::blinkLinkageDisconnection(const double limit)
{
	//transform : 
	//	From
	//	t			t+1		t+2
	//	1(this)-------------2(blinkFuture_)
	//	To
	//	t			t+1		t+2
	//	1(this)-------------NULL
	//	NULL----------------2
	const double DeltaE = blinkFuture_->Likeness() + Likeness() - likenessBlink(blinkFuture_);
	if (!metropolis(DeltaE, limit)) { return 0.; } // deltaLikeness<limit // if the new solution is too bad

	// change the configuration
	blinkFuture_->blinkPast_ = nullptr;
	blinkFuture_ = nullptr;

	return DeltaE;
}

double Node::blinkLinkageSwapBegin(const double limit, Node* node)
{
	//transform : 
	//	From
	//	t			t+1					t+n
	//	1(this)-------------------------2(blinkFuture_)
	//	4---------3(cRand)----...
	//  1=this, 2=blinkFuture_, 3=node, 4=node->past_
	//	To
	//	t			t+1					t+n
	//	4-------------------------------2(blinkFuture_)
	//	1(this)---3(cRand)----...

	const double DeltaE = Likeness(node) + node->past_->likenessBlink(blinkFuture_) - likenessBlink(blinkFuture_) - node->past_->Likeness(node);
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

double Node::blinkLinkageSwapEnd(const double limit, Node* node)
{
	//transform : 
	//	From
	//	t				t+n-1		t+n
	//	1(this)-------------------------2(blinkFuture_)
	//		...------3(cRand)---4
	//  1=this, 2=blinkFuture_, 3=node, 4=node->future_
	//	To
	//	t				t+n-1		t+n
	//	1(this)-------------------------4
	//		...------3(cRand)---2(blinkFuture_)

	const double DeltaE = likenessBlink(node->future_) + node->Likeness(blinkFuture_) - likenessBlink(blinkFuture_) - node->future_->Likeness(node);
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

double Node::blinkLinkageSwapStart(const double limit, Node* node)
{
	//transform : 
	//	From
	//	t			t+1					t+n
	//	1(this)-------------------------2(blinkFuture_)
	//	NULL--------3(cRand2)----...
	//  1=this, 2=blinkFuture_, 3=cRand2->neighbor_, 
	//	To
	//	t			t+1					t+n
	//	NULL----------------------------2(blinkFuture_)
	//	1(this)---3(cRand2)----...

	const double DeltaE = Likeness(node) + blinkFuture_->Likeness() - likenessBlink(blinkFuture_) - node->Likeness();
	if (!metropolis(DeltaE, limit)) { return 0.; } // deltaLikeness<limit // if the new solution is too bad

	// change the configuration
	blinkFuture_->blinkPast_ = nullptr;
	blinkFuture_ = nullptr;
	node->past_ = this;
	future_ = node;

	return DeltaE;
}

double Node::blinkLinkageSwapFinish(const double limit, Node* node)
{
	//transform : 
	//	From
	//	t				t+n-1		t+n
	//	1(this)-------------------------2(blinkFuture_)
	//		...------3(cRand2)---NULL
	//  1=this, 2=blinkFuture_, 3=cRand2->neighbor_
	//	To
	//	t				t+n-1		t+n
	//	1(this)-------------------------NULL
	//		...------3(cRand2)---2(blinkFuture_)
	//  blinkLinkageSwapFinish

	const double DeltaE = Likeness() + node->Likeness(blinkFuture_) - likenessBlink(blinkFuture_) - node->Likeness();
	if (!metropolis(DeltaE, limit)) { return 0.; } // deltaLikeness<limit // if the new solution is too bad

	// change the configuration
	blinkFuture_->blinkPast_ = nullptr;
	blinkFuture_->past_ = node;
	node->future_ = blinkFuture_;
	blinkFuture_ = nullptr;

	return DeltaE;
}
