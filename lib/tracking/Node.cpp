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
	  future_(nullptr), past_(nullptr),
	  fission1Future_(nullptr), fission2Future_(nullptr), fissionPast_(nullptr),
	  fusion1Past_(nullptr), fusion2Past_(nullptr), fusionFuture_(nullptr),
	  trackedIndex_(0), colorIndex_{},
	  neighborsPast_(nullptr), neighborsFuture_(nullptr),
	  neighborsPastCount_(0), neighborsFutureCount_(0) { }

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
double Node::likenessMean(Node* n)
  Compute the likeness from the mean of intensity between this and n.
  New version from 2005/01/12

  We use the hypothesis that
	sigma(mean intensity) proportional to surface^-1/6 simga_noise / Mean in 2D.
*********************************************************************************************/
double Node::likenessMean(const Node* n) const { return !n ? likenessMean_ : getLikeness(mean_ - n->mean_, params_->stdMean_); }

/*********************************************************************************************
double Node::likenessIntensity(Node* n)
  Compute the likeness from the total intensity between this and n.
  New version from 2005/01/12

  We use the hypothesis that sigma(intensity) proportional to simga_noise/Mean.
*********************************************************************************************/
double Node::likenessIntensity(const Node* n) const { return !n ? likenessIntensity_ : getLikeness(intensity_ - n->intensity_, params_->stdIntensity_); }

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
bool Node::isBlinkLinkageAvailable(const Node* node) { return !node->fusion1Past_ && !node->fissionPast_; }

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
