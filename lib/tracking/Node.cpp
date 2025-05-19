/*******************************************************************************************************************
									Node Class
								  Node representing one structure
*******************************************************************************************************************/

/********************************************************************************************************************
Class Node contains the caracteristic off the node, and the chain of the other nodes (time after and time before)
links to the current objet.
********************************************************************************************************************/

#include "Node.hpp"

#include <algorithm>
#include <numbers>

constexpr double WEIGHT_DISTANCE = 0.6;
constexpr double WEIGHT_INTENSITY = 0.2;
constexpr double WEIGHT_MEAN_INTENSITY = 0.2;

//=================================================================
//==================== Constructor, Destructor ====================
//=================================================================
/********************************************************************************************************************
Constructuor(size_t time) of Node class.
The current objets is at the time "time".
********************************************************************************************************************/
Node::Node(const size_t time, Parameters* params)
	: x_(0), y_(0), z_(0), intensity_(0), mean_(0), surf_(0), color_(0), molId_(0),
	  params_(params), time_(time), future_(nullptr), past_(nullptr), trackedIndex_(0), colorIndex_{},
	  distanceToBorder_(0), likeness_(0), likenessDistance_(0), likenessIntensity_(0), likenessMean_(0) { }

/********************************************************************************************************************
Destructuor() of Node class.
Delete all the linkage objets.
********************************************************************************************************************/
Node::~Node()
{
	for (auto& n : neighbors_) { n.clear(); }
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

	for (Node* neighbor : neighbors_[1]) {
		if (!neighbor->past_) {
			const double Val = Likeness(neighbor);
			if (Val > max) {
				max = Val;
				nMax = neighbor;
			}
		}
	}

	if (nMax) {
		nMax->past_ = this;
		future_ = nMax;
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
	double max = -1e10;

	for (const Node* neighbor : neighbors_[1]) {
		const double Val = Likeness(neighbor);
		if (Val > max) {
			max = Val;
			nMax = neighbor;
		}
	}

	if (nMax) {
		intensity += (intensity_ - nMax->intensity_) * (intensity_ - nMax->intensity_);
		mean += (mean_ - nMax->mean_) * (mean_ - nMax->mean_);
		const double Dist = Distance(*nMax);
		distance += Dist * Dist;
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
		if (neighbors_[1].empty()) { return 0.; } //no neighbors
		Node* target = getRandomNeighbor(1);
		// the new random link is the same as future_,try simple deconnection
		if (target == future_) { return simpleLinkageDisconnection(limit); }

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
		for (const Node* neighbor : target->past_->neighbors_[1]) { if (isSimpleLinkageAvailable(neighbor)) { count++; } }

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

bool Node::isSimpleLinkageAvailable(const Node* node) const { return (!node->past_ || node == future_); }

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

double Node::simpleLinkageDoubleSwap(const double limit, Node* node, const size_t count)
{
	//	From
	//	Present		Future
	//	1-----------2
	//	4-----------3(cRand)
	//	NULL--------5(c)

	//find a random one
	const size_t IdRandNeighbor = params_->RandomModulo(count);

	const auto& Neighbors = node->past_->neighbors_[1];
	Node* selected = nullptr;
	size_t iCount = 0;
	for (auto& N : Neighbors) {
		if (isSimpleLinkageAvailable(N)) {
			if (iCount == IdRandNeighbor) {
				selected = N;
				break;
			}
			++iCount;
		}
	}

	//	To
	//	Present		Future
	//	1-----------3(cRand)
	//	4-----------5(c)
	//	NULL--------2
	//if (!c || !isSimpleLinkageAvailable(c->neighbor_)) { return -1; }                             //ERROR
	if (Distance(*node) >= params_->maxSpeed_) { return 0.; }                                 // Distance is not respected
	if (node->past_ && node->past_->Distance(*selected) >= params_->maxSpeed_) { return 0.; } // Distance is not respected

	const double DeltaE = Likeness(node) + node->past_->Likeness(selected) + (!future_ ? 0 : future_->Likeness())
	                      - Likeness(future_) - node->Likeness(node->past_) - selected->Likeness();
	if (!metropolis(DeltaE, limit)) { return 0.; } //deltaLikeness<limit // if the new solution is too bad
	// change the configuration
	//unlink future_

	if (future_) { future_->past_ = nullptr; } //????
	future_ = node;
	future_->past_->future_ = selected;
	future_->past_->future_->past_ = future_->past_;
	future_->past_ = this;

	return DeltaE;
	//return(0.);
}
