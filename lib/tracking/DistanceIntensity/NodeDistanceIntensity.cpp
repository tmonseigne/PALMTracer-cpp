/*******************************************************************************************************************
									NodeDistanceIntensity  Class
									It is derived from the class Nodes
								This class find the neighbor with the best likeness.
								Likeness depend on the distance and the intensity.
*******************************************************************************************************************/


#include <cstdio>
#include <cstdlib>

#include "NodeDistanceIntensity.h"

constexpr double WEIGHT_DISTANCE = 0.6;
constexpr double WEIGHT_INTENSITY = 0.2;
constexpr double WEIGHT_MEAN_INTENSITY = 0.2;
constexpr double MIN_ENERGY_FUSION = 0;
/********************************************************************************************************************
bool metropolis(double deltaE, double temperature) 
Metropolis algorithm. metropolis returns a boolean variable that issues a verdict on whether
to accept a reconfiguration that leads to a change deltaE in the objective function. If deltaE<0,
metropolis = true, while if deltaE>0, metropolis is only true with probability exp(-deltaE/temperature), where
temperature is a temperature determined by the annealing schedule.
********************************************************************************************************************/
bool NodeDistanceIntensity::metropolis(const double deltaE, const double temperature)
{
	if (deltaE >= 0.) { return true; }
	if ((double)(rand() + 1) / (double)RAND_MAX > exp(deltaE / temperature)) { return false; }
	else { return true; }
}


/********************************************************************************************************************
Constructuor(int _time) of NodeDistanceIntensity class.
********************************************************************************************************************/
NodeDistanceIntensity::NodeDistanceIntensity(const int time, Parameters* params): Node(time, params) { }

/********************************************************************************************************************
double getLikeness(double difference, double sigma)
compute the likeness for a difference of value compared to the sigma
likeness=1 if difference=0
likeness=0 if difference=3*sigma
********************************************************************************************************************/
double NodeDistanceIntensity::getLikeness(const double difference, const double sigma) { return 1. - fabs(difference) / (3. * sigma); }


/*********************************************************************************************
double NodeDistanceIntensity::likenessDistance(Node* n)
  Compute the likeness from the distance between this and n.
  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::LikenessDistance(Node* n)
{
	// attention je rajoute le *3 pour les données d'aqi(autoquant)
	if (n == nullptr) { return getLikeness(__min(DistanceToBorder(), parameters_->stdDistance_), parameters_->stdDistance_); }
	return getLikeness(Distance(n), parameters_->stdDistance_);
}


/*********************************************************************************************
double NodeDistanceIntensity::stericEffect(Node* n1, Node* n2)
  Compute the sum of the minimum radius of the nodes n1 and n2
  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::stericEffect(const Node* n1, const Node* n2) { return n1->SphericEffect() + n2->SphericEffect(); }

/*********************************************************************************************
double NodeDistanceIntensity::likenessDistanceFusion(Node* n1, Node* n2)
  Compute the likeness from the distance of a fusion between this and n1 and n2.
  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::LikenessDistanceFusion(Node* n1, Node* n2)
{
	if (n1 == nullptr) { return 0.; }
	if (n2 == nullptr) { return 0.; }

	return getLikeness(fabs(Distance(n1, n2) + n1->Distance(n2) - stericEffect(n1, n2)), parameters_->stdDistance_);
}


/*********************************************************************************************
double NodeDistanceIntensity::likenessMeanIntensity(Node* n)
  Compute the likeness from the mean of intensity between this and n.
  New version from 2005/01/12

  We use the hypothesis that 
	sigma(mean intensity) proportional to surface^-1/6 simga_noise / Mean in 2D.
*********************************************************************************************/
double NodeDistanceIntensity::LikenessMeanIntensity(Node* n)
{
	if (n == nullptr) { return getLikeness(info_->mean_, parameters_->stdMean_ * __max(1., 3.*parameters_->stdDistance_/this->DistanceToBorder())); }

	else { return getLikeness(info_->mean_ - n->info_->mean_, parameters_->stdMean_); }
}

/*********************************************************************************************
double NodeDistanceIntensity::likenessMeanIntensityFusion(Node* n1, Node* n2)
  Compute the likeness from the mean of intensity in fusion between this and n1 and n2 .
  New version from 2005/01/12

*********************************************************************************************/
double NodeDistanceIntensity::LikenessMeanIntensityFusion(Node* n1, Node* n2)
{
	if (n1 == nullptr || n2 == nullptr) { return 0.; } //Error ???
	return getLikeness(info_->mean_ - (n1->info_->intensity_ + n2->info_->intensity_)
	                   / (n1->info_->surf_ + n2->info_->surf_), parameters_->stdMean_);
	//if (Parameters::N_DIM == 2) { return getLikeness(info_->mean_ - (n1->info_->intensity_ + n2->info_->intensity_) / (n1->info_->surf_ + n2->info_->surf_), parameters_->stdMean_);}
	//else { return getLikeness(info_->mean_ - (n1->info_->intensity_ + n2->info_->intensity_) / (n1->info_->GetVolume() + n2->info_->GetVolume()), parameters_->stdMean_);}
}


/*********************************************************************************************
double NodeDistanceIntensity::likenessIntensity(Node* n)
  Compute the likeness from the total intensity between this and n.
  New version from 2005/01/12

  We use the hypothesis that sigma(intensity) proportional to simga_noise/Mean.
*********************************************************************************************/
double NodeDistanceIntensity::LikenessIntensity(Node* n)
{
	//__max(1., 3.*parameters_->stdDistance/this->distanceToBorder())));
	if (n == nullptr) {
		return getLikeness(info_->intensity_, parameters_->stdIntensity_ * __max(1., 3.*parameters_->stdDistance_/this->DistanceToBorder()));
	}
	else { return getLikeness(info_->intensity_ - n->info_->intensity_, parameters_->stdIntensity_); }
}


/*********************************************************************************************
double NodeDistanceIntensity::likenessIntensityFusion(Node* n1, Node* n2
  Compute the likeness for the intensity in a fusion between this and n1 and n2.
  New version from 2005/01/12

  We use the hypothesis that sigma(intensity) proportional to simga_noise/Mean.
*********************************************************************************************/
double NodeDistanceIntensity::LikenessIntensityFusion(Node* n1, Node* n2)
{
	if (n1 == nullptr) { return 0.; } //inutil
	if (n2 == nullptr) { return 0.; } //inutil
	return getLikeness(info_->intensity_ - n1->info_->intensity_ - n2->info_->intensity_, parameters_->stdIntensity_);
}


/*********************************************************************************************
double NodeDistanceIntensity::likeness(Node* n)
  Compute the likeness between this and n.
  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::Likeness(Node* n)
{
	// the /2 is because the sum of the likeness of 2 unconnected nodes has to be compared to the likeness of the connection between those two nodes.
	if (n == nullptr) {
		return (WEIGHT_DISTANCE * LikenessDistance(n) + WEIGHT_INTENSITY * LikenessIntensity(n) + WEIGHT_MEAN_INTENSITY * LikenessMeanIntensity(n)) / 2. -
		       parameters_->costBirthDeath_;
	}


	return WEIGHT_DISTANCE * LikenessDistance(n) + WEIGHT_INTENSITY * LikenessIntensity(n) + WEIGHT_MEAN_INTENSITY * LikenessMeanIntensity(n);
}


/*********************************************************************************************
double NodeDistanceIntensity::likeness(Node* n)
  Compute the likeness between this and n.
  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::LikenessBlink(Node* n)
{
	constexpr double Alpha = 0.7;
	return Likeness(n) * Alpha + (Likeness(nullptr) + n->Likeness(nullptr)) * (1. - Alpha / 2.);
}


/*********************************************************************************************
double NodeDistanceIntensity::likenessFusion(Node* n1, Node* n2)
  Compute the likeness in fusion between this and n1 and n2.
  New version from 2005/01/12
*********************************************************************************************/

double NodeDistanceIntensity::LikenessFusion(Node* n1, Node* n2)
{
	if (n1 == nullptr) { return Likeness(n2); }
	if (n2 == nullptr) { return Likeness(n1); }
	return (WEIGHT_DISTANCE * LikenessDistanceFusion(n1, n2) + WEIGHT_INTENSITY * LikenessIntensityFusion(n1, n2) + WEIGHT_MEAN_INTENSITY *
	        LikenessMeanIntensityFusion(n1, n2) + MIN_ENERGY_FUSION) * 3. / 2.; // the *3/2 is because the sum of 
	// the likeness of 1 unconnected node and 2 linked nodes has to be compared to the likeness of the 
	// fusion connection between those 3 nodes.
}


/*********************************************************************************************
double NodeDistanceIntensity::likenessFission(Node* n1, Node* n2)
  Compute the likeness in fission between this and n1 and n2.
  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::LikenessFission(Node* n1, Node* n2) { return LikenessFusion(n1, n2); }


/*********************************************************************************************
void NodeDistanceIntensity::computeConstantSigmas(double &stdIntensityComputing, 
	double &stdMeanComputing, double &stdDistanceComputing, long &countStd)
  Find in the vicinity of each node the one with the biggest likeness to compute the constant 
  kInt and kMean and kDist
*********************************************************************************************/
void NodeDistanceIntensity::ComputeConstantSigmas(double& stdIntensityComputing, long& countStdInt, double& stdMeanComputing, long& countStdMean,
                                                  double& stdDistanceComputing, long& countStdDist, const bool sigmaClipping, long* hist)
{
	Node* nMax = nullptr;
	double likenessMax = -1e10;

	for (const ChainNeighbor* c = neighborsFuture_; c != nullptr; c = c->next_) {
		const double LikenessVal = Likeness(c->neighbor_);
		if (LikenessVal > likenessMax) {
			likenessMax = LikenessVal;
			nMax = c->neighbor_;
		}
	}


	if (nMax != nullptr) {
		if (hist != nullptr) { hist[(long)Distance(nMax)]++; }
		else {
			if (sigmaClipping) {
				if (fabs(info_->intensity_ - nMax->info_->intensity_) < 3. * parameters_->stdIntensity_) {
					stdIntensityComputing += (info_->intensity_ - nMax->info_->intensity_) * (info_->intensity_ - nMax->info_->intensity_);
					countStdInt++;
				}
				if (fabs(info_->mean_ - nMax->info_->mean_) < 3. * parameters_->stdMean_) {
					stdMeanComputing += (info_->mean_ - nMax->info_->mean_) * (info_->mean_ - nMax->info_->mean_);
					countStdMean++;
				}

				if (Distance(nMax) < 3. * parameters_->stdDistance_) {
					stdDistanceComputing += Distance(nMax) * Distance(nMax);
					countStdDist++;
				}
			}
			else {
				stdIntensityComputing += (info_->intensity_ - nMax->info_->intensity_) * (info_->intensity_ - nMax->info_->intensity_);
				stdMeanComputing += (info_->mean_ - nMax->info_->mean_) * (info_->mean_ - nMax->info_->mean_);
				stdDistanceComputing += Distance(nMax) * Distance(nMax);
				countStdInt++;
				countStdMean++;
				countStdDist++;
			}
		}
	}
}

/*********************************************************************************************
double NodeDistanceIntensity::makeInitializationConnection()
  Make the initial connections. For each strucutres, chose the one with the highest likeness 
  if it is avaliable. If it is not possible no connection.

  It returns the likeness.

  New version from 2005/01/12
*********************************************************************************************/
double NodeDistanceIntensity::MakeInitializationConnection()
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
		return max;
	}

	return Likeness(nullptr);
}


/*********************************************************************************************
double NodeDistanceIntensity::getLikenessOfLocalNode(int lastTime)
  Compute the likeness according to the connections. 

  It returns the likeness.

  New version from 2005/01/13
*********************************************************************************************/
double NodeDistanceIntensity::GetLikenessOfLocalNode(const int lastTime)
{
	double likeness = 0.;

	if (time_ > 0) {
		//to the past
		if (past_ == nullptr && fissionPast_ == nullptr && fusion1Past_ == nullptr) { likeness = Likeness(nullptr); }
		if (fusion1Past_ != nullptr) { likeness = +LikenessFusion(fusion1Past_, fusion2Past_); }
	}

	if (time_ < lastTime - 1) {
		//to the future
		if (future_ == nullptr && fission1Future_ == nullptr && fusionFuture_ == nullptr) { likeness += Likeness(nullptr); }
		if (future_ != nullptr) { likeness += Likeness(future_); }
		if (fission1Future_ != nullptr) { likeness += LikenessFission(fission1Future_, fission2Future_); }
	}
	return likeness;
}

double NodeDistanceIntensity::MakeChange(const double limitLoss)
{
	//choose an operation between:
	//simple linkage
	//strucutre fusion
	//strucutre dissociation

	if (parameters_->allowFusionAndFission_) {
		if (future_ != nullptr) {
			switch (rand() % 4) {
				case 0: return makeChangeFusion(limitLoss);
				case 1: return makeChangeFission(limitLoss);
				default: return makeChangeSimpleLinkage(limitLoss);
			}
		}

		if (Parameters::N_TEMPORAL_WINDOWS > 1 && blinkFuture_ != nullptr) { return makeChangeBlinkLinkage(limitLoss); }
		if (fission1Future_ != nullptr) { return makeChangeFission(limitLoss); }
		if (fusionFuture_ != nullptr) { return makeChangeFusion(limitLoss); }
		if (Parameters::N_TEMPORAL_WINDOWS > 1 && rand() % 2 == 0) { return makeChangeBlinkLinkage(limitLoss); }

		return makeChangeSimpleLinkage(limitLoss);
	}

	if (future_ != nullptr) { return makeChangeSimpleLinkage(limitLoss); }
	if (blinkFuture_ != nullptr) { return makeChangeBlinkLinkage(limitLoss); }
	if (Parameters::N_TEMPORAL_WINDOWS > 2 && rand() % 2 == 0) { return makeChangeBlinkLinkage(limitLoss); }

	return makeChangeSimpleLinkage(limitLoss);
}

double NodeDistanceIntensity::makeChangeFission0(const double limitLoss)
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


	ChainNeighbor* cRand;
	long i;

	if (blinkFuture_ != nullptr) { return -1; } //ERROR


	if (neighborsFutureCount_ < 2) { return 0.; } //no neighbors like 3
	const long IdRandNeighbor = rand() % neighborsFutureCount_;
	//set the node cRand
	for (i = 0, cRand = neighborsFuture_; i < IdRandNeighbor && cRand != nullptr; cRand = cRand->next_, i++) { ; }
	//test if the node is good
	if (cRand == nullptr) { return -1.; }                         //ERROR 
	if (cRand->neighbor_ == future_) { return 0.; }               // nothing to do, the new random link is the same than future_.
	if (cRand->neighbor_->fusion1Past_ != nullptr) { return 0.; } // crand is already implicated in a fusion.
	if (cRand->neighbor_->fissionPast_ != nullptr) { return 0.; } // crand is already implicated in a fission.
	if (cRand->neighbor_->blinkPast_ != nullptr) { return 0.; }   // crand is already implicated in a past blink.

	if (Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected

	const double DeltaE = +LikenessFission(cRand->neighbor_, future_) + (cRand->neighbor_->past_ == nullptr ? 0. : cRand->neighbor_->past_->Likeness(nullptr))
	                      - Likeness(future_) - cRand->neighbor_->Likeness(cRand->neighbor_->past_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
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
}


double NodeDistanceIntensity::makeChangeFission1(const double limitLoss, const ChainNeighbor* cRand)
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

	if (cRand->neighbor_->Distance(fission2Future_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected

	if (cRand->neighbor_->fission1Future_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (cRand->neighbor_->fusionFuture_ != nullptr) { return 0; }   // cRand is already implicated in a fusion.
	if (cRand->neighbor_->blinkFuture_ != nullptr) { return 0.; }   // crand is already implicated in a future blink.

	const double DeltaE = +Likeness(fission1Future_) + cRand->neighbor_->Likeness(fission2Future_) + (
		                      cRand->neighbor_->future_ == nullptr ? 0 : cRand->neighbor_->future_->Likeness(nullptr)) - LikenessFusion(
		                      fission1Future_, fission2Future_) - cRand->neighbor_->Likeness(cRand->neighbor_->future_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		future_ = fission1Future_;
		future_->past_ = this;
		if (cRand->neighbor_->future_ != nullptr) { cRand->neighbor_->future_->past_ = nullptr; }
		cRand->neighbor_->future_ = fission2Future_;
		fission2Future_->past_ = cRand->neighbor_;
		fission1Future_->fissionPast_ = nullptr;
		fission2Future_->fissionPast_ = nullptr;
		fission1Future_ = nullptr;
		fission2Future_ = nullptr;
		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeFission2(const double limitLoss, const ChainNeighbor* cRand)
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

	if (cRand->neighbor_->Distance(fission1Future_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected

	if (cRand->neighbor_->fission2Future_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (cRand->neighbor_->fusionFuture_ != nullptr) { return 0; }   // cRand is already implicated in a fusion.
	if (cRand->neighbor_->blinkFuture_ != nullptr) { return 0.; }   // crand is already implicated in a future blink.

	const double DeltaE = +Likeness(fission2Future_) + cRand->neighbor_->Likeness(fission1Future_) + (
		                      cRand->neighbor_->future_ == nullptr ? 0 : cRand->neighbor_->future_->Likeness(nullptr)) - LikenessFusion(
		                      fission2Future_, fission1Future_) - cRand->neighbor_->Likeness(cRand->neighbor_->future_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		future_ = fission2Future_;
		future_->past_ = this;
		if (cRand->neighbor_->future_ != nullptr) { cRand->neighbor_->future_->past_ = nullptr; }
		cRand->neighbor_->future_ = fission1Future_;
		fission1Future_->past_ = cRand->neighbor_;
		fission2Future_->fissionPast_ = nullptr;
		fission1Future_->fissionPast_ = nullptr;
		fission2Future_ = nullptr;
		fission1Future_ = nullptr;

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeFission3(const double limitLoss)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)

	//possiblilty 3
	//To
	//1(this)---------2(fission1Future_)
	//NULL------------3(fission2Future_)

	const double DeltaE = +Likeness(fission1Future_) + fission2Future_->Likeness(nullptr) - LikenessFission(fission1Future_, fission2Future_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration

		future_ = fission1Future_;
		future_->past_ = this;
		future_->fissionPast_ = nullptr;
		fission2Future_->fissionPast_ = nullptr;
		fission1Future_ = nullptr;
		fission2Future_ = nullptr;

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeFission4(const double limitLoss)
{
	//possiblilty 4
	//To
	//1(this)---------3(fission2Future_)
	//NULL------------2(fission1Future_)
	const double DeltaE = +Likeness(fission2Future_) + fission1Future_->Likeness(nullptr) - LikenessFission(fission1Future_, fission2Future_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		future_ = fission2Future_;
		future_->past_ = this;
		future_->fissionPast_ = nullptr;
		fission1Future_->fissionPast_ = nullptr;
		fission1Future_ = nullptr;
		fission2Future_ = nullptr;
		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeFission5(const double limitLoss, const ChainNeighbor* cRand)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5

	//possibility 5
	//1(this)---------2(fission1Future_)
	//           /----3(fission2Future_)
	//4(cRand)--/-----5

	if (cRand->neighbor_->Distance(fission2Future_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected
	if (cRand->neighbor_->blinkFuture_ != nullptr) { return 0.; }                             // crand is already implicated in a future blink.

	const double DeltaE = +Likeness(fission1Future_) + cRand->neighbor_->LikenessFusion(fission2Future_, cRand->neighbor_->future_) -
	                      LikenessFusion(fission1Future_, fission2Future_) - cRand->neighbor_->Likeness(cRand->neighbor_->future_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration

		future_ = fission1Future_;
		future_->past_ = this;
		cRand->neighbor_->fission1Future_ = fission2Future_;
		cRand->neighbor_->fission1Future_->fissionPast_ = cRand->neighbor_;
		cRand->neighbor_->fission2Future_ = cRand->neighbor_->future_;
		cRand->neighbor_->fission2Future_->fissionPast_ = cRand->neighbor_;

		fission1Future_->fissionPast_ = nullptr;
		fission1Future_ = nullptr;
		fission2Future_ = nullptr;
		cRand->neighbor_->future_->past_ = nullptr;
		cRand->neighbor_->future_ = nullptr;

		return DeltaE;
	}
}


double NodeDistanceIntensity::makeChangeFission8(const double limitLoss, const ChainNeighbor* cRand)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5

	//possibility 8
	//1(this)---------3(fission2Future_)
	//           /----2(fission1Future_)
	//4(cRand)--/-----5

	if (cRand->neighbor_->Distance(fission1Future_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected
	if (cRand->neighbor_->blinkFuture_ != nullptr) { return 0.; }                             // crand is already implicated in a future blink.

	const double DeltaE = +Likeness(fission2Future_) + cRand->neighbor_->LikenessFusion(fission1Future_, cRand->neighbor_->future_) -
	                      LikenessFusion(fission1Future_, fission2Future_) - cRand->neighbor_->Likeness(cRand->neighbor_->future_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration

		future_ = fission2Future_;
		future_->past_ = this;
		cRand->neighbor_->fission2Future_ = fission1Future_;
		cRand->neighbor_->fission2Future_->fissionPast_ = cRand->neighbor_;
		cRand->neighbor_->fission1Future_ = cRand->neighbor_->future_;
		cRand->neighbor_->fission1Future_->fissionPast_ = cRand->neighbor_;

		fission2Future_->fissionPast_ = nullptr;
		fission2Future_ = nullptr;
		fission1Future_ = nullptr;
		cRand->neighbor_->future_->past_ = nullptr;
		cRand->neighbor_->future_ = nullptr;

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeFission9(const double limitLoss, const ChainNeighbor* cRand)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5

	//possibility 8
	//1(this)---------5
	//           /----2(fission1Future_)
	//4(cRand)--/-----3(fission2Future_)
	if (cRand->neighbor_->blinkFuture_ != nullptr) { return 0.; } // crand is already implicated in a future blink.

	if (cRand->neighbor_->future_ != nullptr) {
		if (Distance(cRand->neighbor_->future_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected
	}

	if (cRand->neighbor_->Distance(fission1Future_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected
	if (cRand->neighbor_->Distance(fission2Future_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected

	if (cRand->neighbor_->fission1Future_ != nullptr) { return 0.; }
	if (cRand->neighbor_->fusionFuture_ != nullptr) { return 0.; }

	const double DeltaE = +Likeness(cRand->neighbor_->future_) + cRand->neighbor_->LikenessFusion(fission1Future_, fission2Future_) -
	                      LikenessFusion(fission1Future_, fission2Future_) - cRand->neighbor_->Likeness(cRand->neighbor_->future_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration


		cRand->neighbor_->fission1Future_ = fission1Future_;
		fission1Future_->fissionPast_ = cRand->neighbor_;
		cRand->neighbor_->fission2Future_ = fission2Future_;
		fission2Future_->fissionPast_ = cRand->neighbor_;
		future_ = cRand->neighbor_->future_;
		if (future_ != nullptr) { future_->past_ = this; }
		fission1Future_ = nullptr;
		fission2Future_ = nullptr;
		cRand->neighbor_->future_ = nullptr;
		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeFission6(const double limitLoss, const ChainNeighbor* cRand)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5(cRand->neighbor_->future_)

	//To
	//possibility 6
	//1(this)--\------2(fission1Future_)
	//          \-----5
	//4(cRand)--------3(fission2Future_)

	if (cRand->neighbor_->blinkFuture_ != nullptr) { return 0.; } // crand is already implicated in a future blink.

	if (cRand->neighbor_->future_ != nullptr) {
		if (Distance(cRand->neighbor_->future_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected
	}
	if (cRand->neighbor_->Distance(fission2Future_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected

	const double DeltaE = +LikenessFusion(fission1Future_, cRand->neighbor_->future_) + cRand->neighbor_->Likeness(fission2Future_) -
	                      LikenessFusion(fission1Future_, fission2Future_) - cRand->neighbor_->Likeness(cRand->neighbor_->future_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration

		fission2Future_->past_ = cRand->neighbor_;
		cRand->neighbor_->future_->fissionPast_ = this;
		Node* localN = cRand->neighbor_->future_;
		cRand->neighbor_->future_ = fission2Future_;
		fission2Future_->fissionPast_ = nullptr;
		fission2Future_ = localN;
		fission2Future_->past_ = nullptr;

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeFission7(const double limitLoss, const ChainNeighbor* cRand)
{
	//From
	//1(this)---------2(fission1Future_)
	//          \-----3(fission2Future_)
	//4(cRand)--------5(cRand->neighbor_->future_)

	//To
	//possibility 7
	//1(this)-\-------3(fission2Future_)
	//          \-----5
	//4(cRand)--------2(fission1Future_)

	if (cRand->neighbor_->blinkFuture_ != nullptr) { return 0.; } // crand is already implicated in a future blink.

	if (cRand->neighbor_->future_ != nullptr) {
		if (Distance(cRand->neighbor_->future_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected
	}
	if (cRand->neighbor_->Distance(fission1Future_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected


	const double DeltaE = +LikenessFusion(fission2Future_, cRand->neighbor_->future_) + cRand->neighbor_->Likeness(fission1Future_) -
	                      LikenessFusion(fission2Future_, fission1Future_) - cRand->neighbor_->Likeness(cRand->neighbor_->future_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration

		fission1Future_->past_ = cRand->neighbor_;
		cRand->neighbor_->future_->fissionPast_ = this;
		Node* localN = cRand->neighbor_->future_;
		cRand->neighbor_->future_ = fission1Future_;
		fission1Future_->fissionPast_ = nullptr;
		fission1Future_ = localN;
		fission1Future_->past_ = nullptr;

		return DeltaE;
	}
}


double NodeDistanceIntensity::makeChangeFission(const double limitLoss)
{
	//is the node already makes a Fission with an other structure
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

		return makeChangeFission0(limitLoss);
	}
	else {
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


		switch (rand() % 8) {
			case 0: return makeChangeFission3(limitLoss);
			case 1: return makeChangeFission4(limitLoss);
			default: ;
		}

		ChainNeighbor* cRand;
		long i;


		if (fission2Future_->neighborsPastCount_ < 2) { //no neighbors like "4"
			switch (rand() % 2) {
				case 0: return makeChangeFission3(limitLoss);
				case 1: return makeChangeFission4(limitLoss);
				default: ;
			}
		}

		//try to find a neighbor_ "4" in the neighborhood of "3"
		const long IdRandNeighbor = rand() % fission2Future_->neighborsPastCount_;
		//set the node cRand
		for (i = 0, cRand = fission2Future_->neighborsPast_; i < IdRandNeighbor && cRand != nullptr; cRand = cRand->next_, i++) { ; }
		//test if the node is good

		if (cRand->neighbor_ == this) { return 0; }                   // nothing to do, the new random link is the same than this.
		if (cRand->neighbor_->blinkFuture_ != nullptr) { return 0.; } // crand is already implicated in a future blink.

		if (cRand->neighbor_->future_ != nullptr) {
			switch (rand() % 7) {
				case 0: return makeChangeFission5(limitLoss, cRand);
				case 1: return makeChangeFission6(limitLoss, cRand);
				case 2: return makeChangeFission7(limitLoss, cRand);
				case 3: return makeChangeFission8(limitLoss, cRand);
				default: ;
			}
		}

		switch (rand() % 3) {
			case 0: return makeChangeFission1(limitLoss, cRand);
			case 1: return makeChangeFission2(limitLoss, cRand);
			default: return makeChangeFission9(limitLoss, cRand);
		}
	}
}

double NodeDistanceIntensity::makeChangeFusion0(const double limitLoss)
{
	ChainNeighbor* cRand;
	long i;

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
	const long IdRandNeighbor = rand() % future_->neighborsPastCount_;
	//set the node cRand
	for (i = 0, cRand = future_->neighborsPast_; i < IdRandNeighbor && cRand != nullptr; cRand = cRand->next_, i++) { ; }
	//test if the node is good
	if (cRand == nullptr) { return -1.; }                           //ERROR 
	if (cRand->neighbor_ == this) { return 0; }                     // nothing to do, the new random link is the same than this.
	if (cRand->neighbor_->fission1Future_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (cRand->neighbor_->fusionFuture_ != nullptr) { return 0; }   // cRand is already implicated in a fusion.
	if (cRand->neighbor_->blinkFuture_ != nullptr) { return 0; }    // cRand is already implicated in a blink connection

	if (future_->Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected

	const double DeltaE = +future_->LikenessFusion(cRand->neighbor_, this) + (cRand->neighbor_->future_ != nullptr
		                                                                          ? cRand->neighbor_->future_->Likeness(nullptr)
		                                                                          : 0.) - Likeness(future_) - cRand->neighbor_->Likeness(
		                      cRand->neighbor_->future_);

	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
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
}


double NodeDistanceIntensity::makeChangeFusion1(const double limitLoss, const ChainNeighbor* cRand)
{
	Node* node;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { node = fusionFuture_->fusion2Past_; }
	else { node = fusionFuture_->fusion1Past_; }
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


	if (cRand->neighbor_->fissionPast_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (cRand->neighbor_->fusion1Past_ != nullptr) { return 0; } // cRand is already implicated in a fusion.

	if (cRand->neighbor_->Distance(node) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected


	const double DeltaE = +Likeness(fusionFuture_) + node->Likeness(cRand->neighbor_) + (
		                      cRand->neighbor_->past_ == nullptr ? 0 : cRand->neighbor_->past_->Likeness(nullptr)) - fusionFuture_->LikenessFusion(node, this)
	                      - cRand->neighbor_->Likeness(cRand->neighbor_->past_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		fusionFuture_->past_ = this;
		fusionFuture_->fusion1Past_ = nullptr;
		fusionFuture_->fusion2Past_ = nullptr;
		node->future_ = cRand->neighbor_;
		node->fusionFuture_ = nullptr;
		future_ = fusionFuture_;
		fusionFuture_ = nullptr;

		if (cRand->neighbor_->past_ != nullptr) { cRand->neighbor_->past_->future_ = nullptr; }
		cRand->neighbor_->past_ = node;

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeFusion2(const double limitLoss, const ChainNeighbor* cRand)
{
	Node* node;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { node = fusionFuture_->fusion2Past_; }
	else { node = fusionFuture_->fusion1Past_; }

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

	if (cRand->neighbor_->fissionPast_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (cRand->neighbor_->fusion1Past_ != nullptr) { return 0; } // cRand is already implicated in a fusion.

	if (Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected


	const double DeltaE = +Likeness(cRand->neighbor_) + node->Likeness(fusionFuture_) + (
		                      cRand->neighbor_->past_ == nullptr ? 0 : cRand->neighbor_->past_->Likeness(nullptr)) - fusionFuture_->LikenessFusion(node, this)
	                      - cRand->neighbor_->Likeness(cRand->neighbor_->past_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		//unlink future_

		fusionFuture_->past_ = node;
		fusionFuture_->fusion1Past_ = nullptr;
		fusionFuture_->fusion2Past_ = nullptr;
		node->fusionFuture_ = nullptr;
		node->future_ = fusionFuture_;
		fusionFuture_ = nullptr;
		future_ = cRand->neighbor_;
		if (cRand->neighbor_->past_ != nullptr) { cRand->neighbor_->past_->future_ = nullptr; }
		cRand->neighbor_->past_ = this;

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeFusion3(const double limitLoss)
{
	Node* node;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { node = fusionFuture_->fusion2Past_; }
	else { node = fusionFuture_->fusion1Past_; }
	//From
	//1(this)---------2(fusionFuture_)
	//3(node)----/

	//possiblilty 3
	//To
	//1(this)---------2(fusionFuture_)
	//3(node)---------NULL

	const double DeltaE = +Likeness(fusionFuture_) + node->Likeness(nullptr) - fusionFuture_->LikenessFusion(node, this);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration

		fusionFuture_->past_ = this;
		fusionFuture_->fusion1Past_ = nullptr;
		fusionFuture_->fusion2Past_ = nullptr;
		future_ = fusionFuture_;
		fusionFuture_ = nullptr;
		node->fusionFuture_ = nullptr;
		if (node->future_ != nullptr) { return -1.; } //impossible inutile

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeFusion4(const double limitLoss)
{
	Node* node;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { node = fusionFuture_->fusion2Past_; }
	else { node = fusionFuture_->fusion1Past_; }
	//From
	//1(this)---------2(fusionFuture_)
	//3(node)----/

	//possiblilty 4
	//1(this)---------NULL
	//3(node)---------2(fusionFuture_)
	const double DeltaE = +node->Likeness(fusionFuture_) + Likeness(nullptr) - fusionFuture_->LikenessFusion(node, this);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		fusionFuture_->past_ = node;
		fusionFuture_->fusion1Past_ = nullptr;
		fusionFuture_->fusion2Past_ = nullptr;
		node->future_ = fusionFuture_;
		fusionFuture_ = nullptr;
		node->fusionFuture_ = nullptr;

		if (future_ != nullptr) { return -1.; } //impossible inutile
		return DeltaE;
	}
}


double NodeDistanceIntensity::makeChangeFusion5(const double limitLoss, const ChainNeighbor* cRand)
{
	Node* node;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { node = fusionFuture_->fusion2Past_; }
	else { node = fusionFuture_->fusion1Past_; }

	//From
	//1(this)-----/---2(fusionFuture_)
	//3(node)----/
	//5---------------4(cRand)


	//possibility 5
	//1(this)---------2(fusionFuture_)
	//3(node)----\
	//5-----------\---4(cRand)

	if (node->Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected

	if (cRand->neighbor_->fissionPast_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (cRand->neighbor_->fusion1Past_ != nullptr) { return 0; } // cRand is already implicated in a fusion.

	if (cRand->neighbor_->past_ == nullptr) { return 0; } // not possible to reconnection the fusion.

	const double DeltaE = +Likeness(fusionFuture_) + cRand->neighbor_->LikenessFusion(node, cRand->neighbor_->past_) - fusionFuture_->
	                      LikenessFusion(node, this) - cRand->neighbor_->Likeness(cRand->neighbor_->past_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration

		future_ = fusionFuture_;
		future_->past_ = this;

		node->fusionFuture_ = cRand->neighbor_;
		cRand->neighbor_->fusion1Past_ = node;

		cRand->neighbor_->past_->fusionFuture_ = cRand->neighbor_;
		cRand->neighbor_->fusion2Past_ = cRand->neighbor_->past_;

		future_->fusion1Past_ = nullptr;
		fusionFuture_ = nullptr;

		future_->fusion2Past_ = nullptr;

		cRand->neighbor_->past_->future_ = nullptr;
		cRand->neighbor_->past_ = nullptr;

		return DeltaE;
	}
}


double NodeDistanceIntensity::makeChangeFusion7(const double limitLoss, const ChainNeighbor* cRand)
{
	Node* node;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { node = fusionFuture_->fusion2Past_; }
	else { node = fusionFuture_->fusion1Past_; }

	//From
	//1(this)-----/---2(fusionFuture_)
	//3(node)----/
	//5---------------4(cRand)


	//possibility 7
	//3(node)---------2(fusionFuture_)
	//1(this)----\
	//5-----------\---4(cRand)

	if (cRand->neighbor_->fissionPast_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (cRand->neighbor_->fusion1Past_ != nullptr) { return 0; } // cRand is already implicated in a fusion.

	if (cRand->neighbor_->past_ == nullptr) { return 0; } // not possible to reconnection the fusion.

	if (Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected
	const double DeltaE = +node->Likeness(fusionFuture_) + cRand->neighbor_->LikenessFusion(this, cRand->neighbor_->past_) - fusionFuture_->
	                      LikenessFusion(node, this) - cRand->neighbor_->Likeness(cRand->neighbor_->past_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration

		node->future_ = node->fusionFuture_;
		node->future_->past_ = node;

		fusionFuture_ = cRand->neighbor_;
		cRand->neighbor_->fusion1Past_ = this;

		cRand->neighbor_->past_->fusionFuture_ = cRand->neighbor_;
		cRand->neighbor_->fusion2Past_ = cRand->neighbor_->past_;

		node->future_->fusion1Past_ = nullptr;
		node->fusionFuture_ = nullptr;

		node->future_->fusion2Past_ = nullptr;

		cRand->neighbor_->past_->future_ = nullptr;
		cRand->neighbor_->past_ = nullptr;

		return DeltaE;
	}
}


double NodeDistanceIntensity::makeChangeFusion8(const double limitLoss, const ChainNeighbor* cRand)
{
	Node* node;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { node = fusionFuture_->fusion2Past_; }
	else { node = fusionFuture_->fusion1Past_; }
	//From
	//1(this)---------2(fusionFuture_)
	//3(node)----/
	//5---------------4(cRand)

	//possibility 8
	//5---------------2(fusionFuture_)
	//1(this)----\
	//3(node)-----\---4(cRand)

	if (cRand->neighbor_->fissionPast_ != nullptr) { return 0; } // cRand is already implicated in a fission.
	if (cRand->neighbor_->fusion1Past_ != nullptr) { return 0; } // cRand is already implicated in a fusion.

	if (cRand->neighbor_->past_ == nullptr) { return 0; } // not possible to reconnection the fusion.

	if (Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; }                       //distance is not respected
	if (node->Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; }                 //distance is not respected
	if (cRand->neighbor_->past_->Distance(fusionFuture_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected

	const double DeltaE = +cRand->neighbor_->past_->Likeness(fusionFuture_) + cRand->neighbor_->LikenessFusion(this, node) - fusionFuture_->
	                      LikenessFusion(node, this) - cRand->neighbor_->Likeness(cRand->neighbor_->past_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		fusionFuture_->past_ = cRand->neighbor_->past_;
		fusionFuture_->past_->future_ = fusionFuture_;

		cRand->neighbor_->fusion1Past_ = fusionFuture_->fusion1Past_;
		cRand->neighbor_->fusion2Past_ = fusionFuture_->fusion2Past_;
		cRand->neighbor_->past_ = nullptr;
		//fusionFuture_->past_->future_=NULL;
		node->fusionFuture_ = cRand->neighbor_;
		fusionFuture_->fusion1Past_ = nullptr;
		fusionFuture_->fusion2Past_ = nullptr;
		fusionFuture_ = cRand->neighbor_;

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeFusion6(const double limitLoss, const ChainNeighbor* cRand)
{
	Node* node;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { node = fusionFuture_->fusion2Past_; }
	else { node = fusionFuture_->fusion1Past_; }

	//From
	//1(this)-----/---2(fusionFuture_)
	//3(node)----/
	//5---------------4(cRand)


	//possibility 6
	//1(this)-----/---2(fusionFuture_)
	//5----------/
	//3(node)---------4(cRand)

	if (node->Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected
	if (cRand->neighbor_->past_ != nullptr) {
		if (cRand->neighbor_->past_->Distance(fusionFuture_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected
	}


	const double DeltaE = +fusionFuture_->LikenessFusion(this, cRand->neighbor_->past_) + cRand->neighbor_->Likeness(node) - fusionFuture_->
	                      LikenessFusion(node, this) - cRand->neighbor_->Likeness(cRand->neighbor_->past_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration

		if (fusionFuture_->fusion1Past_ == this) {
			cRand->neighbor_->past_->fusionFuture_ = fusionFuture_;
			cRand->neighbor_->past_->future_ = nullptr;
			fusionFuture_->fusion2Past_ = cRand->neighbor_->past_;
			node->fusionFuture_ = nullptr;
			node->future_ = cRand->neighbor_;
			cRand->neighbor_->past_ = node;
		}
		else {
			cRand->neighbor_->past_->fusionFuture_ = fusionFuture_;
			cRand->neighbor_->past_->future_ = nullptr;
			fusionFuture_->fusion1Past_ = cRand->neighbor_->past_;
			node->fusionFuture_ = nullptr;
			node->future_ = cRand->neighbor_;
			cRand->neighbor_->past_ = node;
		}

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeFusion9(const double limitLoss, const ChainNeighbor* cRand)
{
	Node* node;

	//find "node"
	if (fusionFuture_->fusion1Past_ == this) { node = fusionFuture_->fusion2Past_; }
	else { node = fusionFuture_->fusion1Past_; }
	Node* node2;
	//find "node2"
	if (cRand->neighbor_->fissionPast_->fission1Future_ == cRand->neighbor_) { node2 = cRand->neighbor_->fissionPast_->fission2Future_; }
	else { node2 = cRand->neighbor_->fissionPast_->fission1Future_; }

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

	if (node->Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected


	const double DeltaE = +Likeness(fusionFuture_) + node->Likeness(cRand->neighbor_) + cRand->neighbor_->fissionPast_->Likeness(node2) - fusionFuture_->
	                      LikenessFusion(node, this) - cRand->neighbor_->fissionPast_->LikenessFission(cRand->neighbor_, node2);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss if the new solution is to bad change the configuration

	fusionFuture_->past_ = this;
	future_ = fusionFuture_;
	fusionFuture_->fusion1Past_ = nullptr;
	fusionFuture_->fusion2Past_ = nullptr;
	node->fusionFuture_ = nullptr;
	fusionFuture_ = nullptr;
	node->future_ = cRand->neighbor_;
	cRand->neighbor_->past_ = node;
	cRand->neighbor_->fissionPast_ = nullptr;
	node2->fissionPast_->fission1Future_ = nullptr;
	node2->fissionPast_->fission2Future_ = nullptr;
	node2->fissionPast_->future_ = node2;
	node2->past_ = node2->fissionPast_;
	node2->fissionPast_ = nullptr;

	return DeltaE;
}

double NodeDistanceIntensity::makeChangeFusion(const double limitLoss)
{
	//is the node already making a Fusion with an other structure
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

		return makeChangeFusion0(limitLoss);
	}
	else {
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

		switch (rand() % 8) {
			case 0: return makeChangeFusion3(limitLoss);
			case 1: return makeChangeFusion4(limitLoss);
			default: ;
		}


		ChainNeighbor* cRand;
		long i;

		//get the cRand
		if (neighborsFutureCount_ < 2) { //no neighbors like "4", no cRand
			switch (rand() % 2) {
				case 0: return makeChangeFusion3(limitLoss);
				case 1: return makeChangeFusion4(limitLoss);
				default: ;
			}
		}

		const long IdRandNeighbor = rand() % neighborsFutureCount_;
		//set the node cRand
		for (i = 0, cRand = neighborsFuture_; i < IdRandNeighbor && cRand != nullptr; cRand = cRand->next_, i++) { ; }
		//test if the node is good
		if (cRand == nullptr) { return -1.; } //ERROR

		if (cRand->neighbor_->fissionPast_ != nullptr) { return makeChangeFusion9(limitLoss, cRand); }
		if (cRand->neighbor_->blinkPast_ != nullptr) { return 0; } // cRand is already implicated in a blink connection

		if (cRand->neighbor_->past_ != nullptr) {
			switch (rand() % 7) {
				case 0: return makeChangeFusion5(limitLoss, cRand);
				case 1:
				case 2: return makeChangeFusion6(limitLoss, cRand);
				case 3: return makeChangeFusion7(limitLoss, cRand);
				case 4: return makeChangeFusion8(limitLoss, cRand);
				default: ;
			}
		}

		switch (rand() % 2) {
			case 0: return makeChangeFusion1(limitLoss, cRand);
			default: return makeChangeFusion2(limitLoss, cRand);
		}
	}
}

double NodeDistanceIntensity::makeChangeSimpleLinkage0(const double limitLoss)
{ //simple deconnection

	//From 
	//Present	Future
	//1(this)-----2(future_)

	//To
	//1(this)-----NULL
	//NULL--------2(future_)

	const double DeltaE = +Likeness(nullptr) + future_->Likeness(nullptr) - Likeness(future_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	return 0.;
}


double NodeDistanceIntensity::makeChangeSimpleLinkage1(const double limitLoss, const ChainNeighbor* cRand)
{
	//	From
	//	Present		Future
	//	1(this)-----2
	//	NULL--------3(cRand)

	//	To
	//	Present		Future
	//	1(this)-----3(cRand)
	//	NULL--------2

	if (cRand->neighbor_->blinkPast_ != nullptr) { return 0; } // cRand is already implicated in a blink connection

	if (Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected


	const double DeltaE = +Likeness(cRand->neighbor_) + (future_ != nullptr ? future_->Likeness(nullptr) : 0) - cRand->neighbor_->Likeness(nullptr) -
	                      Likeness(future_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		//unlink future_
		if (future_ != nullptr) { future_->past_ = nullptr; }
		future_ = cRand->neighbor_;
		future_->past_ = this;

		return DeltaE;
	}
	//return(0.);
}

double NodeDistanceIntensity::makeChangeSimpleLinkage2(const double limitLoss, const ChainNeighbor* cRand)
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

	if (Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected

	const double DeltaE = +Likeness(cRand->neighbor_) + (future_ == nullptr ? 0 : future_->Likeness(nullptr)) + (
		                      cRand->neighbor_->past_ == nullptr ? 0 : cRand->neighbor_->past_->Likeness(nullptr)) - Likeness(future_) - cRand->neighbor_->
	                      Likeness(cRand->neighbor_->past_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		//unlink future_


		if (future_ != nullptr) { future_->past_ = nullptr; }
		future_ = cRand->neighbor_;
		future_->past_->future_ = nullptr;
		future_->past_ = this;

		return DeltaE;
	}
	//return(0.);
}

double NodeDistanceIntensity::makeChangeSimpleLinkage3(const double limitLoss, const ChainNeighbor* cRand, const long count)
{
	//	From
	//	Present		Future
	//	1-----------2
	//	4-----------3(cRand)
	//	NULL--------5(c)

	//find a random one
	const long IdRandNeighbor = rand() % count;
	int iCount = 0;

	const ChainNeighbor* c = cRand->neighbor_->past_->neighborsFuture_;
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

	if (Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected
	if (cRand->neighbor_->blinkPast_ != nullptr) { return 0; }               // cRand is already implicated in a blink connection
	if (c->neighbor_->blinkPast_ != nullptr) { return 0; }                   // c is already implicated in a blink connection
	if (cRand->neighbor_->past_ != nullptr) {
		if (cRand->neighbor_->past_->Distance(c->neighbor_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected
	}

	const double DeltaE = +Likeness(cRand->neighbor_) + cRand->neighbor_->past_->Likeness(c->neighbor_) + (future_ == nullptr ? 0 : future_->Likeness(nullptr))
	                      - Likeness(future_) - cRand->neighbor_->Likeness(cRand->neighbor_->past_) - c->neighbor_->Likeness(nullptr);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		//unlink future_

		if (future_ != nullptr) {
			future_->past_ = nullptr; //????
		}
		future_ = cRand->neighbor_;
		future_->past_->future_ = c->neighbor_;
		future_->past_->future_->past_ = future_->past_;
		future_->past_ = this;

		return DeltaE;
	}
	//return(0.);
}

double NodeDistanceIntensity::makeChangeSimpleLinkage4(const double limitLoss, const ChainNeighbor* cRand)
{
	//	From
	//	Present		Future
	//	1(this)-----2
	//	4-----------3(cRand)
	//	To
	//	Present		Future
	//	1-----------3(cRand)
	//	4-----------2


	if (cRand->neighbor_->past_ != nullptr) {
		if (cRand->neighbor_->past_->Distance(future_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected
	}
	if (Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0.; } //distance is not respected

	const double DeltaE = +Likeness(cRand->neighbor_) + future_->Likeness(cRand->neighbor_->past_) - cRand->neighbor_->Likeness(cRand->neighbor_->past_)
	                      - Likeness(future_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration

		cRand->neighbor_->past_->future_ = future_;
		future_->past_ = cRand->neighbor_->past_;
		future_ = cRand->neighbor_;
		cRand->neighbor_->past_ = this;

		return DeltaE;
	}
	//return(0.);
}


double NodeDistanceIntensity::makeChangeBlinkLinkage0(const double limitLoss)
{
	//transform : 
	//	From
	//	t			t+1		t+2
	//	1(this)-------------NULL
	//	NULL----------------2


	//	To
	//	t			t+1		t+2
	//	1(this)-------------2(blinkFuture_)

	const long Dt = rand() % (Parameters::N_TEMPORAL_WINDOWS - 1) + 2;
	long nbNodeAtGoodTime = 0;
	ChainNeighbor* cRand;

	for (cRand = neighborsBlinkFuture_; cRand != nullptr; cRand = cRand->next_) {
		if (cRand->neighbor_->time_ - time_ == Dt && cRand->neighbor_->past_ == nullptr && cRand->neighbor_->blinkPast_ == nullptr && cRand->neighbor_->
		    fissionPast_
		    == nullptr && cRand->neighbor_->fusion1Past_ == nullptr) { nbNodeAtGoodTime++; }
	}
	if (nbNodeAtGoodTime == 0) { return 0; }

	const long IdRandNeighbor = rand() % nbNodeAtGoodTime;
	long i = 0;
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

	const double DeltaE = +LikenessBlink(cRand->neighbor_) - cRand->neighbor_->Likeness(nullptr) - Likeness(nullptr);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } //deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		cRand->neighbor_->blinkPast_ = this;
		blinkFuture_ = cRand->neighbor_;

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeBlinkLinkage1(const double limitLoss)
{
	//transform : 
	//	From
	//	t			t+1		t+2
	//	1(this)-------------2(blinkFuture_)

	//	To
	//	t			t+1		t+2
	//	1(this)-------------NULL
	//	NULL----------------2
	const double DeltaE = +blinkFuture_->Likeness(nullptr) + Likeness(nullptr) - LikenessBlink(blinkFuture_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } // deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration


		blinkFuture_->blinkPast_ = nullptr;
		blinkFuture_ = nullptr;

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeBlinkLinkage2(const double limitLoss, const ChainNeighbor* cRand)
{
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


	const double DeltaE = +Likeness(cRand->neighbor_) + cRand->neighbor_->past_->LikenessBlink(blinkFuture_) - LikenessBlink(blinkFuture_) - cRand->
	                      neighbor_->past_->Likeness(cRand->neighbor_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } // deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		blinkFuture_->blinkPast_ = cRand->neighbor_->past_;
		cRand->neighbor_->past_->blinkFuture_ = blinkFuture_;
		cRand->neighbor_->past_->future_ = nullptr;
		cRand->neighbor_->past_ = this;
		future_ = cRand->neighbor_;
		blinkFuture_ = nullptr;

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeBlinkLinkage3(const double limitLoss, const ChainNeighbor* cRand)
{
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

	const double DeltaE = +LikenessBlink(cRand->neighbor_->future_) + cRand->neighbor_->Likeness(blinkFuture_) - LikenessBlink(blinkFuture_) - cRand->
	                      neighbor_->future_->Likeness(cRand->neighbor_);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } // deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration

		cRand->neighbor_->future_->blinkPast_ = this;
		cRand->neighbor_->future_->past_ = nullptr;
		blinkFuture_->blinkPast_ = nullptr;
		blinkFuture_->past_ = cRand->neighbor_;
		Node* tmp = cRand->neighbor_->future_;
		cRand->neighbor_->future_ = blinkFuture_;
		blinkFuture_ = tmp;

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeBlinkLinkage4(const double limitLoss, const ChainNeighbor* cRand)
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

	const double DeltaE = +Likeness(cRand->neighbor_) + blinkFuture_->Likeness(nullptr) - LikenessBlink(blinkFuture_) - cRand->neighbor_->Likeness(nullptr);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } // deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration
		blinkFuture_->blinkPast_ = nullptr;
		blinkFuture_ = nullptr;
		cRand->neighbor_->past_ = this;
		future_ = cRand->neighbor_;

		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeBlinkLinkage5(const double limitLoss, const ChainNeighbor* cRand)
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


	const double DeltaE = +Likeness(nullptr) + cRand->neighbor_->Likeness(blinkFuture_) - LikenessBlink(blinkFuture_) - cRand->neighbor_->Likeness(nullptr);
	if (!metropolis(DeltaE, limitLoss)) { return 0.; } // deltaLikeness<limitLoss) // if the new solution is to bad
	else {                                             // change the configuration


		blinkFuture_->blinkPast_ = nullptr;
		blinkFuture_->past_ = cRand->neighbor_;
		cRand->neighbor_->future_ = blinkFuture_;
		blinkFuture_ = nullptr;
		return DeltaE;
	}
}

double NodeDistanceIntensity::makeChangeBlinkLinkage(const double limitLoss)
{
	int i;
	if (blinkFuture_ == nullptr) { //nothing connected in the future
		//transform : 
		//	From
		//	t			t+1		t+n
		//	1(this)-------------NULL
		//	NULL----------------2


		//	To
		//	t			t+1		t+n
		//	1(this)-------------2(blinkFuture_)
		return makeChangeBlinkLinkage0(limitLoss);
	}
	else { //blinkFuture_!=NULL
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

		if (rand() % 2) { return makeChangeBlinkLinkage1(limitLoss); }

		//find a potential cRand
		//cRand has to be either just after this or just before blinkFuture_
		if (rand() % 2) { //cRand just after this
			long nbNodeAtGoodTime = 0;
			long nbNodeAtGoodTime2 = 0;
			ChainNeighbor* cRand;

			long idRandNeighbor;

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
			if (nbNodeAtGoodTime == 0 && nbNodeAtGoodTime2 == 0) { return makeChangeBlinkLinkage1(limitLoss); }

			if ((rand() % 2 && nbNodeAtGoodTime != 0) || nbNodeAtGoodTime2 == 0) {
				idRandNeighbor = rand() % nbNodeAtGoodTime;
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

				if (Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0; } //cRand not linked to this!
				if (blinkFuture_->Distance(cRand->neighbor_->past_) >= parameters_->maxSpeed_) { return 0; }
				//blinkFuture_ not linked to cRand->neighbor_->past_
				return makeChangeBlinkLinkage2(limitLoss, cRand);
			}
			else {
				ChainNeighbor* cRand2;
				idRandNeighbor = rand() % nbNodeAtGoodTime2;
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

				if (Distance(cRand2->neighbor_) >= parameters_->maxSpeed_) { return 0; } //cRand2 not linked to this!

				return makeChangeBlinkLinkage4(limitLoss, cRand2);
			}
		}
		else { // cRand is just before blinkFuture_
			long nbNodeAtGoodTime = 0;
			long nbNodeAtGoodTime2 = 0;
			ChainNeighbor* cRand;

			long idRandNeighbor;

			for (cRand = blinkFuture_->neighborsPast_; cRand != nullptr; cRand = cRand->next_) {
				if (blinkFuture_->time_ - cRand->neighbor_->time_ == 1) {
					if (cRand->neighbor_->future_ != nullptr) { nbNodeAtGoodTime++; }
					else {
						if (cRand->neighbor_->fission1Future_ == nullptr && cRand->neighbor_->fusionFuture_ == nullptr && cRand->neighbor_->blinkFuture_ ==
						    nullptr) { nbNodeAtGoodTime2++; }
					}
				}
			}
			if (nbNodeAtGoodTime == 0 && nbNodeAtGoodTime2 == 0) { return makeChangeBlinkLinkage1(limitLoss); }

			if ((rand() % 2 && nbNodeAtGoodTime != 0) || nbNodeAtGoodTime2 == 0) {
				idRandNeighbor = rand() % nbNodeAtGoodTime;
				i = 0;
				bool exitFor = false;
				for (cRand = blinkFuture_->neighborsPast_; !exitFor && cRand != nullptr;) {
					if (blinkFuture_->time_ - cRand->neighbor_->time_ == 1 && cRand->neighbor_->future_ != nullptr) {
						if (i == idRandNeighbor) { exitFor = true; }
						i++;
					}
					if (!exitFor) { cRand = cRand->next_; }
				}
				if (cRand == nullptr) { return -1.; }                                                 //ERROR 
				if (blinkFuture_->Distance(cRand->neighbor_) >= parameters_->maxSpeed_) { return 0; } //cRand not linked to blinkFuture_!
				if (Distance(cRand->neighbor_->future_) >= parameters_->maxSpeed_) { return 0; }      //this  not linked to cRand->neighbor_->future_!


				return makeChangeBlinkLinkage3(limitLoss, cRand);
			}
			else {
				ChainNeighbor* cRand2;
				idRandNeighbor = rand() % nbNodeAtGoodTime2;
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
				if (cRand2 == nullptr) { return -1.; }                                                 //ERROR 
				if (blinkFuture_->Distance(cRand2->neighbor_) >= parameters_->maxSpeed_) { return 0; } //cRand2 not linked to blinkFuture_!
				return makeChangeBlinkLinkage5(limitLoss, cRand2);
			}
		}
	}
}

double NodeDistanceIntensity::makeChangeSimpleLinkage(const double limitLoss)
{
	int i;

	if (future_ == nullptr || rand() % 2 == 0) {
		//cRand is a random node in the vicinity of "this" in the future.

		//compute index of cRand
		if (neighborsFutureCount_ == 0) { return 0.; } //no neighbors
		const long IdRandNeighbor = rand() % neighborsFutureCount_;

		//set the node cRand
		ChainNeighbor* cRand;
		for (i = 0, cRand = neighborsFuture_; i < IdRandNeighbor && cRand != nullptr; cRand = cRand->next_, i++) { ; }
		//test if the node is good
		if (cRand == nullptr) { return -1.; } //ERROR 
		// the new random link is the same than the future_.
		if (cRand->neighbor_ == future_) { return makeChangeSimpleLinkage0(limitLoss); } //try simple deconnection


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

			return makeChangeSimpleLinkage1(limitLoss, cRand);
		}
		else {                                      //there is a conflict, the node cRand->neighbor_ is already connected to a node at the same time_ than this
			if (rand() % 2 && future_ != nullptr) { //just flip the 2 nodes
				//	From
				//	Present		Future
				//	1(this)-----2
				//	4-----------3(cRand)
				//	To
				//	Present		Future
				//	1-----------3(cRand)
				//	4-----------2
				return makeChangeSimpleLinkage4(limitLoss, cRand);
			}
			else { //find a new node to be connected
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
					return makeChangeSimpleLinkage2(limitLoss, cRand);
				}
				else { //there is several available nodes to relink the conflicted node
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
					return makeChangeSimpleLinkage3(limitLoss, cRand, count);
				}
			}
		}
	}
	else { return makeChangeSimpleLinkage0(limitLoss); }
}
