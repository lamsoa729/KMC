/**
 * @file kmc.cpp
 * @brief KMC object and speciallized call functions to implement
 * two stage kmc implementation with code base.
 */
#include "kmc.hpp"

#include <cmath>

/*! \brief Function to update the minimum distance of KMC object to rod
 * (distMinArr_), perpendicular distance from rod axis (distPerpArr_), as well
 * as the point located on the rod from rods center where
 * this minimum distance occurs (muArr_).
 *
 *
 * \param j_rod KMC index of rod used in muArr_ and distPerpArr_
 * \param &rod Rod data structure containing center location, length,
 *        and unit direction vector.
 * \return void, Updates muArr_[j_rod], distMinArr_[j_rod], and
 * distPerpArr_[j_rod] in KMC object.
 */
template <typename TRod>
void KMC<TRod>::UpdateRodDistArr(const int j_rod, const TRod &rod) {
    const double rLen = rod.length; // Vector of rod
    const double *rUVec = rod.direction;
    const double *rCenter = rod.pos;
    double rVec[3], // Rod length vector
        rMinus[3],  // Rod minus end position vector
        rPlus[3],   // Rod plus end position vector
        sepVec[3];  // Vector of rod center to protein
    for (int i = 0; i < 3; ++i) {
        rVec[i] = rLen * rUVec[i];
        rMinus[i] = rCenter[i] - (.5 * rVec[i]);
        rPlus[i] = rCenter[i] + (.5 * rVec[i]);
        sepVec[i] = pos_[i] - rCenter[i];
    }
    // the perpendicular distance & position from protein ends to rod
    double pointMin[3];
    distMinArr_[j_rod] = dist_point_seg(pos_, rMinus, rPlus, pointMin);
    // Closest point of end_pos along rod axis from rod center.
    double mu0 = dot3(sepVec, rUVec);
    muArr_[j_rod] = mu0;
    // Perpendicular distance away from rod axis
    distPerpArr_[j_rod] = sqrt(dot3(sepVec, sepVec) - SQR(mu0));
}

/*! \brief Calculate the probability of a head to bind to surrounding rods.
 *
 * \param rods Array of rod pointers
 * \param &uniqueFlagJ Reference to filter list making sure you do not over
 * count rods
 * \param bindFactor Binding factor of head to rods \return
 * void, Changes prob_tot_ variable of KMC this object
 */
template <typename TRod>
void KMC<TRod>::CalcTotProbsUS(const TRod *const *rods,
                               const std::vector<int> &uniqueFlagJ,
                               const double bindFactor) {
    prob_tot_ = 0;
    for (int j_rod = 0; j_rod < rods_probs_.size(); ++j_rod) {
        if (uniqueFlagJ[j_rod] > 0) {
            rods_probs_[j_rod] = CalcProbUS(j_rod, *(rods[j_rod]), bindFactor);
            prob_tot_ += rods_probs_[j_rod];
        }
    }
}

/*! \brief Calculate the probability of a head to bind to one rod.
 *
 * \param j_rod Index of rod object
 * \param &rod Reference to rod object
 * \param bindFactor Binding factor of head to rod
 * \return Probability of head binding to rod rod
 */
template <typename TRod>
double KMC<TRod>::CalcProbUS(const int j_rod, const TRod &rod,
                             const double bindFactor) {
    // Find and add shortest distance to DistPerp array and the associated
    // locations along rod.
    UpdateRodDistArr(j_rod, rod);
    double distMinSQR = SQR(distMinArr_[j_rod]);
    double r_cutSQR = SQR(r_cutoff_);

    // Find length of line that goes through binding radius
    if (distMinSQR < r_cutSQR) {
        // Closest point of end_pos along rod from rod center.
        double distPerpSQR = SQR(distPerpArr_[j_rod]);
        double tLen = rod.length;
        double mueff = ABS(muArr_[j_rod]);
        double a = sqrt(r_cutSQR - distPerpSQR);
        double max = mueff + a;
        if (max > 0.5 * tLen)
            max = 0.5 * tLen;
        else if (max < -0.5 * tLen)
            max = -0.5 * tLen;

        double min = mueff - a;
        if (min > 0.5 * tLen)
            min = 0.5 * tLen;
        else if (min < -0.5 * tLen)
            min = -0.5 * tLen;

        return bindFactor * (max - min) / (4.0 / 3.0 * M_PI * CUBE(r_cutoff_));
    } else {
        return 0;
    }
}

/*! \brief Calculate the probability of a head unbinding from a rod it
 * is attached too.
 *
 *  Head must be attached to a rod.
 *
 * \param unbindFactor
 * \return void, Changes the prob_tot_ variable of this object.
 */
template <typename TRod>
void KMC<TRod>::CalcProbSU(const double unbindFactor) {
    prob_tot_ = unbindFactor;
}

/*! \brief Calculate the total probability of an unbound head to bind to
 * surrounding rods when other head is attached to another rod.
 *
 *  One head must be bound and its position must be stored in this pos_
 * variable.
 *
 * \param rods Array of surrounding rod pointers
 * \param &uniqueFlagJ Reference to filter list making sure you do not over
 * count rods
 * \param k_spring Spring constant between connected heads
 * \param eqLen Equilibrium length of spring connecting heads
 * \param bindFactor Binding factor of head to rods
 * \return void, Changes tot_prob_ variable of this object
 */
template <typename TRod>
void KMC<TRod>::CalcTotProbsSD(const TRod *const *rods,
                               const std::vector<int> &uniqueFlagJ,
                               const int boundID, const double lambda,
                               const double kappa, const double beta,
                               const double restLen, const double bindFactor) {
    prob_tot_ = 0;
    for (int j_rod = 0; j_rod < rods_probs_.size(); ++j_rod) {
        if (uniqueFlagJ[j_rod] > 0 && rods[j_rod]->gid != boundID) {
            if (LUTablePtr_) {
                // rods_probs_[j_rod] = LUCalcProbSD(j_rod, *(rods[j_rod]),
                // kappa, eqLen, bindFactor);
                rods_probs_[j_rod] =
                    LUCalcProbSD(j_rod, *(rods[j_rod]), bindFactor);
            } else {
                rods_probs_[j_rod] =
                    CalcProbSD(j_rod, *(rods[j_rod]), lambda, kappa, beta,
                               restLen, bindFactor);
            }
            prob_tot_ += rods_probs_[j_rod];
        }
    }
    return;
}

/*! \brief Calculate the probability of a head to bind to one rod when one
 * head is already bound.
 *
 * \param j_rod Index of rod object
 * \param &rod Reference to rod object
 * \param lambda Load sensitivity of unbinding
 * \param kappa Spring constant of protein when stretched
 * \param beta Inverse product of Boltzmann's constant and temperature
 * \param restLen Length of protein when neither compressed or stretched
 * \param bindFactor Binding factor of head to rod
 * \return Probability of head binding to rod rod
 */
template <typename TRod>
double KMC<TRod>::CalcProbSD(const int j_rod, const TRod &rod,
                             const double lambda, const double kappa,
                             const double beta, const double restLen,
                             const double bindFactor) {
    // // Find and add shortest distance to DistPerp array and the associated
    // // locations along rod.
    UpdateRodDistArr(j_rod, rod);

    double distMinSQR = SQR(distMinArr_[j_rod]);
    double r_cutSQR = SQR(r_cutoff_);
    const double D = 2. * rod.radius;

    double result;
    if (r_cutSQR < distMinSQR)
        result = 0;
    else {
        // Non-dimensionalize all exponential factors
        const double distPerp =
            distPerpArr_[j_rod] / D; // Perpendicular distance to rod axis
        const double mu0 = muArr_[j_rod] / D;
        const double M = (1 - lambda) * .5 * kappa * beta * D * D;
        const double ell = restLen / D;
        // Limits of integration
        const double lim0 = -mu0 - 0.5 * (rod.length / D);
        const double lim1 = -mu0 + 0.5 * (rod.length / D);
        lims_[j_rod].first = lim0;
        lims_[j_rod].second = lim1;
        // It's integrating time!
        result = integral(distPerp, lim0, lim1, M, ell);
    }
    return bindFactor * result;
}

/*! \brief Calculate the probability of a head to bind to one rod when
 * one head is already bound.using a pre-made lookup table
 *
 * \param j_rod Index of rod object
 * \param &rod Reference to rod object
 * \param lambda Load sensitivity of unbinding
 * \param kappa Spring constant of protein when stretched
 * \param beta Inverse product of Boltzmann's constant and temperature
 * \param restLen Length of protein when neither compressed or stretched
 * \param bindFactor Binding factor of head to rod
 * \return Probability of head binding to rod rod
 */
template <typename TRod>
double KMC<TRod>::LUCalcProbSD(const int j_rod, const TRod &rod,
                               const double bindFactor) {
    if (LUTablePtr_ == nullptr) {
        std::cerr << " *** Error: Lookup table not initialized ***"
                  << std::endl;
        exit(1);
    }
    // Find and add shortest distance to DistPerp array and the associated
    // locations along rod.
    UpdateRodDistArr(j_rod, rod);

    double result;
    // Bypass probability calculation if protein is too far away
    if (SQR(r_cutoff_) < SQR(distMinArr_[j_rod]))
        result = 0;
    else {
        double mu0 = muArr_[j_rod];
        double distPerp = distPerpArr_[j_rod];   // Perpendicular distance to
                                                 // rod, changes if croslinker
                                                 // is past the rod end point.
        double lim0 = -mu0 - (0.5 * rod.length); // Lower limit of integral
        double lim1 = -mu0 + (0.5 * rod.length); // Upper limit of integral
        lims_[j_rod].first = lim0;
        lims_[j_rod].second = lim1;
        // Get value of integral from end to first bound site.
        //  Since lookup table function is odd you can simply negate the
        //  negative limits.
        // WARNING: LUT stores only half of the table [0,sMax] by symmetry
        // Only the difference between two bounds are needed here.
        // both shifted by a constant
        double term0 = LUTablePtr_->Lookup(distPerp, fabs(lim0)) *
                       ((lim0 < 0) ? -1.0 : 1.0);
        double term1 = LUTablePtr_->Lookup(distPerp, fabs(lim1)) *
                       ((lim1 < 0) ? -1.0 : 1.0);
        result = (term1 - term0);
    }
    return bindFactor * result;
}

/*! \brief Calculate probability of head unbinding from rod
 *
 *
 * \param unbindFactor Unbinding factor for head that is unbinding
 * \return void, Changes the prob_tot_ membr
 */
template <typename TRod>
void KMC<TRod>::CalcProbDS(const double unbindFactor) {
    prob_tot_ = unbindFactor;
    return;
}

/*! \brief Find which MT head will bind to from stored probabilities
 *
 * \param rods Array of rod object pointers
 * \param &bindPos Location along rod relative to rod center that head
 * binds
 * \param roll Uniformally generated random number from 0 to 1 \return
 * Return ID of bound MT
 */
template <typename TRod>
int KMC<TRod>::whichRodBindUS(const TRod *const *rods, double &bindPos,
                              double roll) {
    // assert probabilities are not zero
    double pos_roll = 0.0;

    int i = 0;
    for (auto prob : rods_probs_) {
        if ((pos_roll + prob) > roll) {
            // Use an old random number to get a new uniform random number.
            pos_roll = (roll - pos_roll) / prob;
            break;
        } else {
            pos_roll += prob;
            i++;
        }
    }
    if (i == rods_probs_.size()) {
        // Roll given was too large, CHECK KMC step functions
        return -1;
    }
    const double half_length = .5 * rods[i]->length; // half length of rod
    const double mu = muArr_[i]; // Closest position of particle along rod
                                 // from rod center
    const double distMin = distMinArr_[i];   // Distance from rod
    const double distPerp = distPerpArr_[i]; // Distance from rod
    double bind_range = sqrt(SQR(r_cutoff_) - SQR(distPerp));
    if (std::isnan(bind_range))
        bind_range = 0.0;

    // double diff = sqrt(SQR(distMin) - SQR(distPerp));
    double range_m,
        range_p; // Minus and plus binding limits with respect to rod center
    if ((ABS(mu) - (half_length)) > SMALL) { // Protein is beyond rod
        if (mu < 0) {                        // Protein is beyond the minus end
            range_m = -half_length;
            range_p = ((mu + bind_range) > half_length) ? half_length
                                                        : (mu + bind_range);
        } else { // Protein is beyond the plus end
            range_m = ((mu - bind_range) > half_length) ? half_length
                                                        : (mu - bind_range);
            range_p = half_length;
        }
    } else {
        range_m = ((mu - bind_range) < -half_length) ? -half_length
                                                     : (mu - bind_range);
        // Closest particle can bind to rod plus end
        range_p =
            ((mu + bind_range) > half_length) ? half_length : (mu + bind_range);
    }
    // Find location on rod between plus and minus range
    bindPos = (pos_roll * (range_p - range_m)) + range_m;
    assert(bindPos >= -(half_length + 10e-8) &&
           bindPos <= (half_length + 10e-8));
    return i;
}

/*! \brief Converts 3 uniformaly chosen random numbers into spherical
 *   coordinates where head will be released relative to original bound
 * position.
 *
 * \param R Radius of sphere that protein will unbind to
 * \param rollVec Uniformally generated random number 3 vector all from 0
 * to 1. \return Return new position vector of head once detatched.
 */
template <typename TRod>
void KMC<TRod>::whereUnbindSU(double R, double rollVec[3], double pos[3]) {
    double r = R * std::cbrt(rollVec[0]);
    double costheta = 2. * rollVec[1] - 1.;
    double theta = acos(costheta);
    double phi = 2. * M_PI * rollVec[2];
    pos[0] = r * sin(theta) * cos(phi) + pos_[0];
    pos[1] = r * sin(theta) * sin(phi) + pos_[1];
    pos[2] = r * costheta + pos_[2];
    return;
}

/*! \brief Bind unbound protein head to rod close by.
 *
 *
 * \param roll Uniform random number between 0 and 1.
 * \return Index number of rod that head was bound too.
 */
template <typename TRod>
int KMC<TRod>::whichRodBindSD(double &bindPos, double roll) {
    double pos_roll = 0.0;
    int i = 0; // Index of rods
    for (auto prob : rods_probs_) {
        if ((pos_roll + prob) > roll) { // Found rod to bind to
            // Use an old random number to get a new uniform random number.
            pos_roll = (roll - pos_roll) / prob;
            break;
        } else { // Keep searching for binding rod
            pos_roll += prob;
            i++;
        }
    }
    if (i == rods_probs_.size()) {
        return -1;
    }
    bindPos = RandomBindPosSD(i, pos_roll);
    return i;
}

/*! \brief Where on a rod will a protein bind if already rod to another rod
 *
 * \param j_rod Parameter description
 * \param roll A uniformaly generated random number between 0-1
 * \return Return parameter description
 */
template <typename TRod>
double KMC<TRod>::RandomBindPosSD(int j_rod, double roll) {
    if (LUTablePtr_ == nullptr) {
        std::cerr << " *** Error: Lookup table pointer not initialized ***"
                  << std::endl;
        exit(1);
    }
    // const double D = LUTablePtr_->getRodDiameter();
    // Limits set at the end of the rods
    double lim0 = lims_[j_rod].first;
    double lim1 = lims_[j_rod].second;
    // TODO: Clear up this
    // Lookup table parameter: perpendicular distance away from rod
    const double distPerp = distPerpArr_[j_rod];
    // Range of CDF based on position limits where protein can bind
    //      Since lookup table only stores positive values, take the absolute
    //      value of the position limit and then negate probability if the
    //      position limit was less than 0.
    double pLim0 =
        LUTablePtr_->Lookup(distPerp, fabs(lim0)) * ((lim0 < 0) ? -1.0 : 1.0);
    double pLim1 =
        LUTablePtr_->Lookup(distPerp, fabs(lim1)) * ((lim1 < 0) ? -1.0 : 1.0);
    // Scale random number to be within pLims
    roll = roll * (pLim1 - pLim0) + pLim0;
    // double sbound = 0; // Reset x position parameter
    double preFact = roll < 0 ? -1. : 1.;
    roll *= preFact; // entry in lookup table must be positive
    return preFact * (LUTablePtr_->ReverseLookup(distPerp, roll)) +
           muArr_[j_rod];
}

/*! \brief Where should head be located after unbinding while other head is
 *  is still bound.
 *
 *  At the moment, head just goes to back to location of other bound head.
 * This should not change the behavior of the KMC step presently. Might be
 * changed in the future when crosslinkers are longer.
 *
 * \param boundPos[] Parameter description
 * \param pos[] Parameter description
 * \return Return parameter description
 */
template <typename TRod>
void KMC<TRod>::whereUnbindDS(const double boundPos[], double pos[]) {
    for (int i = 0; i < 3; ++i) {
        pos[i] = boundPos[i];
    }
    return;
}

#include "kmc.tpp"

